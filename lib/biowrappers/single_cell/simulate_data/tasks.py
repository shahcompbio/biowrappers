'''
Created on Oct 26, 2015

@author: Andrew Roth
'''
from __future__ import division

import gzip
import numpy as np
import pandas as pd
import random
import scipy.stats as stats
import yaml

#=======================================================================================================================
# Genotype simulation
#=======================================================================================================================
def simulate_genotypes(config,
                       genotypes_file, 
                       tree_file):
    
    if (config['mutation_rate'] < 0) or (config['mutation_rate'] > 1):
        raise Exception('Mutation rate must be between 0 and 1.') 
    
    if config['seed'] is not None:
        random.seed(config['seed'])
        
        np.random.seed(config['seed'])
    
    normal = Clone(None, [0 for _ in range(config['num_sites'])], 0)

    sim = PhylogenySimulator(normal, 
                             mutation_rate=config['mutation_rate'], 
                             loh_rate=config['loh_rate'],
                             polyclonal=config['polyclonal'])
    
    sim.simulate()
    
    with open(tree_file, 'w') as fh:
        fh.write(sim.root.to_newick_string() + '\n')
    
    clones = sim.clones
    
    if sim.root not in clones:
        clones.append(sim.root)
    
    clones = sorted(clones, key=lambda x: x.name)
    
    genotypes = [x.genotype for x in clones]
    
    index = [x.name for x in clones]
    
    columns = range(0, config['num_sites'])
    
    df = pd.DataFrame(genotypes, index=index, columns=columns)
    
    with gzip.GzipFile(genotypes_file, 'w') as fh:
        df.to_csv(fh, index_label='clone_id', sep='\t')

class Clone(object):
    def __init__(self, ancestor, genotype, name):
        
        self.ancestor = ancestor
        
        self.genotype = genotype
        
        self.name = name
        
        self.children = []
    
    @property
    def is_leaf(self):

        return len(self.children) == 0
    
    @property
    def nodes(self):
        '''
        List of nodes in the tree rooted at self.
        '''
        
        yield self
        
        for descendent in self.descendents:
            yield descendent

    @property
    def descendents(self):
        '''
        List of descendent nodes in the tree rooted at self.
        '''
        
        for child in self.children:
            for node in child.nodes:
                yield node
    
    def _newick_str(self):
        
        if self.is_leaf:
            newick_str = str(self.name)
        
        else:
            newick_str = '(' + ','.join([child._newick_str() for child in self.children]) + ')' + str(self.name)
        
        return newick_str

    def to_newick_string(self):
        return '(' + self._newick_str() + ');'

class PhylogenySimulator(object):
    def __init__(self, root, loh_rate=0.1, mutation_rate=10, polyclonal=False):
        
        self.root = root
        
        self.num_events = len(root.genotype)
        
        self.unmutated_sites = range(self.num_events)
        
        self.clones = [root, ]
        
        self.loh_rate = loh_rate
        
        self.mutation_rate = mutation_rate
        
        self.clone_index = 1
        
        if not polyclonal:
            self.add_child()
            
            self.clones.remove(root)
        
    @property
    def can_add_child(self):
        
        return len(self.unmutated_sites) > 0
    
    def simulate(self):
        while self.can_add_child:
            self.add_child()
    
    def add_child(self):
        
        ancestor = random.choice(self.clones)
        
        while True:
            genotype = self._get_child_genotype(ancestor.genotype)
            
            if genotype != ancestor.genotype:
                break
            
        child = Clone(ancestor, genotype, self.clone_index)
        
        ancestor.children.append(child)
        
        self.clones.append(child)
        
        self.clone_index += 1
        
    def _get_child_genotype(self, ancestor_genotyper):
        
        child_genotype = list(ancestor_genotyper)
        
        print self.num_events * self.mutation_rate
        
        num_mutant_sites = np.random.poisson(self.num_events * self.mutation_rate)
        
        print num_mutant_sites
        
        num_mutant_sites = min(num_mutant_sites, len(self.unmutated_sites))
        
        mutant_sites = random.sample(self.unmutated_sites, num_mutant_sites)
        
        [self.unmutated_sites.remove(x) for x in mutant_sites]
        
        for site in mutant_sites:
            
            child_genotype[site] = 1
            
        for site, value in enumerate(child_genotype):
            
            if value != 1:
                continue
            
            u = random.random()
            
            if u < self.loh_rate:
                child_genotype[site] = 2
        
        return child_genotype
    
#=======================================================================================================================
# Simulate count data
#=======================================================================================================================
def simulate_count_data(config,
                        genotypes_file,
                        counts_file,
                        labels_file,
                        params_file):
    
    if config['seed'] is not None:
        np.random.seed(config['seed'])
    
    G = pd.read_csv(genotypes_file, compression='gzip', index_col='clone_id', sep='\t')
    
    het = pd.read_csv(config['het_positions_file'], compression='gzip', index_col='cell_id', sep='\t')
    
    hom = pd.read_csv(config['hom_positions_file'], compression='gzip', index_col='cell_id', sep='\t')
    
    data = _simulate_count_data(config['doublet_prob'],
                                G.values,
                                het,
                                hom,
                                config['num_data_points'])
    
    background_error_rate = hom.mean().mean()
    
    _write_counts_data(data['X'], counts_file, background_error_rate, config['mean_depth'])
    
    _write_cluster_labels(data, config['num_data_points'], labels_file)
    
    _write_params(data, params_file)
    
def _write_counts_data(X, out_file, background_error_rate, mean_depth):
    X = pd.DataFrame(X, columns=range(X.shape[1]), index=range(X.shape[0]))
    
    X = X.stack(dropna=False).reset_index()
    
    X.columns = 'well_id', 'event_id', 'event_value'
    
    X['well_type'] = 'nucleus'
    
    X['depth'] = np.random.poisson(mean_depth, size=X.shape[0])
    
    X['alt_counts'] = np.random.binomial(X['depth'], X['event_value'])
    
    X['ref_counts'] = X['depth'] - X['alt_counts']
    
    X['ref_p_value'] = stats.binom.sf(X['ref_counts'] - 1,
                                      X['depth'],
                                      background_error_rate)
    
    X['alt_p_value'] = stats.binom.sf(X['alt_counts'] - 1,
                                      X['depth'],
                                      background_error_rate)
    
    X = X[['well_type', 'well_id', 'event_id', 'event_value', 'ref_counts', 'alt_counts', 'ref_p_value', 'alt_p_value']]
    
    with gzip.GzipFile(out_file, 'w') as fh:
        X.to_csv(fh, index=False, sep='\t')

def _write_cluster_labels(data, num_data_points, out_file):
    Z_1 = np.array(data['Z'][1])
    
    Z_1 = pd.DataFrame(Z_1,
                       index=range(num_data_points),
                       columns=['cluster_1', 'cluster_2'])
     
    Z_0 = pd.DataFrame(data['Z'][0],
                       index=range(num_data_points),
                       columns=['cluster'])
    
    Y = data['Y']
    
    Y = pd.DataFrame(Y,
                     index=range(num_data_points),
                     columns=['is_doublet'])
    
    doublet_indices = (Y['is_doublet'] == 1)
    
    labels = Z_1
    
    labels.loc[~doublet_indices, 'cluster_1'] = Z_0.loc[~doublet_indices, 'cluster'].astype(int)
    
    labels.loc[~doublet_indices, 'cluster_2'] = pd.np.nan
    
    labels['is_doublet'] = Y['is_doublet']

    with gzip.GzipFile(out_file, 'w') as fh:
        labels.to_csv(fh, index_label='cell_id', sep='\t')

def _write_params(data, out_file):
    params = {}
    
    for g in ['A', 'AB', 'B']:
        key = '{0}_pos'.format(g)
        
        params[key] = data[key]
    
    params['pi'] = [float(x) for x in data['pi']]
    
    params['d'] = [float(x) for x in data['d']]
    
    with open(out_file, 'w') as fh:
        yaml.dump(params, fh, default_flow_style=False)

def _simulate_count_data(doublet_prob, G, het, hom, num_data_points):
    
    state_map = {
                 0 : [(0, 0)],
                 1 : [(0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1)],
                 2 : [(2, 2)]
                 }
    
    inverse_state_map = _get_inverse_state_map(state_map)
    
    K = G.shape[0]
    
    M = G.shape[1]
    
    N = num_data_points

    d = np.array([1 - doublet_prob, doublet_prob])
    
    e = {
         'A' : [hom.iloc[:, np.random.choice(range(hom.shape[1]))] for _ in range(M)],
         'B' : [1 - hom.iloc[:, np.random.choice(range(hom.shape[1]))] for _ in range(M)],
         'AB' : [het.iloc[:, np.random.choice(range(het.shape[1]))] for _ in range(M)]
         }
    
    pi = np.random.dirichlet(np.ones(K))

    Y = np.random.multinomial(1, d, size=N).argmax(axis=1)

    Z = [np.random.multinomial(1, pi, size=N).argmax(axis=1),
         zip(np.random.multinomial(1, pi, size=N).argmax(axis=1), np.random.multinomial(1, pi, size=N).argmax(axis=1))]

    X = []
    
    for n in range(N):
        z = Z[Y[n]][n]
        
        x_n = []
        
        if np.isscalar(z):
            for m, g in enumerate(G[z]):
                x_n.append(_get_alt_freq(e, g, m))

        else:
            for m, (g_1, g_2) in enumerate(zip(G[z[0]], G[z[1]])):
                g = inverse_state_map[(g_1, g_2)]
                
                x_n.append(_get_alt_freq(e, g, m))
     
        X.append(x_n)
    
    X = np.array(X).astype(float)
    
    return {'d' : d, 
            'pi' : pi, 
            'G' : G, 
            'X' : X, 
            'Y' : Y, 
            'Z' : Z, 
            'A_pos' : [x.name for x in e['A']],
            'AB_pos' : [x.name for x in e['AB']],
            'B_pos' : [x.name for x in e['B']]}

def _get_inverse_state_map(state_map):
    inverse_state_map = {}

    for s in state_map:
        for (u, v) in state_map[s]:
            inverse_state_map[(u, v)] = s
    
    return inverse_state_map

def _get_alt_freq(e, g, m):
    if g == 0:
        geno = 'A'
    
    elif g == 1:
        geno = 'AB'
    
    elif g == 2:
        geno = 'B'
    
    e = e[geno][m]
    
    e = e[~np.isnan(e)]
    
    return np.random.choice(e)

#=======================================================================================================================
# Simulate missing data
#=======================================================================================================================
def simulate_missing_data(in_file, out_file, probability_missing=0, seed=None):
    print in_file
    if seed is not None:
        np.random.seed(seed)
        
    df = pd.read_csv(in_file, compression='gzip', sep='\t')
    
    missing = (np.random.random(df.shape[0]) < probability_missing)
    
    df.loc[missing, ['ref_counts', 'alt_counts']] = 0
    
    df.loc[missing, ['ref_p_value', 'alt_p_value']] = 1
    
    with gzip.GzipFile(out_file, 'w') as fh:
        df.to_csv(fh, sep='\t')