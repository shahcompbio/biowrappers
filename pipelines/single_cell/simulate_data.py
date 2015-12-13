from pipelines.io import make_parent_directory

import os
import pypeliner

from biowrappers.core import PypelinerFile

import biowrappers.single_cell.counts_to_matrix as counts_to_matrix
import biowrappers.single_cell.scg as scg
import biowrappers.single_cell.simulate_data as simulate_data

def main():
    het_positions_file = '/home/andrew/Documents/papers/single_cell_genotyper/supplemental/software/synthetic/het.tsv.gz'
    
    hom_positions_file = '/home/andrew/Documents/papers/single_cell_genotyper/supplemental/software/synthetic/hom.tsv.gz'
    
    template_scg_file = '/home/andrew/version_control/mercurial/biowrappers/lib/biowrappers/single_cell/scg/data/genotyper_template.yaml'
    
    out_dir = '/home/andrew/Desktop/pypeliner_test'
    
    num_counts_files = 2
    
    num_tree_files = 3
    
    tree_config = {'loh_rate' : 0.1, 
                   'mutation_rate' : 0.2, 
                   'num_sites' : 100, 
                   'polyclonal' : False}
    
    counts_config = {'doublet_prob' : 0,
                     'het_positions_file' : het_positions_file,
                     'hom_positions_file' : hom_positions_file, 
                     'mean_depth' : 100,
                     'num_data_points' : 100}
    
    scheduler = pypeliner.scheduler.Scheduler()
    
    sim_files = simulate_data.single_cell_data_simulation_pipeline(
        scheduler,
        counts_config,
        tree_config,
        out_dir,
        missing_prob=0.2,
        num_counts_files=num_counts_files,
        num_tree_files=num_tree_files)
    
    in_file = sim_files['missing_counts']
    
    out_file = PypelinerFile(os.path.join(out_dir, 'input', os.path.basename(in_file.path)), axes=in_file.axes)
    
    make_parent_directory(out_file.path)
    
    counts_to_matrix.convert_counts_to_matrix_pipeline(
        scheduler, 
        in_file, 
        out_file,
        ('event_id',), 
        'snv', 
        None, 
        ('well_id',))
    
    scg.single_cell_genotyper_analysis_pipeline(
        scheduler,
        out_file,
        template_scg_file)
    
    native_spec = '-V -q all.q -l mem_token={mem}G,mem_free={mem}G,h_vmem={mem}G'
    
    with pypeliner.execqueue.LocalJobQueue([simulate_data.tasks]) as exec_queue:
        scheduler.run(exec_queue)
#     
#     tasks = []
#     
#     num_trees = 100
#     
#     num_data_files = 10
#     
#     jm = ClusterJobManager(cleanup_log_files=True, log_dir='log', max_tries=3)
#         
#     try:
#         for clone_seed in range(num_trees):
#             for counts_seed in range(num_data_files):
#                 t = SingleCellDataSimulatorPipeline(job_manager=jm,
#                                                     
#                                                     clone_seed=clone_seed,
#                                                     counts_seed=counts_seed,
#                                                     missing_seed=0,
#                                                     doublet_prob=0,
#                                                     het_positions_file=het_positions_file,
#                                                     hom_positions_file=hom_positions_file,
#                                                     loh_rate=10,
#                                                     missing_prob=40,
#                                                     mutation_rate=20,
#                                                     num_data_points=100,
#                                                     num_sites=100,
#                                                     out_dir='test')
#                 
#                 tasks.append(t)
#         
#         luigi.build(tasks, local_scheduler=False, workers=100)
#     
#     finally:
#         jm.close()

if __name__ == '__main__':
    main()