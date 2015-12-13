'''
Created on Nov 9, 2015

@author: Andrew Roth
'''
import dask.array as da
import numpy as np

from .math import compute_e_log_dirichlet, compute_e_log_q_dirichlet, log_space_normalise_da as log_space_normalise, \
                  log_sum_exp_da as log_sum_exp, safe_multiply


class BatchVariationalInference(object):
    def __init__(self, gamma_prior, kappa_prior, X, gamma=None, kappa=None):
        self.K = len(kappa_prior)
        
        self.N = X.shape[0]
        
        self.S = X.shape[1]
        
        self.D = X.shape[2]
        
        self.gamma_prior = gamma_prior
        
        self.kappa_prior = kappa_prior
        
        if gamma is None:
            self.gamma = np.random.random((self.K, self.S, self.D)) * self.N
        
        else:
            self.gamma = gamma
        
        if kappa is None:
            self.kappa = np.random.random(self.K) * self.N
        
        else:
            self.kappa = kappa
        
        self.X = da.rechunk(X, self._data_chunks)
        
        self.lower_bound = [float('-inf')]
        
        self._debug_lower_bound = [float('-inf')]
                
        self.converged = False
        
        self._init_log_resp()
   
    def _init_log_resp(self):
        shape = (self.N, self.K)

        log_resp = da.random.random(shape, self._resp_chunks)

        log_resp = da.log(log_resp)

        self.log_resp = log_space_normalise(log_resp)

    @property
    def e_log_gamma(self):
        return compute_e_log_dirichlet(self.gamma)
        
    @property
    def e_log_pi(self):
        return compute_e_log_dirichlet(self.kappa)
    
    @property
    def gamma_data_term(self):
        return da.sum(
            self.resp[:, :, np.newaxis, np.newaxis] * self.X[:, np.newaxis, :, :],
            axis=0).compute()

    @property
    def kappa_data_term(self):
        return da.exp(log_sum_exp(self.log_resp, axis=0)).compute()
    
    @property
    def resp(self):
        return da.exp(self.log_resp)
    
    @property
    def _data_chunks(self):
        return (self._row_chunk_size, self.S, self.D)
    
    @property
    def _resp_chunks(self):
        return (self._row_chunk_size, self.K)
    
    @property
    def _row_chunk_size(self):
        return min(self.N, int(1e6) // self.K)
    
    def fit(self, convergence_tolerance=1e-4, debug=False, num_iters=100):
        for i in range(num_iters):
            self._update_resp()
            
            if debug:
                print 'resp', self._diff_lower_bound()
            
            self._update_gamma()
            
            if debug:
                print 'gamma', self._diff_lower_bound()
            
            self._update_kappa()
            
            if debug:
                print 'kappa', self._diff_lower_bound()

            self.lower_bound.append(self._compute_lower_bound())
             
            diff = (self.lower_bound[-1] - self.lower_bound[-2]) / np.abs(self.lower_bound[-1])
             
            print i, self.lower_bound[-1], diff
            
            if abs(diff) < convergence_tolerance:
                print 'Converged'
                
                self.converged = True
                
                break
            
            elif diff < 0:
                print 'Lower bound decreased'
                
                if not debug:
                    self.converged = False
                    
                    break
    
    def _update_gamma(self):
        self.gamma = self.gamma_prior[np.newaxis, :, :] + self.gamma_data_term

    def _update_kappa(self):
        self.kappa = self.kappa_prior + self.kappa_data_term

    def _update_resp(self):
        log_resp = da.sum(
            self.e_log_gamma[np.newaxis, :, :, :] * self.X[:, np.newaxis, :, :],
            axis=(2, 3)
        )
        
        log_resp = self.e_log_pi + log_resp
        
        self.log_resp = log_space_normalise(log_resp)

    def _compute_lower_bound(self):
        return self._compute_e_log_p() - self._compute_e_log_q()

    def _compute_e_log_q(self):
        return sum([
                self._compute_e_log_q_kappa(),
                self._compute_e_log_q_gamma(),
                self._compute_e_log_q_z()
                ])

    def _compute_e_log_q_kappa(self):
        return compute_e_log_q_dirichlet(self.kappa)

    def _compute_e_log_q_gamma(self):
        log_q = 0

        for k in range(self.K):
            for s in range(self.S):
                log_q += compute_e_log_q_dirichlet(self.gamma[k, s, :])

        return log_q
        
    def _compute_e_log_q_z(self):
        return da.sum(self.resp * self.log_resp).compute()

    def _compute_e_log_p(self):
        return sum([
                self._compute_e_log_p_kappa(),
                self._compute_e_log_p_gamma()
                ])

    def _compute_e_log_p_kappa(self):
        return np.sum(safe_multiply(self.e_log_pi, self.kappa_prior + self.kappa_data_term - 1))

    def _compute_e_log_p_gamma(self):
        return np.sum(safe_multiply(self.e_log_gamma, self.gamma_prior + self.gamma_data_term - 1))
    
    def _diff_lower_bound(self):
        self._debug_lower_bound.append(self._compute_lower_bound())
        
        diff = (self._debug_lower_bound[-1] - self._debug_lower_bound[-2]) / np.abs(self._debug_lower_bound[-1])
        
        if diff < 0:
            print 'Bound decreased',
        
        return diff
