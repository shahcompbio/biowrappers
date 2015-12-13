'''
Created on Nov 9, 2015

@author: Andrew Roth
'''
from .math import compute_e_log_dirichlet, log_space_normalise

import numpy as np

class StochaticVariationalInference(object):
    def __init__(
        self,
        gamma_prior,
        kappa_prior,
        X,
        batch_size=1000,
        delay=1,
        forgetting_rate=0.7,
        gamma=None,
        kappa=None):
        
        self.K = len(kappa_prior)
        
        self.N = X.shape[0]
        
        self.S = X.shape[1]
        
        self.D = X.shape[2]
        
        self.batch_size = batch_size
        
        self.delay = delay
        
        self.forgetting_rate = forgetting_rate
        
        self.epoch = 0
        
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
        
        self.X = X
    
    def _init_log_resp(self):
        shape = (self.batch_size, self.K)

        log_resp = np.random.random(shape)

        self.log_resp = log_space_normalise(log_resp)
    
    @property
    def batch_data(self):
        return self._batch_data
    
    @property
    def e_log_gamma(self):
        return compute_e_log_dirichlet(self.gamma)
        
    @property
    def e_log_pi(self):
        return compute_e_log_dirichlet(self.kappa)
    
    @property
    def gamma_data_term(self):
        return np.einsum(
            'nksd, nksd -> nksd',
            self.resp[:, :, np.newaxis, np.newaxis],
            self.batch_data[:, np.newaxis, :, :]
        )
    
    @property
    def kappa_data_term(self):
        return self.resp
    
    @property
    def resp(self):
        return np.exp(self.log_resp)
    
    @property
    def step_size(self):
        return (self.epoch + self.delay) ** (-1 * self.forgetting_rate)
    
    def fit(self, num_iters=100):
        for i in range(num_iters):
            self.epoch = i

            self._update_batch_data()
            
            self._update_resp()
                
            self._update_gamma()

            self._update_kappa()

            print i, len(np.unique(self.resp.argmax(axis=1))), self.step_size
    
    def _update_batch_data(self):
        self._batch_indices = np.random.choice(
            np.arange(self.N),
            replace=False,
            size=self.batch_size
        )
        
        self._batch_indices = np.sort(self._batch_indices)
        
        self._batch_data = np.array(self.X[self._batch_indices])
    
    def _update_gamma(self):
        new_dir = self.gamma_prior[np.newaxis, np.newaxis, :, :] + self.N * self.gamma_data_term
        
        new_dir = np.mean(new_dir, axis=0)
        
        self.gamma = (1 - self.step_size) * self.gamma + self.step_size * new_dir

    def _update_kappa(self):
        new_dir = self.kappa_prior[np.newaxis, :] + self.N * self.kappa_data_term
   
        new_dir = np.mean(new_dir, axis=0)
        
        self.kappa = (1 - self.step_size) * self.kappa + self.step_size * new_dir

    def _update_resp(self):
        log_resp = np.einsum(
            'nksd, nksd -> nk',
            self.e_log_gamma[np.newaxis, :, :, :],
            self.batch_data[:, np.newaxis, :, :])

        log_resp = self.e_log_pi + log_resp
        
        self.log_resp = log_space_normalise(log_resp)
