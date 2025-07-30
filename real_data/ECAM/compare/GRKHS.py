
import os
os.chdir("Z:/experiment/real_data/ECAM/compare/")
import numpy as np
import json
import time
import scipy
import random
from random import choices



import pandas as pd

import sys 
import os

from helper_betad import RKHS_generalized_sketch_microbiome_beta_divergnece


if __name__ == "__main__":
    #np.random.seed(1)
    #random.seed(1)
    np.random.seed(123)
    random.seed(123)
    bias_parameter = 10**(-6)
    s1 = 20   
    s2 = 20
    s3 = 10
    R = 5
    lambda_parameter = 10000
    n_max= 15
    bad_epochs_max = 50
    leanring_rate = 0.01
    iterations_per_epoch = 20
    decrease_lr = False
    gradient_clipping=True
    clipping_threshold = 1
    beta = -0.5 #this means beta = 0.5
    
    with open("microbiome_relative_abundance.json") as f:
        microbiome = json.load(f)
    microbiome = list(microbiome.values())
    result_sketch = RKHS_generalized_sketch_microbiome_beta_divergnece(microbiome, 
                                                                        R = R, 
                                                                        bias_parameter=bias_parameter,
                                                                        bad_epochs_max = bad_epochs_max,
                                                                        lambda_parameter = lambda_parameter, 
                                                                        leanring_rate = leanring_rate, 
                                                                        n_max = n_max, 
                                                                        s1 = s1, s2 = s2, s3 = s3, 
                                                                        test = True, 
                                                                        iterations_per_epoch = iterations_per_epoch,
                                                                        decrease_lr = False, 
                                                                        gradient_clipping = gradient_clipping,
                                                                        clipping_threshold = clipping_threshold,
                                                                        beta = beta, 
                                                                        stop=True,
                                                                        stop_steps=2)
    
    out_A = result_sketch[0][0]
    out_B = result_sketch[0][1]
    
    
    theta_new = result_sketch[0][2]
    
    def Union(lst1, lst2):
        final_list = list(set(lst1) | set(lst2))
        return final_list
    
    def radial_kernel(x,y):
        if np.isscalar(x):
            x = np.array([x])
        if np.isscalar(y):
            y = np.array([y])
        nx = np.shape(x)[0]
        ny = np.shape(y)[0]
        result_matrix = np.zeros((nx,ny))
        for i in range(nx):
          for j in range(ny):
            result_matrix[i,j] = np.exp(-np.abs(x[i]-y[j]))
        return result_matrix
    
    def RKHS_norm_square(theta_tmp):
        return (np.transpose(theta_tmp)@K@theta_tmp)[0][0]
    T = []

    for i in range(len(microbiome)):
        T = Union(T, np.array(microbiome[i][0]))
        
    T.sort() 
    K = radial_kernel(np.array(T), np.array(T))
    
    singular_values = []
    
    for r in range(R):
        theta_tmp = theta_new[r*int(np.shape(theta_new)[0]/R):(r+1)*int(np.shape(theta_new)[0]/R),:]
        theta_norm_tmp = RKHS_norm_square(theta_tmp)** 0.5
        
        singular_values.append(theta_norm_tmp)
        theta_tmp = theta_tmp/theta_norm_tmp
        theta_new[r*len(T):(r+1)*len(T),:] = theta_tmp
    
    out_theta = theta_new.reshape((int(np.shape(theta_new)[0]/R),R))
    
    pd.DataFrame(out_A).to_csv("GRKHS/out_A_R" +str(R)+ ".csv",header=False, index=False)
    pd.DataFrame(out_B).to_csv("GRKHS/out_B_R" +str(R)+ ".csv",header=False, index=False)
    pd.DataFrame(out_theta).to_csv("GRKHS/out_theta_R" +str(R)+ ".csv",header=False, index=False)
    pd.DataFrame(singular_values).to_csv("GRKHS/out_singular_values_R" +str(R)+ ".csv",header=False, index=False)