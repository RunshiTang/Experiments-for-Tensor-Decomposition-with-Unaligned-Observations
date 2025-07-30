# cd "//sscwin/dfsroot/Users/rtang56/Desktop/RKHS/Generalized_RKHS/simulation/poisson2/setting1"

import numpy as np
import json
import time
import scipy
import random
from random import choices



import pandas as pd

np.random.seed(0)
random.seed(0)

import sys 
import os

os.chdir("real_data/ECAM/compare/")


from helpers_v2 import RKHS_microbiome


# true loss -5.315734

if __name__ == "__main__":
    
    
    with open("microbiome.json") as f:
        microbiome = json.load(f)
    
    microbiome = list(microbiome.values())
    
    np.random.seed(0)
    random.seed(0)
    
    s1 = 15
    s2 = 20
    s3 = 10
    R = 10
    n_max = 50
    lambda_parameter = 0.00001
    CG = False
    sketch = True
    stop = True
    
    
    time1 = time.time()
    result_sketch = RKHS_microbiome(microbiome, sketch=True,
                                                        R = R, 
                                                        lambda_parameter = lambda_parameter, 
                                                        stop_steps = 2,
                                                        stop = stop, 
                                                        n_max = n_max, 
                                                        s1 = s1, s2 = s2, s3 = s3, 
                                                        test = True)
    time2 = time.time()
    print(time2 - time1)
    
    out_A = result_sketch[0][0]
    out_B = result_sketch[0][1]
    
    
    theta_new = result_sketch[0][2][:,np.newaxis]
    
    def Union(lst1, lst2):
        final_list = list(set(lst1) | set(lst2))
        return final_list
    def Bernoulli_kernel(x, y):
        k1_x = x-0.5
        k1_y = y-0.5
        k2_x = 0.5*(k1_x**2-1/12)
        k2_y = 0.5*(k1_y**2-1/12)
        xy = abs(np.outer(x,np.ones(len(y))) - np.outer(np.ones(len(x)), y))
        k4_xy = ((xy-0.5)**4-0.5 * (xy-0.5)**2 + 7/240)/24
        kern_xy = np.outer(k1_x, k1_y) + np.outer(k2_x, k2_y) - k4_xy + 1
        return kern_xy
    def RKHS_norm_square(theta_tmp):
        return (np.transpose(theta_tmp)@K@theta_tmp)[0][0]
    T = []

    for i in range(len(microbiome)):
        T = Union(T, np.array(microbiome[i][0]))
        
    T.sort() 
    K = Bernoulli_kernel(np.array(T), np.array(T))
    
    singular_values = []
    
    for r in range(R):
        theta_tmp = theta_new[r*int(np.shape(theta_new)[0]/R):(r+1)*int(np.shape(theta_new)[0]/R),:]
        theta_norm_tmp = RKHS_norm_square(theta_tmp)** 0.5
        
        singular_values.append(theta_norm_tmp)
        theta_tmp = theta_tmp/theta_norm_tmp
        theta_new[r*len(T):(r+1)*len(T),:] = theta_tmp
    
    out_theta = theta_new.reshape((int(np.shape(theta_new)[0]/R),R))
    
    pd.DataFrame(out_A).to_csv("out_A_R" +str(R)+ ".csv",header=False, index=False)
    pd.DataFrame(out_B).to_csv("out_B_R" +str(R)+ ".csv",header=False, index=False)
    pd.DataFrame(out_theta).to_csv("out_theta_R" +str(R)+ ".csv",header=False, index=False)
    pd.DataFrame(singular_values).to_csv("out_singular_values_R" +str(R)+ ".csv",header=False, index=False)
    #s1 100 s2 100 R 168.2128100608838
        
        
        
        




    
