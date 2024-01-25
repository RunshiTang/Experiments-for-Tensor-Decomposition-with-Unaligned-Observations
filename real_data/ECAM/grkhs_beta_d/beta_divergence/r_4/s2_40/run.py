

import os
os.chdir("/home/r/rtang56/RKHS/Generalized_RKHS/real_data/beta_divergence/r_4/s2_40")

import numpy as np
import json
import time
import scipy
import random
from random import choices



import pandas as pd

np.random.seed(1)
random.seed(1)

import sys 
import os
sys.path.append(os.path.abspath("/home/r/rtang56/RKHS/Generalized_RKHS"))
from helper_betad import RKHS_generalized_sketch_microbiome_beta_divergnece
from helper_betad import RKHS_generalized_microbiome_beta_divergnece


if __name__ == "__main__":
    np.random.seed(1)
    random.seed(1)
    bias_parameter = 10**(-6)
    s2 = 40
    s3 = 10
    R = 4
    lambda_parameter = 10000
    n_max= 15
    bad_epochs_max = 50
    leanring_rate = 0.01
    iterations_per_epoch = 10
    decrease_lr = False
    gradient_clipping=True
    clipping_threshold = 1
    beta = -0.5 #this means beta = 3
    
    with open("/home/r/rtang56/RKHS/data/microbiome_relative_abundance.json") as f:
        microbiome = json.load(f)
    microbiome = list(microbiome.values())
    
        
    
    s1_list = [10,20,30]
    time_list_list = []
    acc_list_list = []
    
    for i in range(len(s1_list)): 
        s1 = s1_list[i]
        acc_list = []
        time_list = []
        print(s1)
        for i in range(10):
            print(i)
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
                                                                                beta = beta)
            time_list.append(result_sketch[4])
            acc_list.append(result_sketch[2])
            
        time_list_list += [time_list]
        acc_list_list += [acc_list]
        
        df = pd.DataFrame(time_list)
        df.to_csv('time_s1_'+str(s1) + '.csv')
        
        df = pd.DataFrame(acc_list)
        df.to_csv('acc_s1_'+str(s1) + '.csv')
    
    for i in range(len(s1_list)):
        tmp_array = np.array(time_list_list[i])
        df = pd.DataFrame(tmp_array)
        df.to_csv('time_s1_'+str(s1_list[i]) + '.csv')
        
        tmp_array = np.array(acc_list_list[i])
        df = pd.DataFrame(tmp_array)
        df.to_csv('acc_s1_'+str(s1_list[i]) + '.csv')
    
    
    
    acc_list = []
    time_list = []
    
    for i in range(10):
        print(i)
        result = RKHS_generalized_microbiome_beta_divergnece(microbiome, 
                                                            R = R, 
                                                            bias_parameter=bias_parameter,
                                                            bad_epochs_max = bad_epochs_max,
                                                            lambda_parameter = lambda_parameter, 
                                                            leanring_rate = leanring_rate, 
                                                            n_max = n_max, 
                                                            test = True, 
                                                            iterations_per_epoch = iterations_per_epoch,
                                                            gradient_clipping = gradient_clipping,
                                                            clipping_threshold = clipping_threshold,
                                                            beta = beta)
        time_list.append(result[4])
        acc_list.append(result[2])
    
    tmp_array = np.array(time_list)
    df = pd.DataFrame(tmp_array)
    df.to_csv('time_unsketched.csv')
    
    tmp_array = np.array(acc_list)
    df = pd.DataFrame(tmp_array)
    df.to_csv('acc_unsketched.csv')
    
    
    