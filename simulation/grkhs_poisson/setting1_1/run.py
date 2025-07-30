

import os
os.chdir("Z:/experiment/simulation/grkhs_poisson/setting1")

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
sys.path.append(os.path.abspath(".."))
from helper_counts_v2 import RKHS_generalized_sketch_microbiome_counts
from helper_counts_v2 import RKHS_generalized_microbiome_counts

#-20.78248 true loss

if __name__ == "__main__":
    
    np.random.seed(1)
    random.seed(1)
    bias_parameter = 1**(-10)
    R = 5
    
    n_max = 15
    bad_epochs_max = 50
    leanring_rate = 0.4
    iterations_per_epoch = 10
    decrease_lr = False
    gradient_clipping=True
    clipping_threshold = 0.3
    
    
    with open("../simulated_count8.json") as f:
        microbiome = json.load(f)
        
    time_list_list = []
    acc_list_list = []
    
    lambda_parameter_list = [1000,5000,10000,50000,100000]
    
    for lambda_parameter in lambda_parameter_list: 
        acc_list = []
        time_list = []
        print(lambda_parameter)
        for i in range(10):
            print(i)
            result = RKHS_generalized_microbiome_counts(microbiome, 
                                                R = R, 
                                                bad_epochs_max = bad_epochs_max,
                                                lambda_parameter = lambda_parameter, 
                                                leanring_rate = leanring_rate, 
                                                n_max = n_max, 
                                                test = True, 
                                                iterations_per_epoch = iterations_per_epoch, 
                                                gradient_clipping = gradient_clipping,
                                                clipping_threshold = clipping_threshold)
            time_list.append(result[4])
            acc_list.append(result[2])
            
        time_list_list += [time_list]
        acc_list_list += [acc_list]
        
        df = pd.DataFrame(time_list)
        df.to_csv('time_lambda_'+str(lambda_parameter) + '.csv')
        
        df = pd.DataFrame(acc_list)
        df.to_csv('acc_lambda_'+str(lambda_parameter) + '.csv')
    

    


    
    