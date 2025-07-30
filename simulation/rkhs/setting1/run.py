
import os
os.chdir("simulation/rkhs/setting1/")

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
sys.path.append(os.path.abspath("../"))
from helpers import RKHS_microbiome



if __name__ == "__main__":
    s2 = 40
    s3 = 10
    n_max = 10
    R = 5
    lambda_parameter = 0.0001
    
    with open("../simulated_count7_expected_with_noise.json") as f:
        microbiome = json.load(f)
    
    acc_list = []
    time_list = []
    
    for i in range(10):
        print(i)
        result = RKHS_microbiome(microbiome, R = R, 
                            lambda_parameter = lambda_parameter, n_max = n_max, 
                            CG = False, sketch = False, stop_criterion = False)
        time_list.append(result[4])
        acc_list.append(result[2])
    
    tmp_array = np.array(time_list)
    df = pd.DataFrame(tmp_array)
    df.to_csv('time_unsketched.csv')
    
    tmp_array = np.array(acc_list)
    df = pd.DataFrame(tmp_array)
    df.to_csv('acc_unsketched.csv')
    
    
    
    
    
    
    s1_list = [10,20,30]
    time_list_list = []
    acc_list_list = []
    
    
    for s1 in s1_list:
        acc_list = []
        time_list = []
        print(s1)
        for i in range(10):
            print(i)
            result_sketch = RKHS_microbiome(microbiome, R = R, 
                                lambda_parameter = lambda_parameter, n_max = n_max, 
                                CG = False, sketch = True, 
                                s1 = s1, s2 = s2, s3 = s3, stop_criterion = False)
            time_list.append(result_sketch[4])
            acc_list.append(result_sketch[2])
            
        time_list_list += [time_list]
        acc_list_list += [acc_list]
    
        
    
    for i in range(len(s1_list)):
        tmp_array = np.array(time_list_list[i])
        df = pd.DataFrame(tmp_array)
        df.to_csv('time_s1_'+str(s1_list[i]) + '.csv')
        
        tmp_array = np.array(acc_list_list[i])
        df = pd.DataFrame(tmp_array)
        df.to_csv('acc_s1_'+str(s1_list[i]) + '.csv')








