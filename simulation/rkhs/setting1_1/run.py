
import os
os.chdir("Z:/experiment/simulation/rkhs/setting1_1/")

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
    n_max = 10
    R = 5
    lambda_parameter = 0.0001
    
    with open("../simulated_count7_expected_with_noise.json") as f:
        microbiome = json.load(f)
    

    time_list_list = []
    acc_list_list=[]
    
    lambda_parameter_list = [0.001, 0.0005, 0.0001, 0.00005, 0.00001]
    
    for lambda_parameter in lambda_parameter_list:
        print(lambda_parameter)
        acc_list = []
        time_list = []
        for i in range(10):
            print(i)
            result = RKHS_microbiome(microbiome, R = R, 
                                lambda_parameter = lambda_parameter, n_max = n_max, 
                                CG = False, sketch = False, stop_criterion = False)
            time_list.append(result[4])
            acc_list.append(result[2])
        time_list_list += [time_list]
        acc_list_list += [acc_list]
    
    
    
    for i in range(len(lambda_parameter_list)):
        tmp_array = np.array(time_list_list[i])
        df = pd.DataFrame(tmp_array)
        df.to_csv('time_lambda_'+str(lambda_parameter_list[i]) + '.csv')
        
        tmp_array = np.array(acc_list_list[i])
        df = pd.DataFrame(tmp_array)
        df.to_csv('acc_lambda_'+str(lambda_parameter_list[i]) + '.csv')










