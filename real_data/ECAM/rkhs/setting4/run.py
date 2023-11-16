# cd "//sscwin/dfsroot/Users/rtang56/Desktop/RKHS_new/RKHS/real_data/setting1"

import os
os.chdir("//sscwin/dfsroot/Users/rtang56/Desktop/RKHS_new/RKHS/real_data/setting4")

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
sys.path.append(os.path.abspath("//sscwin/dfsroot/Users/rtang56/Desktop/RKHS_new/RKHS"))
from helpers import RKHS_microbiome



if __name__ == "__main__":
    s2 = 20
    s3 = 10
    n_max = 10
    R = 4
    
    with open("//sscwin/dfsroot/Users/rtang56/Desktop/RKHS/data/microbiome.json") as f:
        microbiome = json.load(f)
        
    microbiome = list(microbiome.values())
    
    acc_list = []
    time_list = []
    
    for i in range(1):
        print(i)
        result = RKHS_microbiome(microbiome, R = R, 
                            lambda_parameter = 0.001, n_max = n_max, 
                            CG = False, sketch = False)
        time_list.append(result[4])
        acc_list.append(result[2])
    
    tmp_array = np.array(time_list)
    df = pd.DataFrame(tmp_array)
    df.to_csv('time_unsketched.csv')
    
    tmp_array = np.array(acc_list)
    df = pd.DataFrame(tmp_array)
    df.to_csv('acc_unsketched.csv')
    
    
    
    
    
    
    s1_list = [2,5,15,30]
    time_list_list = []
    acc_list_list = []
    
    
    for s1 in s1_list:
        acc_list = []
        time_list = []
        print(s1)
        for i in range(1):
            print(i)
            result_sketch = RKHS_microbiome(microbiome, R = 5, 
                                lambda_parameter = 0.001, n_max = n_max, 
                                CG = False, sketch = True, 
                                s1 = s1, s2 = s2, s3 = s3)
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
    




