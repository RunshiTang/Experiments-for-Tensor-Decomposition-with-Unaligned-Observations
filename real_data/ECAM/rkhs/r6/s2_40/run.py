# cd "//sscwin/dfsroot/Users/rtang56/Desktop/RKHS_new/RKHS/real_data/setting1"

import os
os.chdir("/home/r/rtang56/RKHS/RKHS/real_data/r6/s2_40")

import numpy as np
import json
import time
import scipy
import random
from random import choices



import pandas as pd

np.random.seed(123)
random.seed(123)

import sys 
import os
sys.path.append(os.path.abspath("/home/r/rtang56/RKHS/RKHS"))
from helpers import RKHS_microbiome



if __name__ == "__main__":
    s2 = 40
    s3 = 10
    n_max = 10
    R = 6
    
    with open("/home/r/rtang56/RKHS/data/microbiome.json") as f:
        microbiome = json.load(f)
        
    microbiome = list(microbiome.values())
    
    acc_list = []
    time_list = []
    
    for i in range(10):
        print(i)
        result = RKHS_microbiome(microbiome, R = R, 
                            lambda_parameter = 0.0001, n_max = n_max, 
                            CG = False, sketch = False)
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
            result_sketch = RKHS_microbiome(microbiome, R = 5, 
                                lambda_parameter = 0.0001, n_max = n_max, 
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
    














    time_unsketched = pd.read_csv('time_unsketched.csv', sep=",", header=0,index_col=0).to_numpy()
    acc_unsketched = pd.read_csv('acc_unsketched.csv', sep=",", header=0,index_col=0).to_numpy()
    
    time_unsketched_mean = np.mean(time_unsketched, axis=0)
    acc_unsketched_mean = np.mean(acc_unsketched, axis=0)
    
    import matplotlib.pyplot as plt
    color = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    fig, ax = plt.subplots()
    
    
    line1, = ax.plot(time_unsketched_mean, acc_unsketched_mean, 
                     label ="RKHS Decomposition", linestyle='solid', color=color[0], linewidth=0)
    for i in range(np.shape(time_unsketched)[0]):
        ax.plot(time_unsketched[i,:], acc_unsketched[i,:], 
                label ="RKHS Decomposition", 
                linestyle='',
                color=color[0], 
                marker = "o",
                markersize=2, alpha=0.6)
        ax.plot(time_unsketched[i,:], acc_unsketched[i,:], 
                label ="RKHS Decomposition", 
                linestyle='dashed', color=color[0], alpha=0.4)
    line_list=[line1]    
    for j in range(len(s1_list)):
        time_mean = np.mean(time_list_list[j], axis=0)
        acc_mean = np.mean(acc_list_list[j], axis=0)
        line2, = ax.plot(time_mean, acc_mean, 
                         label ="S-RKHS, S1=" + str(s1_list[j]), 
                         linestyle='solid', color=color[j+1], linewidth=0)
        line_list += [line2]
        for i in range(np.shape(time_list_list[j])[0]):
            ax.plot(time_list_list[j][i], acc_list_list[j][i], 
                    label ="S-RKHS, S1=" + str(s1_list[j]), 
                    linestyle='',
                    color=color[j+1],
                    marker = "o",
                    markersize=2, alpha=0.6)
            ax.plot(time_list_list[j][i], acc_list_list[j][i], 
                    label ="S-RKHS, S1=" + str(s1_list[j]), 
                    linestyle='dashed', color=color[j+1], alpha=0.4)
        
    

    plt.xlabel('Time')
    plt.ylabel('Fit')
    # Put a legend below current axis
    #plt.ylim([0.5, 0.85])
    
    plt.legend(handles = line_list, loc='lower right')
    
    plt.savefig('comparison_r' + str(R) + '_s2_'+ str(s2) + '.pdf', format = 'pdf')
    
    
    
    
    
    
    
    
    
    
    
    
    
    fig, ax = plt.subplots()
    
    
    line1, = ax.plot(time_unsketched_mean, acc_unsketched_mean, 
                     label ="RKHS Decomposition", linestyle='solid', color=color[0], linewidth=4)
    for i in range(np.shape(time_unsketched)[0]):
        ax.plot(time_unsketched[i,:], acc_unsketched[i,:], 
                label ="RKHS Decomposition", 
                linestyle='',
                color=color[0], 
                marker = "o",
                markersize=2, alpha=0.6)
        ax.plot(time_unsketched[i,:], acc_unsketched[i,:], 
                label ="RKHS Decomposition", 
                linestyle='dashed', color=color[0], alpha=0.4)
    line_list=[line1]    
    for j in range(len(s1_list)):
        time_mean = np.mean(time_list_list[j], axis=0)
        acc_mean = np.mean(acc_list_list[j], axis=0)
        line2, = ax.plot(time_mean, acc_mean, 
                         label ="S-RKHS, S1=" + str(s1_list[j]), 
                         linestyle='solid', color=color[j+1], linewidth=4)
        line_list += [line2]
        for i in range(np.shape(time_list_list[j])[0]):
            ax.plot(time_list_list[j][i], acc_list_list[j][i], 
                    label ="S-RKHS, S1=" + str(s1_list[j]), 
                    linestyle='',
                    color=color[j+1],
                    marker = "o",
                    markersize=2, alpha=0.6)
            ax.plot(time_list_list[j][i], acc_list_list[j][i], 
                    label ="S-RKHS, S1=" + str(s1_list[j]), 
                    linestyle='dashed', color=color[j+1], alpha=0.4)
        
    

    plt.xlabel('Time')
    plt.ylabel('Fit')
    # Put a legend below current axis
    plt.ylim([0.5, 0.85])
    plt.xlim([0, 12])
    
    plt.legend(handles = line_list, loc='lower right')
    
    plt.savefig('comparison_partial_r' + str(R) + '_s2_'+ str(s2) + '.pdf', format = 'pdf')
    








