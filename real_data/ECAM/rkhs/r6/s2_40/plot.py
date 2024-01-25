


import os
os.chdir("//sscwin/dfsroot/Users/rtang56/Desktop/RKHS_new_new/RKHS/real_data/r6/s2_40")

import pandas as pd
import numpy as np

import glob

acc_s1_csvfiles = []
for file in glob.glob("acc_s1*.csv"):
    acc_s1_csvfiles.append(file)
    
res = [(sub.split('s1_')[1]) for sub in acc_s1_csvfiles] 
res = [int(sub.split('.')[0]) for sub in res]
s1_list = np.array(res)
sort_index = np.argsort(s1_list)
s1_list = s1_list[sort_index]
acc_s1_csvfiles = [acc_s1_csvfiles[i] for i in sort_index]

shape_tmp = (len(acc_s1_csvfiles),) + np.shape(pd.read_csv(acc_s1_csvfiles[0], sep=",", header=0,index_col=0).to_numpy())
acc_list_list = np.zeros(shape_tmp)

for i in range(len(acc_s1_csvfiles)):
    acc_list_list[i,:,:] = pd.read_csv(acc_s1_csvfiles[i], sep=",", header=0,index_col=0).to_numpy()

acc_list_list = acc_list_list

time_s1_csvfiles = []
for file in glob.glob("time_s1*.csv"):
    time_s1_csvfiles.append(file)    
    
res = [(sub.split('s1_')[1]) for sub in time_s1_csvfiles] 
res = [int(sub.split('.')[0]) for sub in res]
s1_list = np.array(res)
sort_index = np.argsort(s1_list)
s1_list = s1_list[sort_index]
time_s1_csvfiles = [time_s1_csvfiles[i] for i in sort_index]

shape_tmp = (len(time_s1_csvfiles),) + np.shape(pd.read_csv(time_s1_csvfiles[0], sep=",", header=0,index_col=0).to_numpy())
time_list_list = np.zeros(shape_tmp)

for i in range(len(time_s1_csvfiles)):
    time_list_list[i,:,:] = pd.read_csv(time_s1_csvfiles[i], sep=",", header=0,index_col=0).to_numpy()    


time_unsketched = pd.read_csv('time_unsketched.csv', sep=",", header=0,index_col=0).to_numpy()
acc_unsketched = pd.read_csv('acc_unsketched.csv', sep=",", header=0,index_col=0).to_numpy()

acc_unsketched = acc_unsketched

time_unsketched_mean = np.mean(time_unsketched, axis=0)
acc_unsketched_mean = np.mean(acc_unsketched, axis=0)

import matplotlib.pyplot as plt
color = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']


fig, ax = plt.subplots(figsize=(3,9/4))

ax.plot(time_unsketched_mean, acc_unsketched_mean, 
                 label ="RKHS", linestyle='solid', color=color[0], linewidth=0)
for i in range(np.shape(time_unsketched)[0]):
    ax.plot(time_unsketched[i,:], acc_unsketched[i,:], 
            label ="RKHS", 
            linestyle='',
            color=color[0], 
            marker = "o",
            markersize=2, alpha=0.6)
    line1, = ax.plot(time_unsketched[i,:], acc_unsketched[i,:], 
            label ="RKHS", 
            linestyle='dashed', color=color[0], alpha=0.4)
line_list=[line1]    

for j in range(len(s1_list)):
    time_mean = np.mean(time_list_list[j], axis=0)
    acc_mean = np.mean(acc_list_list[j], axis=0)
    ax.plot(time_mean, acc_mean, 
                     label ="S-RKHS, S1=" + str(s1_list[j]), 
                     linestyle='solid', color=color[j+1], linewidth=0)
    
    for i in range(np.shape(time_list_list[j])[0]):
        ax.plot(time_list_list[j][i], acc_list_list[j][i], 
                label ="S-RKHS, S1=" + str(s1_list[j]), 
                linestyle='',
                color=color[j+1],
                marker = "o",
                markersize=2, alpha=0.6)
        line2, = ax.plot(time_list_list[j][i], acc_list_list[j][i], 
                label ="S-RKHS, S1=" + str(s1_list[j]), 
                linestyle='dashed', color=color[j+1], alpha=0.4)
    line_list += [line2]


plt.xlabel('Time (sec)')
plt.ylabel('Fit')
# Put a legend below current axis

plt.legend(handles = line_list, loc='lower right')

plt.tight_layout()

#plt.yscale('log')
plt.xlim(-0.4, 10.1)
plt.ylim(0.35, 0.79)
plt.savefig('comparison_r6_s2_40_realdata.pdf', format = 'pdf')

plt.clf()









