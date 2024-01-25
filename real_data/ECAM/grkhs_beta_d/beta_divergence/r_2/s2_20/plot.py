

import os
os.chdir("//sscwin/dfsroot/Users/rtang56/Desktop/RKHS_new_new/Generalized_RKHS/real_data/beta_divergence/r_2/setting1")

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
                 label ="PGD", linestyle='solid', color=color[0], linewidth=0)
for i in range(np.shape(time_unsketched)[0]):
    ax.plot(time_unsketched[i,:], acc_unsketched[i,:], 
            label ="PGD", 
            linestyle='',
            color=color[0], 
            marker = "o",
            markersize=2, alpha=0.6)
    line1, = ax.plot(time_unsketched[i,:], acc_unsketched[i,:], 
            label ="PGD", 
            linestyle='dashed', color=color[0], alpha=0.4)
line_list=[line1]    

for j in range(len(s1_list)):
    time_mean = np.mean(time_list_list[j], axis=0)
    acc_mean = np.mean(acc_list_list[j], axis=0)
    ax.plot(time_mean, acc_mean, 
                     label ="PSGD, S1=" + str(s1_list[j]), 
                     linestyle='solid', color=color[j+1], linewidth=0)
    
    for i in range(np.shape(time_list_list[j])[0]):
        ax.plot(time_list_list[j][i], acc_list_list[j][i], 
                label ="SPGS, S1=" + str(s1_list[j]), 
                linestyle='',
                color=color[j+1],
                marker = "o",
                markersize=2, alpha=0.6)
        line2, = ax.plot(time_list_list[j][i], acc_list_list[j][i], 
                label ="PSGD, S1=" + str(s1_list[j]), 
                linestyle='dashed', color=color[j+1], alpha=0.4)
    line_list += [line2]


plt.xlabel('Time (sec)')
plt.ylabel('Loss')
# Put a legend below current axis

plt.legend(handles = line_list, loc='upper right')

plt.tight_layout()

plt.yscale('log')
plt.ylim(0.1, 10)

plt.savefig('beta_d_r_2_s2_20.pdf', format = 'pdf',bbox_inches='tight')

plt.clf()











ticks = np.arange(len(acc_list_list[0][0]))
 
data = np.concatenate((acc_unsketched[np.newaxis,:], acc_list_list))

lable_list = ["PSGD, S1="+str(i) for i in s1_list]
lable_list = ["PGD"] + lable_list

# create a boxplot for two arrays separately,
# the position specifies the location of the
# particular box in the graph,
# this can be changed as per your wish. Use width
# to specify the width of the plot

plot_list = []

plt.figure(figsize=(8,4))

for i in range(len(lable_list)):
    tmp_plot = plt.boxplot((data[i]),
                           positions=np.array(np.arange(len(data[0][0])))+0.17*(i-len(lable_list)/2),
                           widths=0.15, 
                           flierprops={'marker': 'o', 'markersize': 4, 'markerfacecolor': color[i], "markeredgecolor":color[i]})
    plot_list.append(tmp_plot)
 
# each plot returns a dictionary, use plt.setp()
# function to assign the color code
# for all properties of the box plot of particular group
# use the below function to set color for particular group,
# by iterating over all properties of the box plot
def define_box_properties(plot_name, color_code, label):
    for k, v in plot_name.items():
        plt.setp(plot_name.get(k), color=color_code)
         
    # use plot function to draw a small line to name the legend.
    plt.plot([], c=color_code, label=label)
    plt.legend(loc='upper right')
 
 
# setting colors for each groups

for i in range(len(lable_list)):
    define_box_properties(plot_list[i], color[i],lable_list[i])


# set the x label values
plt.xticks(np.arange(0, len(ticks), 1), ticks+1)

plt.axhline(y = -5.315734, linestyle = ':', color = "#17becf", linewidth = 0.5)
# set the limit for x axis
#plt.xlim(-2, len(ticks)*2)
 
# set the limit for y axis
#plt.ylim(-5.35, -4.6)
#plt.xlim(0.5, 14.5)
plt.xlabel('Epoch')
plt.ylabel('Loss to the Nominal')
plt.yscale('log')
plt.tight_layout()

plt.savefig('loss_iteration.pdf', format = 'pdf')

plt.clf()










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


time_unsketched_mean = np.mean(time_unsketched, axis=0)
acc_unsketched_mean = np.mean(acc_unsketched, axis=0)




plt.figure(figsize=(3,3))

data_concatenated = []
data = np.concatenate((time_unsketched[np.newaxis,:], time_list_list))
for i in range(len(data)):
    for j in range(len(data[i])):
        tmp_v = np.ediff1d(data[i][j], to_begin=data[i][j][0])
        data[i][j] = tmp_v
    data_concatenated.append(np.concatenate(data[i]))
        




lable_list = ["S1="+str(i) for i in s1_list]
lable_list = ["PGD"] + lable_list 

plt.boxplot(data_concatenated, sym = ".", widths = 0.8)

locs, labels = plt.xticks()

plt.xticks(locs, lable_list, rotation=45)

plt.ylabel('Time (sec)')
plt.tight_layout()
plt.savefig('time_box.pdf', format = 'pdf')
plt.clf()






