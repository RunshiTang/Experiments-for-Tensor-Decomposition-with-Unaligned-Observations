# cd "//sscwin/dfsroot/Users/rtang56/Desktop/RKHS/RKHS/simulation/setting1"

import os
os.chdir("simulation/rkhs/setting1_1/")

import pandas as pd
import numpy as np
import re

import glob

acc_lambda_csvfiles = []
for file in glob.glob("acc_lambda_*.csv"):
    acc_lambda_csvfiles.append(file)
    
res = [(sub.split('lambda_')[1]) for sub in acc_lambda_csvfiles] 

res = [float(filename.replace('.csv', '')) for filename in res]

lambda_list = np.array(res)
sort_index = np.argsort(lambda_list)
lambda_list = lambda_list[sort_index]
acc_lambda_csvfiles = [acc_lambda_csvfiles[i] for i in sort_index]

acc_list_list = np.zeros((5,2))

for i in range(len(acc_lambda_csvfiles)):
    tmp_data = pd.read_csv(acc_lambda_csvfiles[i], sep=",", header=0,index_col=0).to_numpy()
    last_column = tmp_data[:, -1]
    acc_list_list[i,0] = np.mean((1-last_column)**2)
    acc_list_list[i,1] = np.std((1-last_column)**2)

standard_form_strings = [f"{num:.10f}".rstrip('0').rstrip('.') for num in lambda_list]
df = pd.DataFrame(acc_list_list.T, columns=standard_form_strings)
df.to_csv('acc_lambda.csv')

