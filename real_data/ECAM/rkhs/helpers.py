import numpy as np
import json
import time
import scipy
import random
from random import choices




def RKHS_microbiome(microbiome, R = 5, 
                    lambda_parameter = 0.001, n_max = 20, 
                    CG = False, sketch = False, 
                    s1 = 30, s2 = 40, s3 = 10, stop_criterion = False):
    
    def conjgrad(A, b, x):
        r = b - np.dot(A, x)
        p = r
        rsold = np.dot(np.transpose(r), r)
        
        for i in range(len(b)):
            Ap = np.dot(A, p)
            alpha = rsold / np.dot(np.transpose(p), Ap)
            x = x + np.dot(alpha, p)
            r = r - np.dot(alpha, Ap)
            rsnew = np.dot(np.transpose(r), r)
            if np.sqrt(rsnew) < 1e-8:
                break
            p = r + (rsnew/rsold)*p
            rsold = rsnew
        return x, i


    def initial_A_B(R, n, p):
        A = []
        for i in range(R):
            tmp_vec = [np.random.uniform(0, 1) for j in range(n)]
            tmp_vec_arr = np.array(tmp_vec)
            A.append(tmp_vec / np.linalg.norm(tmp_vec_arr))
        A = np.transpose(np.array(A))    

        B = []
        for i in range(R):
            tmp_vec = [np.random.uniform(0, 1) for j in range(p)]
            tmp_vec_arr = np.array(tmp_vec)
            B.append(tmp_vec / np.linalg.norm(tmp_vec_arr))
        B = np.transpose(np.array(B))
        return A, B


    def update_Xi(A, B, last_theta = 0, CG = False):
        AB = np.zeros(shape=(p*sum_Ti, R*sum_Ti))
        for i in range(R):
            tmp_v = np.zeros(sum_Ti)
            tmp_sum = 0
            for j in range(n):
                tmp_v[tmp_sum:(tmp_sum+dim_Ti_list[j])] = A[:,i][j]
                tmp_sum += dim_Ti_list[j]
            AB[:,((i)*sum_Ti):((i+1)*sum_Ti)] = np.kron(B[:,i][:, np.newaxis], np.diag(tmp_v))
            
        coef_A = M_bar.transpose() @ AB.transpose() @ AB @ MK_bar + lambda_parameter * np.identity(len(T)*R)
        coef_b = M_bar.transpose() @ AB.transpose() @ y
        
        if CG == False:
            theta1 = np.linalg.solve(coef_A, coef_b)
        else:
            initial_theta = last_theta
            theta1, iterations = conjgrad(coef_A, coef_b, initial_theta)
        
        return theta1

    def update_Xi_sketch(A, B, s1, s2, s3, last_theta = 0, CG = False):
        sample_index_n = choices(range(n), k = s1)
        sample_index_p = choices(range(p), k = s2)
        sample_index_Ti_list = [choices(range(tmp_index), k = s3) for tmp_index in dim_Ti_list]
        sample_T_list = []
        for i in sample_index_n:
            sample_T_list.append([T_list[i][j] for j in sample_index_Ti_list[i]])
            
        sample_T = []
        for i in range(len(sample_T_list)):
            sample_T = Union(sample_T, sample_T_list[i])
        sample_T.sort() 
        
        MK_hat = np.zeros(shape=(s1*s3, len(T)))

        tmp_sum = 0
        for i in range(s1):
            if i == 0: 
                MK_hat[0:s3, :] = Bernoulli_kernel(np.array(sample_T_list[i]), np.array(T))
            else:
                MK_hat[(tmp_sum):(s3 + tmp_sum), :] = Bernoulli_kernel(np.array(sample_T_list[i]), np.array(T))
            tmp_sum += s3

        M_hat = np.zeros(shape=(s1*s3, len(T)))

        tmp_sum = 0
        for i in range(s1):
            tmp_M = np.zeros((s3, len(T)))
            for j in range(s3):
                tmp_index = np.where(T == sample_T_list[i][j])[0]
                tmp_M[j,tmp_index] = 1
            M_hat[(tmp_sum):(s3 + tmp_sum), :] = tmp_M
            tmp_sum += s3

        M_bar_hat = np.kron(np.eye(R,dtype=int),M_hat)
        MK_bar_hat = np.kron(np.eye(R,dtype=int),MK_hat)
        
        y_hat = np.zeros(s2*s1*s3)

        for b_which in range(s2):
            b_index = sample_index_p[b_which]
            for a_which in range(s1):
                a_index = sample_index_n[a_which]
                for xi_which in range(s3):
                    xi_index = sample_index_Ti_list[a_index][xi_which]
                    tmp_index = int(b_which*s1*s3 + a_which*s3 + xi_which)
                    y_hat[tmp_index] = microbiome[a_index][b_index+1][xi_index]
                    
        AB_hat = np.zeros(shape=(s2*s3*s1, R*s1*s3))

        for i in range(R):
            tmp_v = np.zeros(s3*s1)
            tmp_sum = 0
            for j in range(s1):
                tmp_v[tmp_sum:(tmp_sum+s3)] = A[sample_index_n,i][j]
                tmp_sum += s3
            AB_hat[:,(i*s3*s1):((i+1)*s3*s1)] = np.kron(B[sample_index_p,i][:, np.newaxis], np.diag(tmp_v))
                
        coef_A_hat = M_bar_hat.transpose() @ AB_hat.transpose() @ AB_hat @ MK_bar_hat + lambda_parameter * np.identity(len(T)*R)
        coef_b_hat = M_bar_hat.transpose() @ AB_hat.transpose() @ y_hat
            
        if CG == False:
            theta1 = np.linalg.solve(coef_A_hat, coef_b_hat)
        else:
            initial_theta = last_theta
            theta1, iterations = conjgrad(coef_A_hat, coef_b_hat, initial_theta)
        
        return theta1

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

    def update_A(B, theta):
        new_A = np.zeros((n,R))
        theta_tmp = theta.reshape((len(T),R), order = 'F')
        for i in range(n):
            Xi_i = np.array(Bernoulli_kernel(T_list[i], np.array(T))) @ theta_tmp
            coe1_A = scipy.linalg.khatri_rao(B, Xi_i) 
            coe1_b = np.array(microbiome[i][1:]).reshape(p*dim_Ti_list[i])
            solution = np.linalg.lstsq(coe1_A, coe1_b)[0]
            new_A[i,:] = solution
        return new_A

    def update_B(A, theta):
        theta_tmp = theta.reshape((len(T),R), order = 'F')
        
        Xi = MK @ theta_tmp
        
        coe1_A = np.zeros((sum_Ti,R))
        for j in range(R):
            tmp_index_row = 0
            for i in range(n):
                coe1_A[tmp_index_row:(tmp_index_row + dim_Ti_list[i]), j] = A[i, j]*Xi[tmp_index_row:(tmp_index_row + dim_Ti_list[i]), j]
                tmp_index_row += dim_Ti_list[i]
        
        result = np.linalg.lstsq(coe1_A, coe1_b_update_B)
        new_B = np.transpose(result[0])
        
        return new_B
    
    def calculate_fit(A, B, theta):
        AB = np.zeros(shape=(p*sum_Ti, R*sum_Ti))
        for i in range(R):
            tmp_v = np.zeros(sum_Ti)
            tmp_sum = 0
            for j in range(n):
                tmp_v[tmp_sum:(tmp_sum+dim_Ti_list[j])] = A[:,i][j]
                tmp_sum += dim_Ti_list[j]
            AB[:,((i)*sum_Ti):((i+1)*sum_Ti)] = np.kron(B[:,i][:, np.newaxis], np.diag(tmp_v))
            
        fit = 1 - np.linalg.norm(AB @ MK_bar @ theta - y) / np.linalg.norm(y)
        loss = (np.linalg.norm(AB @ MK_bar @ theta - y) )**2 + lambda_parameter * (np.transpose(theta) @ K_bar @ theta)
        return fit, loss
    
    
    T_list = []
    
    for i in range(len(microbiome)):
        T_list.append(np.array(microbiome[i][0])/739)    

    T = []

    for i in range(len(microbiome)):
        T = Union(T, np.array(microbiome[i][0])/739)
        
    T.sort() 

    sum_Ti = 0

    dim_Ti_list = []

    for i in range(len(microbiome)):
        dim_Ti_list.append(len(microbiome[i][0]))
        sum_Ti += len(microbiome[i][0])

    n = len(microbiome)

    p = len(microbiome[1]) - 1
    
    all_T =  [j for i in T_list for j in i]
    
    K = Bernoulli_kernel(np.array(T), np.array(T))
    
    MK = np.zeros(shape=(sum_Ti, len(T)))

    tmp_sum = 0
    for i in range(n):
        if i == 0: 
            MK[0:dim_Ti_list[i], :] = Bernoulli_kernel(np.array(T_list[i]), np.array(T))
        else:
            MK[(tmp_sum):(dim_Ti_list[i] + tmp_sum), :] = Bernoulli_kernel(np.array(T_list[i]), np.array(T))
        tmp_sum += dim_Ti_list[i]

    M = np.zeros(shape=(sum_Ti, len(T)))

    tmp_sum = 0
    for i in range(n):
        tmp_M = np.zeros((dim_Ti_list[i], len(T)))
        for j in range(dim_Ti_list[i]):
            tmp_index = np.where(T == T_list[i][j])[0]
            tmp_M[j,tmp_index] = 1
        M[(tmp_sum):(dim_Ti_list[i] + tmp_sum), :] = tmp_M
        tmp_sum += dim_Ti_list[i]

    M_bar = np.kron(np.eye(R,dtype=int),M)
    MK_bar = np.kron(np.eye(R,dtype=int),MK)
    K_bar = np.kron(np.eye(R,dtype=int),K)
    
    y = np.zeros(p*sum_Ti)

    for b_index in range(p):
        for a_index in range(n):
            for xi_index in range(dim_Ti_list[a_index]):
                tmp_index = int(b_index*sum_Ti + np.array(dim_Ti_list[0:(a_index)]).sum() + xi_index)
                y[tmp_index] = microbiome[a_index][b_index+1][xi_index]
    
    coe1_b_update_B = np.zeros((p,sum_Ti))
    tmp_index_col = 0
    for i in range(n):
        for j in range(dim_Ti_list[i]):
            for h in range(p):
                coe1_b_update_B[h, tmp_index_col + j] = microbiome[i][1+h][j]    
        tmp_index_col += dim_Ti_list[i]
    coe1_b_update_B = np.transpose(coe1_b_update_B)
    
    A, B = initial_A_B(R = R, n = n, p = p)
    
    for r in range(R):
        A[:,r] = A[:,r] / np.linalg.norm( A[:,r])
        B[:,r] = B[:,r] / np.linalg.norm( B[:,r])
    
    theta = 0
    
    last_theta = theta
    if sketch:
        theta = update_Xi_sketch(A, B, s1 = s1, s2 = s2, s3 = s3, 
                                 last_theta = last_theta, CG = CG)
    else:
        theta = update_Xi(A, B, last_theta = last_theta, CG = CG)
    
    fit = []
    loss = []
    
    fit_tmp = calculate_fit(A, B, theta)[0]
    
    fit += [fit_tmp]
    
    loss_tmp = calculate_fit(A, B, theta)[1]
    
    loss += [loss_tmp]
    
    A_list = []
    B_list = []
    theta_list = []
    stopped = False
    time_list = [0]
    for iterator in range(n_max): 
        time1 = time.time()
        A = update_A(B, theta)
        for r in range(R):
            A[:,r] = A[:,r] / np.linalg.norm( A[:,r])

        B = update_B(A, theta)
        for r in range(R):
            B[:,r] = B[:,r] / np.linalg.norm( B[:,r])
            
        last_theta = theta
        if sketch:
            theta = update_Xi_sketch(A, B, s1 = s1, s2 = s2, s3 = s3, 
                                     last_theta = last_theta, CG = CG)
        else:
            theta = update_Xi(A, B, last_theta = last_theta, CG = CG)
        
        A_list += [A]
        B_list += [B]
        theta_list += [theta]
        
        time2 = time.time()
        
        time_list += [time_list[-1] + time2-time1]
        
        fit_tmp = calculate_fit(A, B, theta)[0]
        
        loss_tmp = calculate_fit(A, B, theta)[1]
        
        print("fit: " + str(fit_tmp))
        
        print("loss: " + str(loss_tmp))
        
        fit += [fit_tmp]
        
        loss += [loss_tmp]
        
        final_fit = fit[-1]
        
        
    return [[A, B, theta], final_fit, fit, stopped, time_list, loss]

        
    
    














