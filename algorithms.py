def RKHS(microbiome, R = 5, 
                    lambda_parameter = 0.001, n_max = 20, 
                    CG = False, sketch = False, 
                    stop_steps = 3,
                    stop = False, 
                    s1 = 30, s2 = 40, s3 = 10, stop_criterion = False,
                    test = True):
    
    def conjgrad(A, b, x):
        """
        A function to solve [A]{x} = {b} linear equation system with the 
        conjugate gradient method.
        More at: http://en.wikipedia.org/wiki/Conjugate_gradient_method
        ========== Parameters ==========
        A : matrix 
            A real symmetric positive definite matrix.
        b : vector
            The right hand side (RHS) vector of the system.
        x : vector
            The starting guess for the solution.
        """  
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
        theta = theta.reshape((len(T),R), order = 'F')
        for i in range(n):
            Xi_i = np.array(Bernoulli_kernel(T_list[i], np.array(T))) @ theta
            coe1_A = scipy.linalg.khatri_rao(B, Xi_i) 
            coe1_b = np.array(microbiome[i][1:]).reshape(p*dim_Ti_list[i])
            solution = np.linalg.lstsq(coe1_A, coe1_b)[0]
            new_A[i,:] = solution
        return new_A

    def update_B(A, theta):
        theta = theta.reshape((len(T),R), order = 'F')
        
        Xi = MK @ theta
        
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
        return fit
    
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
    
    fit_tmp = calculate_fit(A, B, theta)
    
    fit += [fit_tmp]
    
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
        
        fit_tmp = calculate_fit(A, B, theta)
        
        fit += [fit_tmp]
        
        final_fit = fit[-1]
        fit_array = np.array(fit)
        
        if stop & (len(fit) > (stop_steps+1)):
            fit_last0 = np.flip(fit_array)[0:(stop_steps)]
            fit_last1 = np.flip(fit_array)[1:(stop_steps+1)]
            if np.all(fit_last0 < fit_last1 + 10**(-4)):
                stopped = True
                return [[A_list[-stop_steps-1], B_list[-stop_steps-1], theta_list[-stop_steps-1]], final_fit, fit, stopped, time_list]
            
        if test:
            print("epoch "+ str(iterator) + " loss:" +str(fit_tmp))
        
        
        
    return [[A, B, theta], final_fit, fit, stopped, time_list]

        
    
    
def RKHS_generalized_sketch_counts(microbiome, R = 5, 
                    lambda_parameter = 1000,
                    leanring_rate = 0.0001, 
                    n_max = 10, 
                    bad_epochs_max = 3,
                    s1 = 30, s2 = 40, s3 = 10, 
                    stop_steps = 3, 
                    iterations_per_epoch = 10,
                    bias_parameter = 0.5,
                    test = False,
                    decrease_lr = False,
                    gradient_clipping = False,
                    clipping_threshold = 0.1):

    def Poisson_loss_f(x,y):
        return_v = x+bias_parameter-y*np.log(x+bias_parameter)
        #return_v[np.isnan(return_v)] = 0
        return return_v

    def Poisson_loss_gradient(x,y):
        return_v = 1-y/(x+bias_parameter)
        return_v[np.isnan(return_v)] = 0
        return return_v
        
    def Bernoulli_kernel(x, y):
        k1_x = x-0.5
        k1_y = y-0.5
        k2_x = 0.5*(k1_x**2-1/12)
        k2_y = 0.5*(k1_y**2-1/12)
        xy = abs(np.outer(x,np.ones(len(y))) - np.outer(np.ones(len(x)), y))
        k4_xy = ((xy-0.5)**4-0.5 * (xy-0.5)**2 + 7/240)/24
        kern_xy = np.outer(k1_x, k1_y) + np.outer(k2_x, k2_y) - k4_xy + 1
        return kern_xy
    
    def radial_kernel(x,y):
        if np.isscalar(x):
            x = np.array([x])
        if np.isscalar(y):
            y = np.array([y])
        nx = np.shape(x)[0]
        ny = np.shape(y)[0]
        result_matrix = np.zeros((nx,ny))
        for i in range(nx):
          for j in range(ny):
            result_matrix[i,j] = np.exp(-np.abs(x[i]-y[j]))
        return result_matrix


    def sampling():
        sample_index_n = choices(range(n), k = s1)
        sample_index_p = choices(range(p), k = s2)
        sample_index_Ti_list = [choices(range(tmp_index), k = s3) for tmp_index in dim_Ti_list]

        sample_T_list = []
        for i in sample_index_n:
            sample_T_list.append([T_list[i][j] for j in sample_index_Ti_list[i]])
         
        y_hat_sketched = np.zeros(p*sum_Ti) + float("nan")
        
        for i in sample_index_n:
            for j in sample_index_p:
                for t in sample_index_Ti_list[i]:
                    tmp_sum = 0
                    for r in range(R):
                        tmp_sum += A[i,r] * B[j,r] * MK[i,:] @ theta[(len(T)*r):(len(T)*(r+1))]
                    y_hat_sketched[j*sum_Ti + int(np.sum(dim_Ti_list[0:i])) + t] = tmp_sum[0]
        
        y_sketched = np.zeros(p*sum_Ti) + float("nan")
        
        for i in sample_index_n:
            for j in sample_index_p:
                for t in sample_index_Ti_list[i]:
                    y_sketched[j*sum_Ti + int(np.sum(dim_Ti_list[0:i])) + t] = y[j*sum_Ti + int(np.sum(dim_Ti_list[0:i])) + t]
        
        MK_hat = np.zeros(shape=(s1*s3, len(T)))

        tmp_sum = 0
        for i in range(s1):
            if i == 0: 
                MK_hat[0:s3, :] = radial_kernel(np.array(sample_T_list[i]), np.array(T))
            else:
                MK_hat[(tmp_sum):(s3 + tmp_sum), :] = radial_kernel(np.array(sample_T_list[i]), np.array(T))
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

        MK_bar_hat = np.kron(np.eye(R,dtype=int),MK_hat)

        y_sketched_for_theta = np.zeros(s2*s1*s3)

        for b_which in range(s2):
            b_index = sample_index_p[b_which]
            for a_which in range(s1):
                a_index = sample_index_n[a_which]
                for xi_which in range(s3):
                    xi_index = sample_index_Ti_list[a_index][xi_which]
                    tmp_index = int(b_which*s1*s3 + a_which*s3 + xi_which)
                    y_sketched_for_theta[tmp_index] = microbiome[a_index][b_index+1][xi_index]

        AB_hat = np.zeros(shape=(s2*s3*s1, R*s1*s3))

        for i in range(R):
            tmp_v = np.zeros(s3*s1)
            tmp_sum = 0
            for j in range(s1):
                tmp_v[tmp_sum:(tmp_sum+s3)] = A[sample_index_n,i][j]
                tmp_sum += s3
            AB_hat[:,(i*s3*s1):((i+1)*s3*s1)] = np.kron(B[sample_index_p,i][:, np.newaxis], np.diag(tmp_v))

        ABMK_bar_hat = AB_hat @ MK_bar_hat
        y_hat_sketched_for_theta = ABMK_bar_hat @ theta            
        
        return y_hat_sketched, y_sketched, y_sketched_for_theta, y_hat_sketched_for_theta, ABMK_bar_hat
                    
    def update_A_sketch():
        grad_A = np.zeros((n, R))
        for i in range(n):
            theta_matrix = theta.reshape((len(T),R))
            Xi_i = np.array(radial_kernel(T_list[i], np.array(T))) @ theta_matrix
            KR_matrix = scipy.linalg.khatri_rao(B, Xi_i)
            tmp_index1 = np.rint(np.sum(dim_Ti_list[0:i])).astype(int)
            tmp_index2 = np.rint(np.sum(dim_Ti_list[0:(i+1)])).astype(int)
            
            y_subset = np.zeros(dim_Ti_list[i]*p)
            y_hat_subset = np.zeros(dim_Ti_list[i]*p)
            
            for j in range(p):
                y_subset[(dim_Ti_list[i]*j) : (dim_Ti_list[i]*(j+1))] = y_sketched[(j*sum_Ti + tmp_index1):(j*sum_Ti + tmp_index2)]
                y_hat_subset[(dim_Ti_list[i]*j) : (dim_Ti_list[i]*(j+1))] = y_hat_sketched[(j*sum_Ti + tmp_index1):(j*sum_Ti + tmp_index2)]
            
            f_grad_A_rowvector = Poisson_loss_gradient(y_hat_subset, y_subset)
            grad_A[i,:] = f_grad_A_rowvector @ KR_matrix /omega_size
        
        #print(np.linalg.norm(grad_A))      
        
        if gradient_clipping: 
            norm_A = np.linalg.norm(grad_A)
            if norm_A>clipping_threshold:
                grad_A = clipping_threshold * grad_A / norm_A
        
        A_new = A - leanring_rate * grad_A * (omega_size / (s1*s2*s3))
        
        return A_new

    def update_B_sketch():
        f_grad_B_matrix = np.zeros((p, sum_Ti))

        for j in range(p):
            f_grad_B_matrix[j,:] = Poisson_loss_gradient(y_hat_sketched[(sum_Ti*j):(sum_Ti*(j+1))], 
                                                    y_sketched[(sum_Ti*j):(sum_Ti*(j+1))])

        theta_matrix = theta.reshape((len(T),R))
        Xi = np.array(radial_kernel(np.array(all_T), np.array(T))) @ theta_matrix
        KR_matrix = np.zeros((sum_Ti, R))

        for r in range(R):
            tmp_index_row = 0
            for i in range(n):
                KR_matrix[tmp_index_row:(tmp_index_row + dim_Ti_list[i]),r] = A[i, r] * Xi[tmp_index_row:(tmp_index_row + dim_Ti_list[i]), r]
                tmp_index_row += dim_Ti_list[i]
        
        grad_B = f_grad_B_matrix @ KR_matrix /omega_size
        
        if gradient_clipping: 
            norm_B = np.linalg.norm(grad_B)
            if norm_B>clipping_threshold:
                grad_B = clipping_threshold * grad_B / norm_B
        
        #print(np.linalg.norm(grad_B))   
        
        B_new = B - leanring_rate * grad_B * (omega_size / (s1*s2*s3))
        
        return B_new
                        
    def update_Xi_sketch():
        grad_theta_sketched = np.transpose(ABMK_bar_hat)@ Poisson_loss_gradient(y_hat_sketched_for_theta, y_sketched_for_theta[:, np.newaxis])/omega_size
        #print(np.linalg.norm(grad_theta_sketched))
        
        if gradient_clipping: 
            norm_Xi = np.linalg.norm(grad_theta_sketched)
            if norm_Xi>clipping_threshold:
                grad_theta_sketched = clipping_threshold * grad_theta_sketched / norm_Xi
        
        theta_new = theta - leanring_rate * grad_theta_sketched * (omega_size / (s1*s2*s3))
        return theta_new
                    
    def get_ABMK_bar():
        AB_new = np.zeros(shape=(p*sum_Ti, R*sum_Ti))
        for i in range(R):
            tmp_v = np.zeros(sum_Ti)
            tmp_sum = 0
            for j in range(n):
                tmp_v[tmp_sum:(tmp_sum+dim_Ti_list[j])] = A[:,i][j]
                tmp_sum += dim_Ti_list[j]
            AB_new[:,((i)*sum_Ti):((i+1)*sum_Ti)] = np.kron(B[:,i][:, np.newaxis], np.diag(tmp_v))

        ABMK_bar = AB_new @ MK_bar
        return ABMK_bar

    def get_y_hat():
        y_hat  = ABMK_bar @ theta
        return y_hat
                    
    def evaluate():
        loss = np.sum(Poisson_loss_f(y_hat.reshape(-1), y))/omega_size
        return loss
    
    def RKHS_norm_square(theta_tmp):
        return (np.transpose(theta_tmp)@K@theta_tmp)[0][0]


    def initial_A_B_theta():
        A = []
        for i in range(R):
            tmp_vec = [np.random.uniform(0, 1) for j in range(n)]
            tmp_vec_arr = np.array(tmp_vec)
            A.append(tmp_vec)
        A = np.transpose(np.array(A))

        B = []
        for i in range(R):
            tmp_vec = [np.random.uniform(0, 1) for j in range(p)]
            tmp_vec_arr = np.array(tmp_vec)
            B.append(tmp_vec)
        B = np.transpose(np.array(B))
        
        theta_last = np.random.uniform(0, 1, len(T)*R)[:, np.newaxis]
        
        return A, B, theta_last

    #prepare
    T_list = []
    for i in range(len(microbiome)):
        T_list.append(np.array(microbiome[i][0])/739)

    def Union(lst1, lst2):
        final_list = list(set(lst1) | set(lst2))
        return final_list

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

    omega_size = sum_Ti * p

    K = radial_kernel(np.array(T), np.array(T))

    MK = np.zeros(shape=(sum_Ti, len(T)))

    tmp_sum = 0
    for i in range(n):
        if i == 0: 
            MK[0:dim_Ti_list[i], :] = radial_kernel(np.array(T_list[i]), np.array(T))
        else:
            MK[(tmp_sum):(dim_Ti_list[i] + tmp_sum), :] = radial_kernel(np.array(T_list[i]), np.array(T))
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

    MK_bar = np.kron(np.eye(R,dtype=int),MK)


    y = np.zeros(p*sum_Ti)

    for b_index in range(p):
        for a_index in range(n):
            for xi_index in range(dim_Ti_list[a_index]):
                tmp_index = int(b_index*sum_Ti + np.array(dim_Ti_list[0:(a_index)]).sum() + xi_index)
                y[tmp_index] = microbiome[a_index][b_index+1][xi_index]
                
     
                
    learning_rate_initial = leanring_rate
    #run 
    
    A, B, theta = initial_A_B_theta()
    
    ABMK_bar = get_ABMK_bar()
    y_hat = get_y_hat()
    
    loss = [evaluate()]

    stopped = False

    time_list = [0]

    bad_epochs = 0
    
    evaluation = np.Inf
    
    for epoch in range(n_max):
        evaluation_old = evaluation
        time1 = time.time()
        for iterator in range(iterations_per_epoch):
            
            if decrease_lr:
                leanring_rate = learning_rate_initial/((iterator+iterations_per_epoch*epoch+1)**2)
            
            
            y_hat_sketched, y_sketched, y_sketched_for_theta, y_hat_sketched_for_theta, ABMK_bar_hat = sampling()
            A_new = update_A_sketch()
            B_new = update_B_sketch()
            theta_new = update_Xi_sketch()
            
            for r in range(R):
                a_norm_tmp = np.linalg.norm(A_new[:,r])
                b_norm_tmp = np.linalg.norm(B_new[:,r])
                theta_tmp = theta_new[r*len(T):(r+1)*len(T),:]
                theta_norm_tmp = RKHS_norm_square(theta_tmp) ** 0.5
                
                norm_prod = a_norm_tmp*b_norm_tmp*theta_norm_tmp
                
                if norm_prod > lambda_parameter:
                    theta_tmp = theta_tmp/theta_norm_tmp * (norm_prod**(1/3))
                    theta_new[r*len(T):(r+1)*len(T),:] = theta_tmp
                    A_new[:,r] = A_new[:,r] / a_norm_tmp * (norm_prod**(1/3))
                    B_new[:,r] = B_new[:,r] / b_norm_tmp * (norm_prod**(1/3))
            
            A_new = np.maximum(A_new, np.zeros(1))
            B_new =  np.maximum(B_new, np.zeros(1))
            theta_new =  np.maximum(theta_new, np.zeros(1))
            
            A = A_new
            B =  B_new
            theta = theta_new
        
            
        time2 = time.time()
        
        time_list += [time2-time1 + time_list[-1]]
        
        ABMK_bar = get_ABMK_bar()
        y_hat = get_y_hat()
        evaluation = evaluate()
        loss += [evaluation]
        
        if evaluation >= evaluation_old + 0.01: 
            bad_epochs += 1
        
        if bad_epochs >= bad_epochs_max:
            stopped = True
            break
        
        if test:
            print("epoch "+ str(epoch) + " loss:" +str(evaluation))
        
        time1 = time.time()
    final_loss = evaluation
    
        
    return [[A, B, theta], final_loss, loss, stopped, time_list]







def RKHS_generalized_counts(microbiome, R = 5, 
                    lambda_parameter = 1000,
                    leanring_rate = 0.0001, 
                    n_max = 10, 
                    bad_epochs_max = 3,
                    iterations_per_epoch = 10,
                    bias_parameter = 0.5,
                    test = False,
                    gradient_clipping = False,
                    clipping_threshold = 0.1):

    
    def Poisson_loss_f(x,y):
        return x+bias_parameter-y*np.log(x+bias_parameter)

    def Poisson_loss_gradient(x,y):
        return_v = 1-y/(x+bias_parameter)
        #return_v[np.isnan(return_v)] = 0
        return return_v
        
    def Bernoulli_kernel(x, y):
        k1_x = x-0.5
        k1_y = y-0.5
        k2_x = 0.5*(k1_x**2-1/12)
        k2_y = 0.5*(k1_y**2-1/12)
        xy = abs(np.outer(x,np.ones(len(y))) - np.outer(np.ones(len(x)), y))
        k4_xy = ((xy-0.5)**4-0.5 * (xy-0.5)**2 + 7/240)/24
        kern_xy = np.outer(k1_x, k1_y) + np.outer(k2_x, k2_y) - k4_xy + 1
        return kern_xy

    def radial_kernel(x,y):
        if np.isscalar(x):
            x = np.array([x])
        if np.isscalar(y):
            y = np.array([y])
        nx = np.shape(x)[0]
        ny = np.shape(y)[0]
        result_matrix = np.zeros((nx,ny))
        for i in range(nx):
          for j in range(ny):
            result_matrix[i,j] = np.exp(-(x[i]-y[j])**2)
        return result_matrix
                    
    def update_A():
        grad_A = np.zeros((n, R))
        for i in range(n):
            theta_matrix = theta.reshape((len(T),R))
            Xi_i = np.array(radial_kernel(T_list[i], np.array(T))) @ theta_matrix
            KR_matrix = scipy.linalg.khatri_rao(B, Xi_i)
            tmp_index1 = np.rint(np.sum(dim_Ti_list[0:i])).astype(int)
            tmp_index2 = np.rint(np.sum(dim_Ti_list[0:(i+1)])).astype(int)
            
            y_subset = np.zeros(dim_Ti_list[i]*p)
            y_hat_subset = np.zeros(dim_Ti_list[i]*p)
            
            for j in range(p):
                y_subset[(dim_Ti_list[i]*j) : (dim_Ti_list[i]*(j+1))] = y[(j*sum_Ti + tmp_index1):(j*sum_Ti + tmp_index2)]
                y_hat_subset[(dim_Ti_list[i]*j) : (dim_Ti_list[i]*(j+1))] = y_hat[(j*sum_Ti + tmp_index1):(j*sum_Ti + tmp_index2)].reshape(-1)
            
            f_grad_A_rowvector = Poisson_loss_gradient(y_hat_subset, y_subset)
            grad_A[i,:] = f_grad_A_rowvector @ KR_matrix /omega_size
        
        if gradient_clipping: 
            norm_A = np.linalg.norm(grad_A)
            if norm_A>clipping_threshold:
                grad_A = clipping_threshold * grad_A / norm_A    
        
        A_new = A - leanring_rate * grad_A
        
        #print("   grad_A:  " + str(np.linalg.norm(grad_A)))     
        
        return A_new

    def update_B():
        f_grad_B_matrix = np.zeros((p, sum_Ti))

        for j in range(p):
            f_grad_B_matrix[j,:] = Poisson_loss_gradient(y_hat[(sum_Ti*j):(sum_Ti*(j+1))].reshape(-1), 
                                                    y[(sum_Ti*j):(sum_Ti*(j+1))])

        theta_matrix = theta.reshape((len(T),R))
        Xi = np.array(radial_kernel(np.array(all_T), np.array(T))) @ theta_matrix
        KR_matrix = np.zeros((sum_Ti, R))

        for r in range(R):
            tmp_index_row = 0
            for i in range(n):
                KR_matrix[tmp_index_row:(tmp_index_row + dim_Ti_list[i]),r] = A[i, r] * Xi[tmp_index_row:(tmp_index_row + dim_Ti_list[i]), r]
                tmp_index_row += dim_Ti_list[i]
                
        grad_B = f_grad_B_matrix @ KR_matrix /omega_size
        
          
        
        if gradient_clipping: 
            norm_B = np.linalg.norm(grad_B)
            if norm_B>clipping_threshold:
                grad_B = clipping_threshold * grad_B / norm_B
        
        #print("   grad_B:  " + str(np.linalg.norm(grad_B)))   
        
        B_new = B - leanring_rate * grad_B
        
        return B_new
                        
    def update_Xi():
        grad_theta = np.transpose(ABMK_bar)@ Poisson_loss_gradient(y_hat, y[:, np.newaxis])/omega_size
        
        if gradient_clipping: 
            norm_Xi = np.linalg.norm(grad_theta)
            if norm_Xi>clipping_threshold:
                grad_theta = clipping_threshold * grad_theta / norm_Xi
        
        #print("   grad_theta:  " + str(np.linalg.norm(grad_theta)))   
        
        theta_new = theta - leanring_rate * grad_theta
        return theta_new               
                    
    def get_ABMK_bar():
        AB_new = np.zeros(shape=(p*sum_Ti, R*sum_Ti))
        for i in range(R):
            tmp_v = np.zeros(sum_Ti)
            tmp_sum = 0
            for j in range(n):
                tmp_v[tmp_sum:(tmp_sum+dim_Ti_list[j])] = A[:,i][j]
                tmp_sum += dim_Ti_list[j]
            AB_new[:,((i)*sum_Ti):((i+1)*sum_Ti)] = np.kron(B[:,i][:, np.newaxis], np.diag(tmp_v))

        ABMK_bar = AB_new @ MK_bar
        return ABMK_bar

    def get_y_hat():
        y_hat  = ABMK_bar @ theta
        return y_hat
                    
    def evaluate():
        loss = np.sum(Poisson_loss_f(y_hat.reshape(-1), y))/omega_size
        return loss
    
    def RKHS_norm_square(theta_tmp):
        return (np.transpose(theta_tmp)@K@theta_tmp)[0][0]


    def initial_A_B_theta():
        A = []
        for i in range(R):
            tmp_vec = [np.random.uniform(0, 1) for j in range(n)]
            tmp_vec_arr = np.array(tmp_vec)
            A.append(tmp_vec)
        A = np.transpose(np.array(A))    

        B = []
        for i in range(R):
            tmp_vec = [np.random.uniform(0, 1) for j in range(p)]
            tmp_vec_arr = np.array(tmp_vec)
            B.append(tmp_vec)
        B = np.transpose(np.array(B))
        
        theta_last = np.random.uniform(0, 1, len(T)*R)[:, np.newaxis]
        
        return A, B, theta_last

    #prepare
    T_list = []
    for i in range(len(microbiome)):
        T_list.append(np.array(microbiome[i][0])/739)

    def Union(lst1, lst2):
        final_list = list(set(lst1) | set(lst2))
        return final_list

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

    omega_size = sum_Ti * p

    K = radial_kernel(np.array(T), np.array(T))

    MK = np.zeros(shape=(sum_Ti, len(T)))

    tmp_sum = 0
    for i in range(n):
        if i == 0: 
            MK[0:dim_Ti_list[i], :] = radial_kernel(np.array(T_list[i]), np.array(T))
        else:
            MK[(tmp_sum):(dim_Ti_list[i] + tmp_sum), :] = radial_kernel(np.array(T_list[i]), np.array(T))
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

    MK_bar = np.kron(np.eye(R,dtype=int),MK)

    y = np.zeros(p*sum_Ti)

    for b_index in range(p):
        for a_index in range(n):
            for xi_index in range(dim_Ti_list[a_index]):
                tmp_index = int(b_index*sum_Ti + np.array(dim_Ti_list[0:(a_index)]).sum() + xi_index)
                y[tmp_index] = microbiome[a_index][b_index+1][xi_index]
                
                
    #run 
    A, B, theta = initial_A_B_theta()
    
    ABMK_bar = get_ABMK_bar()
    y_hat = get_y_hat()
    
    loss = [evaluate()]

    stopped = False

    time_list = [0]

    bad_epochs = 0
    
    evaluation = np.Inf
    
    for epoch in range(n_max):
        evaluation_old = evaluation
        time1 = time.time()
        for iterator in range(iterations_per_epoch):
            ABMK_bar = get_ABMK_bar()
            y_hat = get_y_hat()
            A_new = update_A()
            B_new = update_B()
            theta_new = update_Xi()
            
            for r in range(R):
                a_norm_tmp = np.linalg.norm(A_new[:,r])
                b_norm_tmp = np.linalg.norm(B_new[:,r])
                theta_tmp = theta_new[r*len(T):(r+1)*len(T),:]
                theta_norm_tmp = RKHS_norm_square(theta_tmp) ** 0.5
                
                norm_prod = a_norm_tmp*b_norm_tmp*theta_norm_tmp
                
                if norm_prod > lambda_parameter:
                    theta_tmp = theta_tmp/theta_norm_tmp * (norm_prod**(1/3))
                    theta_new[r*len(T):(r+1)*len(T),:] = theta_tmp
                    A_new[:,r] = A_new[:,r] / a_norm_tmp * (norm_prod**(1/3))
                    B_new[:,r] = B_new[:,r] / b_norm_tmp * (norm_prod**(1/3))
            
            A_new = np.maximum(A_new, np.zeros(1))
            B_new =  np.maximum(B_new, np.zeros(1))
            theta_new =  np.maximum(theta_new, np.zeros(1))
            
            A = A_new
            B =  B_new
            theta = theta_new
            
        time2 = time.time()
        
        time_list += [time2-time1 + time_list[-1]]
        
        ABMK_bar = get_ABMK_bar()
        y_hat = get_y_hat()
        evaluation = evaluate()
        loss += [evaluation]
        
        if evaluation >= evaluation_old + 0.01: 
            bad_epochs += 1
        
        if bad_epochs >= bad_epochs_max:
            stopped = True
            break
        
        if test:
            print("epoch "+ str(epoch) + " loss:" +str(evaluation))
        
        
    final_loss = evaluation
        
    return [[A, B, theta], final_loss, loss, stopped, time_list]



def RKHS_generalized_sketch_beta_divergnece(microbiome, R = 5, beta = 1.5,
                    lambda_parameter = 1000,
                    leanring_rate = 0.0001, 
                    n_max = 10, 
                    bad_epochs_max = 3,
                    s1 = 30, s2 = 40, s3 = 10, 
                    stop_steps = 3, 
                    iterations_per_epoch = 10,
                    bias_parameter = 0.5,
                    test = False,
                    decrease_lr = False,
                    gradient_clipping = False,
                    clipping_threshold = 0.1):

    def Beta_divergnece_loss_f(x,y):
        #x: estimated; y: true
        return_v = ((y+bias_parameter)**(beta+1) + beta * ((x+bias_parameter)** (beta + 1)) - (beta + 1)*(y+bias_parameter)*((x+bias_parameter)**beta))/(beta**2 + beta)
        #return_v[np.isnan(return_v)] = 0
        return return_v

    def Beta_divergnece_loss_gradient(x,y):
        return_v = (x+bias_parameter)**beta - (y+bias_parameter)*((x+bias_parameter)**(beta - 1))
        return_v[np.isnan(return_v)] = 0
        return return_v
        
        
    def Bernoulli_kernel(x, y):
        k1_x = x-0.5
        k1_y = y-0.5
        k2_x = 0.5*(k1_x**2-1/12)
        k2_y = 0.5*(k1_y**2-1/12)
        xy = abs(np.outer(x,np.ones(len(y))) - np.outer(np.ones(len(x)), y))
        k4_xy = ((xy-0.5)**4-0.5 * (xy-0.5)**2 + 7/240)/24
        kern_xy = np.outer(k1_x, k1_y) + np.outer(k2_x, k2_y) - k4_xy + 1
        return kern_xy
    
    def radial_kernel(x,y):
        if np.isscalar(x):
            x = np.array([x])
        if np.isscalar(y):
            y = np.array([y])
        nx = np.shape(x)[0]
        ny = np.shape(y)[0]
        result_matrix = np.zeros((nx,ny))
        for i in range(nx):
          for j in range(ny):
            result_matrix[i,j] = np.exp(-np.abs(x[i]-y[j]))
        return result_matrix


    def sampling():
        sample_index_n = choices(range(n), k = s1)
        sample_index_p = choices(range(p), k = s2)
        sample_index_Ti_list = [choices(range(tmp_index), k = s3) for tmp_index in dim_Ti_list]

        sample_T_list = []
        for i in sample_index_n:
            sample_T_list.append([T_list[i][j] for j in sample_index_Ti_list[i]])
         
        y_hat_sketched = np.zeros(p*sum_Ti) + float("nan")
        
        for i in sample_index_n:
            for j in sample_index_p:
                for t in sample_index_Ti_list[i]:
                    tmp_sum = 0
                    for r in range(R):
                        tmp_sum += A[i,r] * B[j,r] * MK[i,:] @ theta[(len(T)*r):(len(T)*(r+1))]
                    y_hat_sketched[j*sum_Ti + int(np.sum(dim_Ti_list[0:i])) + t] = tmp_sum[0]
        
        y_sketched = np.zeros(p*sum_Ti) + float("nan")
        
        for i in sample_index_n:
            for j in sample_index_p:
                for t in sample_index_Ti_list[i]:
                    y_sketched[j*sum_Ti + int(np.sum(dim_Ti_list[0:i])) + t] = y[j*sum_Ti + int(np.sum(dim_Ti_list[0:i])) + t]
        
        MK_hat = np.zeros(shape=(s1*s3, len(T)))

        tmp_sum = 0
        for i in range(s1):
            if i == 0: 
                MK_hat[0:s3, :] = radial_kernel(np.array(sample_T_list[i]), np.array(T))
            else:
                MK_hat[(tmp_sum):(s3 + tmp_sum), :] = radial_kernel(np.array(sample_T_list[i]), np.array(T))
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

        MK_bar_hat = np.kron(np.eye(R,dtype=int),MK_hat)

        y_sketched_for_theta = np.zeros(s2*s1*s3)

        for b_which in range(s2):
            b_index = sample_index_p[b_which]
            for a_which in range(s1):
                a_index = sample_index_n[a_which]
                for xi_which in range(s3):
                    xi_index = sample_index_Ti_list[a_index][xi_which]
                    tmp_index = int(b_which*s1*s3 + a_which*s3 + xi_which)
                    y_sketched_for_theta[tmp_index] = microbiome[a_index][b_index+1][xi_index]

        AB_hat = np.zeros(shape=(s2*s3*s1, R*s1*s3))

        for i in range(R):
            tmp_v = np.zeros(s3*s1)
            tmp_sum = 0
            for j in range(s1):
                tmp_v[tmp_sum:(tmp_sum+s3)] = A[sample_index_n,i][j]
                tmp_sum += s3
            AB_hat[:,(i*s3*s1):((i+1)*s3*s1)] = np.kron(B[sample_index_p,i][:, np.newaxis], np.diag(tmp_v))

        ABMK_bar_hat = AB_hat @ MK_bar_hat
        y_hat_sketched_for_theta = ABMK_bar_hat @ theta            
        
        return y_hat_sketched, y_sketched, y_sketched_for_theta, y_hat_sketched_for_theta, ABMK_bar_hat
                    
    def update_A_sketch():
        grad_A = np.zeros((n, R))
        for i in range(n):
            theta_matrix = theta.reshape((len(T),R))
            Xi_i = np.array(radial_kernel(T_list[i], np.array(T))) @ theta_matrix
            KR_matrix = scipy.linalg.khatri_rao(B, Xi_i)
            tmp_index1 = np.rint(np.sum(dim_Ti_list[0:i])).astype(int)
            tmp_index2 = np.rint(np.sum(dim_Ti_list[0:(i+1)])).astype(int)
            
            y_subset = np.zeros(dim_Ti_list[i]*p)
            y_hat_subset = np.zeros(dim_Ti_list[i]*p)
            
            for j in range(p):
                y_subset[(dim_Ti_list[i]*j) : (dim_Ti_list[i]*(j+1))] = y_sketched[(j*sum_Ti + tmp_index1):(j*sum_Ti + tmp_index2)]
                y_hat_subset[(dim_Ti_list[i]*j) : (dim_Ti_list[i]*(j+1))] = y_hat_sketched[(j*sum_Ti + tmp_index1):(j*sum_Ti + tmp_index2)]
            
            f_grad_A_rowvector = Beta_divergnece_loss_gradient(y_hat_subset, y_subset)
            grad_A[i,:] = f_grad_A_rowvector @ KR_matrix /omega_size
        
        #print(np.linalg.norm(grad_A))      
        
        if gradient_clipping: 
            norm_A = np.linalg.norm(grad_A)
            if norm_A>clipping_threshold:
                grad_A = clipping_threshold * grad_A / norm_A
        
        A_new = A - leanring_rate * grad_A * (omega_size / (s1*s2*s3))
        
        return A_new

    def update_B_sketch():
        f_grad_B_matrix = np.zeros((p, sum_Ti))

        for j in range(p):
            f_grad_B_matrix[j,:] = Beta_divergnece_loss_gradient(y_hat_sketched[(sum_Ti*j):(sum_Ti*(j+1))], 
                                                    y_sketched[(sum_Ti*j):(sum_Ti*(j+1))])

        theta_matrix = theta.reshape((len(T),R))
        Xi = np.array(radial_kernel(np.array(all_T), np.array(T))) @ theta_matrix
        KR_matrix = np.zeros((sum_Ti, R))

        for r in range(R):
            tmp_index_row = 0
            for i in range(n):
                KR_matrix[tmp_index_row:(tmp_index_row + dim_Ti_list[i]),r] = A[i, r] * Xi[tmp_index_row:(tmp_index_row + dim_Ti_list[i]), r]
                tmp_index_row += dim_Ti_list[i]
        
        grad_B = f_grad_B_matrix @ KR_matrix /omega_size
        
        if gradient_clipping: 
            norm_B = np.linalg.norm(grad_B)
            if norm_B>clipping_threshold:
                grad_B = clipping_threshold * grad_B / norm_B
        
        #print(np.linalg.norm(grad_B))   
        
        B_new = B - leanring_rate * grad_B * (omega_size / (s1*s2*s3))
        
        return B_new
                        
    def update_Xi_sketch():
        grad_theta_sketched = np.transpose(ABMK_bar_hat)@ Beta_divergnece_loss_gradient(y_hat_sketched_for_theta, y_sketched_for_theta[:, np.newaxis])/omega_size
        #print(np.linalg.norm(grad_theta_sketched))
        
        if gradient_clipping: 
            norm_Xi = np.linalg.norm(grad_theta_sketched)
            if norm_Xi>clipping_threshold:
                grad_theta_sketched = clipping_threshold * grad_theta_sketched / norm_Xi
        
        theta_new = theta - leanring_rate * grad_theta_sketched * (omega_size / (s1*s2*s3))
        return theta_new
                    
    def get_ABMK_bar():
        AB_new = np.zeros(shape=(p*sum_Ti, R*sum_Ti))
        for i in range(R):
            tmp_v = np.zeros(sum_Ti)
            tmp_sum = 0
            for j in range(n):
                tmp_v[tmp_sum:(tmp_sum+dim_Ti_list[j])] = A[:,i][j]
                tmp_sum += dim_Ti_list[j]
            AB_new[:,((i)*sum_Ti):((i+1)*sum_Ti)] = np.kron(B[:,i][:, np.newaxis], np.diag(tmp_v))

        ABMK_bar = AB_new @ MK_bar
        return ABMK_bar

    def get_y_hat():
        y_hat  = ABMK_bar @ theta
        return y_hat
                    
    def evaluate():
        loss = np.sum(Beta_divergnece_loss_f(y_hat.reshape(-1), y))/omega_size
        return loss
    
    def RKHS_norm_square(theta_tmp):
        return (np.transpose(theta_tmp)@K@theta_tmp)[0][0]


    def initial_A_B_theta():
        A = []
        for i in range(R):
            tmp_vec = [np.random.uniform(0, 0.01) for j in range(n)]
            tmp_vec_arr = np.array(tmp_vec)
            A.append(tmp_vec)
        A = np.transpose(np.array(A))    

        B = []
        for i in range(R):
            tmp_vec = [np.random.uniform(0, 0.01) for j in range(p)]
            tmp_vec_arr = np.array(tmp_vec)
            B.append(tmp_vec)
        B = np.transpose(np.array(B))
        
        theta_last = np.random.uniform(0, 0.01, len(T)*R)[:, np.newaxis]
        
        return A, B, theta_last

    #prepare
    T_list = []
    for i in range(len(microbiome)):
        T_list.append(np.array(microbiome[i][0])/739)

    def Union(lst1, lst2):
        final_list = list(set(lst1) | set(lst2))
        return final_list

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

    omega_size = sum_Ti * p

    K = radial_kernel(np.array(T), np.array(T))

    MK = np.zeros(shape=(sum_Ti, len(T)))

    tmp_sum = 0
    for i in range(n):
        if i == 0: 
            MK[0:dim_Ti_list[i], :] = radial_kernel(np.array(T_list[i]), np.array(T))
        else:
            MK[(tmp_sum):(dim_Ti_list[i] + tmp_sum), :] = radial_kernel(np.array(T_list[i]), np.array(T))
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

    MK_bar = np.kron(np.eye(R,dtype=int),MK)


    y = np.zeros(p*sum_Ti)

    for b_index in range(p):
        for a_index in range(n):
            for xi_index in range(dim_Ti_list[a_index]):
                tmp_index = int(b_index*sum_Ti + np.array(dim_Ti_list[0:(a_index)]).sum() + xi_index)
                y[tmp_index] = microbiome[a_index][b_index+1][xi_index]
                
     
                
    learning_rate_initial = leanring_rate
    #run 
    
    A, B, theta = initial_A_B_theta()
    
    ABMK_bar = get_ABMK_bar()
    y_hat = get_y_hat()
    
    loss = [evaluate()]

    stopped = False

    time_list = [0]

    bad_epochs = 0
    
    evaluation = np.Inf
    
    for epoch in range(n_max):
        evaluation_old = evaluation
        time1 = time.time()
        for iterator in range(iterations_per_epoch):
            
            if decrease_lr:
                leanring_rate = learning_rate_initial/((iterator+iterations_per_epoch*epoch+1)**2)
            
            
            y_hat_sketched, y_sketched, y_sketched_for_theta, y_hat_sketched_for_theta, ABMK_bar_hat = sampling()
            A_new = update_A_sketch()
            B_new = update_B_sketch()
            theta_new = update_Xi_sketch()
            
            for r in range(R):
                a_norm_tmp = np.linalg.norm(A_new[:,r])
                b_norm_tmp = np.linalg.norm(B_new[:,r])
                theta_tmp = theta_new[r*len(T):(r+1)*len(T),:]
                theta_norm_tmp = RKHS_norm_square(theta_tmp) ** 0.5
                
                norm_prod = a_norm_tmp*b_norm_tmp*theta_norm_tmp
                
                if norm_prod > lambda_parameter:
                    theta_tmp = theta_tmp/theta_norm_tmp * (norm_prod**(1/3))
                    theta_new[r*len(T):(r+1)*len(T),:] = theta_tmp
                    A_new[:,r] = A_new[:,r] / a_norm_tmp * (norm_prod**(1/3))
                    B_new[:,r] = B_new[:,r] / b_norm_tmp * (norm_prod**(1/3))
            
            A_new = np.maximum(A_new, np.zeros(1))
            B_new =  np.maximum(B_new, np.zeros(1))
            theta_new =  np.maximum(theta_new, np.zeros(1))
            
            A = A_new
            B =  B_new
            theta = theta_new
        
            
        time2 = time.time()
        
        time_list += [time2-time1 + time_list[-1]]
        
        ABMK_bar = get_ABMK_bar()
        y_hat = get_y_hat()
        evaluation = evaluate()
        loss += [evaluation]
        
        if evaluation >= evaluation_old + 0.01: 
            bad_epochs += 1
        
        if bad_epochs >= bad_epochs_max:
            stopped = True
            break
        
        if test:
            print("epoch "+ str(epoch) + " loss:" +str(evaluation))
        
        time1 = time.time()
    final_loss = evaluation
    
        
    return [[A, B, theta], final_loss, loss, stopped, time_list, y_hat]







def RKHS_generalized_beta_divergnece(microbiome, R = 5, beta = 0.5,
                    lambda_parameter = 1000,
                    leanring_rate = 0.0001, 
                    n_max = 10, 
                    bad_epochs_max = 3,
                    iterations_per_epoch = 10,
                    bias_parameter = 0.5,
                    test = False,
                    gradient_clipping = False,
                    clipping_threshold = 0.1):

    
    def Beta_divergnece_loss_f(x,y):
        #x: estimated; y: true
        return_v = ((y+bias_parameter)**(beta+1) + beta * ((x+bias_parameter)** (beta + 1)) - (beta + 1)*(y+bias_parameter)*((x+bias_parameter)**beta))/(beta**2 + beta)
        #return_v[np.isnan(return_v)] = 0
        return return_v

    def Beta_divergnece_loss_gradient(x,y):
        return_v = (x+bias_parameter)**beta - (y+bias_parameter)*((x+bias_parameter)**(beta - 1))
        return_v[np.isnan(return_v)] = 0
        return return_v
        
    def Bernoulli_kernel(x, y):
        k1_x = x-0.5
        k1_y = y-0.5
        k2_x = 0.5*(k1_x**2-1/12)
        k2_y = 0.5*(k1_y**2-1/12)
        xy = abs(np.outer(x,np.ones(len(y))) - np.outer(np.ones(len(x)), y))
        k4_xy = ((xy-0.5)**4-0.5 * (xy-0.5)**2 + 7/240)/24
        kern_xy = np.outer(k1_x, k1_y) + np.outer(k2_x, k2_y) - k4_xy + 1
        return kern_xy

    def radial_kernel(x,y):
        if np.isscalar(x):
            x = np.array([x])
        if np.isscalar(y):
            y = np.array([y])
        nx = np.shape(x)[0]
        ny = np.shape(y)[0]
        result_matrix = np.zeros((nx,ny))
        for i in range(nx):
          for j in range(ny):
            result_matrix[i,j] = np.exp(-(x[i]-y[j])**2)
        return result_matrix
                    
    def update_A():
        grad_A = np.zeros((n, R))
        for i in range(n):
            theta_matrix = theta.reshape((len(T),R))
            Xi_i = np.array(radial_kernel(T_list[i], np.array(T))) @ theta_matrix
            KR_matrix = scipy.linalg.khatri_rao(B, Xi_i)
            tmp_index1 = np.rint(np.sum(dim_Ti_list[0:i])).astype(int)
            tmp_index2 = np.rint(np.sum(dim_Ti_list[0:(i+1)])).astype(int)
            
            y_subset = np.zeros(dim_Ti_list[i]*p)
            y_hat_subset = np.zeros(dim_Ti_list[i]*p)
            
            for j in range(p):
                y_subset[(dim_Ti_list[i]*j) : (dim_Ti_list[i]*(j+1))] = y[(j*sum_Ti + tmp_index1):(j*sum_Ti + tmp_index2)]
                y_hat_subset[(dim_Ti_list[i]*j) : (dim_Ti_list[i]*(j+1))] = y_hat[(j*sum_Ti + tmp_index1):(j*sum_Ti + tmp_index2)].reshape(-1)
            
            f_grad_A_rowvector = Beta_divergnece_loss_gradient(y_hat_subset, y_subset)
            grad_A[i,:] = f_grad_A_rowvector @ KR_matrix /omega_size
        
        if gradient_clipping: 
            norm_A = np.linalg.norm(grad_A)
            if norm_A>clipping_threshold:
                grad_A = clipping_threshold * grad_A / norm_A    
        
        A_new = A - leanring_rate * grad_A
        
        #print("   grad_A:  " + str(np.linalg.norm(grad_A)))     
        
        return A_new

    def update_B():
        f_grad_B_matrix = np.zeros((p, sum_Ti))

        for j in range(p):
            f_grad_B_matrix[j,:] = Beta_divergnece_loss_gradient(y_hat[(sum_Ti*j):(sum_Ti*(j+1))].reshape(-1), 
                                                    y[(sum_Ti*j):(sum_Ti*(j+1))])

        theta_matrix = theta.reshape((len(T),R))
        Xi = np.array(radial_kernel(np.array(all_T), np.array(T))) @ theta_matrix
        KR_matrix = np.zeros((sum_Ti, R))

        for r in range(R):
            tmp_index_row = 0
            for i in range(n):
                KR_matrix[tmp_index_row:(tmp_index_row + dim_Ti_list[i]),r] = A[i, r] * Xi[tmp_index_row:(tmp_index_row + dim_Ti_list[i]), r]
                tmp_index_row += dim_Ti_list[i]
                
        grad_B = f_grad_B_matrix @ KR_matrix /omega_size
        
          
        
        if gradient_clipping: 
            norm_B = np.linalg.norm(grad_B)
            if norm_B>clipping_threshold:
                grad_B = clipping_threshold * grad_B / norm_B
        
        #print("   grad_B:  " + str(np.linalg.norm(grad_B)))   
        
        B_new = B - leanring_rate * grad_B
        
        return B_new
                        
    def update_Xi():
        grad_theta = np.transpose(ABMK_bar)@ Beta_divergnece_loss_gradient(y_hat, y[:, np.newaxis])/omega_size
        
        if gradient_clipping: 
            norm_Xi = np.linalg.norm(grad_theta)
            if norm_Xi>clipping_threshold:
                grad_theta = clipping_threshold * grad_theta / norm_Xi
        
        #print("   grad_theta:  " + str(np.linalg.norm(grad_theta)))   
        
        theta_new = theta - leanring_rate * grad_theta
        return theta_new               
                    
    def get_ABMK_bar():
        AB_new = np.zeros(shape=(p*sum_Ti, R*sum_Ti))
        for i in range(R):
            tmp_v = np.zeros(sum_Ti)
            tmp_sum = 0
            for j in range(n):
                tmp_v[tmp_sum:(tmp_sum+dim_Ti_list[j])] = A[:,i][j]
                tmp_sum += dim_Ti_list[j]
            AB_new[:,((i)*sum_Ti):((i+1)*sum_Ti)] = np.kron(B[:,i][:, np.newaxis], np.diag(tmp_v))

        ABMK_bar = AB_new @ MK_bar
        return ABMK_bar

    def get_y_hat():
        y_hat  = ABMK_bar @ theta
        return y_hat
                    
    def evaluate():
        loss = np.sum(Beta_divergnece_loss_f(y_hat.reshape(-1), y))/omega_size
        return loss
    
    def RKHS_norm_square(theta_tmp):
        return (np.transpose(theta_tmp)@K@theta_tmp)[0][0]


    def initial_A_B_theta():
        A = []
        for i in range(R):
            tmp_vec = [np.random.uniform(0, 0.01) for j in range(n)]
            tmp_vec_arr = np.array(tmp_vec)
            A.append(tmp_vec)
        A = np.transpose(np.array(A))    

        B = []
        for i in range(R):
            tmp_vec = [np.random.uniform(0, 0.01) for j in range(p)]
            tmp_vec_arr = np.array(tmp_vec)
            B.append(tmp_vec)
        B = np.transpose(np.array(B))
        
        theta_last = np.random.uniform(0, 0.01, len(T)*R)[:, np.newaxis]
        
        return A, B, theta_last

    #prepare
    T_list = []
    for i in range(len(microbiome)):
        T_list.append(np.array(microbiome[i][0])/739)

    def Union(lst1, lst2):
        final_list = list(set(lst1) | set(lst2))
        return final_list

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

    omega_size = sum_Ti * p

    K = radial_kernel(np.array(T), np.array(T))

    MK = np.zeros(shape=(sum_Ti, len(T)))

    tmp_sum = 0
    for i in range(n):
        if i == 0: 
            MK[0:dim_Ti_list[i], :] = radial_kernel(np.array(T_list[i]), np.array(T))
        else:
            MK[(tmp_sum):(dim_Ti_list[i] + tmp_sum), :] = radial_kernel(np.array(T_list[i]), np.array(T))
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

    MK_bar = np.kron(np.eye(R,dtype=int),MK)

    y = np.zeros(p*sum_Ti)

    for b_index in range(p):
        for a_index in range(n):
            for xi_index in range(dim_Ti_list[a_index]):
                tmp_index = int(b_index*sum_Ti + np.array(dim_Ti_list[0:(a_index)]).sum() + xi_index)
                y[tmp_index] = microbiome[a_index][b_index+1][xi_index]
                
                
    #run 
    A, B, theta = initial_A_B_theta()
    
    ABMK_bar = get_ABMK_bar()
    y_hat = get_y_hat()
    
    #y_hat = np.zeros(p*sum_Ti)
    #evaluate()
    
    loss = [evaluate()]

    stopped = False

    time_list = [0]

    bad_epochs = 0
    
    evaluation = np.Inf
    
    for epoch in range(n_max):
        evaluation_old = evaluation
        time1 = time.time()
        for iterator in range(iterations_per_epoch):
            ABMK_bar = get_ABMK_bar()
            y_hat = get_y_hat()
            A_new = update_A()
            B_new = update_B()
            theta_new = update_Xi()
            
            for r in range(R):
                a_norm_tmp = np.linalg.norm(A_new[:,r])
                b_norm_tmp = np.linalg.norm(B_new[:,r])
                theta_tmp = theta_new[r*len(T):(r+1)*len(T),:]
                theta_norm_tmp = RKHS_norm_square(theta_tmp) ** 0.5
                
                norm_prod = a_norm_tmp*b_norm_tmp*theta_norm_tmp
                
                if norm_prod > lambda_parameter:
                    theta_tmp = theta_tmp/theta_norm_tmp * (norm_prod**(1/3))
                    theta_new[r*len(T):(r+1)*len(T),:] = theta_tmp
                    A_new[:,r] = A_new[:,r] / a_norm_tmp * (norm_prod**(1/3))
                    B_new[:,r] = B_new[:,r] / b_norm_tmp * (norm_prod**(1/3))
            
            A_new = np.maximum(A_new, np.zeros(1))
            B_new =  np.maximum(B_new, np.zeros(1))
            theta_new =  np.maximum(theta_new, np.zeros(1))
            
            A = A_new
            B =  B_new
            theta = theta_new
            
        time2 = time.time()
        
        time_list += [time2-time1 + time_list[-1]]
        
        ABMK_bar = get_ABMK_bar()
        y_hat = get_y_hat()
        evaluation = evaluate()
        loss += [evaluation]
        
        if evaluation >= evaluation_old + 0.01: 
            bad_epochs += 1
        
        if bad_epochs >= bad_epochs_max:
            stopped = True
            break
        
        if test:
            print("epoch "+ str(epoch) + " loss:" +str(evaluation))
        
        
    final_loss = evaluation
        
    return [[A, B, theta], final_loss, loss, stopped, time_list]










