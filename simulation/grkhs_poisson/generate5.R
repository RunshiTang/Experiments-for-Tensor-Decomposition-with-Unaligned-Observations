setwd("C:/Users/97671/Box/HighD/RKHS/code/simulated/poisson")

set.seed(0)

n = 60
p = 51
T_min = 8
T_max = 20
R = 5

dim_Ti = round(runif(n, T_min, T_max))

T_hat = sample(1:739, 251, replace = F)

T_list = list()

for (i in 1:n){
  T_list[[i]] = sample(T_hat, dim_Ti[i], replace = F)
}

all_T = unique(unlist(T_list))

T_size = length(all_T)


radial.kernel = function(x,y, gamma = 1){
  nx = length(x)
  ny = length(y)
  result_matrix = matrix(nrow = nx, ncol = ny)
  for (i in 1:nx) {
    for (j in 1:ny){
      result_matrix[i,j] = exp(-(x[i]-y[j])^2)
    }
  }
  return(result_matrix)
}


Bernoulli.kernel <- function(x, y){
  k1.x = x-0.5
  k1.y = y-0.5
  k2.x = 0.5*(k1.x^2-1/12)
  k2.y = 0.5*(k1.y^2-1/12)
  xy = abs(x %*% t(rep(1,length(y))) - rep(1,length(x)) %*% t(y))
  k4.xy = 1/24 * ((xy-0.5)^4 - 0.5*(xy-0.5)^2 + 7/240)
  kern.xy = k1.x %*% t(k1.y) + k2.x %*% t(k2.y) - k4.xy + 1
  return(kern.xy)
}

K = radial.kernel(all_T/739, all_T/739)
K_bar = diag(x =1,nrow = R) %x% K

a = matrix(nrow = n, ncol = R)
b = matrix(nrow = p, ncol = R)
theta = matrix(nrow = T_size, ncol = R)

for (r in 1:R) {
  a[,r] = runif(n)
  b[,r] = runif(p)
}

tmp_v = runif(T_size)
theta[,1] = tmp_v/(as.vector((tmp_v) %*% K %*% (tmp_v)/100))^(0.5)

tmp_v[1:length(tmp_v)] = 0
tmp_v[1:10] = 1

theta[,2] = tmp_v/(as.vector((tmp_v) %*% K %*% (tmp_v)/100))^(0.5)

tmp_v[1:length(tmp_v)] = 0
tmp_v[20:30] = 1

theta[,3] = tmp_v/(as.vector((tmp_v) %*% K %*% (tmp_v)/100))^(0.5)

tmp_v[1:length(tmp_v)] = 0
tmp_v[50:60] = 1

theta[,4] = tmp_v/(as.vector((tmp_v) %*% K %*% (tmp_v)/100))^(0.5)

tmp_v[1:length(tmp_v)] = 0
tmp_v[90:100] = 1

theta[,5] = tmp_v/(as.vector((tmp_v) %*% K %*% (tmp_v)/100))^(0.5)

simulated_data_expected = list()

for (i in 1:n) {
  data_tmp = matrix(nrow = p+1, ncol = dim_Ti[i])
  data_tmp[1,] = T_list[[i]]
  for (j in 1:(p)){
    for (t in 1:dim_Ti[i]) {
      lambda = 0
      for (r in 1:R) {
        lambda = lambda + a[i,r] * b[j,r] * radial.kernel(T_list[[i]][t]/739, all_T/739) %*% theta[,r]
      }
      data_tmp[j+1,t] = lambda+1
    }
  }
  simulated_data_expected[[i]] = data_tmp
}

simulated_data = list()

Poisson_loss_f = function(x,y){
  return(x-y*log(x))
}


loss = 0

for (i in 1:n) {
  data_tmp = matrix(nrow = p+1, ncol = dim_Ti[i])
  data_tmp[1,] = T_list[[i]]
  for (j in 1:(p)){
    for (t in 1:dim_Ti[i]) {
      lambda = simulated_data_expected[[i]][j+1,t]
      data_tmp[j+1,t] = rpois(1, lambda)
      loss = loss + Poisson_loss_f(simulated_data_expected[[i]][j+1,t], data_tmp[j+1,t])
    }
  }
  simulated_data[[i]] = data_tmp
}

loss/(sum(dim_Ti)*p)





simulated_xi = function(x, r = 1){
  return((radial.kernel(x, all_T/739) %*% theta[,r]))
}
simulated_xi(0, r = 1)

x_v = 0:100/100

plot(x_v, simulated_xi(x_v, r = 1), type = "l", lty = 1)
plot(x_v, simulated_xi(x_v, r = 2), type = "l", lty = 1)



library(jsonlite)
write_json(simulated_data, "simulated_count5.json")
write_json(simulated_data_expected, "simulated_count5_expected.json", digits = 12)
