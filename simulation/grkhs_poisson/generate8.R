setwd("simulation/grkhs_poisson")

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
for (r in 1:R) {
  a[,r] = runif(n)
  b[,r] = runif(p)
}


set.seed(0)
theta = matrix(0, R, 10)
for (i in 1:10){
  theta[,i] = runif(R, -1/i, 1/i)
}


xi = function(x, r, theta){
  v = matrix(0, length(x), 10)
  for (j in 1:length(x)) {
    for (i in 1:10) {
      if (i %% 2 == 0){
        v[j,i] = theta[r,i] * sqrt(2) * cos((2-1)*pi*x[j])
      }else{
        v[j,i] = theta[r,i] * sqrt(2) * sin((2-1)*pi*x[j])
      }
    }
  }
  return(rowSums(v) + 1)
}

simulated_data_expected = list()

for (i in 1:n) {
  data_tmp = matrix(nrow = p+1, ncol = dim_Ti[i])
  data_tmp[1,] = T_list[[i]]
  for (j in 1:(p)){
    for (t in 1:dim_Ti[i]) {
      lambda = 0
      for (r in 1:R) {
        lambda = lambda + 10 * sqrt(r) * a[i,r] * b[j,r] * 
          xi(T_list[[i]][t]/739, r = r, theta = theta)
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



x_v = 0:100/100

library(ggplot2)

library(latex2exp)
la = c(TeX('$\\xi_1$'), 
       TeX('$\\xi_2$'),
       TeX('$\\xi_3$'),
       TeX('$\\xi_4$'),
       TeX('$\\xi_5$'))


ggplot() + 
  geom_line(aes(x = x_v, y = xi(x_v, r = 1, theta = theta), color = "Xi 1"))+
  geom_line(aes(x = x_v, y = xi(x_v, r = 2, theta = theta), color = "Xi 2"))+
  geom_line(aes(x = x_v, y = xi(x_v, r = 3, theta = theta), color = "Xi 3"))+
  geom_line(aes(x = x_v, y = xi(x_v, r = 4, theta = theta), color = "Xi 4"))+
  geom_line(aes(x = x_v, y = xi(x_v, r = 5, theta = theta), color = "Xi 5"))+
  scale_color_manual(labels = la, values = 1:5)+
  labs(x = "x",
       y = TeX('$\\xi(x)$'),
       color = NULL)
ggsave("xi_poisson.pdf", width = 4,height = 3)

library(jsonlite)
write_json(simulated_data, "simulated_count8.json")
write_json(simulated_data_expected, "simulated_count8_expected.json", digits = 12)
