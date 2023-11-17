library(gtools)

format_tfpca <- function(taxon_table, time_point, subjectID, threshold=0.8, pseudo_count=0.5){
  # format data table into a list as input for tfsvd(), 
  # read counts all have 1/2 added, before being normalized into proportions and log transformation
  # check length of subID and time_point
  if (length(subjectID)!=nrow(taxon_table)) 
    stop('length of subjectID does not match taxon_table!')
  if (length(time_point)!=nrow(taxon_table)) 
    stop('length of time_point does not match taxon_table!')
  # keep taxon that has non-zeros in >=1-threshold samples
  taxon_table <- taxon_table[,colMeans(taxon_table==0)<threshold]
  taxon_table <- taxon_table+pseudo_count
  taxon_table <- t(log(taxon_table/rowSums(taxon_table)))
  taxon_table <- rbind(time_point, taxon_table)
  rownames(taxon_table)[1] <- 'time_point'
  subID <- unique(subjectID)
  nsub <- length(subID)
  
  # construct list of data matrices, each element representing one subject
  datlist <- vector("list", length = nsub)
  names(datlist) <- subID
  
  # Each slice represents an individual (unequal sized matrix).
  for (i in 1:nsub){
    # print(i)
    datlist[[i]] <- taxon_table[, subjectID==subID[i]]
    datlist[[i]] <- datlist[[i]][,order(datlist[[i]][1,])]
    datlist[[i]] <- datlist[[i]][,!duplicated(datlist[[i]][1,])]
    colnames(datlist[[i]]) <- NULL
  }
  return(datlist)
}



# RKHS functional kernel calculations
#########################
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

interpolate <- function(x){
  n = length(x)
  i = 1
  if(is.na(x[i])){
    while(is.na(x[i])) i = i + 1
    x[1:(i-1)] = x[i]
  }
  i = n
  if(is.na(x[i])){
    while (is.na(x[i])) i = i - 1
    x[(i+1):n] = x[i]
  }
  for (i in 1:n) {
    if (is.na(x[i])){
      imin = imax = i
      while(is.na(x[imin])) imin = imin - 1
      while(is.na(x[imax])) imax = imax + 1
      x[i] = (x[imax]-x[imin])/(imax-imin)*(i-imin) + x[imin]
    }
  }
  return(x)
}


#' Power iteration for estimating single singular component. The functional domains
#' are assumed to be from 0 to 1.
#' @param f.tensor An order-d tabular tensors.
#' @param f.grid A list recording the grids for all functional modes.
#' @param init A list of initializations for each mode. 
#' @param resolution A list of number and NULL indicating the desired functional resolutions.
#' @param penalty The tuning parameter in RKHS regularized empirical risk minimization.
#' @param t_max The maximum number for power iterations.
FTSVD.iteration <- function(f.tensor, f.grid, init, resolution, penalty, t_max = 10){
  # Check the input validity.
  if(class(f.tensor) != "Tensor") stop("invalid input: f.tensor should be of tensor type.")
  p = dim(f.tensor)
  d = length(p)
  if(d != length(f.grid)) stop("The order of f.tensor should be the same as the length of f.grid.")
  if(d != length(init)) stop("The order of f.tensor should be the same as the length of init.")
  for (k in 1:d) {
    if (is.null(f.grid) & length(f.grid[k]) != p[k]) stop("The sizes of tensor and grid does not match.")
    if (is.null(f.grid) & length(init[k]) != p[k]) stop("The sizes of tensor and initializations does not match.")
  }
  
  est = init
  
  # prepare the kernel and its inverse matrices.
  Kers = list()
  Kers.inverse = list()
  for(k in 1:d){
    if (is.null(f.grid[[k]])){
      Kers = c(Kers, list(NULL))
      Kers.inverse = c(Kers.inverse, list(NULL))
    }
    else{
      Kernel.matrix = Bernoulli.kernel(f.grid[[k]], f.grid[[k]])
      Kers = c(Kers, list(Kernel.matrix))
      temp = eigen(Kernel.matrix, symmetric = TRUE)
      temp$value[temp$value < 0] = 0
      temp$value = 1/(temp$value + penalty)
      Kers.inverse = c(Kers.inverse, list((temp$vector)%*%(t(temp$vector)*temp$value)))
    }
  }
  
  loadings = list()
  for (k in 1:d) {
    loadings = c(loadings, list(t(as.matrix(init[[k]]))))
  }
  
  # power iteration
  for (t in 1:t_max) {
    for (k in 1:length(p)) {
      # tabular mode
      if (is.null(f.grid[[k]])){

        uk = ttl(f.tensor, loadings[-k], (1:d)[-k])
        uk = k_unfold(uk, k)@data
        uk = uk / norm(uk, type='F')
        loadings[[k]] = t(uk)
        
        if (t==t_max){
          est[[k]] = uk
        }
      }
      # functional mode
      else {
        uk = ttl(f.tensor, loadings[-k], (1:d)[-k])
        uk = k_unfold(uk, k)@data
        beta = Kers.inverse[[k]] %*% uk
        xi.n = Kers[[k]] %*% beta
        loadings[[k]] = t(xi.n / norm(xi.n, type='F'))
        
        if (t==t_max){
          uniform.grid = seq(0,1,length.out = resolution[[k]])
          K = Bernoulli.kernel(uniform.grid, f.grid[[k]])
          xi = K %*% beta
          est[[k]] = xi / norm(xi, type = 'F')
          
        }
      }
    }
  }
  lambda = ttl(f.tensor, loadings, (1:d))@data
  for (k in 1:d){
    loadings[[k]] = t(loadings[[k]])
  }
  est.tensor = ttl(as.tensor(lambda), loadings, 1:d)
  return(c(list(lambda),est,list(est.tensor)))
}


#' Rank-1 initialization for tabular-functional tensor data
#' @param f.tensor An order-d tabular tensor.
#' @param f.grid A list recording the grids for all functional modes.
#' @return A list of loading vectors for all d modes.
FTSVD.spectral <- function(f.tensor, f.grid, penalty){
  # Check the input validity.
  if(class(f.tensor) != "Tensor") stop("invalid input: f.tensor should be of tensor type.")
  p = dim(f.tensor)
  d = length(p)
  if(d != length(f.grid)) stop("The order of f.tensor should be the same as the length of f.grid.")
  
  init = list()
  
  # prepare the kernel and its inverse matrices.
  Kers = list()
  Kers.inverse = list()
  for(k in 1:d){
    if (is.null(f.grid[[k]])){
      Kers = c(Kers, list(NULL))
      Kers.inverse = c(Kers.inverse, list(NULL))
    }
    else{
      Kernel.matrix = Bernoulli.kernel(f.grid[[k]], f.grid[[k]])
      Kers = c(Kers, list(Kernel.matrix))
      temp = eigen(Kernel.matrix, symmetric = TRUE)
      temp$value[temp$value < 0] = 0
      temp$value = 1/(temp$value + penalty)
      Kers.inverse = c(Kers.inverse, list((temp$vector)%*%(t(temp$vector)*temp$value)))
    }
  }
  
  for (k in 1:d) {
    # tabular mode
    if (is.null(f.grid[[k]])){
      uk = svd(k_unfold(f.tensor, k)@data)$u[,1]
      init = c(init, list(uk))
    }
    # functional mode: ERM
    else{
      temp = svd(k_unfold(f.tensor, k)@data)
      uk = temp$d[1] * temp$u[,1]
      beta = Kers.inverse[[k]] %*% uk
      xi = Kers[[k]] %*% beta
      xi = as.vector( xi / norm(xi, type = 'F'))
      init = c(init, list(xi))
    }
  }
  return(init)
}



#' Main function for functional tensor svd.
#' @param f.tensor An order-d tabular tensor.
#' @param f.grid A list recording the grids for all functional modes.
#' @param rank Tensor rank.
#' @param resolution Resolution level for estimated functions.
#' @param penalty Tuning parameter in regularized ERM.
#' @param L Number of random initializations.
FTSVD <- function(f.tensor, f.grid, rank, resolution, penalty, L = 20){
  # Check the input validity.
  if(class(f.tensor) != "Tensor") stop("invalid input: f.tensor should be of tensor type.")
  p = dim(f.tensor)
  d = length(p)
  if(d != length(f.grid)) stop("The order of f.tensor should be the same as the length of f.grid.")
  
  if(missing("resolution")) {
    resolution = rep(list(NULL), d)
    for (k in 1:d) {
      if (!is.null(f.grid[[k]])){
        resolution[[k]] = 101
      }
    }
  }
  
  if(length(resolution)!=d) stop("The length of resolution does not match tensor orders.")
  for (k in 1:d) {
    if(xor(is.null(f.grid[[k]]),is.null(resolution[[k]])))
      stop("Functional mode indicators does not match for f.grid and resolution.")
  }
  
  if (rank==1){
    # spectral initialization + power iteration for rank-1 estimation
    init = FTSVD.spectral(f.tensor, f.grid, penalty)
    results = FTSVD.iteration(f.tensor, f.grid, init, resolution, penalty)
    return(results)
  }
  else{
    # obtain the iterated result from spectral initializations
    
    #a.init = svd(k_unfold(f.tensor,1)@data)$u[,1:rank]
      
    lambda = rep(0, rank)
    loadings = rep(list(NULL),d)
    for (s in 1:rank) {
      a = svd(k_unfold(f.tensor,1)@data)$u[,1]
      b = svd(k_unfold(ttm(f.tensor,t(a),1),2)@data)$u[,1]
        
        # prepare the kernel and its inverse matrices.
      Kers = list()
      Kers.inverse = list()
      for(k in 1:d){
        if (is.null(f.grid[[k]])){
          Kers = c(Kers, list(NULL))
          Kers.inverse = c(Kers.inverse, list(NULL))
        }
        else{
          Kernel.matrix = Bernoulli.kernel(f.grid[[k]], f.grid[[k]])
          Kers = c(Kers, list(Kernel.matrix))
          temp = eigen(Kernel.matrix, symmetric = TRUE)
          temp$value[temp$value < 0] = 0
          temp$value = 1/(temp$value + penalty)
          Kers.inverse = c(Kers.inverse, list((temp$vector)%*%(t(temp$vector)*temp$value)))
        }
      }
      beta = Kers.inverse[[k]] %*% k_unfold(ttl(f.tensor,c(list(t(a)),list(t(b))),c(1,2)),3)@data
      xi = Kers[[k]] %*% beta
      xi = xi / norm(xi, type = 'F')
      s.init = c(list(a),list(b),list(xi))
      
      s.init = FTSVD.spectral(f.tensor, f.grid, penalty)
      s.est = FTSVD.iteration(f.tensor, f.grid, s.init, resolution, penalty)
      lambda[s] = s.est[[1]]
      for (k in 1:(length(s.est)-2)) {
        loadings[[k]] = cbind(loadings[[k]], s.est[[k+1]])
      }
        
      f.tensor = f.tensor - s.est[[k+2]]
    }
    return(results = c(list(lambda), loadings))
    
  }
}

CP <- function(f.tensor, f.grid, rank, resolution = 101, L = 20){
  # Check the input validity.
  if(class(f.tensor) != "Tensor") stop("invalid input: f.tensor should be of tensor type.")
  p = dim(f.tensor)
  d = length(p)
  if(d != length(f.grid)) stop("The order of f.tensor should be the same as the length of f.grid.")
  
  if(missing("resolution")) {
    resolution = rep(list(NULL), d)
    for (k in 1:d) {
      if (!is.null(f.grid[[k]])){
        resolution[[k]] = 101
      }
    }
  }
  
  Kers = list()
  Kers.inverse = list()
  for(k in 1:d){
    if (is.null(f.grid[[k]])){
      Kers = c(Kers, list(NULL))
      Kers.inverse = c(Kers.inverse, list(NULL))
    }
    else{
      Kernel.matrix = Bernoulli.kernel(f.grid[[k]], f.grid[[k]])
      Kers = c(Kers, list(Kernel.matrix))
      temp = eigen(Kernel.matrix, symmetric = TRUE)
      temp$value[temp$value < 0] = 0
      temp$value = 1/(temp$value + 1e-10)
      Kers.inverse = c(Kers.inverse, list((temp$vector)%*%(t(temp$vector)*temp$value)))
    }
  }
  
  if (rank==1){
    # spectral initialization
    init = NULL
    for (k in 1:d) {
      uk = svd(k_unfold(f.tensor, k)@data)$u[,1]
      init = c(init, list(uk))
    }

    est = init
    # prepare the kernel and its inverse matrices.
    
    loadings = list()
    for (k in 1:d) {
      loadings = c(loadings, list(t(as.matrix(init[[k]]))))
    }
    
    # power iteration
    for (t in 1:10) {
      for (k in 1:length(p)) {
        # tabular mode
        uk = ttl(f.tensor, loadings[-k], (1:d)[-k])
        uk = k_unfold(uk, k)@data
        uk = uk / norm(uk, type='F')
        loadings[[k]] = t(uk)
          
        if (t==10){
          est[[k]] = uk
        }
      }
    }
    lambda = ttl(f.tensor, loadings, (1:d))@data
    for (k in 1:d) {
      if(!is.null(f.grid[[k]])){
        uniform.grid = seq(0,1,length.out = resolution[[k]])
        temp = Bernoulli.kernel(uniform.grid, f.grid[[k]]) %*% Kers.inverse[[k]] %*% est[[k]]
        temp = temp / norm(temp, type='F')
        est[[k]] = temp
      }
    }
    return(c(list(lambda),est))
  }
  else{
    lambda.set = rep(0,L)
    loadings.set = rep(list(NULL),L)
    for (l in 1:L) {
      l.init = NULL
      Y = f.tensor
      for (k in 1:d) {
          if (k==1){
            # random initialization
            uk = as.matrix(rnorm(dim(Y)[1]))
            uk = uk / norm(uk,type = 'F')
          }
          else{
            uk = svd(k_unfold(Y, k)@data)$u[,1]
          }
          l.init = c(l.init, list(uk))
          Y = ttm(Y, t(uk), k)
      }
      loadings = list()
      for (k in 1:d) {
        loadings = c(loadings, list(t(as.matrix(l.init[[k]]))))
      }
      # power iteration
      l.est = l.init
      for (t in 1:10) {
        for (k in 1:length(p)) {
          # tabular mode
          uk = ttl(f.tensor, loadings[-k], (1:d)[-k])
          uk = k_unfold(uk, k)@data
          uk = uk / norm(uk, type='F')
          loadings[[k]] = t(uk)
          
          if (t==10){
            l.est[[k]] = uk
          }
        }
      }
      
      lambda.set[l] = ttl(f.tensor, loadings, (1:d))@data
      loadings.set[[l]] = l.est
    }
    
    lambda = rep(0, rank)
    loadings = rep(list(NULL),d)
    # print(lambda.set)
    # clustering
    for (s in 1:rank) {
      l.max = which.max(abs(lambda.set))
      lambda[s] = lambda.set[l.max]
      for (k in 1:d) {
        loadings[[k]] = cbind(loadings[[k]], loadings.set[[l.max]][[k]])
      }
      for (l in 1:L) {
        for (k in 1:d) {
          if (is.null(f.grid[[k]]) & abs(t(loadings.set[[l]][[k]]) %*% loadings.set[[l.max]][[k]]) > 0.5)
            lambda.set[l] = 0
        }
      }
    }
    
    # transform the functional mode to functional.
    for (k in 1:d) {
      if(!is.null(f.grid[[k]])){
        uniform.grid = seq(0,1,length.out = resolution[[k]])
        for (s in 1:r) {
          functional_loading = NULL
          temp = Bernoulli.kernel(uniform.grid, f.grid[[k]]) %*% Kers.inverse[[k]] %*% loadings[[k]][,s]
          temp = temp / norm(temp, type='F')
          functional_loading = cbind(functional_loading, temp)
        }
        loadings[[k]] = functional_loading
      }
    }
    
    
    return(results = c(list(lambda), loadings))
    
  }
}


xi.generate <- function(tn,Q){
  xi = rep(0,101)
  xin = rep(0, n)
  for (q in 1:Q) {
    alpha = (runif(1) - 0.5) / q
    if (q==1) {
      xi = xi + alpha
      xin = xin + alpha
    }
    else {
      xi = xi + alpha * sqrt(2) * cos((q-1) * pi * seq(0,1,length.out = 101))
      xin = xin + alpha * sqrt(2) * cos((q-1) * pi * tn)
    }
  }
  xi = as.matrix(xi); xin = as.matrix(xin)
  xin = xin / norm(xi, type = 'F')
  xi = xi / norm(xi, type = 'F')
  return(c(list(xin),list(xi)))
}


dist.measure <- function(A,B){
  r = ncol(A)
  perm = permutations(r, r, 1:r)
  loss = 1
  for (i in 1:nrow(perm)) {
    B.temp = as.matrix(B[, perm[i,]])
    dist.max = 0
    for (s in 1:r) {
      dist.max = dist.max + sqrt(1 - abs(t(A[,s])%*%B.temp[,s])^2)
    }
    dist.max = dist.max / r
    loss = min(loss, dist.max) 
  }
  return(loss)
}
