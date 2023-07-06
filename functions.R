library(matrixcalc) 
library(mvtnorm)
library(expm)

#function to simulate p A matrices of dimension m
#elements of each A matrix are randomly generated using independent normal distributions
A_matrix_simulation=function(m,p){
    A_mats=list()
    for(j in 1:p){
      A_mats[[j]]=matrix(rnorm(m^2,0,1),nrow=m,ncol=m)
    }
  return(A_mats)
}

#function to map from the A matrices to the partial autocorrelation matrices P
A_to_P=function(A,p,m){
  P=list()
  #repeat mapping for each matrix
  for(j in 1:p){
    curr_B=sqrtm(diag(rep(1,m))+A[[j]]%*%t(A[[j]]))
    P[[j]]=solve(curr_B)%*%A[[j]]
  }
  return(P)
}

#function to map from partial autocorrelation matrices P to phi matrices, using 
#mapping algorithm in Heaps (2023)
P_to_phi=function(Sigma,P,m,p){
  
  if(m ==1 ) {
    Sigma = matrix(Sigma, 1, 1)
    for(i in 1:p) {
      P[[i]] = matrix(P[[i]], 1, 1)
    }
  }
  
  Sigma_splusone=Sigma
  S_splusone=sqrtm(Sigma_splusone)
  Sigma_s=Sigma_splusone
  
  for(s in (p-1):0){
    Sigma_splusone=Sigma_s
    B_splusone_inv=sqrtm(diag(rep(1,m))-P[[s+1]]%*%t(P[[s+1]]))
    S_s=solve(B_splusone_inv)%*%sqrtm(Re(B_splusone_inv%*%Sigma_splusone%*%B_splusone_inv))%*%solve(B_splusone_inv)
    Sigma_s=S_s%*%t(S_s)
  }
  Gamma=list()
  Gamma[[1]]=Sigma_s
  
  Sigma_0=Gamma[[1]]
  Sigma_0_star=Gamma[[1]]
  S_0=sqrtm(Re(Sigma_0))
  S_0_star=S_0
  
  phisplusone=array(S_0%*%P[[1]]%*%solve(S_0_star),dim=c(m,m,1))
  phisplusone_star=array(S_0_star%*%t(P[[1]])%*%solve(S_0),dim=c(m,m,1))
  Sigma_splusone=Sigma_0-phisplusone[,,1]%*%Sigma_0_star%*%t(phisplusone[,,1])
  Sigma_splusone_star=Sigma_0_star-phisplusone_star[,,1]%*%Sigma_0%*%t(phisplusone_star[,,1])
  S_splusone=sqrtm(Re(Sigma_splusone))
  S_splusone_star=sqrtm(Re(Sigma_splusone_star))
  Gamma[[2]]=t(phisplusone[,,1]%*%Sigma_0_star)
  
  if(p>1){
    
    
    for(splusone in 2:p){
      phis=phisplusone
      phis_star=phisplusone_star
      phisplusone=array(dim=c(m,m,splusone))
      phisplusone_star=array(dim=c(m,m,splusone))
      Sigma_s=Sigma_splusone
      Sigma_s_star=Sigma_splusone_star
      S_s=S_splusone
      S_s_star=S_splusone_star
      
      phisplusone[,,splusone]=S_s%*%P[[splusone]]%*%solve(S_s_star)
      phisplusone_star[,,splusone]=S_s_star%*%t(P[[splusone]])%*%solve(S_s)
      
      for(i in 1:(splusone-1)){
        phisplusone[,,i]=phis[,,i]-phisplusone[,,splusone]%*%phis_star[,,splusone-i]
        phisplusone_star[,,i]=phis_star[,,i]-phisplusone_star[,,splusone]%*%phis[,,splusone-i]
      }
      Sigma_splusone=Sigma_s-phisplusone[,,splusone]%*%Sigma_s_star%*%t(phisplusone[,,splusone])
      Sigma_splusone_star=Sigma_s_star-phisplusone_star[,,splusone]%*%Sigma_s%*%t(phisplusone_star[,,splusone])
      S_splusone=sqrtm(Re(Sigma_splusone))
      S_splusone_star=sqrtm(Sigma_splusone_star)
      
      Gamma_temp=phisplusone[,,splusone]%*%Sigma_s_star
      for(i in 1:(splusone-1)){
        Gamma_temp=Gamma_temp+phis[,,i]%*%t(Gamma[[splusone-i+1]])
      }
      Gamma[[splusone+1]]=t(Gamma_temp)
    }
    
  }
  phi = list()
  for(i in 1:p){
    phi[[i]] = phisplusone[,,i]
  }
  list(phi=phi,Gamma=Gamma)
}


#check that the autogressive coefficients are in the stationary region
check_stationarity = function(phi) { 
  # Convert m-variate AR(p) model to mp-variate AR(1) model then
  # check all eigenvalues of mp x mp coefficient matrix (the 
  # companion matrix) are less than 1 in modulus.
  if(!is.list(phi)) stop("phi must be a list.")
  p = length(phi)
  for(i in 1:p) phi[[i]] = as.matrix(phi[[i]])
  m = nrow(phi[[1]])
  companion_mat = matrix(0, p * m, p * m)
  for(i in 1:p) {
    companion_mat[1:m, (i-1)*m + 1:m] = phi[[i]]
    if(i>1) {
      companion_mat[(i-1)*m + 1:m, (i-2)*m + 1:m] = diag(m)
    }
  }
  all(Mod(eigen(companion_mat)$values)<1)
}

#simulate a VAR process froma autoregressive coefficients
simulate_VAR = function(N, phi, Sigma, mu=rep(0, nrow(Sigma)), Gamma) { 
  # Checks on inputs
  if(!is.list(phi)) stop("phi must be a list")
  p = length(phi)
  if(is.matrix(phi[[1]]) & is.matrix(Sigma)) {
    m = nrow(Sigma)
    if(ncol(Sigma) != m) stop("Sigma must be square.")
    if(!isSymmetric(Sigma)) stop("Sigma must be symmetric.")
    if(!is.positive.definite(Sigma)) stop("Sigma must be positive definite.")
    for(i in 1:p) {
      if(nrow(phi[[i]])!=m || ncol(phi[[i]])!=m) stop("phis must be m x m matrices.")
    }
    if(length(mu)!=m) stop("mu must be an m-vector")
  }
  else {
    m = 1
    Sigma = matrix(Sigma, 1, 1)
    for(i in 1:p) phi[[i]] = matrix(phi[[i]], 1, 1)
  }
  # Check parameters lie in stationary region
  if(!check_stationarity(phi)) stop("Parameters do not lie in stationary region.")
  #param = phi2a(Sigma, phi, acf=TRUE, sym=TRUE)
  # Construct stationary distribution
  G = matrix(NA, p*m, p*m)
  for(i in 1:p) {
    for(j in 1:p) {
      if(i<=j) G[((i-1)*m+1):(i*m), ((j-1)*m+1):(j*m)] = Gamma[[j-i+1]]
      else G[((i-1)*m+1):(i*m), ((j-1)*m+1):(j*m)] = t(Gamma[[i-j+1]])
    }
  }
  y = matrix(NA, N, m)
  # Initialise in stationary distribution
  mut_init = numeric(p*m)
  mut_rest = matrix(0, N-p, m)
  for(t in 1:p) mut_init[((t-1)*m+1):(t*m)] = mu
  y[1:p,] = matrix(rmvnorm(1, mut_init, G), p, m, byrow=TRUE)
  # Simulate rest of process
  for(t in (p+1):N) {
    mut_rest[t-p,] = mu
    for(i in 1:p) {
      mut_rest[t-p,] = mut_rest[t-p,] + phi[[i]] %*% (y[t-i,] - mu)
    }
    y[t,] = rmvnorm(1, mut_rest[t-p,], Sigma)
  }
  return(list(y=y, mu=mu, phi=phi, Sigma=Sigma))
}

#generate a stationary VAR process and return model parameters and time series
generate_data_set = function(N, m, p, Sigma, mu) {
  A = A_matrix_simulation(m, p)
  P = A_to_P(A, p, m)
  phiGamma = P_to_phi(Sigma, P, m, p)
  phi = phiGamma[[1]]
  Gamma = phiGamma[[2]]

  y = simulate_VAR(N, phi, Sigma, mu, Gamma)
  y$A = A
  y$P = P
  return(y)
}

#calculate the threshold for the truncation criteria
calculate_threshold = function(m, N, beta) {
  q = qnorm((beta^(1/m^2)+1)/2)/sqrt(N)
  return(q)
}

##Function to apply truncation criterion when given MCMC output
eff_order = function(P, pmax, threshold, singval=FALSE) {
  iters = nrow(P)
  m = sqrt(ncol(P) / pmax)
  nonzero = matrix(NA, iters, pmax)
  for(iter in 1:iters) {
    Parr = array(P[iter,], c(pmax, m, m))
    if(!singval) nonzero[iter,] = apply(Parr, 1, function(mat) !all(abs(mat)<threshold))
    else nonzero[iter,] = apply(Parr, 1, function(mat) svd(mat)$d[1]>=threshold)
  }
  pstar = rowSums(nonzero)
  posterior_mass = tabulate(pstar+1, nbins=pmax+1) / nrow(P)
  names(posterior_mass) = paste("Pr(p*=", 0:pmax, "| ...)", sep="")
  return(list(nonzero=nonzero, samples=pstar, pmf=posterior_mass))
}
