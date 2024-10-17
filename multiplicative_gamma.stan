
functions{
  //function to calculate the symmetric matrix-square-root of a matrix
  matrix sqrtm(matrix A) { 
    int m = rows(A);
    tuple(matrix[m, m], vector[m]) eigen = eigendecompose_sym(A);
    matrix[m, m] eprod = diag_post_multiply(eigen.1, sqrt(sqrt(eigen.2)));
    return tcrossprod(eprod);
  }
  
  //function to map from the transformed partial autocorrelation
  //matrices A to the partial autocorrelation matrices P
  matrix A_to_P(matrix A){
    int m = cols(A);
    matrix[m,m] C;
    matrix[m,m] B;
    matrix[m,m] P;
    C=diag_matrix(rep_vector(1,m))+tcrossprod(A);
    B=sqrtm(C);
    P=mdivide_left_spd(B,A);
    return P;
  }
  
  //function to map from the partial autocorrelation matrices
  //P to the autogressive matrices phi
  //also outputs the theoretical autocovariances
  array[,] matrix P_to_AR(array[] matrix P,matrix Sigma){
    int p = size(P);
    int m = cols(Sigma);
    array[p+1] matrix[m,m] Sigma_s;
    array[p+1] matrix[m,m] S_s;
    array[p+1] matrix[m,m] Sigma_star_s;
    array[p+1] matrix[m,m] S_star_s;
    array[p+1] matrix[m,m] Gamma;
    array[2,p] matrix[m,m] phiGamma;
    matrix[m,m] Y;
    matrix[m,m] sqY;
    array[p] matrix[m,m] phi_splusone;
    array[p] matrix[m,m] phi_splusone_star;
    // initialise
    Sigma_s[p+1] = Sigma;
    S_s[p+1]=sqrtm(Sigma_s[p+1]);


    //find gamma_0 using mapping algorithm
    for(s in 0:(p-1)){
      Y=diag_matrix(rep_vector(1,m))-tcrossprod(P[p-s]);
      sqY=sqrtm(Y);
      S_s[p-s] = mdivide_right_spd(mdivide_left_spd(sqY, 
                              sqrtm(quad_form_sym(Sigma_s[p-s+1], sqY))), 
                                                                    sqY);
      Sigma_s[p-s] = tcrossprod(S_s[p-s]);
    }
    Gamma[1]=Sigma_s[1];
    Sigma_star_s[1]=Gamma[1];
    S_star_s[1]=S_s[1];
    phi_splusone[1]=mdivide_right(S_s[1]*P[1],S_star_s[1]);
    phi_splusone_star[1]=mdivide_right(S_star_s[1]*P[1]',S_s[1]);
    Sigma_s[2]=Sigma_s[1]-quad_form_sym(Sigma_star_s[1],phi_splusone[1]');
    Sigma_star_s[2]=Sigma_star_s[1]-
                          quad_form_sym(Sigma_s[1],phi_splusone_star[1]');
    S_s[2]=sqrtm(Sigma_s[2]);
    S_star_s[2]=sqrtm(Sigma_star_s[2]);
    Gamma[2]=(phi_splusone[1]*Sigma_star_s[1])';
    //return initialisation if p=1
    if(p==1){
      
      phiGamma[1,p]=phi_splusone[1];
      phiGamma[2,p]=Gamma[1];
      return phiGamma;
    }
    //proceed with mapping algorithm if p>1
    else{
      //store theoretical autocovariance for i=0,...,p-1
      array[p] matrix[m,m] store;
      array[p] matrix[m,m] store_star;
      
      store[1]=phi_splusone[1];
      store_star[1]=phi_splusone_star[1];
      
      //rest of algorithm
      for(splusone in 2:p){
        //phi_s in algorithm
        array[splusone-1] matrix[m,m] phi_s;
        array[splusone-1] matrix[m,m] phi_s_star;
        //used to store gamma in reverse order
        matrix[m,m] Gamma_temp;
        matrix[m,m] matrix_sum;
        //Initialise phi_s
        phi_s=store[1:(splusone-1)];
        phi_s_star=store_star[1:(splusone-1)];
        //from algorithm
        phi_splusone[splusone]=mdivide_right(S_s[splusone]*P[splusone],
                                                   S_star_s[splusone]);
        phi_splusone_star[splusone]=mdivide_right(S_star_s[splusone]*
                                           P[splusone]',S_s[splusone]);
        for(i in 1:(splusone-1)){
          //From algorithm
          phi_splusone[i]=phi_s[i]-
                          phi_splusone[splusone]*phi_s_star[splusone-i];
          phi_splusone_star[i]=phi_s_star[i]-
                          phi_splusone_star[splusone]*phi_s[splusone-i];
        }
        //stores phi_{s+1}
        store[1:splusone]=phi_splusone[1:splusone];
        store_star[1:splusone]=phi_splusone_star[1:splusone];
        //from algorithm
        Sigma_s[splusone+1]=Sigma_s[splusone]-
           quad_form_sym(Sigma_star_s[splusone],phi_splusone[splusone]');
        S_s[splusone+1]=sqrtm(Sigma_s[splusone+1]);
        Sigma_star_s[splusone+1]=Sigma_star_s[splusone]-
           quad_form_sym(Sigma_s[splusone],phi_splusone_star[splusone]');
        S_star_s[splusone+1]=sqrtm(Sigma_star_s[splusone+1]);
        //reverse s gamma values
        matrix_sum = phi_splusone[splusone]*Sigma_star_s[splusone];
        for(i in 1:(splusone-1)){
          Gamma_temp=Gamma[splusone-i+1];
          matrix_sum += phi_s[i]*Gamma_temp';
        }
        //from algorithm
        
        Gamma[splusone+1]=matrix_sum';
      }
      //output is vector containing p autoregressive parameters followed
      //by Gamma 0,...,p-1
      phiGamma[1,]=phi_splusone;
      phiGamma[2,]=Gamma[1:p];
      return phiGamma;
    }
  }

  
}


data {
  int<lower=1> m; //Dimension of VAR process
  int<lower=0> N; //Length of time series
  int<lower=1> p; //Maximum order
  int<lower=0> n_miss; //Number of missing data points
  array[N] vector[m] y;    //Time series
  array[n_miss] int<lower=1> ind_miss; //Indices of missing data points
  
  //Hyperparameters in inverse Wishart prior for Sigma
  int<lower=m+3> df; // Degrees of freedom 
                     //(limit ensures variance is finite)
  real<lower=0> scale_diag;  // Diagonal element in scale matrix
  real<lower=-scale_diag/(m-1)> scale_offdiag; /* Off-diagonal element in 
                                                  scale matrix */  
  //Hyperparameters in multiplicative gamma process prior for 
  //elements of A matrices
  real<lower=0> a;
  real<lower=0> a1;
  real<lower=0> a2;
}



transformed data {
  vector[m] mu;
  matrix[m, m] scale_mat;  // Scale-matrix in prior for Sigma
  mu = rep_vector(0.0, m); //Zero mean of VAR process
  for(i in 1:m) {
    for(j in 1:m) {
      if(i==j) scale_mat[i, j] = scale_diag;
      else scale_mat[i, j] = scale_offdiag;
    }
  }
}
 
 
parameters {
  array[p] matrix[m,m] A; //The A matrices
  cov_matrix[m] Sigma; //Error variance, Sigma
  array[n_miss] vector[m] y_miss; //Missing data points
  //Parameters in multiplicative gamma process prior for the elements 
  //of the A matrices
  vector<lower=0>[p] delta;
  array[p] matrix<lower=0>[m,m] lambda;
}



transformed parameters {
  array[p] matrix[m,m] P; //Partial autocorrelation matrices, P_s
  array[p] matrix[m,m] phi; //The phi_i
  array[2,p] matrix[m,m] phiGamma; //place to store output 
                                   //of mapping function
  array[p] matrix[m,m] Gamma; //Autocovariances
  cov_matrix[m*p] G; //Stationary variance of (y_1, ..., y_p)
  vector<lower=0>[p] tau; //Parameter in multiplicative gamma process prior
  vector[p*m] y_init; //Initial values
  vector[p*m] mu_long; //mp-vector with mu repeated p times
  array[N] vector[m] y_complete; //Data with missing data inferred
  y_complete = y;
  if(n_miss>0) {
    for(i in 1:n_miss){
      y_complete[ind_miss[i]] = y_miss[i];
    }
  }
  
  for(i in 1:p){
    y_init[(1+(i-1)*m):i*m]=y_complete[i];
  }
  for(i in 1:p){
    P[i]=A_to_P(A[i]);
  }
  for(i in 1:p){
    mu_long[(1+(i-1)*m):i*m] = mu;
  }
  phiGamma=P_to_AR(P,Sigma);
  phi=phiGamma[1,];
  Gamma=phiGamma[2,];
  
  for(i in 1:p){
      for(j in 1:p){
        if(i<=j){
          G[(1+m*(i-1)):m*i,(1+m*(j-1)):m*j]=Gamma[j-i+1];
        }
        else{
          G[(1+m*(i-1)):m*i,(1+m*(j-1)):m*j]=Gamma[i-j+1]';
        }
        
      }
    }
    {
      vector[p] logtau;
      logtau[1]=log(delta[1]);
      for(j in 2:p){
        logtau[j]=logtau[j-1]+log(delta[j]);
      }
      tau = exp(logtau);
    }
  
}



model {
  y_init~multi_normal(mu_long,G); //initial distribution

  //for distribution when t>p
  for(t in (p+1):N){
    vector[m] AR_mean; //conditional mean of y_t
    array[p] vector[m] ys; //to store previous p values in time series
    for(i in 1:p){
      ys[i]=y_complete[t-i];
    }
    AR_mean=mu;
    for(i in 1:p){
      AR_mean=AR_mean+phi[i]*(ys[i] - mu);
    }
    y_complete[t]~multi_normal(AR_mean,Sigma); //Distribution when t>p
  }
  
  // Multiplicative gamma process prior for A matrices
  for(s in 1:p){
    for(i in 1:m){
      for(j in 1:m){
        lambda[s,i,j]~gamma(a/2,a/2);
        A[s,i,j]~normal(0,1/sqrt(lambda[s,i,j]*tau[s]));
      }
    }
    
  }
  delta[1]~gamma(a1,1);
  for(i in 2:p){
    delta[i]~gamma(a2,1);
  }
  //Inverse Wishart prior for Sigma
  Sigma ~ inv_wishart(df,scale_mat);
}