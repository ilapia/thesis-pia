 data {
    int<lower=1> N;
    int<lower=1> Dx;
    matrix[N,Dx] x;
    int<lower=1> Ds;
    matrix[N,Ds] s;
    int<lower=0> y[N];
    vector[N] V;
  }
  transformed data {
    vector[N] mu;
    matrix[N, N] Dist_spatial;
    matrix[N, N] Sigma_lin;
    real s2_lin;
    s2_lin = 10;
    for (i in 1:N)
      mu[i] = 0;
    // off-diagonal elements
    for (i in 1:(N-1)) {
      for (j in (i+1):N) {
        Dist_spatial[i, j] = pow(dot_self(s[i] - s[j]),0.5)  ;    
        Sigma_lin[i, j] = s2_lin * dot_product(x[i],x[j]);     // linear covariance function
        
        // Fill in the other half
        Dist_spatial[j, i] = Dist_spatial[i, j];
        Sigma_lin[j, i] = Sigma_lin[i, j];
      }
    }
    // diagonal elements
    for (k in 1:N){
      Dist_spatial[k, k] = 0;
      Sigma_lin[k, k] = s2_lin * dot_product(x[k],x[k]) + 1e-6;   // add also some jitter
    }
  }
  parameters {
    real<lower=0> l;
    real<lower=0> s2_matern;
    vector[N] z;
    real<lower=0> r;
  }
  transformed parameters {
    matrix[N, N] Sigma;
    matrix[N, N] L;
    real<lower=0> inv_l;
  
    inv_l = inv(l);
  Sigma = s2_matern*exp(-inv_l*Dist_spatial ) + Sigma_lin; // Exponential   
   //     Sigma = s2_matern*(1 + pow(3,0.5)*inv_l*Dist_spatial).*exp(-pow(3,0.5)*inv_l*Dist_spatial)  + Sigma_lin ;  //Matern
    L = cholesky_decompose(Sigma);
  }
  model {
    vector[N] ff;
    
       // A weakly informative prior for magnitude
    s2_matern ~ student_t(4, 0, 1);
  
    // A weakly informative prior for l, that shrinks to 0 cor. among locations with more than 50km dist
    inv(l) ~ student_t(4, 0, 1);
    target += -2*log(l);   // The Jacobian part for the 1/l ~ Student_t
  
// A weakly informative prior for overdispersion
    r ~ gamma(2, 0.1 );
    
    z ~ normal(0, 1);
    ff = L*z;
  
    for (n in 1:N) {
      y[n] ~ neg_binomial_2(V[n]*exp(ff[n]), r);
    }
      
    
  }
  generated quantities {
    vector[N] f;
    // derived quantity (transform)
    f = L*z;
  }
