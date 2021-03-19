data {
  int<lower=0> n_K;
  int<lower=0> K_N;
  int<lower=0> K_M;
  
  int<lower=0> n_Z;
  int<lower=0> Z_X;
  int<lower=0> Z_Y;
}

parameters {
  real<lower=0,upper=1> rho;
  real<lower=0> M_K;
  real<lower=0> M_Z;
}

transformed parameters {
  real<lower=0> lambda_N0;
  real<lower=0> lambda_M0;
  real<lower=0> lambda_N1;
  real<lower=0> lambda_M1;
  real<lower=0> lambda_X;
  real<lower=0> lambda_Y;
  real<lower=0,upper=1> p_N0;
  real<lower=0,upper=1> p_M0;
  real<lower=0,upper=1> p_N1;
  real<lower=0,upper=1> p_M1;
  real<lower=0,upper=1> p_X;
  real<lower=0,upper=1> p_Y;
  
  vector[2] lp_N;
  vector[2] lp_M;
  
  // Concentrations of N and M if foetus heterozygous
  lambda_N0 = M_K / 2;
  lambda_M0 = M_K / 2;
  
  // Concentrations of N and M if foetus homozygous
  lambda_N1 = M_K * (1 + rho) / 2;
  lambda_M1 = M_K * (1 - rho) / 2;
  
  // Concentrations of X and Y
  lambda_X = M_Z * (2 - rho) / 2;
  lambda_Y = M_Z * rho / 2;
  
  // Translate into proportions of droplets (Poisson correction)
  p_N0 = 1 - exp(-lambda_N0);
  p_M0 = 1 - exp(-lambda_M0);

  p_M1 = 1 - exp(-lambda_M1);
  p_N1 = 1 - exp(-lambda_N1);
  
  p_X = 1 - exp(-lambda_X);
  p_Y = 1 - exp(-lambda_Y);
  
  lp_N[1] = log(0.5) + binomial_lpmf(K_N | n_K, p_N0);
  lp_N[2] = log(0.5) + binomial_lpmf(K_N | n_K, p_N1);
  
  lp_M[1] = log(0.5) + binomial_lpmf(K_M | n_K, p_M0);
  lp_M[2] = log(0.5) + binomial_lpmf(K_M | n_K, p_M1);

}

model {

  // rho ~ uniform(0, 1);  // Do not need to include as already implied
  M_K ~ gamma(1e-4, 1e-4);
  M_Z ~ gamma(1e-4, 1e-4);
  Z_X ~ binomial(n_Z, p_X);
  Z_Y ~ binomial(n_Z, p_Y);
  
  target += log_sum_exp(lp_N);
  target += log_sum_exp(lp_M);

}

generated quantities {
  simplex[2] pG;
  pG = softmax(lp_N + lp_M);
}

