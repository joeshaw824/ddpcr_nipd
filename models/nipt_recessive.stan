// Model from Tristan Snowsill. Variables renamed to be consistent with other stan models.

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
  real<lower=0> lambda_N_NN;
  real<lower=0> lambda_N_NM;
  real<lower=0> lambda_N_MM;
  real<lower=0> lambda_M_NN;
  real<lower=0> lambda_M_NM;
  real<lower=0> lambda_M_MM;
  real<lower=0> lambda_X;
  real<lower=0> lambda_Y;
  real<lower=0,upper=1> p_N_NN;
  real<lower=0,upper=1> p_N_NM;
  real<lower=0,upper=1> p_N_MM;
  real<lower=0,upper=1> p_M_NN;
  real<lower=0,upper=1> p_M_NM;
  real<lower=0,upper=1> p_M_MM;
  real<lower=0,upper=1> p_X;
  real<lower=0,upper=1> p_Y;
  
  vector[3] lp_N;
  vector[3] lp_M;
  
  // Concentrations of A and S if foetus homozygous AA
  lambda_N_NN = M_K * (1 + rho) / 2;
  lambda_M_NN = M_K * (1 - rho) / 2;
  
  // Concentrations of A and S if foetus heterozygous AS
  lambda_N_NM = M_K / 2;
  lambda_M_NM = M_K / 2;
  
  // Concentrations of A and S if foetus homozygous SS
  lambda_N_MM = M_K * (1 - rho) / 2;
  lambda_M_MM = M_K * (1 + rho) / 2;
  
  // Concentrations of X and Y
  lambda_X = M_Z * (2 - rho) / 2;
  lambda_Y = M_Z * rho / 2;
  
  // Translate into proportions of droplets (Poisson correction)
  p_N_NN = 1 - exp(-lambda_N_NN);
  p_M_NN = 1 - exp(-lambda_M_NN);
  
  p_N_NM = 1 - exp(-lambda_N_NM);
  p_M_NM = 1 - exp(-lambda_M_NM);
  
  p_N_MM = 1 - exp(-lambda_N_MM);
  p_M_MM = 1 - exp(-lambda_M_MM);

  
  p_X = 1 - exp(-lambda_X);
  p_Y = 1 - exp(-lambda_Y);
  
  lp_N[1] = log(0.25) + binomial_lpmf(K_N | n_K, p_N_NN);
  lp_N[2] = log(0.50) + binomial_lpmf(K_N | n_K, p_N_NM);
  lp_N[3] = log(0.25) + binomial_lpmf(K_N | n_K, p_N_MM);
  
  lp_M[1] = log(0.25) + binomial_lpmf(K_M | n_K, p_M_NN);
  lp_M[2] = log(0.50) + binomial_lpmf(K_M | n_K, p_M_NM);
  lp_M[3] = log(0.25) + binomial_lpmf(K_M | n_K, p_M_MM);

}

model {
  M_K ~ gamma(1e-4, 1e-4);
  M_Z ~ gamma(1e-4, 1e-4);
  Z_X ~ binomial(n_Z, p_X);
  Z_Y ~ binomial(n_Z, p_Y);
  
  target += log_sum_exp(lp_N);
  target += log_sum_exp(lp_M);
}

generated quantities {
  simplex[3] pG;
  pG = softmax(lp_N + lp_M);
}