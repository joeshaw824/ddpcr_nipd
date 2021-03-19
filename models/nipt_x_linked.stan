data {
  // X-linked assay
  int<lower=0> n_K;  // Number of accepted droplets
  int<lower=0> K_N;  // Positive droplets for normal allele
  int<lower=0> K_M;  // Positive droplets for mutant allele
  
  // Foetal fraction assay
  int<lower=0> n_Z;  // Number of accepted droplets
  int<lower=0> Z_X;  // Positive droplets for maternal allele
  int<lower=0> Z_Y;  // Positive droplets for paternal allele
}

parameters {
  real<lower=0,upper=1> rho;  // Foetal fraction (0-1)
  real<lower=0> M_K;          // DNA concentration in X-linked assay
  real<lower=0> M_Z;          // DNA concentration in foetal fraction assay
}

transformed parameters {
  // DNA concentration estimates
  real<lower=0> lambda_N_G0;  // Concentration of normal allele if foetus unaffected
  real<lower=0> lambda_N_G1;  // Concentration of normal allele if foetus affected
  real<lower=0> lambda_M_G0;  // Concentration of mutant allele if foetus unaffected
  real<lower=0> lambda_M_G1;  // Concentration of mutant allele if foetus affected
  real<lower=0> lambda_X;     // Concentration of maternal allele
  real<lower=0> lambda_Y;     // Concentration of paternal allele
  
  // Proportion of droplets positive estimates
  real<lower=0,upper=1> p_N_G0;
  real<lower=0,upper=1> p_N_G1;
  real<lower=0,upper=1> p_M_G0;
  real<lower=0,upper=1> p_M_G1;
  real<lower=0,upper=1> p_X;
  real<lower=0,upper=1> p_Y;
  
  // Marginalised log likelihoods
  // (see https://mc-stan.org/docs/2_24/stan-users-guide/latent-discrete-chapter.html)
  vector[2] lp_N;
  vector[2] lp_M;
  
  // Concentrations of N and M if foetus unaffected (XY)
  lambda_N_G0 = M_K / 2;
  lambda_M_G0 = M_K * (1 - rho) / 2;
  
  // Concentrations of N and M if foetus affected (xY)
  lambda_N_G1 = M_K * (1 - rho) / 2;
  lambda_M_G1 = M_K / 2;
  
  // Concentrations of X and Y
  lambda_X = M_Z * (2 - rho) / 2;
  lambda_Y = M_Z * rho / 2;
  
  // Translate into proportions of droplets (Poisson correction)
  p_N_G0 = 1 - exp(-lambda_N_G0);
  p_M_G0 = 1 - exp(-lambda_M_G0);
  
  p_N_G1 = 1 - exp(-lambda_N_G1);
  p_M_G1 = 1 - exp(-lambda_M_G1);

  
  p_X = 1 - exp(-lambda_X);
  p_Y = 1 - exp(-lambda_Y);
  
  // Marginalise log likelihoods
  lp_N[1] = log(0.5) + binomial_lpmf(K_N | n_K, p_N_G0);
  lp_N[2] = log(0.5) + binomial_lpmf(K_N | n_K, p_N_G1);
  
  lp_M[1] = log(0.5) + binomial_lpmf(K_M | n_K, p_M_G0);
  lp_M[2] = log(0.5) + binomial_lpmf(K_M | n_K, p_M_G1);
  
}

model {

  // PRIOR DISTRIBUTIONS
  // rho ~ uniform(0, 1);  // Do not need as flat prior is default
  M_K ~ gamma(1e-4, 1e-4);
  M_Z ~ gamma(1e-4, 1e-4);
  
  // LIKELIHOOD FOR FOETAL FRACTION ASSAY
  Z_X ~ binomial(n_Z, p_X);
  Z_Y ~ binomial(n_Z, p_Y);
  
  // LIKELIHOOD FOR X-LINKED ASSAY
  target += log_sum_exp(lp_N);
  target += log_sum_exp(lp_M);

}

generated quantities {
  // Calculate probability of each genotype (unaffected and affected)
  simplex[2] pG;
  pG = softmax(lp_N + lp_M);
}