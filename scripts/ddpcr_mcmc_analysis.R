#############################################################
## MonteCarlo Markov Chain ddPCR Analysis
## March 2021
## Joseph.Shaw@gosh.nhs.uk 
## This script is compiled from scripts and modesl supplied by 
## Tristan Snowsill (Exeter). Exeter use a threshold of
## 0.95 for classifying genotypes. This script can be used with 
## the ddPCR data prepared in the ddpcr_nipd script, named 
## ddpcr_data_mcmc
#############################################################

#############################################################
# Load Resources and Models
#############################################################

# Data structure
# n_K	= number of droplets tested for variant assay
# K_M	= number of droplets positive for variant (mutant) allele
# K_N	= number of droplets positive for normal (reference) allele
# n_Z	= number of droplets tested for fetal fraction assay
# Z_X	= number of droplets positive for maternal homozygous allele
# Z_Y	= number of droplets positive for paternal allele

# Load cmdstanr
library(cmdstanr)
library(tidyverse)

# Load functions
source("functions/ddpcr_nipd_functions.R")

# Compile the models

dominant_model <- cmdstan_model("models/nipt_dominant.stan")

x_linked_model <- cmdstan_model("models/nipt_x_linked.stan")

recessive_model <- cmdstan_model("models/nipt_recessive.stan")

# Intiialise the chains

initialise_chains_dominant <- function() list(rho = runif(1, 0.1, 0.5), 
                                              M_K = runif(1, 0.1, 0.5), 
                                              M_Z = runif(1, 0.1, 0.5))

initialise_chains_xlinked <- function() list(rho = rbeta(1, 4, 32),
                               M_K = abs(rnorm(1, sd = 0.05)),
                               M_Z = abs(rnorm(1, sd = 0.05)))


initialise_chains_recessive <- function() list(rho = runif(1, 0.1, 0.5),
                                     M_K = runif(1, 0.1, 0.5),
                                     M_Z = runif(1, 0.1, 0.5))

# Specify which assays are bi-allelic
biallelic_assays <- c("ADAR c.2997 G>T", "RNASEH2C c.205C>T",
                      "FGFR3 c.1138G>A", "ADA c.556G>A", "PMM2 c.691G>A", "HBB c.20A>T")

#############################################################
# Dominant Conditions Analysis
#############################################################

# Add on fits and extract the probabilities that the fetus is heterozygous (pG[1])
# and homozygous reference (pG[2]).
dominant_with_fits <- ddpcr_data_mcmc %>%
  filter(Inheritance == "autosomal" & !variant_assay %in% biallelic_assays) %>%
  nest(data = n_K:Z_Y) %>%
  mutate(
    data    = map(data, as.list),
    fit     = map(data, ~ dominant_model$sample(data = .,
                                                init = initialise_chains_dominant,
                                                step_size = 0.2,
                                                parallel_chains = parallel::detectCores())),
    results = map(fit,  ~ setNames(pivot_wider(.$summary(c("pG", "rho"), "mean"),
                                               names_from = "variable",
                                               values_from = "mean"),
                                   c("p_G1", "p_G2", "rho_est")))
  ) %>%
  unnest_wider(results)
  
dominant_mcmc_calls <- dominant_with_fits %>%
  select(-c(data, fit)) %>%
  dplyr::rename(fetal_fraction = rho_est) %>%
  mutate(mcmc_prediction = case_when(
    p_G1 > 0.95 ~"heterozygous",
    p_G2 > 0.95 ~"homozygous reference",
    p_G1 < 0.95 & p_G2 < 0.95 ~"no call"))

#############################################################
# X-linked Conditions Analysis
#############################################################

# Add on fits and extract the probabilities that the fetus is hemizygous reference (p_G0)
# and hemizygous variant (p_G1).

x_linked_with_fits <- ddpcr_data_mcmc %>%
  filter(Inheritance == "x_linked") %>%
  nest(data = n_K:Z_Y) %>%
  mutate(
    data    = map(data, as.list),
    fit     = map(data, ~ x_linked_model$sample(data = .,
                                                init = initialise_chains_xlinked,
                                                step_size = 0.2,
                                                parallel_chains = parallel::detectCores())),
    results = map(fit,  ~ setNames(pivot_wider(.$summary(c("pG", "rho"), "mean"),
                                               names_from = "variable",
                                               values_from = "mean"),
                                   c("p_G0", "p_G1", "rho_est")))
  ) %>%
  unnest_wider(results)

x_linked_mcmc_calls <- x_linked_with_fits %>%
  select(-c(data, fit)) %>%
  dplyr::rename(fetal_fraction = rho_est) %>%
  mutate(mcmc_prediction = case_when(
    p_G0 > 0.95 ~"hemizygous reference",
    p_G1 > 0.95 ~"hemizygous variant",
    p_G0 < 0.95 & p_G1 < 0.95 ~"no call"))


#############################################################
# Rare Recessive Conditions Analysis
#############################################################

# Extract the probabilities that the fetus is homozygous reference (pG[1]),
# heterozygous (pG[2]) and homozygous variant (pG[3]).

recessive_with_fits <- ddpcr_data_mcmc %>%
  filter(variant_assay %in% biallelic_assays & variant_assay != "HBB c.20A>T") %>%
  nest(data = n_K:Z_Y) %>%
  mutate(
    data    = map(data, as.list),
    fit     = map(data, ~ recessive_model$sample(data = .,
                                                init = initialise_chains_recessive,
                                                step_size = 0.2,
                                                parallel_chains = parallel::detectCores())),
    results = map(fit,  ~ setNames(pivot_wider(.$summary(c("pG", "rho"), "mean"),
                                               names_from = "variable",
                                               values_from = "mean"),
                                   c("p_G1", "p_G2", "p_G3", "rho_est")))
  ) %>%
  unnest_wider(results)

recessive_mcmc_calls <- recessive_with_fits %>%
  select(-c(data, fit)) %>%
  dplyr::rename(fetal_fraction = rho_est) %>%
  mutate(mcmc_prediction = case_when(
    p_G1 > 0.95 ~"homozygous reference",
    p_G2 > 0.95 ~"heterozygous",
    p_G3 > 0.95 ~"homozygous variant",
    p_G1 < 0.95 & p_G2 < 0.95 & p_G2 < 0.95 ~"no call"
  ))

#############################################################
# Sickle Cell Disease Analysis
#############################################################

# Extract the probabilities that the fetus is homozygous reference (pG[1]),
# heterozygous (pG[2]) and homozygous variant (pG[3]).

sickle_with_fits <- ddpcr_data_mcmc %>%
  filter(variant_assay == "HBB c.20A>T") %>%
  nest(data = n_K:Z_Y) %>%
  mutate(
    data    = map(data, as.list),
    fit     = map(data, ~ recessive_model$sample(data = .,
                                                 init = initialise_chains_recessive,
                                                 step_size = 0.2,
                                                 parallel_chains = parallel::detectCores())),
    results = map(fit,  ~ setNames(pivot_wider(.$summary(c("pG", "rho"), "mean"),
                                               names_from = "variable",
                                               values_from = "mean"),
                                   c("p_G1", "p_G2", "p_G3", "rho_est")))
  ) %>%
  unnest_wider(results)

sickle_mcmc_calls <- sickle_with_fits %>%
  select(-c(data, fit)) %>%
  dplyr::rename(fetal_fraction = rho_est) %>%
  mutate(mcmc_prediction = case_when(
    p_G1 > 0.95 ~"homozygous reference",
    p_G2 > 0.95 ~"heterozygous",
    p_G3 > 0.95 ~"homozygous variant",
    p_G1 < 0.95 & p_G2 < 0.95 & p_G2 < 0.95 ~"no call"
  ))

#############################################################
# Compare to SPRT
#############################################################

# Compare to the SPRT results

x_linked_comparison <- left_join(
  # First table
  x_linked_mcmc_calls,
  # Second table
  bespoke_cohort_blinded %>%
    mutate(r_number = as.character(r_number))%>%
    filter(Inheritance == "x_linked") %>%
    select(r_number,Likelihood_ratio, SPRT_prediction), 
  # Join by
  by = "r_number")

dominant_comparison <- left_join(
  # First table
  dominant_mcmc_calls,
  # Second table
  bespoke_cohort_blinded %>%
    mutate(r_number = as.character(r_number))%>%
    filter(Inheritance == "autosomal") %>%
    select(r_number, Likelihood_ratio, SPRT_prediction), 
  # Join by
  by = "r_number")

rare_recessive_comparison <- left_join(
  # First table
  recessive_mcmc_calls,
  # Second table
  bespoke_cohort_blinded %>%
    mutate(r_number = as.character(r_number))%>%
    filter(Inheritance == "autosomal") %>%
    select(r_number, Likelihood_ratio, SPRT_prediction), 
  # Join by
  by = "r_number")

scd_comparison <- left_join(
  # First table
  sickle_mcmc_calls,
  # Second table
  sickle_cell_blinded %>%
    mutate(r_number = as.character(r_number))%>%
    select(r_number, Likelihood_ratio, SPRT_prediction), 
  # Join by
  by = "r_number")

# Format the columns the same and bind together
mcmc_vs_sprt <- rbind(x_linked_comparison %>%
  mutate(p_G2 = "") %>%
  mutate(p_G3 = "") %>%
  select(r_number, Inheritance, variant_assay, fetal_fraction, 
         p_G0, p_G1, p_G2,  p_G3, Likelihood_ratio, mcmc_prediction,  
         SPRT_prediction),
  dominant_comparison %>%
  mutate(p_G0 = "") %>%
  mutate(p_G3 = "") %>%
  select(r_number, Inheritance, variant_assay, fetal_fraction, 
         p_G0, p_G1, p_G2,  p_G3, Likelihood_ratio, mcmc_prediction,  
         SPRT_prediction),
  rare_recessive_comparison %>%
  mutate(p_G0 = "") %>%
  select(r_number, Inheritance, variant_assay, fetal_fraction, 
         p_G0, p_G1, p_G2,  p_G3, Likelihood_ratio, mcmc_prediction,  
         SPRT_prediction),
  scd_comparison %>%
  mutate(p_G0 = "") %>%
  select(r_number, Inheritance, variant_assay, fetal_fraction, 
         p_G0, p_G1, p_G2,  p_G3, Likelihood_ratio, mcmc_prediction,  
         SPRT_prediction)) %>%
  mutate(concordant = ifelse(mcmc_prediction == SPRT_prediction,
                             "yes", "no"))

mcmc_vs_sprt_outcomes <- left_join(
  mcmc_vs_sprt,
  RAPID_biobank %>%
    mutate(r_number = as.character(r_number)) %>%
    select(r_number, confirmed_diagnosis, mutation_genetic_info_fetus),
  by = "r_number"
)

view(mcmc_vs_sprt_outcomes)

#############################################################
# 4 - Output csvs
#############################################################

current_time <- Sys.time()

write.csv(mcmc_vs_sprt_outcomes, 
          file = paste0("analysis_outputs/mcmc_vs_sprt_outcomes", 
                        format(current_time, "%Y%m%d_%H%M%S"), ".csv"),
          row.names = FALSE)
