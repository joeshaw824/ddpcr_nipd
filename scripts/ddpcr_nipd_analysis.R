#############################################################
## ddPCR for Non Invasive Prenatal Diagnosis (NIPD)
## January 2021
## Joseph.Shaw@gosh.nhs.uk
## Analysis script for prediction of fetal genotypes from 
## cfDNA testing using ddPCR for maternally-inherited
## variants. This includes modules for analysis of ddPCR data
## using the sequential probability ratio test (SPRT) (Lo et al,
## 2008; PMID:  17664418), and via MonteCarlo Markov Chain 
## (MCMC) analysis (Caswell et al, 2020; PMID:  32533152).
## The MCMC analysis section  is compiled from scripts and
## models supplied by Tristan Snowsill (Exeter).
#############################################################

#############################################################
# Load libraries and resources
#############################################################

## Load necessary packages
library(tidyverse)
# cmdstanr required for running stan models
library(cmdstanr)
# bayesplot and posterior required for MCMC diagnostic modelling
library(bayesplot)
library(posterior)

## Set working directory
setwd("W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/ddPCR_R_Analysis/ddpcr_nipd")

# Load resources
controls <- readr::read_csv("resources/controls.csv")
ddpcr_target_panel <- readr::read_csv("resources/ddpcr_target_panel.csv") 

# Source functions
source("functions/ddPCR_nipd_functions.R")
source("functions/RAPID_biobank.R")

#############################################################
# Read in cfDNA data and wrangle into shape
#############################################################

## Load in all the csv files exported from QuantaSoft.

dataPath <- "data/ddpcr_data/"

ddpcr_files <- list.files(dataPath)

#Empty data frame
ddpcr_data <- data.frame()

# Read and collate each worksheet csv
for (dataFile in ddpcr_files){
  tmp_dat <- read_csv(paste0(dataPath,dataFile), col_names = TRUE)
  # Add a worksheet identifier to make things easier later on.
  tmp_dat_ws <- tmp_dat %>%
    mutate(Worksheet = dataFile) %>%
    mutate(Worksheet_well = paste(Worksheet, Well, sep = "_"))
  ddpcr_data <-rbind(ddpcr_data, tmp_dat_ws)
  rm(tmp_dat_ws)
  rm(tmp_dat)
}

# Remove single wells and controls
ddpcr_data_merged_samples <- ddpcr_data %>%
  # Remove single wells
  drop_na("MergedWells") %>%
  # Remove controls
  filter(!(Sample %in% controls$Sample))

# Reshape the data frame to sum all values by Target
ddpcr_data_reshape <- ddpcr_data_merged_samples %>% 
  group_by(Sample, Target) %>% 
  summarise(Positives = sum(Positives),
            AcceptedDroplets = sum(AcceptedDroplets),
            .groups="drop")

ddpcr_with_target <- ddpcr_data_reshape %>% 
  left_join(ddpcr_target_panel %>%
              select(Target, Target_category), by = "Target")

# Use pivot_wider to get one row per sample
pivotted_ddpcr <- ddpcr_with_target %>% 
  pivot_wider(id_cols = Sample,
              names_from = Target_category,
              values_from = c(AcceptedDroplets, Positives))

# Add on assay and inheritance pattern to the table.
ref_table_var <- left_join(
  # First table
  ddpcr_with_target %>%
    filter(Target_category == "variant"),
  # Second table
  ddpcr_target_panel %>%
    select(Assay, Inheritance_chromosomal, Inheritance_pattern, Target) %>%
    # Rename the assay column to avoid clash with fetal fraction assay
    rename(variant_assay = Assay),
  # Join by
  by = "Target")
  
ref_table_ff <- left_join(
  # First table
  ddpcr_with_target %>%
    filter(Target_category == "ff_allele1"),
  # Second table
  ddpcr_target_panel %>%
    select(Assay, Target) %>%
    rename(ff_assay = Assay),
  # Join by
  by = "Target")

ddpcr_data_tbl <- pivotted_ddpcr %>%
  left_join(ref_table_var %>%
              select(Sample, Inheritance_chromosomal, Inheritance_pattern, variant_assay), by = "Sample",
            .groups="drop") %>% 
  left_join(ref_table_ff %>%
              select(Sample, ff_assay), by = "Sample",
            .groups="drop") %>%
  
  # Remove any samples which haven't had both assays performed.
  filter(!is.na(Positives_variant) & !is.na(Positives_ff_allele1)) %>%
  
  # Remove duplicate columns and rename to be compatible with functions
  select(-c(AcceptedDroplets_ff_allele2, AcceptedDroplets_reference)) %>%
  rename(AcceptedDroplets_FetalFrac = AcceptedDroplets_ff_allele1) %>%
  rename(AcceptedDroplets_Variant_assay = AcceptedDroplets_variant) %>%
  
  # Determine the maternal and paternally-inherited alleles for the fetal fraction calculation
  mutate(Positives_maternal = pmax(Positives_ff_allele1, Positives_ff_allele2)) %>%
  mutate(Positives_paternal = pmin(Positives_ff_allele1, Positives_ff_allele2))

#############################################################
# Wrangle gDNA data into shape
#############################################################

# Get single well controls only without NTC
ddpcr_controls <- ddpcr_data %>%
  filter(Sample %in% controls$Sample & Sample != "NTC" & !is.na(CopiesPer20uLWell))

ddpcr_controls_with_target <- ddpcr_controls %>% 
  left_join(ddpcr_target_panel %>%
              select(Target, Target_category), by = "Target")

# No need to pivot to get one row per sample because we want the control
# data as individual wells.
# Use pivot_wider to get one row per well for each well tested with 
# a variant assay
pivotted_controls_var <- ddpcr_controls_with_target %>%
  filter(Target_category %in% c("variant", "reference")) %>%
  pivot_wider(id_cols = c(Worksheet_well, Sample),
              names_from = Target_category,
              values_from = c(AcceptedDroplets, Positives)) %>%
  select(-AcceptedDroplets_reference)

# Use pivot_wider to get one row per well for each well tested with 
# a fetal fraction assay
pivotted_controls_ff <- ddpcr_controls_with_target %>%
  filter(Target_category %in% c("ff_allele2", "ff_allele1")) %>%
  pivot_wider(id_cols = c(Worksheet_well, Sample),
              names_from = Target_category,
              values_from = c(AcceptedDroplets, Positives)) %>%
  select(-AcceptedDroplets_ff_allele2)

# Add on assay and inheritance pattern to the variant control table.
control_table_var <- left_join(
  # First table
  ddpcr_controls_with_target %>%
    filter(Target_category == "variant") %>%
    select(Worksheet_well, Sample, Target),
  # Second table
  ddpcr_target_panel %>%
    select(Assay, Inheritance_chromosomal, Inheritance_pattern, Target) %>%
    # Rename the assay column to avoid clash with fetal fraction assay
    rename(variant_assay = Assay),
  # Join by
  by = "Target")

# Add on assay to the fetal fraction control table.
control_table_ff <- left_join(
  # First table
  ddpcr_controls_with_target %>%
    filter(Target_category == "ff_allele1") %>%
    select(Worksheet_well, Sample, Target),
  # Second table
  ddpcr_target_panel %>%
    select(Assay, Target) %>%
    rename(ff_assay = Assay),
  # Join by
  by = "Target")

# Get the control information for the control wells tested with a variant assay
ddpcr_control_tbl_var <- pivotted_controls_var %>%
  left_join(control_table_var %>%
              select(Worksheet_well, Inheritance_chromosomal, Inheritance_pattern, variant_assay), by = "Worksheet_well",
            .groups="drop") %>%
  # Remove duplicate columns and rename to be compatible with functions
  rename(AcceptedDroplets_Variant_assay = AcceptedDroplets_variant) %>%
  mutate(Molecules_variant = Poisson_correct(AcceptedDroplets_Variant_assay,Positives_variant)) %>%
  mutate(Molecules_reference = Poisson_correct(AcceptedDroplets_Variant_assay,Positives_reference))

# Get the control information for the control wells tested with a fetal fraction assay
ddpcr_control_tbl_ff <- pivotted_controls_ff %>%
  left_join(control_table_ff %>%
              select(Worksheet_well, ff_assay), by = "Worksheet_well",
            .groups="drop") %>%
  # Remove duplicate columns and rename to be compatible with functions
  rename(AcceptedDroplets_FetalFrac = AcceptedDroplets_ff_allele1) %>%
  mutate(Positives_maternal = pmax(Positives_ff_allele1, Positives_ff_allele2)) %>%
  mutate(Positives_paternal = pmin(Positives_ff_allele1, Positives_ff_allele2)) %>%
  mutate(Molecules_maternal = Poisson_correct(AcceptedDroplets_FetalFrac,Positives_maternal)) %>%
  mutate(Molecules_paternal = Poisson_correct(AcceptedDroplets_FetalFrac,Positives_paternal))

#############################################################
# cfDNA SPRT analysis
#############################################################

# Set likelihood ratio threshold
LR_threshold <- 250

ddpcr_data_analysed <- ddpcr_data_tbl %>%
  # Rename sample to r_number to allow merge with RAPID Biobank data in next step.
  # Have to specify dplyr for rename.
  dplyr::rename(r_number = Sample) %>%
  
  # Perform Poisson correction to determine the total number of molecules detected for each species.
  
  mutate(Molecules_variant = Poisson_correct(AcceptedDroplets_Variant_assay,Positives_variant)) %>%
  mutate(Molecules_reference = Poisson_correct(AcceptedDroplets_Variant_assay,Positives_reference)) %>%
  mutate(Molecules_maternal = Poisson_correct(AcceptedDroplets_FetalFrac,Positives_maternal)) %>%
  mutate(Molecules_paternal = Poisson_correct(AcceptedDroplets_FetalFrac,Positives_paternal)) %>%
  
  # Find the difference between the numbers of molecules for each allele
  
  mutate(Molecules_difference = abs(Molecules_reference - Molecules_variant)) %>%
  
  # Find the total number of molecules for each assay
  
  mutate(Molecules_variant_assay = Molecules_variant + Molecules_reference) %>%
  
  mutate(Molecules_ff_assay = Molecules_maternal + Molecules_paternal) %>%
  
  mutate(Molecules_total = Molecules_variant_assay + Molecules_ff_assay) %>%
  
  # Calculate the fetal fraction
  mutate(Fetal_fraction = calc_ff(Molecules_maternal, Molecules_paternal)) %>%
  
  # Convert to a percentage in case I want to print as a table later on (percentages are easier to look at than decimals)
  mutate(Fetal_fraction_percent = Fetal_fraction*100) %>%
  
  # Calculate the fractional abundance of each allele
  
  mutate(Reference_fraction = Molecules_reference / Molecules_variant_assay) %>%
  mutate(Variant_fraction = Molecules_variant / Molecules_variant_assay) %>%
  mutate(Variant_fraction_percent  = Variant_fraction * 100) %>%
  
  # Calculate the fraction of the overrepresented allele
  
  mutate(Over_represented_fraction = pmax(Reference_fraction, Variant_fraction)) %>%
  mutate(Over_represented_fraction_percent = (pmax(Reference_fraction, Variant_fraction))*100) %>%
  
  # Calculate copies per droplet (cpd) for each allele
  mutate(Cpd_variant = Molecules_variant / AcceptedDroplets_Variant_assay) %>%
  mutate(Cpd_reference = Molecules_reference / AcceptedDroplets_Variant_assay) %>%
  mutate(Cpd_paternal = Molecules_paternal / AcceptedDroplets_FetalFrac) %>%
  mutate(Cpd_maternal = Molecules_maternal / AcceptedDroplets_FetalFrac) %>%
  
  # Calculate the 95% confidence intervals of the numbers of molecules for each allele
  mutate(Molecules_variant_max = Poisson_max(Cpd_variant, AcceptedDroplets_Variant_assay)) %>%
  mutate(Molecules_variant_min = Poisson_min(Cpd_variant, AcceptedDroplets_Variant_assay)) %>%
  mutate(Molecules_reference_max = Poisson_max(Cpd_reference, AcceptedDroplets_Variant_assay)) %>%
  mutate(Molecules_reference_min = Poisson_min(Cpd_reference, AcceptedDroplets_Variant_assay)) %>%
  mutate(Molecules_paternal_max = Poisson_max(Cpd_paternal, AcceptedDroplets_FetalFrac)) %>%
  mutate(Molecules_paternal_min = Poisson_min(Cpd_paternal, AcceptedDroplets_FetalFrac)) %>%
  mutate(Molecules_maternal_max = Poisson_max(Cpd_maternal, AcceptedDroplets_FetalFrac)) %>%
  mutate(Molecules_maternal_min = Poisson_min(Cpd_maternal, AcceptedDroplets_FetalFrac)) %>%
  
  # Calculate the 95% confidence intervals of the fractional abundances
  
  mutate(Variant_fraction_max_percent = (Poisson_fraction_max(Molecules_variant_max, Molecules_reference))*100) %>%
  mutate(Variant_fraction_min_percent = (Poisson_fraction_min(Molecules_variant_min, Molecules_reference))*100) %>%
  mutate(Reference_fraction_max_percent = (Poisson_fraction_max(Molecules_reference_max, Molecules_variant))*100) %>%
  mutate(Reference_fraction_min_percent = (Poisson_fraction_min(Molecules_reference_min, Molecules_variant))*100) %>%
  mutate(Over_represented_fraction_max_percent = pmax(Variant_fraction_max_percent, Reference_fraction_max_percent)) %>%
  mutate(Over_represented_fraction_min_percent = pmax(Variant_fraction_min_percent, Reference_fraction_min_percent)) %>%
  
  # When calculating for the fetal fraction, the paternal fraction must be multiplied by 2.
  mutate(Fetal_fraction_max_percent = 200* (Poisson_fraction_max(Molecules_paternal_max, Molecules_maternal))) %>%
  mutate(Fetal_fraction_min_percent = 200* (Poisson_fraction_min(Molecules_paternal_min, Molecules_maternal))) %>%
  
  # Perform SPRT and return the likelihood ratio
  mutate(Likelihood_ratio = case_when(
    Inheritance_chromosomal == "x_linked" ~ calc_LR_X_linked(Fetal_fraction, Over_represented_fraction, Molecules_variant_assay),
    Inheritance_chromosomal == "autosomal" ~ calc_LR_autosomal(Fetal_fraction, Over_represented_fraction, Molecules_variant_assay))) %>%
  
  # Classify based on likelihood ratio threshold supplied
  mutate(SPRT_prediction = case_when(
    Inheritance_chromosomal == "autosomal" & Likelihood_ratio > LR_threshold & Over_represented_fraction == Reference_fraction ~ "homozygous reference",
    Inheritance_chromosomal == "autosomal" & Likelihood_ratio > LR_threshold & Over_represented_fraction == Variant_fraction ~ "homozygous variant",
    Inheritance_chromosomal == "autosomal" & Likelihood_ratio < (1/LR_threshold) ~ "heterozygous",
    Inheritance_chromosomal == "autosomal" & Likelihood_ratio < LR_threshold & Likelihood_ratio > (1/LR_threshold) ~ "no call",
    Inheritance_chromosomal == "x_linked" & Likelihood_ratio > LR_threshold & Over_represented_fraction == Reference_fraction ~ "hemizygous reference",
    Inheritance_chromosomal == "x_linked" & Likelihood_ratio > LR_threshold & Over_represented_fraction == Variant_fraction ~ "hemizygous variant",
    Inheritance_chromosomal == "x_linked" & Likelihood_ratio < LR_threshold ~ "no call")) %>%
  
  mutate(Call = ifelse(SPRT_prediction == "no call", "no call", "call"))


#############################################################
# Compare SPRT results to RAPID Biobank results
#############################################################

ddpcr_nipd_unblinded <- left_join(
  ddpcr_data_analysed,
  RAPID_biobank %>%
    # Change r_number to a character to match ddpcr_data_analysed
    mutate(r_number = as.character(r_number)) %>%
    filter(r_number %in% ddpcr_data_analysed$r_number) %>%
    select(r_number, study_id, gestation_weeks, gestation_days, Gestation_total_weeks, 
           date_of_blood_sample, vacutainer, mutation_genetic_info_fetus, Partner_sample_available),
  by = "r_number")

#############################################################
# Sickle cell disease cohort
#############################################################

# Ammend the fetal genotype prediction for samples which failed QC steps
samples_failing_qc <- c(13262, 20763, 20810)

sickle_cell_unblinded <- ddpcr_nipd_unblinded %>%
  filter(variant_assay == "HBB c.20A>T") %>%
  # Add in a column for the overall prediction based on the data.
  mutate(overall_prediction = ifelse(r_number %in% samples_failing_qc, "no call", SPRT_prediction))

# Ammend the dataframe for the sample from a twin pregnancy (20915).

sickle_cell_unblinded$SPRT_prediction[sickle_cell_unblinded$r_number == 20915] <- "twin pregnancy"
sickle_cell_unblinded$overall_prediction[sickle_cell_unblinded$r_number == 20915] <- "twin pregnancy"

nrow(sickle_cell_unblinded %>%
       filter(overall_prediction != "no call"))

#############################################################
# Perform MCMC analysis
#############################################################
#########################
# Prepare data
#########################

# n_K	= number of droplets tested for variant assay
# K_M	= number of droplets positive for variant (mutant) allele
# K_N	= number of droplets positive for normal (reference) allele
# n_Z	= number of droplets tested for fetal fraction assay
# Z_X	= number of droplets positive for maternal homozygous allele
# Z_Y	= number of droplets positive for paternal allele

ddpcr_data_mcmc <- ddpcr_data_tbl %>%
  dplyr::rename(r_number = Sample,
                n_K = AcceptedDroplets_Variant_assay,
                K_M = Positives_variant,
                K_N = Positives_reference,
                n_Z = AcceptedDroplets_FetalFrac,
                Z_X = Positives_maternal,
                Z_Y = Positives_paternal) %>%
  arrange(Inheritance_chromosomal, Inheritance_pattern, variant_assay) %>%
  select(r_number, Inheritance_chromosomal, Inheritance_pattern, variant_assay, n_K, n_Z, 
         K_N, K_M, n_Z, Z_X, Z_Y)



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

#########################
# Dominant conditions
#########################

# Add on fits and extract the probabilities that the fetus is heterozygous (pG[1])
# and homozygous reference (pG[2]).
dominant_with_fits <- ddpcr_data_mcmc %>%
  filter(Inheritance_chromosomal == "autosomal" & !variant_assay %in% biallelic_assays) %>%
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

#########################
# X-linked Conditions
#########################

# Add on fits and extract the probabilities that the fetus is hemizygous reference (p_G0)
# and hemizygous variant (p_G1).

x_linked_with_fits <- ddpcr_data_mcmc %>%
  filter(Inheritance_chromosomal == "x_linked") %>%
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

#########################
# Rare recessive conditions
#########################

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

#########################
# Sickle cell disease
#########################

# Extract the probabilities that the fetus is homozygous reference (pG[1]),
# heterozygous (pG[2]) and homozygous variant (pG[3]).

sickle_with_fits <- ddpcr_data_mcmc %>%
  filter(variant_assay == "HBB c.20A>T") %>%
  filter(r_number %in% c(30113, 20611)) %>%
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
# Compare MCMC to SPRT
#############################################################

# Compare to the SPRT results

x_linked_comparison <- left_join(
  # First table
  x_linked_mcmc_calls,
  # Second table
  ddpcr_nipd_unblinded %>%
    mutate(r_number = as.character(r_number))%>%
    filter(Inheritance_chromosomal == "x_linked") %>%
    select(r_number,Likelihood_ratio, SPRT_prediction), 
  # Join by
  by = "r_number")

dominant_comparison <- left_join(
  # First table
  dominant_mcmc_calls,
  # Second table
  ddpcr_nipd_unblinded %>%
    mutate(r_number = as.character(r_number))%>%
    filter(Inheritance_chromosomal == "autosomal") %>%
    select(r_number, Likelihood_ratio, SPRT_prediction), 
  # Join by
  by = "r_number")

rare_recessive_comparison <- left_join(
  # First table
  recessive_mcmc_calls,
  # Second table
  ddpcr_nipd_unblinded %>%
    mutate(r_number = as.character(r_number))%>%
    filter(Inheritance_chromosomal == "autosomal") %>%
    select(r_number, Likelihood_ratio, SPRT_prediction), 
  # Join by
  by = "r_number")

scd_comparison <- left_join(
  # First table
  sickle_mcmc_calls,
  # Second table
  ddpcr_nipd_unblinded %>%
    filter(variant_assay == "HBB c.20A>T") %>%
    mutate(r_number = as.character(r_number))%>%
    select(r_number, Likelihood_ratio, SPRT_prediction), 
  # Join by
  by = "r_number")

# Format the columns the same and bind together
mcmc_vs_sprt <- rbind(x_linked_comparison %>%
                        mutate(p_G2 = "") %>%
                        mutate(p_G3 = "") %>%
                        select(r_number, Inheritance_chromosomal, variant_assay, fetal_fraction, 
                               p_G0, p_G1, p_G2,  p_G3, Likelihood_ratio, mcmc_prediction,  
                               SPRT_prediction),
                      dominant_comparison %>%
                        mutate(p_G0 = "") %>%
                        mutate(p_G3 = "") %>%
                        select(r_number, Inheritance_chromosomal, variant_assay, fetal_fraction, 
                               p_G0, p_G1, p_G2,  p_G3, Likelihood_ratio, mcmc_prediction,  
                               SPRT_prediction),
                      rare_recessive_comparison %>%
                        mutate(p_G0 = "") %>%
                        select(r_number, Inheritance_chromosomal, variant_assay, fetal_fraction, 
                               p_G0, p_G1, p_G2,  p_G3, Likelihood_ratio, mcmc_prediction,  
                               SPRT_prediction),
                      scd_comparison %>%
                        mutate(p_G0 = "") %>%
                        select(r_number, Inheritance_chromosomal, variant_assay, fetal_fraction, 
                               p_G0, p_G1, p_G2,  p_G3, Likelihood_ratio, mcmc_prediction,  
                               SPRT_prediction)) %>%
  mutate(concordant = ifelse(mcmc_prediction == SPRT_prediction,
                             "yes", "no"))

#############################################################
# ROC curve analysis
############################################################

# This part is for plotting ROC curves for the sickle cell 
# disease data.

mcmc_vs_sprt_scd <- left_join(
  scd_comparison,
  RAPID_biobank %>%
    mutate(r_number = as.character(r_number)) %>%
    select(r_number, mutation_genetic_info_fetus),
  by = "r_number") %>%
  # Remove twin pregnancy, sample without outcome and HbAC case
  filter(!r_number %in% c(30230, 20915, 17004)) %>%
  # Convert invasive results to binary outcomes
  mutate(unbalanced = case_when(
    mutation_genetic_info_fetus %in% c("HbSS", "HbAA") ~"TRUE",
    mutation_genetic_info_fetus == "HbAS" ~"FALSE")) %>%
  # Convert "unbalanced" column to Boolean vector
  mutate(unbalanced = as.logical(unbalanced)) %>%
  # Convert the MCMC calls to a binary outcome
  mutate(mcmc_hom_call = pmax(p_G1, p_G3))

sprt_roc <- mcmc_vs_sprt_scd %>%
  arrange(desc(Likelihood_ratio)) %>%
  mutate(true_positive_rate = cumsum(unbalanced)/sum(unbalanced)) %>%
  mutate(false_positive_rate = cumsum(!unbalanced)/sum(!unbalanced)) %>%
  mutate(analysis_type = "sprt")%>%
  select(analysis_type, true_positive_rate, false_positive_rate)

mcmc_roc <- mcmc_vs_sprt_scd %>%
  arrange(desc(mcmc_hom_call)) %>%
  mutate(true_positive_rate = cumsum(unbalanced)/sum(unbalanced)) %>%
  mutate(false_positive_rate = cumsum(!unbalanced)/sum(!unbalanced)) %>%
  mutate(analysis_type = "mcmc") %>%
  select(analysis_type, true_positive_rate, false_positive_rate)

total_roc <- rbind(sprt_roc, mcmc_roc)

ggplot(total_roc, aes(x = false_positive_rate, y = true_positive_rate))+
  geom_line(size = 2, alpha = 0.2)+
  facet_wrap(~analysis_type)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(linetype = "dashed")+
  ylim(0, 1)+
  xlim(0,1)+
  labs(x = "False positive rate", y = "True positive rate", title = "Sickle cell disease ddPCR cohort ROC curve")

############################################################
# MCMC diagnostic plotting
#############################################################

# Use this web address: https://mc-stan.org/cmdstanr/reference/cmdstanr-package.html
# Looks like the URL that Tristan sent uses stanr, whereas we/he use cmdstanr.
# You can still use bayesplot with cmdstanr, but you'll need to follow the 
# cmdstanr specific tutorial in the link above.

# This section is not finished yet.


## Use one sample only.
ddpcr_data_mcmc %>%
  filter(r_number == 19261)

# Data for one sample
data_19261 <- list(n_K = 254853, n_Z = 240083, K_N = 2159, K_M = 2282, Z_X = 4352, Z_Y = 457)
data_19711 <- list(n_K = 361743, n_Z = 372037, K_N =   1129, K_M = 1009, Z_X = 2494, Z_Y = 116)
data_12945 <- list(n_K = 367794, n_Z = 336160, K_N =   1728 , K_M = 1798, Z_X = 8123, Z_Y = 281)

data_all <- dominant_with_fits$data

# Run the single sample according to the code upstream
fit_19261 <- dominant_model$sample(
  data = data_19261,
  init = initialise_chains_dominant,
  step_size = 0.2,
  parallel_chains = parallel::detectCores())

fit_19711 <- dominant_model$sample(
  data = data_19711,
  init = initialise_chains_dominant,
  step_size = 0.2,
  parallel_chains = parallel::detectCores())

fit_12945 <- x_linked_model$sample(
  data = data_12945,
  init = initialise_chains_xlinked,
  step_size = 0.2,
  parallel_chains = parallel::detectCores())

# Run all dominant samples
fit_dominant_all <- dominant_model$sample(
  data = data_all,
  init = initialise_chains_dominant,
  step_size = 0.2,
  parallel_chains = parallel::detectCores())

# In order to plot some graphs it looks like you have to convert to a 
# stanfit object
stanfit_12945 <- rstan::read_stan_csv(fit_12945$output_files())
stanfit_19711 <- rstan::read_stan_csv(fit_19711$output_files())

lp_stanfit_12945 <- log_posterior(stanfit_12945)

posterior_12945 <- as.array(stanfit_12945)
posterior_19711 <- as.array(stanfit_19711)

np_stanfit_12945 <- nuts_params(stanfit_12945)
np_stanfit_19711 <- nuts_params(stanfit_19711)

mcmc_parcoord(stanfit_12945, np = np_stanfit_12945)

plot_1 <- mcmc_pairs(posterior_19711, np = np_stanfit_19711, pars = c("rho","M_K","M_Z"))
plot_2 <- mcmc_pairs(posterior_19261, np = np_stanfit_19261, pars = c("rho","M_K","M_Z"))

mcmc_trace(posterior_12945, pars = "rho", np = stanfit_12945)

mcmc_trace(posterior_cp, pars = "tau", np = np_cp)

x_linked_test <- ddpcr_data_mcmc %>%
  filter(Inheritance_chromosomal == "x_linked") %>%
  nest(data = n_K:Z_Y) %>%
  mutate(data = map(data, as.list))

x_linked_fit_test <- map(x_linked_test$data, ~ x_linked_model$sample(data = .,
                                                init = initialise_chains_xlinked,
                                                step_size = 0.2,
                                                parallel_chains = parallel::detectCores()))

# In order to plot some graphs it looks like you have to convert to a 
# stanfit object
stanfit_x_linked <- rstan::read_stan_csv(x_linked_fit_test$output_files())

#############################################################
# Sickle cell gDNA analysis
#############################################################

# Worksheets containing 20ng 19RG-220G0190 per replicate
worksheets_20ng <- c("20-4046.csv", "20-3824.csv", "20-3756_FAM_HEX.csv", "20-4321.csv", 
                     "20-4317.csv", "20-4227.csv")

het_DNA_ids <- read.csv("resources/het_DNA_ids.csv")

mat_het_gDNA <- ddpcr_data %>%
  filter(Sample == "19RG-220G0190") %>%
  filter(is.na(MergedWells)) %>%
  filter(Worksheet %in% worksheets_20ng) %>%
  filter(Target == "HbS") %>%
  select(Worksheet_well, Positives, FractionalAbundance, 
         PoissonFractionalAbundanceMax, PoissonFractionalAbundanceMin)

pat_cfDNA <- ddpcr_data %>%
  filter(Sample %in% c(30139, 30130)) %>%
  filter(!is.na(MergedWells)) %>%
  filter(Target == "HbS") %>%
  select(Worksheet_well, Positives, FractionalAbundance, 
         PoissonFractionalAbundanceMax, PoissonFractionalAbundanceMin)

het_DNA <- rbind(mat_het_gDNA, pat_cfDNA)

het_gDNA_arranged <- left_join(het_DNA, het_DNA_ids, by = "Worksheet_well") %>%
  select(Control_id, FractionalAbundance, PoissonFractionalAbundanceMax, PoissonFractionalAbundanceMin)

colnames(het_gDNA_arranged) <- c("Identifier", "Variant_fraction_percent", "Variant_fraction_max_percent", 
                                 "Variant_fraction_min_percent")

mat_cfDNA_het <- sickle_cell_unblinded %>%
  filter(invasive_result == "HbAS") %>%
  arrange(overall_prediction) %>%
  select(Identifier, Variant_fraction_percent, Variant_fraction_max_percent, Variant_fraction_min_percent)

all_het_DNA <- rbind(het_gDNA_arranged, mat_cfDNA_het)

het_DNA_plot <- ggplot(het_gDNA_arranged, aes(x = Identifier, y = Variant_fraction_percent))+
  geom_point(size = 4)+
  geom_pointrange(aes(ymin = Variant_fraction_min_percent, ymax = Variant_fraction_max_percent))+
  geom_hline(aes(yintercept=50), linetype="dashed", size = 1)+
  geom_hline(aes(yintercept=51.2), linetype="dashed", colour = "grey", size = 1)+
  geom_hline(aes(yintercept=48.9), linetype="dashed", colour = "grey", size = 1)+
  ylim(40, 60)+
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        plot.title = element_text(size=20),
        legend.position = c(0.85, 0.85), legend.background = element_rect(fill="white"), 
        legend.title = element_text(size= 18), legend.text = element_text(size= 15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "", y = "Variant fractional abundance (%)")+
  theme(axis.text.x = element_text(angle = 90))

#############################################################
# Sickle cell disease limit of detection study
#############################################################

LOD_data <- read_csv("data/20-1557_LOD.csv", col_names = TRUE)

LOD_data_longer <- LOD_data %>%
  mutate(unique_identifier = paste(Input_molecules, Sample)) %>%
  pivot_wider(id_cols = c(unique_identifier, Sample, Input_molecules, Mass_molecules, fetal_fraction),
              names_from = Target,
              values_from = c(AcceptedDroplets, Positives, FractionalAbundance, PoissonFractionalAbundanceMax, 
                              PoissonFractionalAbundanceMin)) %>%
  select(-c(AcceptedDroplets_HbS, FractionalAbundance_HbA, PoissonFractionalAbundanceMax_HbA, PoissonFractionalAbundanceMin_HbA)) %>%
  rename(AcceptedDroplets_Variant_assay = AcceptedDroplets_HbA) %>%
  mutate(HbS_molecules = Poisson_correct(AcceptedDroplets_Variant_assay, Positives_HbS)) %>%
  mutate(HbA_molecules = Poisson_correct(AcceptedDroplets_Variant_assay, Positives_HbA)) %>%
  mutate(total_DNA_molecules = HbS_molecules + HbA_molecules) %>%
  mutate(Reference_fraction = HbA_molecules / total_DNA_molecules) %>%
  mutate(Variant_fraction = HbS_molecules / total_DNA_molecules) %>%
  mutate(Over_represented_fraction = pmax(Reference_fraction, Variant_fraction)) %>%
  mutate(Likelihood_ratio = calc_LR_autosomal(fetal_fraction, Over_represented_fraction, total_DNA_molecules)) %>%
  mutate(SPRT_prediction = case_when(
    Likelihood_ratio > 250 & Over_represented_fraction == Reference_fraction ~ "homozygous reference",
    Likelihood_ratio > 250 & Over_represented_fraction == Variant_fraction ~ "homozygous variant",
    Likelihood_ratio < (1/250) ~ "heterozygous",
    Likelihood_ratio < 250 & Likelihood_ratio > (1/250) ~ "no call"))

LOD_data_even_longer <- LOD_data_longer %>%
  select(Sample, Mass_molecules, Input_molecules, HbS_molecules, HbA_molecules, AcceptedDroplets_Variant_assay) %>%
  pivot_longer(cols = c(HbS_molecules, HbA_molecules), names_to = "Target", values_to = "molecules") %>%
  # Order the factors
  mutate(Sample = factor(Sample, levels = c("AA 12%", "AA 10%", "AA 8%", "AA 6%", "AA 4%","AA 2%",
                                            "0%", "SS 2%", "SS 4%", "SS 6%", 
                                            "SS 8%", "SS 10%", "SS 12%"))) %>%
  mutate(Mass_molecules = factor(Mass_molecules, levels = c("9.9ng (~3,000 molecules)", "19.8ng (~6,000 molecules)",
                                                            "29.7ng (~9,000 molecules)", "39.6ng (~12,000 molecules)"))) %>%
  
  # Add in the 95% Poisson confidence intervals
  
  # Calculate copies per droplet (cpd) for each allele
  mutate(Cpd = molecules / AcceptedDroplets_Variant_assay) %>%
  
  # Calculate the 95% confidence intervals of the numbers of molecules for each allele
  mutate(Molecules_max = Poisson_max(Cpd, AcceptedDroplets_Variant_assay)) %>%
  mutate(Molecules_min = Poisson_min(Cpd, AcceptedDroplets_Variant_assay)) %>%
  mutate(Percentage = case_when(
    Sample %in% c("SS 2%", "AA 2%") ~"2%",
    Sample %in% c("SS 4%", "AA 4%") ~"4%",
    Sample %in% c("SS 6%", "AA 6%") ~"6%",
    Sample %in% c("SS 8%", "AA 8%") ~"8%",
    Sample %in% c("SS 10%", "AA 10%") ~"10%",
    Sample %in% c("SS 12%", "AA 12%") ~"12%",
    Sample %in% c("0") ~"0%"))

# Plot the results for presentation
ggplot(LOD_data_even_longer %>%
         filter(Input_molecules %in% c(12000)), aes(x = Sample, y = molecules, fill = Target))+
  geom_col(width = 0.5, position = position_dodge(width =0.5), colour = "black", alpha = 0.6)+
  geom_errorbar(aes(ymin = Molecules_min, ymax = Molecules_max), width = .2, position=position_dodge(width=0.5))+
  scale_fill_manual(values=c("#3366FF", "#FF0000"), labels= c("Reference", "Variant"))+
  theme_bw()+
  theme(axis.text=element_text(size=10), axis.title = element_text(size=14), 
        plot.title = element_text(size=20), legend.position = "bottom", legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "", y = "Molecules detected by ddPCR")+
  ylim(0, 8000)+
  scale_x_discrete(labels=c("12%","10%","8%","6%","4%","2%","0%",
                            "2%","4%","6%","8%","10%","12%"))

# A different plot - not as nice to look at.
ggplot(LOD_data_longer, aes(x = total_DNA_molecules, y = FractionalAbundance_HbS, colour = Sample))+
  geom_point(size = 4)+
  geom_pointrange(aes(ymin = PoissonFractionalAbundanceMin_HbS, ymax = PoissonFractionalAbundanceMax_HbS))+
  geom_hline(yintercept=50, linetype="dashed", size = 1)


#############################################################
# cfDNA genomic equivalents per ml plasma (GE/ml)
#############################################################

# cfDNA concentration in plasma depends on the extraction volume,
# number of replicates and input volume, all of which was changed
# for different samples.

# 12585: ZFXY tested on 20-4420, Test9 and 19-3581, F8 tested on
# 20-4515 and Test7
# 16319: ZFXY tested on 19-3581 and Test9
# 16881: ZFXY tested on 20-4420, Test9 and 19-3581; F8 tested on
# 20-4515 and Test8
# 18164: ZFXY tested on Test9 and 20-3156
# 18385: ZFXY tested on 20-4420 and Test27; ATP7A tested on 
# 20-4468 and Test26
# 18891: ZFXY tested on 20-4126 and Test27; FOXP3 tested on
# 20-4126 and Test22

# Specify the sample_worksheet combinations to exclude
rows_to_exclude <- c("12585_Test9.csv", "12585_Test7.csv", "12585_19-3581.csv",
                     "16319_Test9.csv", "16881_Test9.csv", "16881_19-3581.csv",
                     "16881_Test8.csv", "18385_Test26.csv",
                     "18164_Test9.csv", "18385_Test27.csv", "18891_Test27.csv",
                     "18891_Test22.csv")

# This tibble can then replace the ddpcr_data_reshape tibble from earlier.
# This should get all the bespoke ddPCR cohort with one worksheet each for the
# variant and fetal fraction assays, so that the GE/ml can be calculated.

for_GE_calc <- ddpcr_data_merged_samples %>%
  filter(Sample %in% bespoke_cohort_analysed$r_number) %>%
  # Change concentration to a numeric
  mutate(Concentration = as.numeric(Concentration))%>%
  mutate(PoissonConfMax = as.numeric(PoissonConfMax))%>%
  mutate(PoissonCNVMin = as.numeric(PoissonCNVMin))%>%
  mutate(replicates = str_count(MergedWells, ",")+1) %>%
  mutate(Sample_worksheet = paste(Sample, Worksheet, sep = "_"))%>%
  # Filter out the repeat samples
  filter(!Sample_worksheet %in% rows_to_exclude) %>%
  arrange(Sample)

# Which samples/worksheets had odd elution or extraction volumes?
# 13625 on Test23 (variant assay) would have had a 2ml extraction
# 16319 on Test8 (variant assay) would have had a 2ml extraction
# 16468 on Test24 (variant assay) would have had a 2ml extraction
# 18891 on 20-4126 had 9ul input into each replicate for the variant assay
# 20817 would have had a 2ml extraction on both assays

typeof(for_GE_calc$Positives)

# Add on the target category
for_GE_calc_with_target <- for_GE_calc %>% 
  left_join(ddpcr_target_panel %>%
              select(Target, Target_category), by = "Target") %>%
  select(Well, Sample, Target, Target_category, Concentration, PoissonConfMax, PoissonConfMin,
         Positives, AcceptedDroplets, Sample_worksheet, replicates) %>%
  arrange(Sample) %>%
  pivot_wider(id_cols = c(Sample),
              names_from = Target_category,
              values_from = c(Concentration, PoissonConfMax, PoissonConfMin, Positives,
                              AcceptedDroplets, replicates)) %>%
  mutate(ff_assay_replicate_input_ul = ifelse(Sample == 18891, 10, 5)) %>%
  mutate(variant_assay_replicate_input_ul = ifelse(Sample == 18891, 9, 5)) %>%
  mutate(ff_extraction_input_ml = ifelse(Sample %in% c(13625, 16319, 16468, 20817), 2, 4)) %>%
  mutate(variant_extraction_input_ml = ifelse(Sample %in% c(13625, 16319, 16468, 20817), 2, 4)) %>%
  
  # Determine which targets are the fetals signal and the shared maternal allele
  mutate(Positives_maternal = pmax(Positives_ff_allele1, Positives_ff_allele2)) %>%
  mutate(Positives_paternal = pmin(Positives_ff_allele1, Positives_ff_allele2)) %>%
  mutate(Concentration_maternal = pmax(Concentration_ff_allele1, Concentration_ff_allele2)) %>%
  mutate(Concentration_paternal = pmin(Concentration_ff_allele1, Concentration_ff_allele2)) %>%
  mutate(PoissonConfMax_maternal = pmax(PoissonConfMax_ff_allele1, PoissonConfMax_ff_allele2)) %>%
  mutate(PoissonConfMax_maternal = pmax(PoissonConfMax_ff_allele1, PoissonConfMax_ff_allele2)) %>%
  mutate(PoissonConfMin_maternal = pmax(PoissonConfMin_ff_allele1, PoissonConfMin_ff_allele2)) %>%
  mutate(PoissonConfMax_paternal = pmin(PoissonConfMax_ff_allele1, PoissonConfMax_ff_allele2)) %>%
  mutate(PoissonConfMin_paternal = pmin(PoissonConfMin_ff_allele1, PoissonConfMin_ff_allele2)) %>%
  
  # Get rid of duplicated columns 
  select(-c("AcceptedDroplets_ff_allele2", "AcceptedDroplets_reference", "replicates_ff_allele2",
            "replicates_reference")) %>%
  
  # Calculate the total cfDNA detected by each assay
  mutate(Positives_ff_assay = Positives_maternal + Positives_paternal) %>%
  mutate(Positives_variant_assay = Positives_reference + Positives_variant) %>%
  
  # Calculate concentration in copies per ul for each assay. Each droplet is 0.00085ul (0.85nl)
  mutate(Concentration_ff_assay = -log((AcceptedDroplets_ff_allele1 - Positives_ff_assay)/ 
                                         AcceptedDroplets_ff_allele1) / 0.00085) %>%
  mutate(Concentration_variant_assay = -log((AcceptedDroplets_variant - Positives_variant_assay)/ 
                                         AcceptedDroplets_variant) / 0.00085) %>%
  
  mutate(Concentration_ff_assay_max = Concentration_ff_assay + 1.96*(sqrt(Concentration_ff_assay / AcceptedDroplets_ff_allele1*0.00085))) %>%
  mutate(Concentration_ff_assay_min = Concentration_ff_assay - 1.96*(sqrt(Concentration_ff_assay / AcceptedDroplets_ff_allele1*0.00085))) %>%
  mutate(Concentration_variant_assay_max = Concentration_variant_assay + 1.96*(sqrt(Concentration_variant_assay / AcceptedDroplets_variant*0.00085))) %>%
  mutate(Concentration_variant_assay_min = Concentration_variant_assay - 1.96*(sqrt(Concentration_variant_assay / AcceptedDroplets_variant*0.00085))) %>%
  
  # Calculate the copies per ul of the cfDNA extraction elution
  mutate(cfDNA_elution_concentration_ff_assay = Concentration_ff_assay * (replicates_ff_allele1*22) /
           (replicates_ff_allele1*ff_assay_replicate_input_ul)) %>%
  mutate(cfDNA_GE_ml_plasma_ff_assay = (cfDNA_elution_concentration_ff_assay * 60) / ff_extraction_input_ml) %>%
  
  mutate(cfDNA_elution_concentration_ff_assay_max = Concentration_ff_assay_max * (replicates_ff_allele1*22) /
           (replicates_ff_allele1*ff_assay_replicate_input_ul)) %>%
  mutate(cfDNA_GE_ml_plasma_ff_assay_max = (cfDNA_elution_concentration_ff_assay_max * 60) / ff_extraction_input_ml) %>%
  
  mutate(cfDNA_elution_concentration_ff_assay_min = Concentration_ff_assay_min * (replicates_ff_allele1*22) /
           (replicates_ff_allele1*ff_assay_replicate_input_ul)) %>%
  mutate(cfDNA_GE_ml_plasma_ff_assay_min = (cfDNA_elution_concentration_ff_assay_min * 60) / ff_extraction_input_ml) %>%
  
  mutate(cfDNA_elution_concentration_variant_assay = Concentration_variant_assay * (replicates_variant*22) /
           (replicates_variant*variant_assay_replicate_input_ul)) %>%
  mutate(cfDNA_GE_ml_plasma_variant_assay = (cfDNA_elution_concentration_variant_assay * 60) / variant_extraction_input_ml) %>%
  
  mutate(cfDNA_elution_concentration_variant_assay_max = Concentration_variant_assay_max * (replicates_variant*22) /
           (replicates_variant*variant_assay_replicate_input_ul)) %>%
  mutate(cfDNA_GE_ml_plasma_variant_assay_max = (cfDNA_elution_concentration_variant_assay_max * 60) / variant_extraction_input_ml) %>%
  
  mutate(cfDNA_elution_concentration_variant_assay_min = Concentration_variant_assay_min * (replicates_variant*22) /
           (replicates_variant*variant_assay_replicate_input_ul)) %>%
  mutate(cfDNA_GE_ml_plasma_variant_assay_min = (cfDNA_elution_concentration_variant_assay_min * 60) / variant_extraction_input_ml) %>%
  
  mutate(cffDNA_elution_concentration_ff_assay = (Concentration_paternal * (replicates_ff_allele1*22)) /
           (replicates_ff_allele1*ff_assay_replicate_input_ul)) %>%
  
  mutate(cffDNA_elution_concentration_ff_assay_max = PoissonConfMax_paternal * (replicates_ff_allele1*22) / 
           (replicates_ff_allele1*ff_assay_replicate_input_ul)) %>%
  
  mutate(cffDNA_elution_concentration_ff_assay_min = PoissonConfMin_paternal * (replicates_ff_allele1*22) / 
           (replicates_ff_allele1*ff_assay_replicate_input_ul)) %>%
  
  mutate(cffDNA_GE_ml_plasma_ff_assay = (cffDNA_elution_concentration_ff_assay * 60) / ff_extraction_input_ml) %>%
  
  mutate(cffDNA_GE_ml_plasma_ff_assay_max = (cffDNA_elution_concentration_ff_assay_max * 60) / ff_extraction_input_ml) %>%
  
  mutate(cffDNA_GE_ml_plasma_ff_assay_min = (cffDNA_elution_concentration_ff_assay_min * 60) / ff_extraction_input_ml)
  
cfDNA_conc_gestations <- left_join(
  bespoke_cohort_unblinded, 
  for_GE_calc_with_target %>%
    mutate(r_number = as.numeric(Sample)),
  by = "r_number")

ggplot(cfDNA_conc_gestations, aes(x = cfDNA_GE_ml_plasma_variant_assay, y = cfDNA_GE_ml_plasma_ff_assay, 
                                  shape = Call))+
  geom_point(size = 5, alpha = 0.6)+
  scale_shape_manual(values = c(19, 1))+
  geom_errorbarh(aes(xmin = cfDNA_GE_ml_plasma_variant_assay_min, xmax = cfDNA_GE_ml_plasma_variant_assay_max)) +
  geom_errorbar(aes(ymin = cfDNA_GE_ml_plasma_ff_assay_min, ymax = cfDNA_GE_ml_plasma_ff_assay_max))+
  ylim(0, 7000)+
  xlim(0, 7000)+
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none")+
  geom_abline(linetype = "dashed", alpha = 0.5)+
  labs(x = "cfDNA concentration from variant/reference ddPCR (GE/ml plasma)",
       y = "cfDNA concentration from fetal fraction ddPCR (GE/ml plasma)")


ggplot(cfDNA_conc_gestations, aes(x = Gestation_total_weeks, y = cffDNA_GE_ml_plasma_ff_assay,
                                  shape = Call))+
  geom_point(size = 3)+
  scale_shape_manual(values = c(19, 1))+
  geom_errorbar(aes(ymin = cffDNA_GE_ml_plasma_ff_assay_min, ymax = cffDNA_GE_ml_plasma_ff_assay_max))+
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none")+
  xlim(0, 40)+
  labs(x = "Gestational Age (weeks)", y = "Fetal specific allele (GE/ml plasma)")

# Perform a t test to see if the fetal fraction and variant fraction assays amplify cfDNA equally
t.test(cfDNA_conc_gestations$cfDNA_GE_ml_plasma_variant_assay, cfDNA_conc_gestations$cfDNA_GE_ml_plasma_ff_assay,
       paired = TRUE)

#############################################################
# Sickle cell disease twin case (20915)
#############################################################

twin_sample <- sickle_cell_blinded %>%
  filter(r_number %in% c(20915, 13751))

colnames(sickle_cell_blinded)

sample_20915 <- sickle_cell_blinded %>%
  filter(r_number == 20915)

sample_13751 <- sickle_cell_blinded %>%
  filter(r_number == 13751)

ff_error_20915 <- sample_20915$Fetal_fraction_max_percent - sample_20915$Fetal_fraction_min_percent
ff_20915 <- sample_20915$Fetal_fraction_percent

ff_error_13751 <- sample_13751$Fetal_fraction_max_percent - sample_13751$Fetal_fraction_min_percent
ff_13751 <- sample_13751$Fetal_fraction_percent

# Plot the results of singleton vs twin pregnancy
ggplot(twin_sample, aes(x = r_number, y = Variant_fraction_percent))+
  geom_point()+
  geom_errorbar(aes(ymin = Variant_fraction_min_percent, 
                    ymax = Variant_fraction_max_percent), width =0.1)+
  theme_bw()+
  
  ## Twin sample classification areas
  
  # Region for classification as both heterozygous
  geom_rect(aes(xmin = 1.5, xmax = Inf, ymin = 50+(ff_error_20915/2), ymax = 50-(ff_error_20915/2)),
            alpha = 0.1,
            fill = "red")+
  # Region for classification as one heterozygous, one homozygous variant
  geom_rect(aes(xmin = 1.5, xmax = Inf, ymin = 50+(ff_20915/4) -(ff_error_20915/2), 
            ymax = 50+(ff_20915/4) +(ff_error_20915/2)),
            alpha = 0.1,
            fill = "red")+
  # Region for classification as both homozygous variant
  geom_rect(aes(xmin = 1.5, xmax = Inf, ymin = 50+(ff_20915/2) -(ff_error_20915/2), 
                ymax = 50+(ff_20915/2) +(ff_error_20915/2)),
            alpha = 0.1,
            fill = "red")+
  # Region for classification as one heterozygous, one homozygous reference
  geom_rect(aes(xmin = 1.5, xmax = Inf, ymin = 50-(ff_20915/4) -(ff_error_20915/2), 
                ymax = 50-(ff_20915/4) +(ff_error_20915/2)),
            alpha = 0.1,
            fill = "red")+
  # Region for classification as both homozygous reference
  geom_rect(aes(xmin = 1.5, xmax = Inf, ymin = 50-(ff_20915/2) -(ff_error_20915/2), 
                ymax = 50-(ff_20915/2) +(ff_error_20915/2)),
            alpha = 0.1,
            fill = "red")+
  
  ## Singleton sample classification areas
  ## Twin sample classification areas
  
  # Region for classification as both heterozygous
  geom_rect(aes(xmin = -Inf, xmax = 1.5, ymin = 50+(ff_error_13751/2), ymax = 50-(ff_error_13751/2)),
            alpha = 0.1,
            fill = "red")+
  # Region for classification as homozygous variant
  geom_rect(aes(xmin = -Inf, xmax = 1.5, ymin = 50+(ff_13751/2) -(ff_error_13751/2), 
                ymax = 50+(ff_13751/2) +(ff_error_13751/2)),
            alpha = 0.1,
            fill = "red")+
  # Region for classification as homozygous reference
  geom_rect(aes(xmin = -Inf, xmax = 1.5, ymin = 50-(ff_13751/2) -(ff_error_13751/2), 
                ymax = 50-(ff_13751/2) +(ff_error_13751/2)),
            alpha = 0.1,
            fill = "red")+
  labs(x = "", y = "Variant fraction (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(aes(xintercept=1.5), linetype="dashed", size = 1)+
  
  # Add on labels
  annotate(geom="text", x=1.65, y=58, label="Twin 1: HbSS
Twin 2: HbSS")+
  annotate(geom="text", x=1.65, y=54, label="Twin 1: HbAS
Twin 2: HbSS")+
  annotate(geom="text", x=1.7, y=50, label="Twin 1: HbAS or HbSS
Twin 2: HbAS or HbAA")+
  annotate(geom="text", x=1.65, y=46, label="Twin 1: HbAS
Twin 2: HbAA")+
  annotate(geom="text", x=1.65, y=42, label="Twin 1: HbAA
Twin 2: HbAA") +
  annotate(geom="text", x=0.55, y=46, label="Fetus HbAA")+
  annotate(geom="text", x=0.55, y=50, label="Fetus HbAS")+
  annotate(geom="text", x=0.55, y=54, label="Fetus HbSS")+
  annotate(geom="text", x=2.1, y=60, label="Sample 20915: dichorionic diamniotic twin 
pregnancy sample with fetal fraction of 16.1%", size = 4, fontface = "bold")+
  annotate(geom="text", x=1, y=60, label="Sample 13751: singleton pregnancy 
sample with fetal fraction of 9.0%", size = 4, fontface = "bold")

# Set the likelihood ratio
LR <- 250

sample_20915_twin_SPRT <- sample_20915 %>%
  mutate(q0 = 0.5) %>%
  # Prediction for one homozygous fetus, one heterozygous
  mutate(q1 = 0.5+(Fetal_fraction/4)) %>%
  # Prediction for both fetuses homozygous
  mutate(q2 = 0.5+(Fetal_fraction/2)) %>%
  mutate(Delta_1 = (1- q1)/(1-q0)) %>%
  mutate(Delta_2 = (1- q2)/(1-q1)) %>%
  mutate(Gamma_1 = ((q1 * (1-q0))/ (q0*(1-q1)))) %>%
  mutate(Gamma_2 = ((q2 * (1-q1))/ (q1*(1-q2)))) %>%
  # Likelihood ratio for one homozygous
  mutate(LR_1 = exp((((Over_represented_fraction*log(Gamma_1)) + log(Delta_1))*Molecules_variant_assay))) %>%
  # Likelihood ratio for both homozygous
  mutate(LR_2 = exp((((Over_represented_fraction*log(Gamma_2)) + log(Delta_2))*Molecules_variant_assay))) %>%
  # Calculate the threshold variant fractions for each classification
  mutate(threshold_AS_SS_lower = (((log(LR)/Molecules_variant_assay) - log(Delta_1)) / log(Gamma_1))*100) %>%
  mutate(threshold_SS_SS = (((log(LR)/Molecules_variant_assay) - log(Delta_2)) / log(Gamma_2))*100) %>%
  mutate(threshold_AS_AS_upper = (((log(1/LR)/Molecules_variant_assay) - log(Delta_1)) / log(Gamma_1))*100) %>%
  mutate(threshold_AS_SS_upper = (((log(1/LR)/Molecules_variant_assay) - log(Delta_2)) / log(Gamma_2))*100) %>%
  mutate(threshold_AA_AA = 50-(threshold_SS_SS-50)) %>%
  mutate(threshold_AS_AA_upper = 50- (threshold_AS_SS_lower-50)) %>%
  mutate(threshold_AS_AS_lower = 50-(threshold_AS_AS_upper-50)) %>%
  mutate(threshold_AS_AA_lower = 50-(threshold_AS_SS_upper-50))


ggplot(sample_20915_twin_SPRT, aes(x = Identifier, y = Over_represented_fraction_percent))+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = Over_represented_fraction_min_percent, 
                    ymax = Over_represented_fraction_max_percent), width =0)+
  # Region for classification as both homozygous variant
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = threshold_SS_SS, ymax = Inf),
            alpha = 0.1,
            fill = "grey")+
  # Region for classification as one fetus homozygous variant
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = threshold_AS_SS_lower, 
                ymax = threshold_AS_SS_upper),
            alpha = 0.1,
            fill = "grey")+
  # Region for classification heterozygous
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = threshold_AS_AS_lower, 
                ymax = threshold_AS_AS_upper),
            alpha = 0.1,
            fill = "grey")+
  # Region for classification as one fetus homozygous reference
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = threshold_AS_AA_lower, 
                ymax = threshold_AS_AA_upper),
            alpha = 0.1,
            fill = "grey")+
  # Region for classification as both homozygous reference
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = threshold_AA_AA),
            alpha = 0.1,
            fill = "grey")+
  ylim(42, 58)+
  labs(x = "", y = "")+
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18), 
        plot.title = element_text(size=20),
        legend.position = "bottom", legend.background = element_rect(fill="white"), 
        legend.title = element_blank(), legend.text = element_text(size= 15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  # Add on labels
  annotate(geom="text", x=0.65, y=58, label="HbSS:HbSS", size = 4)+
  annotate(geom="text", x=0.65, y=54, label="HbAS:HbSS", size = 4)+
  annotate(geom="text", x=0.65, y=50, label="HbAS:HbAS
or
HbAA:HbSS", size = 4)+
  annotate(geom="text", x=0.65, y=46, label="HbAS:HbAA", size = 4)+
  annotate(geom="text", x=0.65, y=42, label="HbAA:HbAA",size = 4)


# Plot for presentation
ggplot(sample_20915_twin_SPRT, aes(x = Identifier, y = Over_represented_fraction_percent))+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = Over_represented_fraction_min_percent, 
                    ymax = Over_represented_fraction_max_percent), width =0)+
  # Region for classification as both homozygous variant
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = threshold_SS_SS, ymax = Inf),
            alpha = 0.2,
            fill = "blue")+
  # Region for classification as one fetus homozygous variant
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = threshold_AS_SS_lower, 
                ymax = threshold_AS_SS_upper),
            alpha = 0.2,
            fill = "blue")+
  # Region for classification heterozygous
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = threshold_AS_AS_lower, 
                ymax = threshold_AS_AS_upper),
            alpha = 0.2,
            fill = "blue")+
  # Region for classification as one fetus homozygous reference
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = threshold_AS_AA_lower, 
                ymax = threshold_AS_AA_upper),
            alpha = 0.2,
            fill = "blue")+
  # Region for classification as both homozygous reference
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = threshold_AA_AA),
            alpha = 0.2,
            fill = "blue")+
  ylim(42, 58)+
  labs(x = "", y = "Variant fraction (%)")+
  theme_bw()+
  theme(axis.text.y=element_text(size=15), axis.text.x=element_blank(),
        axis.title = element_text(size=18), 
        plot.title = element_text(size=20),
        legend.position = "bottom", legend.background = element_rect(fill="white"), 
        legend.title = element_blank(), legend.text = element_text(size= 15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


## Empty plot
ggplot(sample_20915_twin_SPRT %>%
         filter(Fetal_fraction_percent > 100), aes(x = Identifier, y = Over_represented_fraction_percent))+
  geom_point(size = 3)+
  ylim(42, 58)+
  labs(x = "", y = "Variant fraction (%)")+
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18), 
        plot.title = element_text(size=20),
        legend.position = "bottom", legend.background = element_rect(fill="white"), 
        legend.title = element_blank(), legend.text = element_text(size= 15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#############################################################
# Enzyme cut site analysis
#############################################################

# Specify DNASE1L3 end motifs from Serpas et al PMID: 30593563 in both forward and reverse
# complement.
dnase1l3_end_motifs <- c("CCCA", "CCAG", "CCTG", "CCAA", "CCCT", "CCAT", 
                         "TGGG", "CTGG", "CAGG", "TTGG", "AGGG", "ATGG") 

var_amplicons <- c("IDS_c182_189del_var", "ABCD1_c3GA_var", 
                   "TCOF1_c3611_C_A_v2_var", "COL4A5_c1295GAv2_var")

top_CC_assays <- c("IDS_c182_189del_var", "IDS_c182_189del_ref", "ABCD1_c3GA_var", 
                   "ABCD1_c3GA_ref", "TCOF1_c3611_C_A_v2_var", "TCOF1_c3611_C_A_v2_ref", 
                   "COL4A5_c1295GAv2_var", "COL4A5_c1295GAv2_ref")

# Read in amplicons
amplicons <- read_excel("W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/ddPCR Assay Design/ddPCR_Assay_Ordering.xlsx",
                        sheet = "amplicons")

reference_amplicons <- amplicons %>%
  mutate(dnase1l3_CCCA = str_count(amplicon_reference, "CCCA")) %>%
  mutate(dnase1l3_CCAG = str_count(amplicon_reference, "CCAG")) %>%
  mutate(dnase1l3_CCTG = str_count(amplicon_reference, "CCTG")) %>%
  mutate(dnase1l3_CCAA = str_count(amplicon_reference, "CCAA")) %>%
  mutate(dnase1l3_CCCT = str_count(amplicon_reference, "CCCT")) %>%
  mutate(dnase1l3_CCAT = str_count(amplicon_reference, "CCAT")) %>%
  mutate(dnase1l3_TGGG = str_count(amplicon_reference, "TGGG")) %>%
  mutate(dnase1l3_CTGG = str_count(amplicon_reference, "CTGG")) %>%
  mutate(dnase1l3_CAGG = str_count(amplicon_reference, "CAGG")) %>%
  mutate(dnase1l3_TTGG = str_count(amplicon_reference, "TTGG")) %>%
  mutate(dnase1l3_AGGG = str_count(amplicon_reference, "AGGG")) %>%
  mutate(dnase1l3_ATGG = str_count(amplicon_reference, "ATGG"))

dnase1l3_sum <- rowSums(reference_amplicons[,3:14])

amplicons_motifs <- cbind(reference_amplicons, dnase1l3_sum)

ggplot(amplicons_motifs %>%
         filter(!assay_name %in% var_amplicons) %>%
         mutate(assay_name = fct_reorder(assay_name, dnase1l3_sum)), 
       aes(x = dnase1l3_sum, y = assay_name))+
  geom_col()+
  theme_bw()+
  labs( x = "Number of DNAse1L3 CC motifs in target amplicon", y = "ddPCR assay",
        title = "Frequency of CC motifs in target reference amplicons")

ggplot(amplicons_motifs %>%
         filter(assay_name %in% top_CC_assays) %>%
         arrange(assay_name), aes(x = dnase1l3_sum, y = assay_name))+
  geom_col()+
  theme_bw()+
  labs( x = "Number of DNAse1L3 CC motifs in target amplicon", y = "ddPCR amplicon",
        title = "Frequency of CC motifs in reference and variant amplicons")



#############################################################
# Individual relative mutation dosage plots
#############################################################

# This is the function for plotting relative mutation dosage (rmd) 
# results for ddPCR, including parental gDNA controls. The amount of 
# information on this plot can be modified to suit user preference.

plot_rmd_graph <- function(cfdna_sample, maternal, paternal){
  # Get the sample variant fraction information
  variant_cfdna_sample <- ddpcr_data_analysed %>%
    filter(r_number == cfdna_sample) %>%
    select(r_number, Molecules_variant, variant_assay,
           Molecules_reference, AcceptedDroplets_Variant_assay) %>%
    pivot_longer(
      cols = !c(r_number, variant_assay, AcceptedDroplets_Variant_assay),
      names_to = "Target",
      values_to = "count"
    ) %>%
    rename(assay = variant_assay, AcceptedDroplets = AcceptedDroplets_Variant_assay) %>%
    mutate(identity = "cfDNA") %>%
    select(r_number, assay, AcceptedDroplets, identity, Target, count)
  
  # Get the sample fetal fraction information
  ff_cfdna_sample <- ddpcr_data_analysed %>%
    filter(r_number == cfdna_sample) %>%
    select(r_number, ff_assay, AcceptedDroplets_FetalFrac, Molecules_maternal,
           Molecules_paternal) %>%
    pivot_longer(
      cols = !c(r_number, ff_assay, AcceptedDroplets_FetalFrac),
      names_to = "Target",
      values_to = "count"
    ) %>%
    rename(assay = ff_assay, AcceptedDroplets = AcceptedDroplets_FetalFrac) %>%
    mutate(identity = "cfDNA") %>%
    select(r_number, assay, AcceptedDroplets, identity, Target, count)
  
  # Get the parental control information
  ddpcr_control_tbl_ff_no_id <- pivotted_controls_ff %>%
    left_join(control_table_ff %>%
                select(Worksheet_well, ff_assay), by = "Worksheet_well",
              .groups="drop") %>%
    # Remove duplicate columns and rename to be compatible with functions
    rename(AcceptedDroplets_FetalFrac = AcceptedDroplets_ff_allele1)
  
  # Specify whether maternal or paternal
  ddpcr_control_tbl_ff_id <- ddpcr_control_tbl_ff_no_id %>%
    left_join(controls %>%
                select(Sample, identity), by = "Sample",
              .groups="drop") %>%
    filter(Sample == maternal | Sample == paternal)
  
  # Calculate the molecules for each target type (maternal or paternal)
  ddpcr_control_tbl_ff_id_molecules <- ddpcr_control_tbl_ff_id %>%
    mutate(Positives_maternal = case_when(
      identity == "paternal gDNA" ~pmin(Positives_ff_allele1, Positives_ff_allele2),
      identity == "maternal gDNA" ~ pmax(Positives_ff_allele1, Positives_ff_allele2))) %>%
    mutate(Positives_paternal = case_when(
      identity == "paternal gDNA" ~pmax(Positives_ff_allele1, Positives_ff_allele2),
      identity == "maternal gDNA" ~ pmin(Positives_ff_allele1, Positives_ff_allele2))) %>% 
    mutate(Molecules_maternal = Poisson_correct(AcceptedDroplets_FetalFrac,Positives_maternal)) %>%
    mutate(Molecules_paternal = Poisson_correct(AcceptedDroplets_FetalFrac,Positives_paternal))
  
  ddpcr_control_tbl_var_id <- ddpcr_control_tbl_var %>%
    left_join(controls %>%
                select(Sample, identity), by = "Sample",
              .groups="drop") %>%
    filter(Sample == maternal | Sample == paternal)
  
  # Get the parental control variant fraction information
  control_variant_all <- ddpcr_control_tbl_var_id %>%
    select(Worksheet_well, Molecules_variant, variant_assay, AcceptedDroplets_Variant_assay,
           Molecules_reference, identity) %>%
    pivot_longer(
      cols = !c(Worksheet_well, variant_assay, AcceptedDroplets_Variant_assay, identity),
      names_to = "Target",
      values_to = "count"
    ) %>%
    rename(r_number = Worksheet_well, assay = variant_assay,
           AcceptedDroplets = AcceptedDroplets_Variant_assay) %>%
    select(r_number, assay, AcceptedDroplets, identity, Target, count)
  
  # Get the parental control ff information
  control_ff_all <- ddpcr_control_tbl_ff_id_molecules %>%
    select(Worksheet_well, ff_assay, AcceptedDroplets_FetalFrac, Molecules_maternal,
           Molecules_paternal, identity) %>%
    pivot_longer(
      cols = !c(Worksheet_well, ff_assay, AcceptedDroplets_FetalFrac, identity),
      names_to = "Target",
      values_to = "count") %>%
    rename(r_number = Worksheet_well, assay = ff_assay, 
           AcceptedDroplets = AcceptedDroplets_FetalFrac) %>%
    select(r_number, assay, AcceptedDroplets, identity, Target, count)
  
  
  # Split the parental data and select only one well from each
  # If no paternal sample is available, this should just return and
  # empty tibble.
  
  maternal_ff <- head(control_ff_all %>%
                        filter(identity == "maternal gDNA") %>%
                        arrange(r_number), n = 2)
  
  paternal_ff <- head(control_ff_all %>%
                        filter(identity == "paternal gDNA") %>%
                        arrange(r_number), n = 2)
  
  maternal_var <- head(control_variant_all %>%
                         filter(identity == "maternal gDNA") %>%
                         arrange(r_number), n = 2)
  
  paternal_var <- head(control_variant_all %>%
                         filter(identity == "paternal gDNA") %>%
                         arrange(r_number), n = 2)
  
  # Bind the tables
  plot_table_cfdna_sample <- rbind(maternal_ff, paternal_ff, maternal_var, paternal_var,
                                   ff_cfdna_sample, variant_cfdna_sample) %>%
    mutate(r_number = as.character(r_number))%>%
    mutate(id = paste(r_number, assay)) %>%
    mutate(Cpd = count/AcceptedDroplets) %>%
    mutate(count_max = Poisson_max(Cpd, AcceptedDroplets)) %>%
    mutate(count_min = Poisson_min(Cpd, AcceptedDroplets)) %>%
    mutate(new_target = case_when(
      Target == "Molecules_maternal" ~"Maternal SNP",
      Target == "Molecules_paternal" ~"Paternal SNP",
      Target == "Molecules_variant" ~"Variant allele",
      Target == "Molecules_reference" ~"Reference allele")) %>%
    mutate(assay_type = case_when(
      new_target %in% c("Paternal SNP", "Maternal SNP") ~"fetal_fraction",
      new_target %in% c("Variant allele", "Reference allele") ~"variant_fraction")) %>%
    mutate(case = cfdna_sample)
  
  # Make them factors to define the order in the plot
  plot_table_cfdna_sample$new_target <- factor(plot_table_cfdna_sample$new_target, levels = c("Maternal SNP",
                                                                                              "Paternal SNP",
                                                                                              "Reference allele",
                                                                                              "Variant allele"))
  
  plot_table_cfdna_sample$identity <- factor(plot_table_cfdna_sample$identity, levels = c("maternal gDNA",
                                                                                          "paternal gDNA", "cfDNA"))
  
  plot_table_cfdna_sample$assay_type <- factor(plot_table_cfdna_sample$assay_type, levels = c("fetal_fraction",
                                                                                              "variant_fraction"))
  
  plot_table_cfdna_sample <- plot_table_cfdna_sample%>%
    arrange(assay_type, identity, new_target)
  
  
  # Add the information for the plot subtitle
  sprt_result <- as.character(mcmc_vs_sprt %>% 
                                filter(r_number == cfdna_sample) %>%
                                select(SPRT_prediction))
  
  mcmc_result <- as.character(mcmc_vs_sprt %>%
                                filter(r_number == cfdna_sample) %>%
                                select(mcmc_prediction))
  
  # Fetal fraction result
  ff_result <- as.character(round(ddpcr_nipd_unblinded %>% 
                                    filter(r_number == cfdna_sample) %>%
                                    select(Fetal_fraction_percent), digits = 1))
  
  # Variant fraction result
  vf_result <- as.character(round(ddpcr_nipd_unblinded %>% 
                                    filter(r_number == cfdna_sample) %>%
                                    select(Variant_fraction_percent), digits = 1))
  
  sprt_result_subtitle <- paste0("SPRT: ", sprt_result, "   ", "MCMC: ", mcmc_result, 
                                 "   ", "Fetal fraction: ", ff_result, "%", 
                                 "   ", "Variant fraction: ", vf_result, "%")
  
  # Plot the results
  rmd_plot <- ggplot(plot_table_cfdna_sample, aes(x = identity, y = count, fill = new_target))+
    geom_col(position = position_dodge(width = 0.9), colour="black", alpha = 0.6)+
    
    scale_fill_manual(values = c("#99FFFF", "#FFCC99", "#3366FF", "#FF0000"))+
    geom_errorbar(aes(ymin = count_min, ymax = count_max, width = 0.3), position = position_dodge(width = 0.9))+
    theme_bw()+
    theme(axis.text=element_text(size=18), axis.title = element_text(size=18),
          legend.position = "bottom", legend.title = element_blank(),
          legend.text = element_text(size= 18),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x = "", y = "Molecules", title = paste("Sample:", cfdna_sample, " ", "Variant assay:", 
                                                variant_cfdna_sample$assay, "  ",
                                                "Fetal fraction assay:", ff_cfdna_sample$assay),
         subtitle = paste(sprt_result_subtitle))
  
  return(rmd_plot)
}

# Examples
#plot_rmd_graph(30113, "21RG-126G0134", "")
#plot_rmd_graph(20611, "21RG-126G0140", "")
#plot_rmd_graph(13519, "20RG-307G0060", "20RG-307G0062")

#############################################################
# Plotting graphs of sample cohorts
#############################################################

# Plot the results of fetal fraction versus variant fraction, including inconclusive results.
ggplot(sickle_cell_unblinded, aes(x = Fetal_fraction_percent, y = Variant_fraction_percent))+
  geom_point(aes(shape = overall_prediction, alpha = overall_prediction), size = 4)+
  # Change the shape of the shapes displayed
  scale_shape_manual(values = c(15, 17, 18, 19, 19))+
  scale_colour_manual(values = c("#FFFFFF"))+
  geom_errorbarh(aes(xmin = Fetal_fraction_min_percent, xmax = Fetal_fraction_max_percent,
                     alpha = overall_prediction)) +
  geom_errorbar(aes(ymin = Variant_fraction_min_percent, ymax = Variant_fraction_max_percent,
                    alpha = overall_prediction))+
  scale_alpha_manual(values = c(0.8,0.8,0.8,0.2, 0.8))+
  fifty_percent_line +
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18), 
        plot.title = element_text(size=20),
        legend.position = "bottom", legend.background = element_rect(fill="white"), 
        legend.title = element_blank(), legend.text = element_text(size= 15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Fetal fraction (%)", y = "Variant fraction (%)")+
  ylim(42.5, 57.5)+
  xlim(0, 22)+
  geom_text(data=subset(sickle_cell_blinded, Identifier %in% c("HBB-57", "HBB-62", "HBB-64")),
            aes(label=Identifier, vjust = 1.5, hjust = 1.2))

# Plot the results without calls or legend
ggplot(sickle_cell_unblinded, aes(x = Fetal_fraction_percent, y = Variant_fraction_percent))+
  geom_point(size = 8, alpha = 0.5)+
  # Change the shape of the shapes displayed
  geom_errorbarh(aes(xmin = Fetal_fraction_min_percent, xmax = Fetal_fraction_max_percent)) +
  geom_errorbar(aes(ymin = Variant_fraction_min_percent, ymax = Variant_fraction_max_percent))+
  fifty_percent_line +
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18), 
        plot.title = element_text(size=20),
        legend.position = "none", legend.background = element_rect(fill="white"), 
        legend.title = element_blank(), legend.text = element_text(size= 15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Fetal fraction (%)", y = "Variant fraction (%)")+
  ylim(42.5, 57.5)+
  xlim(0, 22)


# Empty plot without results
ggplot(sickle_cell_unblinded %>%
         filter(Fetal_fraction_percent > 100), aes(x = Fetal_fraction_percent, y = Variant_fraction_percent))+
  geom_point(size = 8, alpha = 0.5)+
  # Change the shape of the shapes displayed
  geom_errorbarh(aes(xmin = Fetal_fraction_min_percent, xmax = Fetal_fraction_max_percent)) +
  geom_errorbar(aes(ymin = Variant_fraction_min_percent, ymax = Variant_fraction_max_percent))+
  fifty_percent_line +
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18), 
        plot.title = element_text(size=20),
        legend.position = "none", legend.background = element_rect(fill="white"), 
        legend.title = element_blank(), legend.text = element_text(size= 15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Fetal fraction (%)", y = "Variant fraction (%)")+
  ylim(42.5, 57.5)+
  xlim(0, 22)


# Plot the results without calls or legend or confidence intervals
ggplot(sickle_cell_unblinded, aes(x = Fetal_fraction_percent, y = Variant_fraction_percent))+
  geom_point(size = 8, alpha = 0.5)+
  # Change the shape of the shapes displayed
  fifty_percent_line +
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18), 
        plot.title = element_text(size=20),
        legend.position = "none", legend.background = element_rect(fill="white"), 
        legend.title = element_blank(), legend.text = element_text(size= 15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Fetal fraction (%)", y = "Variant fraction (%)")+
  ylim(42.5, 57.5)+
  xlim(0, 22)


# Plot the SCD results in colour without error bars
ggplot(sickle_cell_unblinded, aes(x = Fetal_fraction_percent, y = Variant_fraction_percent))+
  # Add colours for the classifications. Also put them in the right order for the legend: homozygous
  # variant at the top
  # 0000FF is dark blue
  # 999999 is grey
  # 000000 is black
  # 99CCFF is light blue
  # Factor order is alphabetical: het, hom normal, hom variant, no call
  scale_fill_manual(values=c("#000000", "#99CCFF", "#999999", "#0000FF","#FF0000"),
                    breaks=c("homozygous variant", "heterozygous", "no call", "homozygous reference", "twin pregnancy"))+
  geom_point(size = 8, aes(fill = overall_prediction), colour="black", pch=21, alpha = 0.8)+
  scale_colour_manual(values = c("#FFFFFF"))+
  fifty_percent_line +
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18), 
        plot.title = element_text(size=20),
        legend.background = element_rect(fill="white"),
        legend.position = c(0.85, 0.85), 
        legend.title = element_blank(), legend.text = element_text(size= 15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Fetal fraction (%)", y = "Variant fraction (%)")+
  ylim(42.5, 57.5)+
  xlim(0, 22)

# Plot the SCD results for ESHG 2021 Poster
ggplot(sickle_cell_unblinded %>%
         filter(!r_number %in% c(13262)), aes(x = Fetal_fraction_percent, y = Variant_fraction_percent))+
  # Add colours for the classifications. Also put them in the right order for the legend: homozygous
  # variant at the top
  # 0000FF is dark blue
  # 999999 is grey
  # 000000 is black
  # 99CCFF is light blue
  # Factor order is alphabetical: het, hom normal, hom variant, no call
  scale_fill_manual(values=c("#000000", "#99CCFF", "#999999", "#0000FF","#FF0000"),
                    breaks=c("homozygous variant", "heterozygous", "no call", "homozygous reference", "twin pregnancy"))+
  geom_point(size = 8, aes(fill = overall_prediction), colour="black", pch=21, alpha = 0.8)+
  scale_colour_manual(values = c("#FFFFFF"))+
  fifty_percent_line +
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18), 
        plot.title = element_text(size=20),
        legend.background = element_rect(colour = 'black', fill = 'white', 
                                         size = 1),
        legend.position = c(0.8, 0.15), 
        legend.title = element_blank(), legend.text = element_text(size= 18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Fetal fraction (%)", y = "Variant fraction (%)")+
  xlim(0, 25)



cfDNA_for_graph <- sickle_cell_unblinded %>%
  filter(overall_prediction != "no call") %>%
  arrange(overall_prediction, Identifier) %>%
  mutate(Identifier=factor(Identifier, levels=Identifier)) %>%
  select(Identifier, Variant_fraction_percent, Variant_fraction_max_percent, Variant_fraction_min_percent)

all_DNA_for_graph <- rbind(het_gDNA_arranged, cfDNA_for_graph) %>%
  mutate(Identifier=factor(Identifier, levels=Identifier))

ggplot(sickle_cell_unblinded %>%
         filter(overall_prediction != "no call") %>%
         arrange(overall_prediction, Identifier) %>%
         mutate(Identifier=factor(Identifier, levels=Identifier)), aes(x = Identifier, y = Variant_fraction_percent))+
  geom_point(stat="identity")+
  geom_errorbar(aes(ymin = Variant_fraction_min_percent, ymax = Variant_fraction_max_percent))+
  fifty_percent_line +
  ddPCR_plot_theme+  
  labs(x = "", y = "Variant fraction (%)")

variant_fraction_plot <- ggplot(all_DNA_for_graph, aes(x = Identifier, y = Variant_fraction_percent))+
  geom_point(stat="identity")+
  geom_errorbar(aes(ymin = Variant_fraction_min_percent, ymax = Variant_fraction_max_percent))+
  geom_hline(aes(yintercept=50), linetype="dashed", size = 1)+
  geom_vline(aes(xintercept=23.5), linetype="dashed", colour = "grey", size = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 15), axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Samples", y = "Variant fraction (%)")+
  ylim(42, 58)

variant_fraction_plot + annotate(geom="text", x=10, y=57, label="HbAS gDNA and paternal cfDNA") +
  annotate(geom="text", x=35, y=57, label="cfDNA from pregnant HbAS carriers")

# Plot the results of all samples with HbAS fetuses by fetal and variant fraction.
ggplot(sickle_cell_unblinded %>%
         filter(invasive_result == "HbAS"), aes(x = Fetal_fraction_percent, y = Variant_fraction_percent,
                                                colour = overall_prediction))+
  geom_point()+
  geom_errorbar(aes(ymin = Variant_fraction_min_percent, ymax = Variant_fraction_max_percent))+
  geom_errorbarh(aes(xmin = Fetal_fraction_min_percent, xmax = Fetal_fraction_max_percent))+
  fifty_percent_line +
  ddPCR_plot_theme+  
  labs(x = "Total DNA molecules", y = "Variant fraction (%)")+
  ylim(47, 53)

# Plot the results of all samples with HbAS fetuses by variant fraction and total number of molecules
ggplot(sickle_cell_unblinded %>%
         filter(invasive_result == "HbAS"), aes(x = Molecules_variant_assay, y = Variant_fraction_percent,
                                                colour = overall_prediction))+
  geom_point()+
  geom_errorbar(aes(ymin = Variant_fraction_min_percent, ymax = Variant_fraction_max_percent))+
  fifty_percent_line +
  ddPCR_plot_theme+  
  labs(x = "Total DNA molecules", y = "Variant fraction (%)")+
  ylim(47, 53)

ggplot(bespoke_cohort_analysed, aes(x = Fetal_fraction_percent, y = Variant_fraction_percent))+
  geom_point(size = 4)+
  geom_errorbarh(aes(xmin = Fetal_fraction_min_percent, xmax = Fetal_fraction_max_percent)) +
  geom_errorbar(aes(ymin = Variant_fraction_min_percent, ymax = Variant_fraction_max_percent))+
  fifty_percent_line +
  ddPCR_plot_theme+
  labs(x = "Fetal fraction (%)", y = "Variant fraction (%)")

# Plot graph of variant imbalance vs fetal fraction.
autosomal_cohort <- ddpcr_data_tbl %>%
  filter(Inheritance == "autosomal")

autosomal_cohort_analysed <- calc_SPRT(calc_conf_intervals(calc_molecules(autosomal_cohort)), 8)

autosomal_cohort_imbalance <- autosomal_cohort_analysed %>%
  filter(SPRT_prediction != "heterozygous")

ggplot(autosomal_cohort_imbalance %>%
         filter(Sample != 20915), aes(x = Fetal_fraction_percent, y = Over_represented_fraction_percent))+
  geom_point(aes(shape = variant_assay))+
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18), 
        plot.title = element_text(size=20),
        legend.position = "bottom", legend.background = element_rect(fill="white"), 
        legend.title = element_text(size= 18), legend.text = element_text(size= 15))+
  geom_errorbarh(aes(xmin = Fetal_fraction_min_percent, xmax = Fetal_fraction_max_percent)) +
  geom_errorbar(aes(ymin = Over_represented_fraction_min_percent, ymax = Over_represented_fraction_max_percent))+
  ylim(48, 60)+
  xlim(0, 20)+
  geom_abline(intercept = 50, linetype="dashed", colour = "grey", size = 1, slope = 0.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Fetal fraction %", y = "Over represented fraction %",
       title = "Over-represented fraction versus fetal fraction
       in ddPCR cases with homozygous fetal genotypes")+
  geom_smooth(se = FALSE, method = "lm")

#############################################################
# Export csvs with time stamps
#############################################################

current_time <- Sys.time()

write.csv(ddpcr_nipd_unblinded, 
          file = paste0("analysis_outputs/ddpcr_nipd_unblinded", 
                        format(current_time, "%Y%m%d_%H%M%S"), ".csv"),
          row.names = FALSE)

write.csv(mcmc_vs_sprt_outcomes, 
          file = paste0("analysis_outputs/mcmc_vs_sprt_outcomes", 
                        format(current_time, "%Y%m%d_%H%M%S"), ".csv"),
          row.names = FALSE)
