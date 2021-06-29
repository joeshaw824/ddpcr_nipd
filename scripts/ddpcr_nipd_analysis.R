#############################################################
## ddPCR for Non Invasive Prenatal Diagnosis (NIPD)
## January 2021
## Joseph.Shaw@gosh.nhs.uk
## This is an analysis script for the prediction of fetal genotypes from 
## cfDNA testing using ddPCR for maternally-inherited
## variants. This includes modules for analysis of ddPCR data
## using the sequential probability ratio test (SPRT) (Lo et al,
## 2008; PMID:  17664418), and via MonteCarlo Markov Chain 
## (MCMC) analysis (Caswell et al, 2020; PMID:  32533152).
## The MCMC analysis section  is compiled from scripts and
## models supplied by Tristan Snowsill (Exeter).
#############################################################

## Naming convention
# For consistency, variables are named as: target_category_qualifier
# Target examples: variant, reference, paternal, maternal, fetal, difference, 
#                   major_allele, minor_allele, ff_assay, vf_assay, total
# Category examples: molecules, fraction, percent, positives
# Qualifier examples: max, min

# Example: fetal_percent_max

#############################################################
# Load libraries and resources
#############################################################

## Load necessary packages
library(tidyverse)
# cmdstanr required for running stan models
library(cmdstanr)
# bayesplot and posterior can be added for MCMC diagnostic modelling
#library(bayesplot)
#library(posterior)

## Set working directory
setwd("W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/ddPCR_R_Analysis/ddpcr_nipd")

# Load resources
controls <- readr::read_csv("resources/controls.csv")
ddpcr_target_panel <- readr::read_csv("resources/ddpcr_target_panel.csv") 

# Source functions
source("functions/ddPCR_nipd_functions.R")
source("functions/RAPID_biobank.R")

#############################################################
# Read in ddPCR data 
#############################################################

# The data wrangling in this section is due to collating 
# data acquired over 3 years with various lab workflows.

# Load in all the csv files exported from QuantaSoft.

dataPath <- "data/ddpcr_data/"

ddpcr_files <- list.files(dataPath)

#Empty data frame
ddpcr_data <- data.frame()

# Read and collate each worksheet csv
for (dataFile in ddpcr_files){
  tmp_dat <- read_csv(paste0(dataPath,dataFile), col_names = TRUE)
  # Add a worksheet identifier and a unique identifier for each
  # well on each worksheet
  tmp_dat_ws <- tmp_dat %>%
    mutate(Worksheet = dataFile,
           Worksheet_well = paste(Worksheet, Well, sep = "_"))
  ddpcr_data <-rbind(ddpcr_data, tmp_dat_ws)
  rm(tmp_dat_ws)
  rm(tmp_dat)
}

#########################
# Molecule calculation functions
#########################

# These functions calculate lots of stuff
# This is to prevent massive duplication of code

var_ref_calculations <- function(data_input) {
  
  # The data_input needs these variables for this function to work:
  # AcceptedDroplets_Variant_assay, Positives_variant, Positives_reference
  
  stopifnot(c("AcceptedDroplets_Variant_assay", "Positives_variant", 
              "Positives_reference")
            %in% colnames(data_input))
  
  data_output <- data_input %>%
  
  mutate(
  # Perform Poisson correction to determine the total number of molecules 
  # detected for each target.
  variant_molecules = Poisson_correct(AcceptedDroplets_Variant_assay,
                                      Positives_variant), 
  reference_molecules = Poisson_correct(AcceptedDroplets_Variant_assay,
                                        Positives_reference),   
  
  # Calculate the total molecules for the variant fraction (vf) assay
  vf_assay_molecules = variant_molecules + reference_molecules,
  
  # Calculate the fractional abundance of each allele
  reference_fraction = reference_molecules / vf_assay_molecules,
  variant_fraction = variant_molecules / vf_assay_molecules,
  variant_percent  = variant_fraction * 100,
  reference_percent  = reference_fraction * 100,
  
  # Select the major and minor alleles
  major_allele = case_when(
    # If both fractions are exactly 50, then the variant is arbitrarily 
    # chosen as the major allele
    variant_percent >= reference_percent ~ "variant allele",
    reference_percent > variant_percent ~ "reference allele"),
  
  minor_allele = case_when(
    major_allele == "variant allele" ~ "reference allele",
    major_allele == "reference allele" ~ "variant allele"),
  
  # Choose the percentage of the major and minor alleles
  major_allele_percent = case_when(
    major_allele == "variant allele" ~variant_percent,
    major_allele == "reference allele" ~reference_percent),
  
  # Choose the molecules for each allele
  major_allele_molecules = case_when(
    major_allele == "variant allele" ~variant_molecules,
    major_allele == "reference allele" ~reference_molecules),
  
  minor_allele_molecules = case_when(
    minor_allele == "variant allele" ~variant_molecules,
    minor_allele == "reference allele" ~reference_molecules),
  
  # Find the difference between the numbers of molecules for each allele
  difference_molecules = abs(major_allele_molecules - 
                               minor_allele_molecules),
  
  # Calculate the 95% confidence intervals of the numbers of molecules
  # for each allele
  variant_molecules_max = Poisson_max((variant_molecules/ 
                                         AcceptedDroplets_Variant_assay), 
                                      AcceptedDroplets_Variant_assay),
  variant_molecules_min = Poisson_min((variant_molecules/
                                         AcceptedDroplets_Variant_assay),
                                      AcceptedDroplets_Variant_assay),
  reference_molecules_max = Poisson_max((reference_molecules/
                                           AcceptedDroplets_Variant_assay),
                                        AcceptedDroplets_Variant_assay),
  reference_molecules_min = Poisson_min((reference_molecules/
                                           AcceptedDroplets_Variant_assay), 
                                        AcceptedDroplets_Variant_assay),
  
  major_allele_molecules_max = case_when(
    major_allele == "variant allele" ~variant_molecules_max,
    major_allele == "reference allele" ~reference_molecules_max),
  
  major_allele_molecules_min = case_when(
    major_allele == "variant allele" ~variant_molecules_min,
    major_allele == "reference allele" ~reference_molecules_min),
  
  minor_allele_molecules_max = case_when(
    minor_allele == "variant allele" ~variant_molecules_max,
    minor_allele == "reference allele" ~reference_molecules_max),
  
  minor_allele_molecules_min = case_when(
    minor_allele == "variant allele" ~variant_molecules_min,
    minor_allele == "reference allele" ~reference_molecules_min),
  
  difference_molecules_max = major_allele_molecules_max - 
    minor_allele_molecules_min,
  difference_molecules_min = major_allele_molecules_min - 
    minor_allele_molecules_max,
  
  # Calculate the 95% confidence intervals of the fractional abundances
  variant_percent_max = (Poisson_fraction_max(
    variant_molecules_max, reference_molecules))*100,
  variant_percent_min = (Poisson_fraction_min(
    variant_molecules_min, reference_molecules))*100,
  reference_percent_max = (Poisson_fraction_max(
    reference_molecules_max, variant_molecules))*100,
  reference_percent_min = (Poisson_fraction_min(
    reference_molecules_min, variant_molecules))*100,
  
  major_allele_percent_max = case_when(
    major_allele == "variant allele" ~variant_percent_max,
    major_allele == "reference allele" ~reference_percent_max),
  
  major_allele_percent_min = case_when(
    major_allele == "variant allele" ~variant_percent_min,
    major_allele == "reference allele" ~reference_percent_min),
  
  # Calculate the 95% confidence intervals of the numbers of molecules 
  # detected by each assay
  vf_assay_molecules_max = Poisson_max((
    vf_assay_molecules/AcceptedDroplets_Variant_assay), 
    AcceptedDroplets_Variant_assay),
  vf_assay_molecules_min = Poisson_min((
    vf_assay_molecules/AcceptedDroplets_Variant_assay),
    AcceptedDroplets_Variant_assay))
  
  return(data_output)
}

ff_calculations <- function(data_input) {
  
  # The data_input needs these variables for this function to work:
  # maternal_molecules, paternal_molecules, AcceptedDroplets_FetalFrac
  
  stopifnot(c("maternal_positives", "paternal_positives",
              "AcceptedDroplets_FetalFrac")
            %in% colnames(data_input))
  
  data_output <- data_input %>%
    mutate(
    
    # Calculate the same variables for the fetal fraction assay
    
    maternal_molecules = Poisson_correct(
      AcceptedDroplets_FetalFrac,maternal_positives),
    paternal_molecules = Poisson_correct(
      AcceptedDroplets_FetalFrac,paternal_positives),
    
    # Calculate the fetal fraction
    fetal_fraction = calc_ff(maternal_molecules, paternal_molecules),
    fetal_percent = fetal_fraction*100,
    
    ff_assay_molecules = maternal_molecules + paternal_molecules,
    
    paternal_molecules_max = Poisson_max((
      paternal_molecules / AcceptedDroplets_FetalFrac), 
      AcceptedDroplets_FetalFrac),
    paternal_molecules_min = Poisson_min((
      paternal_molecules / AcceptedDroplets_FetalFrac), 
      AcceptedDroplets_FetalFrac),
    
    maternal_molecules_max = Poisson_max((
      maternal_molecules / AcceptedDroplets_FetalFrac), 
      AcceptedDroplets_FetalFrac),
    maternal_molecules_min = Poisson_min((
      maternal_molecules / AcceptedDroplets_FetalFrac), 
      AcceptedDroplets_FetalFrac),
    
    ff_assay_molecules_max = Poisson_max((
      ff_assay_molecules/AcceptedDroplets_FetalFrac), 
      AcceptedDroplets_FetalFrac),
    ff_assay_molecules_min = Poisson_min((
      ff_assay_molecules/AcceptedDroplets_FetalFrac), 
      AcceptedDroplets_FetalFrac),
    
    # When calculating for the fetal fraction, the paternal fraction 
    # must be multiplied by 2.
    fetal_percent_max = 200* (Poisson_fraction_max(
      paternal_molecules_max, maternal_molecules)),
    fetal_percent_min = 200* (Poisson_fraction_min(
      paternal_molecules_min, maternal_molecules)))
  
  return(data_output)
}

#########################
# Wrangle cfDNA data into shape
#########################

# Remove single wells and controls
ddpcr_data_merged_samples <- ddpcr_data %>%
  # Remove single wells and controls
  filter(substr(Well, 1, 1) == "M" & !(Sample %in% controls$Sample)) %>%
  # Count the number of wells which have been merged together,
  # which is the number of commas in "MergedWells" string plus one
  mutate(num_wells = str_count(MergedWells, ",")+1)

# Reshape the data frame to sum all values by Target.
# Often cfDNA from the same sample was tested across
# multiple worksheets and then must be summed together
# for the analysis.
ddpcr_data_reshape <- ddpcr_data_merged_samples %>% 
  group_by(Sample, Target) %>% 
  summarise(Positives = sum(Positives),
            AcceptedDroplets = sum(AcceptedDroplets),
            num_wells = sum(num_wells),
            .groups="drop")

# Add on the category for each ddPCR target
# Data model: 4 categories - variant, reference, ff_allele1 and ff_allele2
ddpcr_with_target <- ddpcr_data_reshape %>% 
  left_join(ddpcr_target_panel %>%
              select(Target, Target_category), by = "Target")

# Use pivot_wider to get one row per sample
pivotted_ddpcr <- ddpcr_with_target %>% 
  pivot_wider(id_cols = Sample,
              names_from = Target_category,
              values_from = c(AcceptedDroplets, Positives, num_wells))

# Add on assay and inheritance pattern to the table.
# Do this first by adding the information to the "variant" rows, then to
# the "ff_allele1" rows.
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
              select(Sample, Inheritance_chromosomal, Inheritance_pattern, 
                     variant_assay), by = "Sample",
            .groups="drop") %>% 
  left_join(ref_table_ff %>%
              select(Sample, ff_assay), by = "Sample",
            .groups="drop") %>%
  
  # Remove any samples which haven't had both assays performed.
  # This removes several samples which didn't have both assays performed
  # for various reasons.
  filter(!is.na(Positives_variant) & !is.na(Positives_ff_allele1)) %>%
  
  # Remove duplicate columns and rename to be compatible with functions
  select(-c(AcceptedDroplets_ff_allele2, AcceptedDroplets_reference, 
            num_wells_reference,
            num_wells_ff_allele2)) %>%
  rename(AcceptedDroplets_FetalFrac = AcceptedDroplets_ff_allele1,
         AcceptedDroplets_Variant_assay = AcceptedDroplets_variant,
         ff_assay_num_wells = num_wells_ff_allele1,
         vf_assay_num_wells = num_wells_variant) %>%
  
  # Determine the maternal and paternally-inherited alleles for the 
  # fetal fraction calculation
  mutate(maternal_positives = pmax(Positives_ff_allele1, Positives_ff_allele2),
         paternal_positives = pmin(Positives_ff_allele1, Positives_ff_allele2))

#########################
# Wrangle gDNA data into shape
#########################

# Get single well controls only without NTCs
ddpcr_controls <- ddpcr_data %>%
  filter(Sample %in% controls$Sample & Sample != "NTC" & 
           substr(Well, 1, 1) != "M")

ddpcr_controls_with_target <- ddpcr_controls %>% 
  left_join(ddpcr_target_panel %>%
              select(Target, Target_category), by = "Target")

# No need to pivot to get one row per sample because we want the control
# data as individual wells. Use pivot_wider to get one row per well for 
# each well tested with a variant assay
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
              select(Worksheet_well, Inheritance_chromosomal, 
                     Inheritance_pattern, variant_assay), 
            by = "Worksheet_well",
            .groups="drop") %>%
  
  # Remove duplicate columns and rename to be compatible with functions
  rename(AcceptedDroplets_Variant_assay = AcceptedDroplets_variant)

# Get the control information for the control wells tested with a 
# fetal fraction assay
ddpcr_control_tbl_ff <- pivotted_controls_ff %>%
  left_join(control_table_ff %>%
              select(Worksheet_well, ff_assay), by = "Worksheet_well",
            .groups="drop") %>%
  # Remove duplicate columns and rename to be compatible with functions
  rename(AcceptedDroplets_FetalFrac = AcceptedDroplets_ff_allele1)

#############################################################
# cfDNA SPRT analysis
#############################################################

# Set likelihood ratio threshold
LR_threshold <- 250

ddpcr_sprt_analysed <- ff_calculations(
  var_ref_calculations(ddpcr_data_tbl)) %>% 
  
  # Rename sample to r_number to allow merge with RAPID Biobank data.
  rename(r_number = Sample) %>%
  
  mutate(
    total_molecules = vf_assay_molecules + ff_assay_molecules,

  # Perform SPRT and return the likelihood ratio
  Likelihood_ratio = case_when(
     Inheritance_chromosomal == "x_linked" ~ 
       calc_LR_X_linked(fetal_fraction, (major_allele_percent/100), 
                        vf_assay_molecules),
     Inheritance_chromosomal == "autosomal" ~ 
       calc_LR_autosomal(fetal_fraction, (major_allele_percent/100), 
                         vf_assay_molecules)),
   
  # Classify based on likelihood ratio threshold supplied
  # Fetal genotype predictions are named consistently as 
  # "inconclusive", "heterozygous", "homozygous/hemizygous reference/variant"
  SPRT_prediction = case_when(
      
      Inheritance_chromosomal == "autosomal" &
      Likelihood_ratio > LR_threshold &
      major_allele == "reference allele" 
          ~ "homozygous reference",
      
      Inheritance_chromosomal == "autosomal" &
      Likelihood_ratio > LR_threshold &
      major_allele == "variant allele" 
          ~ "homozygous variant",
      
      Inheritance_chromosomal == "autosomal" &
      Likelihood_ratio < (1/LR_threshold)
          ~ "heterozygous",
      
      Inheritance_chromosomal == "autosomal" &
      Likelihood_ratio < LR_threshold &
      Likelihood_ratio > (1/LR_threshold) 
          ~ "inconclusive",
      
      Inheritance_chromosomal == "x_linked" &
      Likelihood_ratio > LR_threshold &
      major_allele == "reference allele" 
          ~ "hemizygous reference",
      
      Inheritance_chromosomal == "x_linked" &
      Likelihood_ratio > LR_threshold &
      major_allele == "variant allele" 
          ~ "hemizygous variant",
      
      Inheritance_chromosomal == "x_linked" &
      Likelihood_ratio < LR_threshold 
          ~ "inconclusive"))

#############################################################
# cfDNA MCMC analysis
#############################################################

# Prepare ddPCR data for MCMC

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
                Z_X = maternal_positives,
                Z_Y = paternal_positives) %>%
  arrange(Inheritance_chromosomal, Inheritance_pattern, variant_assay) %>%
  select(r_number, Inheritance_chromosomal, Inheritance_pattern, 
         variant_assay, n_K, n_Z, K_N, K_M, n_Z, Z_X, Z_Y)

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

# Set probability threshold for accepting fetal genotype predictions
mcmc_threshold <- 0.95

# Analyse the entire ddPCR cohort using Tristan's MCMC pipeline.
# (11/06/2021: this takes approximately 13 minutes)
# This generates the probabilities for the fetal genotype based on 
# the inheritance pattern of the variant.

# Autosomal dominant:
# p_G1: probability fetus is heterozygous
# p_G2: probability fetus is homozygous reference

# Autosomal recessive:
# p_G1: probability fetus is homozygous reference
# p_G2: probability fetus is heterozygous
# p_G3: probability fetus is homozygous variant

# X-linked:
# p_G0: probability fetus is hemizygous reference
# p_G1: probability fetus is hemizygous variant

ddpcr_with_fits <- ddpcr_data_mcmc %>%
  nest(data = n_K:Z_Y) %>%
  mutate(
    
  data = map(data, as.list),
         
  fit = case_when(
  Inheritance_chromosomal == "autosomal" &
  Inheritance_pattern == "dominant" ~map(
    data, ~ dominant_model$sample(data = .,
                                  init = initialise_chains_dominant,
                                  step_size = 0.2,
                                  parallel_chains = parallel::detectCores())),
   
  Inheritance_chromosomal == "autosomal" & 
  Inheritance_pattern == "recessive" ~ map(
    data, ~ recessive_model$sample(data = .,
                                  init = initialise_chains_recessive,
                                  step_size = 0.2,
                                  parallel_chains = parallel::detectCores())),
   
  Inheritance_chromosomal == "x_linked" ~map(
    data, ~ x_linked_model$sample(data = .,
                                 init = initialise_chains_xlinked,
                                 step_size = 0.2,
                                 parallel_chains = parallel::detectCores()))),
 
 results = case_when(
    Inheritance_chromosomal == "autosomal" &
    Inheritance_pattern == "dominant" ~map(
      fit,  ~ setNames(pivot_wider(.$summary(c("pG", "rho"), "mean"),
                                  names_from = "variable",
                                  values_from = "mean"),
                       c("p_G1", "p_G2", "rho_est"))),
   
   Inheritance_chromosomal == "autosomal" & 
    Inheritance_pattern == "recessive" ~map(
      fit,  ~ setNames(pivot_wider(.$summary(c("pG", "rho"), "mean"),
                                 names_from = "variable",
                                 values_from = "mean"),
                     c("p_G1", "p_G2", "p_G3", "rho_est"))),
   
   Inheritance_chromosomal == "x_linked" ~map(
      fit,  ~ setNames(pivot_wider(.$summary(c("pG", "rho"), "mean"),
                                  names_from = "variable",
                                  values_from = "mean"),
                      c("p_G0", "p_G1", "rho_est"))))) %>%
   
   unnest_wider(results)

# Add on fetal genotype predictions based on the previously set threshold.

ddpcr_mcmc_analysed <- ddpcr_with_fits %>%
  select(-c(data, fit)) %>%
  rename(fetal_fraction = rho_est) %>%
  mutate(mcmc_prediction = case_when(
    
    # Dominant predictions
    Inheritance_chromosomal == "autosomal" & 
    Inheritance_pattern == "dominant" & 
    p_G1 > mcmc_threshold ~"heterozygous",
    
    Inheritance_chromosomal == "autosomal" &
    Inheritance_pattern == "dominant" & 
    p_G2 > mcmc_threshold ~"homozygous reference",
    
    Inheritance_chromosomal == "autosomal" &
    Inheritance_pattern == "dominant" & 
    p_G1 < mcmc_threshold & 
    p_G2 < mcmc_threshold ~"inconclusive",
    
    # Recessive predictions
    Inheritance_chromosomal == "autosomal" &
    Inheritance_pattern == "recessive" & 
    p_G1 > mcmc_threshold ~"homozygous reference",
    
    Inheritance_chromosomal == "autosomal" &
    Inheritance_pattern == "recessive" &
    p_G2 > mcmc_threshold ~"heterozygous",
    # ok
    
    Inheritance_chromosomal == "autosomal" &
    Inheritance_pattern == "recessive" &
    p_G3 > mcmc_threshold ~"homozygous variant",
    
    Inheritance_chromosomal == "autosomal" &
    Inheritance_pattern == "recessive" & 
    p_G1 < mcmc_threshold & 
    p_G2 < mcmc_threshold ~"inconclusive",
    
    # X-linked predictions
    Inheritance_chromosomal == "x_linked" &
    p_G0 > mcmc_threshold ~"hemizygous reference",
    
    Inheritance_chromosomal == "x_linked" &
    p_G1 > mcmc_threshold ~"hemizygous variant",
    
    Inheritance_chromosomal == "x_linked" &
    p_G0 < mcmc_threshold &
    p_G1 < mcmc_threshold ~"inconclusive"))

#############################################################
# Collation of results
#############################################################

# This stage joins the SPRT and MCMC results together, and then compares
# them to the diagnostic results from RAPID Biobank.

# Samples to exlude:
# 13262 - this sample had contamination
# 17004 - this sample was actually HbAC
# 20915 - this sample was from a twin pregnancy
samples_to_exclude <- c(13262, 17004, 20915)

ddpcr_analysed <- left_join(
  ddpcr_sprt_analysed,
  ddpcr_mcmc_analysed %>%
    select(r_number, p_G0, p_G1, p_G2, p_G3, mcmc_prediction),
  by = "r_number") %>%
  filter(!r_number %in% samples_to_exclude)

ddpcr_nipd_unblinded <- left_join(
  ddpcr_analysed,
  RAPID_biobank %>%
    # Change r_number to a character to match ddpcr_analysed
    mutate(r_number = as.character(r_number)) %>%
    filter(r_number %in% ddpcr_analysed$r_number) %>%
    select(r_number, study_id, gestation_weeks, gestation_days, Gestation_total_weeks, gestation_character, 
           date_of_blood_sample, vacutainer, mutation_genetic_info_fetus, Partner_sample_available, 
           original_plasma_vol, tubes_plasma_current, report_acquired),
  by = "r_number")

# Rerunning the mcmc pipeline takes a while, so this is 
# a checkpoint for later.
current_time <- Sys.time()

write.csv(ddpcr_nipd_unblinded, 
          file = paste0("analysis_outputs/ddpcr_nipd_unblinded", 
                        format(current_time, "%Y%m%d_%H%M%S"), ".csv"),
          row.names = FALSE)

#############################################################
# Sickle cell disease analysis for paper
#############################################################

# The "secondary cohort" refers to the most recent phase of the 
# ddPCR sickle cell disease project, when all samples were extracted 
# using a 6ml protocol and the lab workflow was finalised.

secondary_cohort <- c("14182", "19868", "20238", "20611", 
                     "20874", "30063", "30068", "30113", "30142", 
                     "30206", "30228", "30230", "30078", "30065", 
                    "13402", "20939", "30215", "30203")

scd_ddpcr <- ddpcr_nipd_unblinded %>%
  filter(variant_assay == "HBB c.20A>T") %>%
  
  # Create an "algorithm prediction" based on the results of the two
  # statistical analyses and the technical data (see paper draft v8)
  
  mutate(algorithm_prediction = case_when(
    (SPRT_prediction == "heterozygous" | mcmc_prediction == "heterozygous") &
    variant_percent > 48.9 &
    variant_percent < 51.1 &
    difference_molecules < 200 
          ~"heterozygous",
    
    (SPRT_prediction == "homozygous reference" | 
       mcmc_prediction == "homozygous reference") & 
    variant_percent < 48.9 &
    difference_molecules > 200
          ~"homozygous reference",
    
    (SPRT_prediction == "homozygous variant" | 
       mcmc_prediction == "homozygous variant") &
    variant_percent > 51.1 & 
    difference_molecules > 200 
          ~"homozygous variant",
    
          TRUE ~"inconclusive"),
    
    # Create an "clinical prediction" of whether the fetus will
    # be affected or not.
    
    clinical_prediction = case_when(
      algorithm_prediction == "heterozygous" | 
        algorithm_prediction == "homozygous reference" 
                ~ "unaffected",
      algorithm_prediction == "homozygous variant" 
                ~ "affected",
      algorithm_prediction == "inconclusive" 
                ~"inconclusive"),
    
    clinical_outcome = case_when(
      mutation_genetic_info_fetus == "HbAA" | 
        mutation_genetic_info_fetus == "HbAS" ~"unaffected",
      mutation_genetic_info_fetus == "HbSS" ~"affected"),
      
      cohort = ifelse(r_number %in% secondary_cohort, "secondary", "primary"))

#########################
# gDNA and cfDNA plots
#########################

# The parental gDNA controls were tested at ~20ng per ddPCR well. 
# In order to have an equivalent dataset for comparison with cfDNA data, 
# I have down-sampled the genomic controls by selecting wells that bring 
# the total molecules tested to different levels for each sample
# (~6000 or ~9000 molecules)

wells_for_9000 <- c("21-1116.csv_C02", "21-1116.csv_D02", "21-1116.csv_G02", 
                    "21-1116.csv_H02", "21-1116.csv_H03", "21-1116.csv_A04", 
                    "21-1116.csv_D04", "21-1116.csv_E04", "21-1116.csv_E05", 
                    "21-1116.csv_F05", "21-1116.csv_A06", "21-1116.csv_B06",
                    "21-1116.csv_C06", "21-1413.csv_G03", "21-1413.csv_H03", 
                    "21-1413.csv_A04", "21-1413.csv_B04", "21-1413.csv_C04", 
                    "21-1413.csv_D04", "21-1413.csv_E06", "21-1413.csv_F06", 
                    "21-1413.csv_E06", "21-1413.csv_H06", "21-1413.csv_A07",
                    "21-1413.csv_B07", "21-1413.csv_A09", "21-1413.csv_B09", 
                    "21-1413.csv_D09", "21-1413.csv_E09", "21-1413.csv_F09", 
                    "21-1227.csv_F02", "21-1227.csv_G02", "21-1227.csv_A03", 
                    "21-1227.csv_B03", "21-1227.csv_C03", "21-1227.csv_H04",
                    "21-1227.csv_A05", "21-1227.csv_C05", "21-1227.csv_D05", 
                    "21-1227.csv_B07", "21-1227.csv_D07", "21-1227.csv_E07", 
                    "21-1227.csv_F07", "21-1227.csv_D09", "21-1227.csv_E09", 
                    "21-1227.csv_G09", "21-1227.csv_H09", "21-1705.csv_G04",
                    "21-1705.csv_H04", "21-1705.csv_B05", "21-1705.csv_C05", 
                    "21-1863.csv_G07", "21-1863.csv_H07", "21-1863.csv_A08", 
                    "21-1863.csv_B08", "21-1863.csv_C08", "21-1863.csv_E08", 
                    "21-1863.csv_F08", "21-1863.csv_H08", "21-1863.csv_A09",
                    "21-1946.csv_E08", "21-1946.csv_F08", "21-1946.csv_C09",
                    "21-1946.csv_D09")

wells_for_6000 <- c("21-1116.csv_C02", "21-1116.csv_G02", "21-1116.csv_H03", 
                    "21-1116.csv_D04", "21-1116.csv_E05", "21-1116.csv_F05", 
                    "21-1116.csv_A06", "21-1116.csv_B06", "21-1413.csv_G03",
                    "21-1413.csv_H03", "21-1413.csv_B04", "21-1413.csv_E06",
                    "21-1413.csv_F06", "21-1413.csv_E06", "21-1413.csv_H06", 
                    "21-1413.csv_A07", "21-1413.csv_A09", "21-1413.csv_D09",
                    "21-1413.csv_F09", "21-1227.csv_F02", "21-1227.csv_G02",
                    "21-1227.csv_A03", "21-1227.csv_B03", "21-1227.csv_H04",
                    "21-1227.csv_A05", "21-1227.csv_C05", "21-1227.csv_D05", 
                    "21-1227.csv_B07", "21-1227.csv_C07", "21-1227.csv_E07",
                    "21-1227.csv_F07", "21-1227.csv_D09", "21-1227.csv_E09", 
                    "21-1227.csv_G09", "21-1705.csv_G04", "21-1705.csv_B05", 
                    "21-1863.csv_G07", "21-1863.csv_A08", "21-1863.csv_B08", 
                    "21-1863.csv_C08", "21-1863.csv_E08", "21-1863.csv_H08",
                    "21-1946.csv_E08", "21-1946.csv_C09", "21-1946.csv_F09", 
                    "21-1946.csv_G10")

# Perform the var_ref_calculations for the parental controls in the 
# downsampled dataset

parents_downsampled <- var_ref_calculations(
  
              ddpcr_data %>%
              filter(Worksheet_well %in% wells_for_6000) %>%
              mutate(worksheet_sample = paste0(Worksheet,"_",Sample))%>%
              group_by(worksheet_sample, Target) %>% 
              summarise(Positives = sum(Positives),
                        AcceptedDroplets = sum(AcceptedDroplets),
                        .groups="drop") %>%
              pivot_wider(id_cols = worksheet_sample,
                          names_from = Target,
                          values_from = 
                            c(AcceptedDroplets, Positives)) %>%
              rename(
                AcceptedDroplets_Variant_assay =AcceptedDroplets_HbS,
                Positives_variant = Positives_HbA,
                Positives_reference = Positives_HbS) %>%
              mutate(
                sample_type = "gDNA",
                genotype = "het gDNA")) %>%
  
              rename(
                r_number = worksheet_sample,
                algorithm_prediction = genotype) %>%
              select(r_number, algorithm_prediction, 
                     vf_assay_molecules, sample_type, 
                     variant_molecules, reference_molecules,
                     variant_percent, variant_percent_max, 
                     variant_percent_min, difference_molecules, 
                     difference_molecules_max, 
                     difference_molecules_min)
  
cfDNA_compare <- scd_ddpcr %>%
  mutate(sample_type = "cfDNA") %>%
  select(r_number, algorithm_prediction, vf_assay_molecules, sample_type, 
         variant_molecules, reference_molecules,
         variant_percent, variant_percent_max, variant_percent_min,
         difference_molecules, difference_molecules_max, 
         difference_molecules_min)

# Bind cfDNA and genomic DNA results together
cfDNA_gDNA_bind <- rbind(cfDNA_compare, parents_downsampled) %>%
  mutate(x_axis_label = case_when(
    algorithm_prediction == "heterozygous" ~"cfDNA (heterozygous)",
    algorithm_prediction == "homozygous reference" ~"cfDNA (homozygous reference)",
    algorithm_prediction == "homozygous variant" ~"cfDNA (homozygous variant)",
    algorithm_prediction == "het gDNA" ~"gDNA (heterozgous)",
    algorithm_prediction == "inconclusive" ~ "inconclusive"),
    
    balanced_or_not = case_when(
      algorithm_prediction == "heterozygous" ~"het cfDNA",
      algorithm_prediction == "homozygous reference" ~"unbalanced",
      algorithm_prediction == "homozygous variant" ~"unbalanced",
      algorithm_prediction == "het gDNA" ~"het gDNA",
      algorithm_prediction == "inconclusive" ~ "inconclusive"))

# Downsampling plot
downsampling_plot <- ggplot(cfDNA_gDNA_bind, aes(x = 
                                       sample_type, y = vf_assay_molecules))+
  geom_boxplot()+
  theme_bw()+
  no_gridlines+
  ylim(0, 35000)+
  labs(y = "Total molecules detected by HBB c.20A>T ddPCR", 
       x = "",
       title = "Molecules detected in cfDNA and down-sampled gDNA datasets")

ggsave(filename = "plots/downsampling_plot.png", 
       plot =  downsampling_plot, dpi = 300)

# Make the two new variables factors to define the order in subsequent plots
cfDNA_gDNA_bind$x_axis_label <- factor(
  cfDNA_gDNA_bind$x_axis_label, levels = c("cfDNA (homozygous reference)",
                                            "gDNA (heterozgous)", 
                                            "cfDNA (heterozygous)", 
                                            "cfDNA (homozygous variant)",
                                            "inconclusive"))

cfDNA_gDNA_bind$balanced_or_not <- factor(
  cfDNA_gDNA_bind$balanced_or_not, levels = c("het gDNA", "het cfDNA",
                                              "unbalanced", "inconclusive"))

# Organise data for variant fraction plot
variant_fraction_data <- cfDNA_gDNA_bind %>%
  arrange(x_axis_label, variant_percent) %>%
  # Fix the order of samples
  mutate(row_number = as.numeric(row.names(cfDNA_gDNA_bind)))

# Variant fraction plot
variant_fraction_plot <- ggplot(variant_fraction_data, 
                                aes(x = row_number, 
                                    y = variant_percent, 
                                    colour = x_axis_label))+
  geom_errorbar(colour = "black", 
                aes(ymin = variant_percent_min, ymax = variant_percent_max))+
  # "#0000FF" = dark blue; "#FFFFFF" = white; "#99CCFF" = light blue
  # "#000000" = black; "#999999" = grey
  scale_fill_manual(values=c("#0000FF", "#FFFFFF", 
                             "#99CCFF", "#000000", "#999999"))+
  geom_point(size = 4, aes(fill = x_axis_label), colour = "black", pch=21)+
  theme_bw()+
  theme(axis.text.x = element_blank())+
  ylim(42, 58)+
  no_gridlines+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(), 
        legend.position = "bottom")+
  geom_hline(yintercept = 48.9, linetype = "dashed")+
  geom_hline(yintercept = 51, linetype = "dashed")+
  labs(x = "", y = "Variant fraction (%)",
       title = "Variant fraction in gDNA and cfDNA ddPCR results")

ggsave(filename = "plots/variant_fraction_plot.png", 
       plot = variant_fraction_plot, dpi = 300)

# Organise data for molecules of difference plot
molecules_difference_data <- cfDNA_gDNA_bind %>%
  arrange(balanced_or_not, difference_molecules) %>%
  mutate(row_number = as.numeric(row.names(cfDNA_gDNA_bind)))

# Change factor level for plot legend
molecules_difference_data$x_axis_label <- factor(
  molecules_difference_data$x_axis_label, levels = c(
                                            "gDNA (heterozgous)", 
                                            "cfDNA (heterozygous)",
                                            "cfDNA (homozygous reference)",
                                            "cfDNA (homozygous variant)",
                                            "inconclusive"))

# Molecules of difference plot
molecules_difference_plot <- ggplot(molecules_difference_data, 
                                    aes(x = row_number, 
                                        y = difference_molecules))+
  geom_errorbar(colour = "black", aes(ymin = difference_molecules_min, 
                                      ymax = difference_molecules_max))+
  scale_fill_manual(values=
                      c("#FFFFFF", "#000000", "#99CCFF", "#FF0000", "#999999"))+
  geom_point(size = 4, aes(fill = x_axis_label), colour = "black", pch=21)+
  theme_bw()+
  no_gridlines+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),legend.position = "bottom")+
  labs(x = "", y = "Molecules of difference",
       title = "Molecules of difference in gDNA and cfDNA ddPCR results")+
  geom_hline(yintercept = 200, linetype = "dashed")+
  scale_y_continuous(breaks = c(0, 200, 1000, 2000, 3000))

ggsave(filename = "plots/molecules_difference_plot.png", 
       plot = molecules_difference_plot, dpi = 300)

# This function is for selecting different sections of the
# cfDNA_gDNA_bind table and generating a plot.
# Input is the algorithm_prediction defined above

plot_cohort <- function(cohort) {
  
  stopifnot(cohort %in% c("heterozygous", "homozygous reference", 
                          "homozygous variant", "het gDNA",
                          "inconclusive"))
  
  cohort_data <- cfDNA_gDNA_bind %>%
  filter(algorithm_prediction == cohort) %>%
  select(r_number, algorithm_prediction, variant_molecules, 
         reference_molecules) %>%
  arrange(algorithm_prediction, variant_molecules) %>%
  mutate(row_number = seq(1, nrow(cfDNA_gDNA_bind %>%
                                filter(algorithm_prediction == cohort))))%>%
  pivot_longer(
    cols = c(variant_molecules, reference_molecules),
    names_to = "target",
    values_to = "count")
  
  cohort_plot <- ggplot(cohort_data, aes(x = row_number, y = count, 
                                         fill = target))+
    geom_col(width = 0.7, position = "dodge", colour="black", alpha = 0.4)+
    scale_fill_manual(values = c("#3366FF", "#FF0000"))+
    xlim(0, 43)+
    ylim(0, 17000)+
    theme_bw()+
    labs(y = "Molecules detected")+
    theme(axis.text.x=element_blank(), axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14),
          legend.position = "bottom", legend.title = element_blank(),
          legend.text = element_text(size= 12),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  return(cohort_plot)
}

plot_cohort("heterozygous")+
  labs(title = "cfDNA heterozygous predictions")

plot_cohort("homozygous reference")+
  labs(title = "cfDNA homozygous reference predictions")

plot_cohort("homozygous variant")+
  labs(title = "cfDNA homozygous variant predictions")

plot_cohort("het gDNA")+
  labs(title = "Heterozygous gDNA")

plot_cohort("inconclusive")+
  labs(title = "Inconclusive results")

#########################
# Cohort plots
#########################
no_gridlines <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Plot the results of fetal fraction versus variant fraction, including inconclusive results.
ggplot(scd_ddpcr, aes(x =fetal_percent, y = variant_percent))+
  # Add colours for the classifications. Also put them in the right order for the legend: homozygous
  # variant at the top
  # 0000FF is dark blue
  # 999999 is grey
  # 000000 is black
  # 99CCFF is light blue
  # Factor order is alphabetical: het, hom normal, hom variant, no call
  
  scale_fill_manual(values=c("#000000", "#99CCFF", "#999999", "#0000FF"),
                    breaks=c("homozygous variant", "heterozygous", "inconclusive", "homozygous reference"))+
  geom_point(size = 8, aes(fill = algorithm_prediction), colour="black", pch=21, alpha = 0.8)+
  scale_colour_manual(values = c("#FFFFFF"))+
  geom_errorbarh(alpha = 0.2, aes(xmin = fetal_percent_min, xmax = fetal_percent_max)) +
  geom_errorbar(alpha = 0.2, aes(ymin = variant_percent_min, ymax = variant_percent_max))+
  fifty_percent_line +
  theme_bw()+
  ylim(42, 58)+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18), 
        plot.title = element_text(size=20),
        legend.position = "bottom", legend.background = element_rect(fill="white"), 
        legend.title = element_blank(), legend.text = element_text(size= 15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Fetal fraction (%)", y = "Variant fraction (%)")

ggplot(scd_ddpcr, aes(x = Gestation_total_weeks, y = Fetal_fraction_percent))+
  geom_point(size = 3, alpha = 0.4)+
  geom_errorbar(alpha = 0.2, aes(ymin = Fetal_fraction_min_percent, ymax = Fetal_fraction_max_percent))+  
  theme_bw()+
  xlim(0, 40)

count(scd_ddpcr, algorithm_prediction, mutation_genetic_info_fetus)

# Primary and secondary cohort plot
primary_secondary_cohort_plot <- ggplot(scd_ddpcr, aes(x = cohort, y = Molecules_variant_assay))+
  geom_boxplot()+
  theme_bw()+
  ylim(0, 35000)+
  no_gridlines+
  labs(y = "Total molecules detected by HBB c.20A>T ddPCR", x = "",
       title = "Molecules detected by HBB c.20A>T ddPCR in primary and secondary cfDNA cohorts")

ggsave(filename = "plots/primary_secondary_cohort_plot.png", plot =  primary_secondary_cohort_plot, dpi = 300)

#########################
# RMD plot function
#########################

# This is the function for plotting relative mutation dosage (rmd) 
# results for ddPCR, including parental gDNA controls. The amount of 
# information on this plot can be modified to suit user preference.

plot_rmd_graph <- function(cfdna_sample, maternal, paternal){
  # Get the sample variant fraction information
  variant_cfdna_sample <- ddpcr_analysed %>%
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
  ff_cfdna_sample <- ddpcr_analysed %>%
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
    mutate(paternal_molecules = Poisson_correct(AcceptedDroplets_FetalFrac,Positives_paternal))
  
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
           paternal_molecules, identity) %>%
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
  
  
  # Add the information for the plot title and subtitle
  sprt_result <- as.character(ddpcr_analysed %>% 
                                filter(r_number == cfdna_sample) %>%
                                select(SPRT_prediction))
  
  mcmc_result <- as.character(ddpcr_analysed %>%
                                filter(r_number == cfdna_sample) %>%
                                select(mcmc_prediction))
  
  
  sample_id <- as.character(ddpcr_nipd_unblinded %>%
                              filter(r_number == cfdna_sample) %>%
                              select(r_number))
  
  # Fetal fraction result
  ff_result <- as.character(round(ddpcr_analysed %>% 
                                    filter(r_number == cfdna_sample) %>%
                                    select(Fetal_fraction_percent), digits = 1))
  
  # Variant fraction result
  vf_result <- as.character(round(ddpcr_analysed %>% 
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
          legend.text = element_text(size= 14),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x = "", y = "Molecules", title = paste("Sample:", sample_id, " ", "Variant assay:", 
                                                variant_cfdna_sample$assay, "  ",
                                                "Fetal fraction assay:", ff_cfdna_sample$assay),
         subtitle = paste(sprt_result_subtitle))
  
  return(rmd_plot)
}

#########################
# Individual plots
#########################

# Plot individual graphs for each sickle cell disease case
rmd_13402 <- plot_rmd_graph(13402, "21RG-120G0077", "")
rmd_14182 <- plot_rmd_graph(14182, "21RG-062G0108", "21RG-062G0111")
rmd_19868 <- plot_rmd_graph(19868, "21RG-083G0126", "21RG-083G0132")
rmd_30113 <- plot_rmd_graph(30113, "21RG-126G0134", "")
rmd_20238 <- plot_rmd_graph(20238, "21RG-103G0120", "21RG-103G0122")
rmd_20611 <- plot_rmd_graph(20611, "21RG-126G0140", "")
rmd_20874 <- plot_rmd_graph(20874, "21RG-103G0118", "21RG-103G0119")
rmd_20939 <- plot_rmd_graph(20939, "21RG-126G0131", "")
rmd_30063 <- plot_rmd_graph(30063, "21RG-103G0115", "21RG-103G0117")
rmd_30068 <- plot_rmd_graph(30068, "21RG-120G0099", "21RG-120G0103")
rmd_30113 <- plot_rmd_graph(30113, "21RG-126G0134", "")
rmd_30142 <- plot_rmd_graph(30142, "21RG-083G0112",	"21RG-083G0120")
rmd_30206 <- plot_rmd_graph(30206, "21RG-120G0092", "21RG-120G0097")
rmd_30228 <- plot_rmd_graph(30228, "21RG-103G0061", "21RG-103G0062")
rmd_30230 <- plot_rmd_graph(30230, "21RG-103G0120", "21RG-103G0122")
rmd_30065 <- plot_rmd_graph(30065, "21RG-126G0126", "")
rmd_30078 <- plot_rmd_graph(30078, "21RG-126G0124", "")
rmd_30203 <- plot_rmd_graph(30203, "21RG-138G0155", "")
rmd_30215 <- plot_rmd_graph(30215, "21RG-138G0159", "")

# Save plots in a composite pdf
pdf(file = "plots/scd_cohort_individual_plots.pdf", onefile = TRUE, 
    width=10, height=6)
rmd_13402
rmd_14182
rmd_19868
rmd_20238
rmd_20611
rmd_20874
rmd_20939
rmd_30063
rmd_30065
rmd_30068
rmd_30078
rmd_30113
rmd_30142
rmd_30203
rmd_30206
rmd_30215
rmd_30228
rmd_30230
dev.off()

# Save plots individually using ggsave
ggsave(filename = "plots/rmd_14182.png", plot =  rmd_14182, dpi = 300)
ggsave(filename = "plots/rmd_30113.png", plot =  rmd_30113, dpi = 300)
ggsave(filename = "plots/rmd_19868.png", plot =  rmd_19868, dpi = 300)
ggsave(filename = "plots/rmd_20238.png", plot =  rmd_20238, dpi = 300)
ggsave(filename = "plots/rmd_20611.png", plot =  rmd_20611, dpi = 300)
ggsave(filename = "plots/rmd_20874.png", plot =  rmd_20874, dpi = 300)
ggsave(filename = "plots/rmd_30063.png", plot =  rmd_30063, dpi = 300)
ggsave(filename = "plots/rmd_30068.png", plot =  rmd_30068, dpi = 300)
ggsave(filename = "plots/rmd_30113.png", plot =  rmd_30113, dpi = 300)
ggsave(filename = "plots/rmd_30142.png", plot =  rmd_30142, dpi = 300)
ggsave(filename = "plots/rmd_30206.png", plot =  rmd_30206, dpi = 300)
ggsave(filename = "plots/rmd_30228.png", plot =  rmd_30228, dpi = 300)
ggsave(filename = "plots/rmd_30230.png", plot =  rmd_30230, dpi = 300)
ggsave(filename = "plots/rmd_30078.png", plot =  rmd_30078, dpi = 300)
ggsave(filename = "plots/rmd_30065.png", plot =  rmd_30065, dpi = 300)
ggsave(filename = "plots/rmd_13402.png", plot =  rmd_13402, dpi = 300)
ggsave(filename = "plots/rmd_20939.png", plot =  rmd_20939, dpi = 300)

## Replace the NAs with ""
## Could have a loop function to plot the graphs for the whole cohort
# for i in table {
# plot <- plot_rmd_graph(r_number, mat, pat)



# Plot individual graphs for each bespoke case
# X linked cases
plot_rmd_graph(10280, "18G10515", "") # Needs ZFXY for control
plot_rmd_graph(11928, "18G10519", "") # Needs ZFXY for control
plot_rmd_graph(12585, "18G05794", "") # Needs ZFXY for control
plot_rmd_graph(12945, "18G10521", "") # Needs ZFXY for control
plot_rmd_graph(13625, "18G10690", "") # Needs ZFXY for control and variant assay at appropriate dilution
plot_rmd_graph(13965, "18G10691", "")
plot_rmd_graph(14247, "18G10693", "") # Needs ZFXY for control
plot_rmd_graph(14917, "18G10694", "") # Needs ZFXY for control
plot_rmd_graph(16319, "18G05795", "") # Needs ZFXY for control and variant assay at appropriate dilution
plot_rmd_graph(16468, "18G10769", "") # Needs ZFXY for control and variant assay at appropriate dilution
plot_rmd_graph(16881, "18G05796", "") # Needs ZFXY for control
plot_rmd_graph(17667, "21RG-027G0004", "21RG-027G0006") # Need to pull dad's data through
plot_rmd_graph(18164, "18G05797", "") # Needs ZFXY for control
plot_rmd_graph(18385, "18G10771", "") # Needs ZFXY for control
plot_rmd_graph(18891, "18G10772", "")
plot_rmd_graph(19397, "18G10773", "") # Needs ZFXY for control
plot_rmd_graph(19611, "21RG-027G0070", "21RG-027G0075")
plot_rmd_graph(20817, "18G10522", "") # Needs ZFXY for control and variant assay at appropriate dilution
plot_rmd_graph(20980, "21RG-027G0010", "21RG-027G0014")
plot_rmd_graph(30030, "21RG-027G0017", "")

# Dominant cases
plot_rmd_graph(12990, "20RG-148G0079", "")
plot_rmd_graph(13519, "20RG-307G0060", "20RG-307G0062")
plot_rmd_graph(14116, "20RG-148G0076", "") 
plot_rmd_graph(14142, "21RG-047G0089", "21RG-047G0093")
plot_rmd_graph(14491, "20RG-148G0081", "")
plot_rmd_graph(14522, "20RG-148G0075", "")
plot_rmd_graph(19261, "20RG-307G0064", "20RG-307G0070")
plot_rmd_graph(19711, "20RG-148G0080", "")

# Recessive cases
plot_rmd_graph(17531, "20RG-336G0056", "") # Need SNP data for control
plot_rmd_graph(17841, "21RG-047G0111", "21RG-047G0112")
plot_rmd_graph(19102, "20RG-336G0066", "") # Need SNP data for control

# How many cases have plasma remaining?
ddpcr_nipd_unblinded %>%
  filter(variant_assay != "HBB c.20A>T") %>%
  select(r_number, variant_assay, tubes_plasma_current)

#########################
# cfDNA Amplifiability
#########################

scd_qubit <- read_excel("W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/ddPCR Worksheets/scd_labwork_plan.xlsx", sheet = "scd_labwork_plan")

phase3_results <- left_join(ddpcr_nipd_unblinded %>%
                  filter(r_number %in% phase3_samples),
                scd_qubit %>%
                  # Comvert r_number to character to allow join
                  mutate(r_number = as.character(r_number)) %>%
                  select(r_number, cfDNA_Qubit, cfDNA_vol),
                by = "r_number") %>%
  
                  # Calculate the molecules per ml plasma input
                  mutate(
                    ff_locus_molecules_ul = Molecules_ff_assay / (num_wells_ff_assay*5),
                    ff_locus_molecules_ul_max = Molecules_ff_assay_max / (num_wells_ff_assay*5),
                    ff_locus_molecules_ul_min = Molecules_ff_assay_min / (num_wells_ff_assay*5),
                    ff_locus_molecules_ml = as.integer((ff_locus_molecules_ul * cfDNA_vol) / original_plasma_vol),
                    ff_locus_molecules_ml_max = as.integer((ff_locus_molecules_ul_max * cfDNA_vol) / original_plasma_vol),
                    ff_locus_molecules_ml_min = as.integer((ff_locus_molecules_ul_min * cfDNA_vol) / original_plasma_vol),
                    fetal_specific_molecules_ul = Molecules_paternal / (num_wells_ff_assay*5),
                    fetal_specific_molecules_ul_max = Molecules_paternal_max / (num_wells_ff_assay*5),
                    fetal_specific_molecules_ul_min = Molecules_paternal_min / (num_wells_ff_assay*5),
                    fetal_specific_molecules_ml = as.integer((fetal_specific_molecules_ul * cfDNA_vol) / original_plasma_vol),
                    fetal_specific_molecules_ml_max = as.integer((fetal_specific_molecules_ul_max * cfDNA_vol) / original_plasma_vol),
                    fetal_specific_molecules_ml_min = as.integer((fetal_specific_molecules_ul_min * cfDNA_vol) / original_plasma_vol),
                    
                    # Now calculate how many molecules we'd expect from the Qubit for the HbASv3 assay
                    # All wells had 5ul cfDNA input, 1ng is 1000 pg and 1 haploid human genome weighs 3.3pg
                    expected_copies_from_Qubit = ((cfDNA_Qubit * num_wells_variant_assay *5)*1000)/3.3,
                    
                    # Calculate the total cfDNA extracted
                    total_cfDNA_ng = cfDNA_Qubit * cfDNA_vol,
                    
                    # Amplification ratio is the observed over expected counts
                    amplifcation_ratio = (Molecules_variant_assay / expected_copies_from_Qubit)*100)

# Compare molecules detected by ddPCR to molecules expected from Qubit
# Error bars are too small to see so aren't plotted
amplifiability_plot <- ggplot(phase3_results, aes(x = expected_copies_from_Qubit, Molecules_variant_assay))+
  geom_point(size = 3)+
  theme_bw()+
  xlim(0, 51000)+
  ylim(0, 51000)+
  geom_abline(linetype = "dashed")+
  no_gridlines+
  labs(x = "Expected molecules from Qubit reading", y = "Molecules detected by ddPCR",
       title = "Expected versus observed molecular counts for the HBB c.20A>T ddPCR assay")

ggsave(filename = "plots/amplifiability_plot.png", plot =  amplifiability_plot, dpi = 300)

# Plot the molecules of total cfDNA and cffDNA per millilitre plasma
total_cfDNA_plot <- ggplot(phase3_results, aes(x = Gestation_total_weeks, y = ff_locus_molecules_ml))+
  geom_errorbar(aes(ymin = ff_locus_molecules_ml_min, ymax = ff_locus_molecules_ml_max))+
  geom_point()+
  theme_bw()+
  xlim(0, 40)+
  no_gridlines+
  labs(x = "Gestation (weeks)", y = "Total haploid genome equivalents per millilitre plasma (GE/ml)",
       tile = "Total cfDNA by gestation")

ggsave(filename = "plots/total_cfDNA_plot.png", plot =  total_cfDNA_plot, dpi = 300)

# Plot the molecules of fetal-specific cffDNA per millilitre plasma
cffDNA_plot <- ggplot(phase3_results, aes(x = Gestation_total_weeks, y = fetal_specific_molecules_ml))+
  geom_point()+
  geom_errorbar(aes(ymin = fetal_specific_molecules_ml_min, ymax = fetal_specific_molecules_ml_max))+
  theme_bw()+
  xlim(0, 40)+
  scale_y_continuous(breaks = seq(0, 300, by = 50))+
  no_gridlines+
  labs(x = "Gestation (weeks)", y = "Fetal-specific haploid genome equivalents per millilitre plasma (GE/ml)",
       title = "Fetal-specific cfDNA by gestation")

ggsave(filename = "plots/cffDNA_plot.png", plot = cffDNA_plot, dpi = 300)

#########################
# Collating SNP genotypes
#########################

snp_panel <- read_csv("W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/ddPCR SNP Panel/Camunas_Soler_panel_47_GnomAD_frequencies.csv")

dataPath <- "data/ddPCR_SNP_genotyping/"

ddpcr_files <- list.files(dataPath)

snp_panel_24 <- c(unique(snp_panel %>%
                           # Selct only the first 24 SNPs
                           filter(GOSH_ID_ds < 25) %>%
                           select(dbSNP)))
#Empty data frame
SNP_data <- data.frame()

# Read and collate each worksheet csv
for (dataFile in ddpcr_files){
  tmp_dat <- read_csv(paste0(dataPath,dataFile), col_names = TRUE)
  SNP_data <-rbind(SNP_data, tmp_dat)
  rm(tmp_dat)
}

view(ddpcr_nipd_unblinded %>%
  select(r_number, variant_assay, tubes_plasma_current))

write.csv(ddpcr_nipd_unblinded, "analysis_outputs/ddpcr_nipd_unblinded.csv", row.names = FALSE)

nrow(scd)
(17/85)*100


#########################
# Compare pre-amplification and non-pre-amplification results
#########################

get_minor_fractions <- function(dataframe){
  new_dataframe <- dataframe %>%
    dplyr::rename(r_number = Sample,
                  fraction_a = FractionalAbundance,
                  fractionmax_a = PoissonFractionalAbundanceMax,
                  fractionmin_a = PoissonFractionalAbundanceMin) %>%
    
    mutate(fraction_b = 100 - fraction_a,
           fractionmax_b = 100 - fractionmin_a,
           fractionmin_b = 100 - fractionmax_a,
           
           minor_fraction = pmin(fraction_a, fraction_b),
           minor_fractionmax = pmin(fractionmax_a, fractionmax_b),
           minor_fractionmin = pmin(fractionmin_a, fractionmin_b)) %>%
    
    left_join(ddpcr_target_panel %>%
                select(Target, Assay),
              by = "Target") %>%
    
    select(c(r_number, Target, Assay, minor_fraction, minor_fractionmax, minor_fractionmin))
  
  return(new_dataframe)
}

pre_amplification_fractions <- get_minor_fractions(SNP_data %>%
                                                     filter(Sample != "NTC" & TargetType == "Ch1Unknown")) %>%
  # Rename the columns to allow graph plotting
  dplyr::rename(preamp_minor_fraction = minor_fraction,
                preamp_minor_fractionmax = minor_fractionmax,
                preamp_minor_fractionmin = minor_fractionmin)

non_preamp_fractions <- get_minor_fractions(ddpcr_data_merged_samples %>%
                                              filter(Sample %in% phase3_samples & TargetType == "Ch1Unknown"))


ff_comparison <- inner_join(non_preamp_fractions, pre_amplification_fractions, 
           by = c("r_number", "Assay"))

# Plot graph for paper supplementary information.
preamp_comparison_plot <- ggplot(ff_comparison, aes(x = preamp_minor_fraction, y = minor_fraction))+
  geom_point(size = 3)+
  geom_errorbarh(aes(xmin = preamp_minor_fractionmin, xmax = preamp_minor_fractionmax,
                     height = 0.1), alpha = 0.5)+
  geom_errorbar(aes(ymin = minor_fractionmin, ymax = minor_fractionmax,
                    width = 0.1), alpha = 0.5)+
  labs(x = "Fetal-specific fraction - pre-amplified cfDNA (%)",
       y = "Fetal-specific fraction - original cfDNA (%)",
       title = "Fetal-specific fraction in cfDNA measured by ddPCR with and without pre-amplification")+
  geom_abline(linetype = "dashed")+
  theme_bw()+
  xlim(0, 8)+
  ylim(0, 8)+
  no_gridlines+
  # Add on the Pearson correlation coefficient
  annotate(geom="text", x=7, y=1, label= paste0("r = ", round(cor(ff_comparison$minor_fraction, 
                                                                  ff_comparison$preamp_minor_fraction, 
                                                                  method = "pearson"), digits = 2), size = 6))

ggsave(filename = "plots/preamp_comparison_plot.png", plot =  preamp_comparison_plot, dpi = 300)

#########################
# Sickle cell disease limit of detection study
#########################

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
  mutate(molecules_difference = pmax(HbS_molecules, HbA_molecules) - pmin(HbS_molecules, HbA_molecules)) %>%
  mutate(total_DNA_molecules = HbS_molecules + HbA_molecules) %>%
  mutate(Reference_fraction = HbA_molecules / total_DNA_molecules) %>%
  mutate(Variant_fraction = HbS_molecules / total_DNA_molecules) %>%
  mutate(Variant_fraction_percent = round(Variant_fraction*100, 1)) %>%
  mutate(Over_represented_fraction = pmax(Reference_fraction, Variant_fraction)) %>%
  mutate(Likelihood_ratio = calc_LR_autosomal(fetal_fraction, Over_represented_fraction, total_DNA_molecules)) %>%
  mutate(SPRT_prediction = case_when(
    Likelihood_ratio > 250 & Over_represented_fraction == Reference_fraction ~ "homozygous reference",
    Likelihood_ratio > 250 & Over_represented_fraction == Variant_fraction ~ "homozygous variant",
    Likelihood_ratio < (1/250) ~ "heterozygous",
    Likelihood_ratio < 250 & Likelihood_ratio > (1/250) ~ "no call"))

view(LOD_data_longer %>%
       select(unique_identifier, molecules_difference))

LOD_data_even_longer <- LOD_data_longer %>%
  select(Sample, Mass_molecules, Input_molecules, HbS_molecules, HbA_molecules, AcceptedDroplets_Variant_assay, Variant_fraction_percent) %>%
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
    Sample %in% c("0") ~"0%")) %>%
  mutate(label_y = ifelse(Target == "HbS_molecules", NA, paste0(Variant_fraction_percent, "%")))

view(LOD_data_longer)
         
# Plot the results for paper
lod_plot <- ggplot(LOD_data_even_longer %>%
                     filter(Input_molecules %in% c(6000)), aes(x = Sample, y = molecules, fill = Target))+
  geom_col(width = 0.5, position = position_dodge(width =0.5), colour = "black", alpha = 0.6)+
  geom_errorbar(aes(ymin = Molecules_min, ymax = Molecules_max), width = .2, position=position_dodge(width=0.5))+
  scale_fill_manual(values=c("#3366FF", "#FF0000"), labels= c("Reference", "Variant"))+
  theme_bw()+
  theme(axis.text=element_text(size=10), axis.title = element_text(size=14), 
        plot.title = element_text(size=20), legend.position = "bottom", legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "", y = "Molecules detected by ddPCR")+
  ylim(0, 5000)+
  scale_x_discrete(labels=c("12%","10%","8%","6%","4%","2%","0%",
                            "2%","4%","6%","8%","10%","12%"))+
  geom_text(aes(y = 4500, label = label_y), fontface = "bold")
  

scd_ddpcr %>%
  filter(Fetal_fraction_percent < 4 & algorithm_prediction != "inconclusive") %>%
  select(r_number, algorithm_prediction, Fetal_fraction_percent, Molecules_variant_assay)




ggsave(filename = "plots/lod_plot.png", plot =  lod_plot, dpi = 300)

LOD_data_longer %>%
  filter(Input_molecules %in% c(6000)) %>%
  select(unique_identifier, Variant_fraction)
  
           
         
# Say why it is ok for us to go below 4% fetal fraction        
         
         
         

# A different plot - not as nice to look at.
ggplot(LOD_data_longer, aes(x = total_DNA_molecules, y = FractionalAbundance_HbS, colour = Sample))+
  geom_point(size = 4)+
  geom_pointrange(aes(ymin = PoissonFractionalAbundanceMin_HbS, ymax = PoissonFractionalAbundanceMax_HbS))+
  geom_hline(yintercept=50, linetype="dashed", size = 1)

#############################################################
# ROC curve analysis
############################################################

# This part is for plotting ROC curves for the sickle cell 
# disease data.

mcmc_vs_sprt_scd <- ddpcr_nipd_unblinded %>%
  filter(variant_assay == "HBB c.20A>T" & !r_number %in% samples_to_exclude) %>%
  # Convert invasive results to binary outcomes
  mutate(unbalanced = case_when(
    mutation_genetic_info_fetus %in% c("HbSS", "HbAA") ~"TRUE",
    mutation_genetic_info_fetus == "HbAS" ~"FALSE"),
    
    # Convert "unbalanced" column to Boolean vector
    unbalanced = as.logical(unbalanced),
    # Convert the MCMC calls to a binary outcome
    mcmc_hom_call = pmax(p_G1, p_G3)) %>%
  select(r_number, Likelihood_ratio, mcmc_hom_call, unbalanced) %>%
  filter(!is.na(unbalanced))

sprt_roc <- mcmc_vs_sprt_scd %>%
  arrange(desc(Likelihood_ratio)) %>%
  mutate(sprt_true_positive_rate = cumsum(unbalanced)/sum(unbalanced)) %>%
  mutate(sprt_false_positive_rate = cumsum(!unbalanced)/sum(!unbalanced)) %>%
  select(r_number, Likelihood_ratio, sprt_true_positive_rate, sprt_false_positive_rate, unbalanced)

mcmc_roc <- mcmc_vs_sprt_scd %>%
  arrange(desc(mcmc_hom_call)) %>%
  mutate(mcmc_true_positive_rate = cumsum(unbalanced)/sum(unbalanced)) %>%
  mutate(mcmc_false_positive_rate = cumsum(!unbalanced)/sum(!unbalanced)) %>%
  select(r_number, mcmc_hom_call, mcmc_true_positive_rate, mcmc_false_positive_rate)

total_roc <- sprt_roc %>%
  left_join(mcmc_roc,
            by = "r_number") %>%
  select(r_number, Likelihood_ratio, mcmc_hom_call, sprt_true_positive_rate,
         mcmc_true_positive_rate, sprt_false_positive_rate, mcmc_false_positive_rate, unbalanced)


ggplot(total_roc, aes(x = unbalanced, y = Likelihood_ratio))+
  geom_jitter()+
  scale_y_log10()+
  geom_hline(yintercept = 250, linetype = "dashed")

ggplot(total_roc, aes(x = unbalanced, y = mcmc_hom_call))+
  geom_jitter()

write.csv(sprt_roc, "analysis_outputs/sprt_roc.csv", row.names = FALSE)

write.csv(mcmc_roc, "analysis_outputs/mcmc_roc.csv", row.names = FALSE)


ggplot(total_roc, aes(x = false_positive_rate, y = true_positive_rate))+
  geom_line(size = 2, aes(colour = analysis_type))+
  theme_bw()+
  facet_wrap(~analysis_type)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(linetype = "dashed")+
  ylim(0, 1)+
  xlim(0,1)+
  labs(x = "False positive rate", y = "True positive rate", title = "Sickle cell disease ddPCR cohort ROC curve")

view(total_roc)

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
# Sickle cell gDNA analysis (out of date)
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
# cfDNA genomic equivalents per ml plasma (GE/ml) (out of date)
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



# Plotting graphs of sample cohorts
#############################################################
# Sample cohort plots
#############################################################


# Plot the results of fetal fraction versus variant fraction, including inconclusive results.
ggplot(ddpcr_nipd_unblinded %>%
         filter(variant_assay == "HBB c.20A>T"), aes(x = Fetal_fraction_percent, y = Variant_fraction_percent))+
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
ggplot(ddpcr_nipd_unblinded %>%
         filter(variant_assay == "HBB c.20A>T"), aes(x = Fetal_fraction_percent, y = Variant_fraction_percent))+
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
  ylim(40, 60)+
  labs(x = "Fetal fraction (%)", y = "Variant fraction (%)")


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

write.csv(sickle_cell_unblinded, 
          file = paste0("analysis_outputs/sickle_cell_unblinded", 
                        format(current_time, "%Y%m%d_%H%M%S"), ".csv"),
          row.names = FALSE)

write.csv(ddpcr_nipd_unblinded %>%
            filter(variant_assay != "HBB c.20A>T"), 
          file = paste0("analysis_outputs/bespoke_cohort_unblinded", 
                        format(current_time, "%Y%m%d_%H%M%S"), ".csv"),
          row.names = FALSE)

write.csv(mcmc_vs_sprt_outcomes, 
          file = paste0("analysis_outputs/mcmc_vs_sprt_outcomes", 
                        format(current_time, "%Y%m%d_%H%M%S"), ".csv"),
          row.names = FALSE)


## Additional code (non essential)

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






