################################################################################
## ddPCR Analysis for All Samples
## joseph.shaw3@nhs.net
## This analysis script performs 3 different analyses on the entire cohort
## of samples tested using ddPCR: the sequential probability ratio test (SPRT)
## (Lo et al, 2008; PMID:  17664418), MonteCarlo Markov Chain 
## (MCMC) analysis (Caswell et al, 2020; PMID:  32533152) and z score
## analysis.
################################################################################

# DISCLAIMER: This software is intended for research purposes only, and is not 
# validated for clinical use. No guarantee is provided for third party usage.

# Clear environment
rm(list=ls())

##################################################
# Load packages, functions and data
##################################################

## Packages
library(tidyverse)

# janitor for cleaning header names
library(janitor)

# stringr for counting string patterns
library(stringr)

# epiR for sensitivity calculations
library(epiR)

# ggpubr for plot collation
library(ggpubr)

# cmdstanr required for running stan models
library(cmdstanr)

# pROC for plotting ROC curves
library(pROC)

# Working directory
setwd("W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/ddPCR_R_Analysis/ddpcr_nipd")

# Source functions
source("functions/ddPCR_nipd_functions.R")

# Sample outcome information (generated on 23/10/2022 from the RAPID
# biobank, with study ID and maternal date of birth removed.)
sample_details <- read.csv("resources/sample_details.csv") %>%
  mutate(r_number = as.character(r_number))

# Load ddPCR SNP panel (from Camunas Soler et al 2018)
ddpcr_snp_panel <- read.csv("resources/camunas_soler_snp_panel.csv")

# Assay gene information
gene_info <- read.csv("resources/vf_assay_gene_information.csv")

# Load sample wells for drawing RMD plots
sample_wells <- read.csv("resources/sample_wells.csv")

##################################################
# LOAD DATA
##################################################
#########################
# Load resources
#########################

# Control gDNA information
controls <- readr::read_csv("resources/controls.csv")

# Target panel
ddpcr_target_panel <- readr::read_csv("resources/ddpcr_target_panel.csv") 

# Number and type of plasma extractions
plasma_extractions <- read.csv("resources/cfdna_plasma_extraction_volumes.csv") %>%
  group_by(r_number) %>%
  summarise(plasma_volume_ml = (sum(tubes_removed))*2) 

plasma_replicates <- read.csv("resources/cfdna_plasma_extraction_volumes.csv") %>%
  group_by(r_number) %>%
  summarise(extraction_replicates = n())

# Type of invasive sampling
invasive_sampling <- read.csv("resources/confirmation_testing.csv")

#########################
# Read in ddPCR data 
#########################

# Load in all the csv files exported from QuantaSoft.

dataPath <- "data/ddpcr_data/"

ddpcr_files <- list.files(dataPath)

#Empty data frame
ddpcr_data <- data.frame()

# Read and collate each worksheet csv
for (dataFile in ddpcr_files){
  tmp_dat <- readr::read_csv(paste0(dataPath,dataFile), col_names = TRUE,
                             show_col_types = FALSE)
  # Add a worksheet identifier and a unique identifier for each
  # well on each worksheet, and for each sample
  tmp_dat_ws <- tmp_dat %>%
    janitor::clean_names() %>%
    mutate(worksheet = dataFile,
           worksheet_well = paste(worksheet, well, sep = "_"),
           worksheet_well_sample = 
             paste(worksheet, well, sample, sep = "_")) %>%
    dplyr::rename(droplets = accepted_droplets)
  ddpcr_data <-rbind(ddpcr_data, tmp_dat_ws)
  rm(tmp_dat_ws)
  rm(tmp_dat)
}

# Check all files loaded in
stopifnot(length(setdiff(ddpcr_files, ddpcr_data$worksheet)) == 0)

ddpcr_data <- ddpcr_data %>%
  # Remove limit of detection data
  filter(worksheet != "20-1557.csv")

#########################
# Wrangle cfDNA data into shape
#########################

# The data wrangling in this section is due to collating 
# data acquired over 3 years with various lab workflows.

ddpcr_data_merged_samples <- ddpcr_data %>%
  # Remove single wells and controls
  filter(substr(well, 1, 1) == "M" & !(sample %in% controls$sample)) %>%
  # Count the number of wells which have been merged together,
  # which is the number of commas in "merged_wells" string plus one
  mutate(num_wells = stringr::str_count(merged_wells, ",")+1)

# Reshape the data frame to sum all values by target.
# Often cfDNA from the same sample was tested across
# multiple worksheets and then must be summed together
# for the analysis.
ddpcr_data_reshape <- ddpcr_data_merged_samples %>% 
  group_by(sample, target) %>% 
  summarise(positives = sum(positives),
            droplets = sum(droplets),
            num_wells = sum(num_wells),
            .groups="drop")

# Add on the category for each ddPCR target
# Data model: 4 categories - variant, reference, ff_allele1 and ff_allele2
ddpcr_with_target <- ddpcr_data_reshape %>% 
  left_join(ddpcr_target_panel %>%
              select(target, target_category), 
            by = "target")

# Use pivot_wider to get one row per sample.
# There should be data for 151 cfDNA samples in total, 24 do not have both 
# variant fraction and fetal fraction assays run, which leaves the 
# 127 included in the paper.
pivotted_ddpcr <- ddpcr_with_target %>% 
  pivot_wider(id_cols = sample,
              names_from = target_category,
              values_from = c(droplets, positives, num_wells),
              # Use names_glue to keep new columns names with naming
              # convention
              names_glue = "{target_category}_{.value}")

# Add on assay and inheritance pattern to the table.
# Do this first by adding the information to the "variant" rows, then to
# the "ff_allele1" rows.
ref_table_var <- left_join(
  # First table
  ddpcr_with_target %>%
    filter(target_category == "variant"),
  # Second table
  ddpcr_target_panel %>%
    select(assay, inheritance_chromosomal, inheritance_pattern, target) %>%
    # Rename the assay column to avoid clash with fetal fraction assay
    dplyr::rename(vf_assay = assay),
  # Join by
  by = "target")

ref_table_ff <- left_join(
  # First table
  ddpcr_with_target %>%
    filter(target_category == "ff_allele1"),
  # Second table
  ddpcr_target_panel %>%
    select(assay, target) %>%
    dplyr::rename(ff_assay = assay),
  # Join by
  by = "target")

cfdna_ddpcr_data <- pivotted_ddpcr %>%
  left_join(ref_table_var %>%
              select(sample, inheritance_chromosomal, inheritance_pattern, 
                     vf_assay), by = "sample",
            .groups="drop") %>% 
  left_join(ref_table_ff %>%
              select(sample, ff_assay), by = "sample",
            .groups="drop") %>%
  
  # Remove any samples which haven't had both assays performed.
  # This removes 24 samples which didn't have both assays performed
  # for various reasons.
  filter(!is.na(variant_positives) & !is.na(ff_allele1_positives)) %>%
  
  # Remove duplicate columns. Reference and variant targets have the same
  # number of droplets because they originate from the same ddPCR well.
  select(-c(ff_allele2_droplets, reference_droplets, 
            reference_num_wells,
            ff_allele2_num_wells)) %>%
  # Rename to be compatible with functions
  dplyr::rename(ff_assay_droplets = ff_allele1_droplets,
                vf_assay_droplets = variant_droplets,
                ff_assay_num_wells = ff_allele1_num_wells,
                vf_assay_num_wells = variant_num_wells) %>%
  
  # Determine the maternal and paternally-inherited alleles for the 
  # fetal fraction assay for the fetal fraction calculation.
  # This step assumes that the paternally-inherited allele will be 
  # lower (fetal cfDNA only) than the maternally-inherited allele (fetal and
  # maternal cfDNA)
  mutate(maternal_positives = pmax(ff_allele1_positives, ff_allele2_positives),
         paternal_positives = pmin(ff_allele1_positives, ff_allele2_positives))

#########################
# Perform Poisson calculations for cfDNA data
#########################

cfdna_ddpcr_data_molecules <- ff_calculations(
  var_ref_calculations(cfdna_ddpcr_data)) %>% 
  
  # Rename sample to r_number to allow merge sample outcome data.
  rename(r_number = sample) %>%
  
  # Add information about extractions
  left_join(plasma_extractions %>%
              mutate(r_number = as.character(r_number)), by = "r_number") %>%
  left_join(plasma_replicates %>%
              mutate(r_number = as.character(r_number)), by = "r_number") %>%
  left_join(invasive_sampling %>%
              mutate(r_number = as.character(r_number)), by = "r_number") %>%
  mutate(
    total_molecules = vf_assay_molecules + ff_assay_molecules,
    totalGE_ml_plasma = total_molecules / plasma_volume_ml)

#########################
# Wrangle gDNA data into shape
#########################

# Check that controls csv does not contain duplicate sample names, which
# will effect later pivot_wider steps
stopifnot(duplicated(controls$sample) == FALSE)

# Get single well and merged well control data without NTCs
ddpcr_controls <- ddpcr_data %>%
  filter(sample %in% controls$sample & sample != "NTC")

# Add on target category 
ddpcr_controls_with_target <- ddpcr_controls %>% 
  left_join(ddpcr_target_panel %>%
              select(target, target_category, assay), by = "target") %>%
  # Add on sample identity as "maternal", "paternal" or "generic"
  left_join(controls %>%
              select(sample,	identity),
            by = "sample")

# Use pivot_wider to get one row per well for 
# each well tested with a variant assay, and then calculate
# molecular counts using var_ref_calculations

parent_gDNA_var_ref <- var_ref_calculations(
  
  ddpcr_controls_with_target %>%
    filter(target_category %in% c("variant", "reference")) %>%
    select(worksheet_well_sample, sample, identity, assay, 
           target_category, droplets, positives) %>%
    pivot_wider(id_cols = c(worksheet_well_sample, sample, identity, assay),
                names_from = target_category,
                values_from = c(droplets, positives),
                names_glue = "{target_category}_{.value}") %>%
    select(-reference_droplets) %>%
    dplyr::rename(vf_assay_droplets = variant_droplets,
                  vf_assay = assay)
)

# Use pivot_wider to get one row per well for each well tested with 
# a fetal fraction assay, and then calculate
# molecular counts using ff_calculations

parent_gDNA_ff <- ff_calculations(
  
  ddpcr_controls_with_target %>%
    filter(target_category %in% c("ff_allele2", "ff_allele1") &
             identity %in% c("maternal gDNA", "paternal gDNA")) %>%
    select(worksheet_well_sample, sample, identity, assay,
           target_category, droplets, positives) %>%
    pivot_wider(id_cols = c(worksheet_well_sample, sample, identity, assay),
                names_from = target_category,
                values_from = c(droplets, positives),
                names_glue = "{target_category}_{.value}") %>%
    dplyr::rename(ff_assay_droplets = ff_allele1_droplets,
                  ff_assay = assay) %>%
    
    # Need to create the maternal_positives and paternal_positives columns.
    # "Maternal positives" refers to the positive droplets for the homozygous
    # maternal allele.
    mutate(
      maternal_positives = case_when(
        identity == "paternal gDNA" ~pmin(ff_allele1_positives, 
                                          ff_allele2_positives),
        identity == "maternal gDNA" ~ pmax(ff_allele1_positives, 
                                           ff_allele2_positives)),
      # "Paternal positives" refers to the positive droplets for the 
      # homozygous paternal allele.
      paternal_positives = case_when(
        identity == "paternal gDNA" ~pmax(ff_allele1_positives, 
                                          ff_allele2_positives),
        identity == "maternal gDNA" ~ pmin(ff_allele1_positives, 
                                           ff_allele2_positives))) %>%
    select(-c(ff_allele2_droplets, ff_allele2_positives, 
              ff_allele1_positives)))

##################################################
# ANALYSIS
##################################################
#########################
# Heterozygous gDNA cohort
#########################

# Upper limit of number of cfDNA molecules per case
cfdna_molecules_max <- max(cfdna_ddpcr_data_molecules$vf_assay_molecules)

# Lower limit for analysis (gDNA variation is high below this)
vf_assay_molecules_limit <- 2000

# Remove duplicates in ddPCR target panel
ddpcr_targets <- ddpcr_target_panel[!duplicated(ddpcr_target_panel$assay),]

het_controls <- controls %>%
  filter(identity %in% c("maternal gDNA", "paternal gDNA") &
           zygosity == "heterozygous")

# Select heterozygous gDNA data
het_gdna <- parent_gDNA_var_ref %>%
  filter(sample %in% het_controls$sample &
           # Remove empty well 
           reference_positives != 0) %>%
  # New column to distinguish control gDNA
  mutate(sample_type = "gDNA") %>%
  dplyr::rename(r_number = sample) %>%
  # Add inheritance patterns
  left_join(ddpcr_targets %>%
              dplyr::rename(vf_assay = assay) %>%
              select(vf_assay, inheritance_chromosomal, inheritance_pattern),
            by = "vf_assay")

het_gdna_range <- het_gdna %>%
  # Filter the gDNA dataset to a similar range to the cfDNA dataset
  filter(vf_assay_molecules < (cfdna_molecules_max +1000) &
           vf_assay_molecules > vf_assay_molecules_limit)

# vp is variant percent
het_gDNA_mean_vp <- mean(het_gdna_range$variant_percent)            

stopifnot(het_gDNA_mean_vp > 49 & het_gDNA_mean_vp < 51)

# sd is standard deviation
het_gDNA_sd_vp <- sd(het_gdna_range$variant_percent)

# cv is coefficient of variation
het_gDNA_cv_vp <- (het_gDNA_sd_vp/het_gDNA_mean_vp) * 100

#########################
# SPRT analysis - cfDNA
#########################

all_samples_sprt <- cfdna_ddpcr_data_molecules %>%
  mutate(
    
    # Perform SPRT and return the likelihood ratio
    likelihood_ratio = case_when(
      inheritance_chromosomal == "x_linked" ~ 
        calc_lr_x_linked(fetal_fraction, (major_allele_percent/100), 
                         vf_assay_molecules),
      inheritance_chromosomal == "autosomal" ~ 
        calc_lr_autosomal(fetal_fraction, (major_allele_percent/100), 
                          vf_assay_molecules)))

# Predict genotypes with a likelihood ratio threshold of 8
all_samples_sprt <- predict_sprt_genotypes(all_samples_sprt, 8)

#########################
# MCMC analysis - cfDNA
#########################

# Prepare ddPCR data for MCMC

ddpcr_data_mcmc <- cfdna_ddpcr_data %>%
  dplyr::rename(r_number = sample,
                n_K = vf_assay_droplets,
                K_M = variant_positives,
                K_N = reference_positives,
                n_Z = ff_assay_droplets,
                Z_X = maternal_positives,
                Z_Y = paternal_positives) %>%
  arrange(inheritance_chromosomal, inheritance_pattern, vf_assay) %>%
  select(r_number, inheritance_chromosomal, inheritance_pattern, 
         vf_assay, n_K, n_Z, K_N, K_M, n_Z, Z_X, Z_Y)

# 23/10/2022 - 17 minutes 19 seconds for 127 cfDNA samples
# 30/10/2022 - 21 minutes 20 seconds for 127 cfDNA samples
# 07/12/2022 - 20 minutes 4 seconds for 127 cfDNA samples
all_samples_mcmc <- run_mcmc(ddpcr_data_mcmc, 0.95)

#########################
# Z score analysis - cfDNA
#########################

zscore_imbalance_threshold <- 3
zscore_balance_threshold <- 2

all_samples_zscore <- cfdna_ddpcr_data_molecules %>%
  mutate(zscore = (variant_percent - het_gDNA_mean_vp) / het_gDNA_sd_vp,
         
         zscore_prediction = case_when(
           
           # Predict fetal genotypes for X-linked cases
           inheritance_chromosomal == "x_linked" &
             zscore > zscore_imbalance_threshold ~"hemizygous variant",
           
           inheritance_chromosomal == "x_linked" &
             zscore < -zscore_imbalance_threshold ~"hemizygous reference",
           
           inheritance_chromosomal == "x_linked" &
             zscore > -zscore_imbalance_threshold &
             zscore < zscore_imbalance_threshold ~"inconclusive",
           
           # Predict fetal genotypes for autosomal cases
           inheritance_chromosomal == "autosomal" &
             zscore > zscore_imbalance_threshold ~"homozygous variant",
           
           inheritance_chromosomal == "autosomal" &
             zscore < -zscore_imbalance_threshold ~"homozygous reference",
           
           inheritance_chromosomal == "autosomal" &
             zscore < zscore_balance_threshold &
             zscore > -zscore_balance_threshold ~"heterozygous",
           
           inheritance_chromosomal == "autosomal" &
             zscore < zscore_imbalance_threshold &
             zscore > zscore_balance_threshold ~"inconclusive",
           
           inheritance_chromosomal == "autosomal" &
             zscore > -zscore_imbalance_threshold &
             zscore < -zscore_balance_threshold ~ "inconclusive"))

#########################
# Binary predictions
#########################

# Convert SPRT predictions into a binary format, for assistance with
# sensitivity calculations later on
all_samples_sprt <- binary_predictions(df = all_samples_sprt, 
                                       prediction = quo(sprt_prediction)) %>%
  dplyr::rename(sprt_prediction_binary_pre_qc = binary_call)

all_samples_mcmc <- binary_predictions(
  df = all_samples_mcmc,
  prediction = quo(mcmc_prediction)) %>%
  dplyr::rename(mcmc_prediction_binary_pre_qc = binary_call)

all_samples_zscore <- binary_predictions(
  df = all_samples_zscore,
  prediction = quo(zscore_prediction)) %>%
  dplyr::rename(zscore_prediction_binary_pre_qc = binary_call)

#########################
# Collate cfDNA analyses and quality filtering
#########################

fetal_fraction_threshold <- 4

all_samples_blinded <- left_join(
  all_samples_sprt,
  all_samples_mcmc %>%
    select(r_number, p_G0, p_G1, p_G2, p_G3, mcmc_prediction,
           mcmc_prediction_binary_pre_qc),
  by = "r_number") %>%
  left_join(
    all_samples_zscore %>%
      select(r_number, zscore, zscore_prediction,
             zscore_prediction_binary_pre_qc),
    by = "r_number") %>%
  # Quality filters
  mutate(quality_filter = case_when(
    fetal_percent < fetal_fraction_threshold | 
      vf_assay_molecules < vf_assay_molecules_limit ~"fail",
    fetal_percent > fetal_fraction_threshold & 
      vf_assay_molecules > vf_assay_molecules_limit ~"pass")) %>%
  # Rename prediction variables 
  dplyr::rename(
    sprt_prediction_pre_qc = sprt_prediction,
    mcmc_prediction_pre_qc = mcmc_prediction,
    zscore_prediction_pre_qc = zscore_prediction) %>%
  mutate(
    # Amend predictions following quality control (qc)
    sprt_prediction = case_when(
      quality_filter == "fail" ~"inconclusive",
      quality_filter == "pass" ~sprt_prediction_pre_qc),
    mcmc_prediction  = case_when(
      quality_filter == "fail" ~"inconclusive",
      quality_filter == "pass" ~mcmc_prediction_pre_qc),
    zscore_prediction  = case_when(
      quality_filter == "fail" ~"inconclusive",
      quality_filter == "pass" ~zscore_prediction_pre_qc),
    # Amend binary predictions following qc
    sprt_prediction_binary = case_when(
      quality_filter == "fail" ~"inconclusive",
      quality_filter == "pass" ~sprt_prediction_binary_pre_qc),
    mcmc_prediction_binary  = case_when(
      quality_filter == "fail" ~"inconclusive",
      quality_filter == "pass" ~mcmc_prediction_binary_pre_qc),
    zscore_prediction_binary  = case_when(
      quality_filter == "fail" ~"inconclusive",
      quality_filter == "pass" ~zscore_prediction_binary_pre_qc),
    
    cohort = ifelse(vf_assay == "HBB c.20A>T", "sickle cell disease",
                       "bespoke design"),
       cohort = factor(cohort, levels = c("bespoke design",
                                          "sickle cell disease")),
    
    # Add on abbreviation of inheritance patterns
    inheritance_abbreviation = case_when(
      inheritance_chromosomal == "autosomal" &
        inheritance_pattern == "recessive" ~"AR",
      inheritance_chromosomal == "autosomal" &
        inheritance_pattern == "dominant" ~"AD",
      inheritance_chromosomal == "x_linked" &
        inheritance_pattern == "dominant" ~"XLD",
      inheritance_chromosomal == "x_linked" &
        inheritance_pattern == "recessive" ~"XLR"),
      
      # Factorise for ordering
      inheritance_abbreviation = factor(inheritance_abbreviation,
                                        levels = c("AR", "AD", 
                                                   "XLR", "XLD"))) %>%
  # Bind to gene information
  left_join(gene_info,
            by = "vf_assay") %>%
  # Order by cohort and inheritance pattern
  arrange(cohort, inheritance_abbreviation, vf_assay)

export_timestamp(all_samples_blinded)

#########################
# Compare analysis predictions against invasive outcomes
#########################

# The "all_samples_unblinded" table is the central table for creating plots
# and tables for the paper. It should contain all the information needed
# about each cfDNA case.

# Samples to exclude
# 13262 - contamination
# 17004 - father was HbAC
# 20915 - twin pregnancy
samples_to_exclude <- c("13262", "20915", "17004")

all_samples_unblinded <- all_samples_blinded %>%
  # Add cfDNA number for referencing within the paper
  mutate(sample_id = 
           paste0("cfDNA-", as.character(row.names(all_samples_blinded))),
         
         # Specify the method used to determine the fetal fraction
         ff_determination = case_when(
           inheritance_chromosomal == "x_linked" ~"ZFXY",
           inheritance_chromosomal == "autosomal" & ff_assay %in%
             ddpcr_snp_panel$dbSNP ~ "ddPCR SNP panel - workflow 2",
           TRUE ~"NGS SNP panel - workflow 1")) %>%
  filter(!r_number %in% samples_to_exclude) %>%
  left_join(sample_details, 
          by = "r_number") %>%
  # Change the nomenclature of the outcomes for sickle cell disease
  mutate(fetal_genotype = case_when(
    mutation_genetic_info_fetus == "HbSS" ~"homozygous variant",
    mutation_genetic_info_fetus == "HbAS" ~"heterozygous",
    mutation_genetic_info_fetus == "HbAA" ~"homozygous reference",
    TRUE ~mutation_genetic_info_fetus),
    
    # Compare predictions to invasive testing
    sprt_outcome = case_when(
      sprt_prediction == "inconclusive" ~"inconclusive",
      sprt_prediction == fetal_genotype ~"correct",
      TRUE ~"incorrect"),
    
    mcmc_outcome = case_when(
      mcmc_prediction == "inconclusive" ~"inconclusive",
      mcmc_prediction == fetal_genotype ~"correct",
      TRUE ~"incorrect"),
    
    zscore_outcome = case_when(
      zscore_prediction == "inconclusive" ~"inconclusive",
      zscore_prediction == fetal_genotype ~"correct",
      TRUE ~"incorrect")) %>%
  select(-mutation_genetic_info_fetus) %>%
  
  mutate(
    # Factorise for plotting
         fetal_genotype = factor(fetal_genotype,
                                     levels = c("hemizygous variant",
                                                "homozygous variant",
                                                "heterozygous",
                                                "homozygous reference",
                                                "hemizygous reference")),
         zscore_outcome = factor(zscore_outcome,
                                 levels = c("correct", "incorrect", 
                                            "inconclusive")),
         mcmc_outcome = factor(mcmc_outcome,
                               levels = c("correct", "incorrect", 
                                          "inconclusive")),
         sprt_outcome = factor(sprt_outcome,
                               levels = c("correct", "incorrect", 
                                          "inconclusive")),
         
         sprt_prediction = factor(sprt_prediction, levels = 
                                    c("hemizygous variant",
                                      "homozygous variant",
                                      "heterozygous",
                                      "homozygous reference",
                                      "hemizygous reference",
                                      "inconclusive")),
         
         mcmc_prediction = factor(mcmc_prediction, levels = 
                                    c("hemizygous variant",
                                      "homozygous variant",
                                      "heterozygous",
                                      "homozygous reference",
                                      "hemizygous reference",
                                      "inconclusive")),
         
         zscore_prediction = factor(zscore_prediction, levels = 
                                    c("hemizygous variant",
                                      "homozygous variant",
                                      "heterozygous",
                                      "homozygous reference",
                                      "hemizygous reference",
                                      "inconclusive"))) %>%
  arrange(family_number)

# Families/patients
# There should be 6 patients ("families") with more than one sample present
# in the dataset, with the following R numbers:
# 17841	and 18735	
# 20764	and 14992	
# 20788	and 14182	
# 20813, 14743	and 30112	
# 30173	and 14639	
# 30230	and 20238	

export_timestamp(all_samples_unblinded)

##################################################
# TABLES
##################################################
#########################
# Summary of all cfDNA results: Supplementary Data
#########################

# The supplementary table is a streamlined version of the 
# "all_samples_unblinded" table, with some columns renamed and others removed
# for brevity.

supplementary_table <- all_samples_unblinded %>%
  select(
    # Case information (R number removed)
    sample_id, family_number, cohort, inheritance_abbreviation, condition, 
    gestation_character, partner_sample_available,
    # Assay information
    vf_assay, ff_determination, ff_assay,
    variant_percent, fetal_percent, vf_assay_molecules,
    # Quality control
    quality_filter, 
    # Predictions before quality control
    sprt_prediction_pre_qc, mcmc_prediction_pre_qc,
    zscore_prediction_pre_qc,
    # Predictions after quality control
    sprt_prediction, mcmc_prediction, zscore_prediction,
    # Outcomes
    fetal_genotype,
    sprt_outcome, mcmc_outcome, zscore_outcome,
    # Sampling information
    vacutainer, hours_to_first_spin, days_to_storage, diagnostic_sampling,
    plasma_volume_ml, extraction_replicates,
    # ddpcr information,
    vf_assay_num_wells, ff_assay_num_wells, variant_positives, 
    reference_positives, vf_assay_droplets, 
    maternal_positives, paternal_positives,
    ff_assay_droplets, variant_molecules,
    reference_molecules, maternal_molecules, paternal_molecules,
    # Confidence intervals
    variant_percent_max, variant_percent_min, 
    # Analysis information
    likelihood_ratio, p_G0, p_G1, p_G2, p_G3, zscore) %>%
  dplyr::rename(
    inheritance = inheritance_abbreviation,
    gestation = gestation_character,
    variant_fraction_assay = vf_assay,
    fetal_fraction_determination = ff_determination,
    fetal_fraction_assay =  ff_assay,
    `variant_fraction_%` = variant_percent,
    `fetal_fraction_%` = fetal_percent,
    `variant_fraction_%_max` = variant_percent_max, 
    `variant_fraction_%_min` = variant_percent_min, 
    variant_fraction_assay_molecules = vf_assay_molecules,
    variant_positive_droplets = variant_positives,
    reference_positive_droplets = reference_positives,
    maternal_marker_positive_droplets = maternal_positives,
    paternal_marker_positive_droplets = paternal_positives,
    maternal_marker_molecules = maternal_molecules,
    paternal_marker_molecules = paternal_molecules,
    sprt_likelihood_ratio = likelihood_ratio)

export_timestamp(supplementary_table)

#########################
# Bespoke cohort individual results table
#########################

bespoke_individual_results_table <- all_samples_unblinded %>%
  filter(vf_assay != "HBB c.20A>T") %>%
  select(inheritance_abbreviation, sample_id, condition,
         gene, variant_dna, fetal_genotype, sprt_prediction, 
         mcmc_prediction, zscore_prediction) %>%
  
  # Remove factor levels for renaming
  mutate(
    fetal_genotype = as.character(fetal_genotype),
    sprt_prediction = as.character(sprt_prediction),
    mcmc_prediction = as.character(mcmc_prediction),
    zscore_prediction = as.character(zscore_prediction)) %>%
  
  dplyr::rename(
    Inheritance = inheritance_abbreviation,
    `Sample number` = sample_id,
    Condition = condition,
    Gene = gene,
    DNA = variant_dna,
    `Fetal genotype` = fetal_genotype,
    SPRT = sprt_prediction,
    MCMC = mcmc_prediction,
    `Z score` = zscore_prediction)

# Rename with shortened names to allow easier presentation in the paper
bespoke_individual_results_table[bespoke_individual_results_table == "homozygous reference"] <- "hom ref"
bespoke_individual_results_table[bespoke_individual_results_table == "homozygous variant"] <- "hom var"
bespoke_individual_results_table[bespoke_individual_results_table == "heterozygous"] <- "het"
bespoke_individual_results_table[bespoke_individual_results_table == "hemizygous reference"] <- "hemi ref"
bespoke_individual_results_table[bespoke_individual_results_table == "hemizygous variant"] <- "hemi var"

export_timestamp(bespoke_individual_results_table)

#########################
# Incorrect results table
#########################

incorrect_results_table <- all_samples_unblinded %>%
  filter(sprt_outcome == "incorrect" |
           mcmc_outcome == "incorrect" |
           zscore_outcome == "incorrect") %>%
  mutate(
    fetal_percent = round(fetal_percent, 1),
    variant_percent = round(variant_percent, 1)) %>%
  select(sample_id, vf_assay, fetal_percent, variant_percent,
         vf_assay_molecules, fetal_genotype, sprt_prediction, 
         mcmc_prediction, zscore_prediction) %>%
  
  mutate(
    fetal_genotype = as.character(fetal_genotype),
    sprt_prediction = as.character(sprt_prediction),
    mcmc_prediction = as.character(mcmc_prediction),
    zscore_prediction = as.character(zscore_prediction)) %>%
  
  dplyr::rename(
    `Sample number` = sample_id,
    Variant = vf_assay,
    `Fetal fraction (%)` = fetal_percent,
    `Variant fraction (%)` = variant_percent,
    `Molecules detected` = vf_assay_molecules,
    `Fetal genotype` = fetal_genotype,
    SPRT = sprt_prediction,
    MCMC = mcmc_prediction,
    `Z score` = zscore_prediction)

incorrect_results_table[incorrect_results_table == "homozygous reference"] <- "hom ref"
incorrect_results_table[incorrect_results_table == "homozygous variant"] <- "hom var"
incorrect_results_table[incorrect_results_table == "heterozygous"] <- "het"
incorrect_results_table[incorrect_results_table == "hemizygous reference"] <- "hemi ref"
incorrect_results_table[incorrect_results_table == "hemizygous variant"] <- "hemi var"

export_timestamp(incorrect_results_table)

#########################
# Sickle cell results table
#########################

get_scd_table <- function(analysis_type) {
  
  genotype_count <- nrow(all_samples_unblinded %>%
                           filter(cohort == "sickle cell disease"))
  
  count_incon_excluded <- nrow(all_samples_unblinded %>%
                                 filter(cohort == "sickle cell disease"
                                        & !!analysis_type != "inconclusive"))
  
  table_1 <- all_samples_unblinded %>%
    filter(cohort == "sickle cell disease") %>%
    group_by(!!analysis_type, fetal_genotype) %>%
    summarise(total = n()) %>%
    pivot_wider(id_cols = !!analysis_type,
                names_from = fetal_genotype,
                values_from = total) %>%
    replace(is.na(.), 0) %>%
    mutate(total_count = `homozygous variant` + heterozygous +
             `homozygous reference`,
           percent = round((total_count/genotype_count)*100,0),
           percent_2 = round((total_count/count_incon_excluded)*100,0)) %>%
    dplyr::rename(HbSS = `homozygous variant`,
                  HbAS = heterozygous,
                  HbAA = `homozygous reference`)
  
  return(table_1)
  
}

scd_results_table <- cbind(get_scd_table(quo(sprt_outcome)), 
      get_scd_table(quo(mcmc_outcome)),
      get_scd_table(quo(zscore_outcome))) %>%
  select(-c(mcmc_outcome, zscore_outcome)) %>%
  dplyr::rename(Outcome = sprt_outcome)

export_timestamp(scd_results_table)

#########################
# Bespoke results summary table
#########################

get_bespoke_table <- function(analysis_type) {
  
  genotype_count <- nrow(all_samples_unblinded %>%
                           filter(cohort == "bespoke design"))
  
  count_incon_excluded <- nrow(all_samples_unblinded %>%
                                 filter(cohort == "bespoke design"
                                        & !!analysis_type != "inconclusive"))
  
  table_1 <- all_samples_unblinded %>%
    filter(cohort == "bespoke design") %>%
    group_by(!!analysis_type) %>%
    summarise(total = n()) %>%
    mutate(percent = round((total/genotype_count)*100,0),
           percent_2 = round((total/count_incon_excluded)*100,0))
  
  return(table_1)
  
}

bespoke_table <- cbind(get_bespoke_table(quo(sprt_outcome)),
                       get_bespoke_table(quo(mcmc_outcome)),
                       get_bespoke_table(quo(zscore_outcome)))

export_timestamp(bespoke_table)

#########################
# Whole cohort summary table
#########################

get_summary_table <- function(analysis_type) {
  
  genotype_count <- nrow(all_samples_unblinded)
  
  count_incon_excluded <- nrow(all_samples_unblinded %>%
                                 filter(!!analysis_type != "inconclusive"))
  
  table_1 <- all_samples_unblinded %>%
    group_by(!!analysis_type) %>%
    summarise(total = n()) %>%
    mutate(percent = round((total/genotype_count)*100,0),
           percent_2 = round((total/count_incon_excluded)*100,0))
  
}

summary_table <- cbind(get_summary_table(quo(sprt_outcome)),
                       get_summary_table(quo(mcmc_outcome)),
                       get_summary_table(quo(zscore_outcome)))

export_timestamp(summary_table)

##################################################
# PLOTS
##################################################
#########################
# HBB c.20A>T limit of detection study
#########################

lod_data <- read_csv("data/ddpcr_data/20-1557.csv", col_names = TRUE,
                     show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  dplyr::rename(droplets = accepted_droplets)

# There are 4 replicates per spike-in. Merging them together gives
# outputs of 3000, 6000, 9000 and 12000 genome equivalents.

# Wells for 3000 GE dataset
wells_3000_GE <- c("A02", "A03", "E03", "A04", "E04", "A05", "E05", 
                   "A06", "E06", "A07", "E07", "A08", "E08")

wells_6000_GE <- c(wells_3000_GE, "B02", "B03", "F03", "B04", "F04", 
                   "B05", "F05", "B06", "F06", "B07", "F07", "B08", "F08")

wells_9000_GE <- c(wells_6000_GE, "C02", "C03", "G03", "C04", "G04", 
                   "C05", "G05", "C06", "G06", "C07", "G07", "C08", "G08")

wells_12000_GE <- c(wells_9000_GE, "D02", "D03", "H03", "D04", "H04", 
                    "D05", "H05", "D06", "H06", "D07", "H07", "D08", "H08")

# Merge values together with a function

merge_lod_wells <- function(input_wells, GE_level){
  
  output <- lod_data %>% 
    select(well, sample, target, positives, droplets) %>%
    filter(well %in% input_wells) %>%
    group_by(sample, target) %>% 
    summarise(positives = sum(positives),
              droplets = sum(droplets),
              .groups="drop") %>%
    pivot_wider(id_cols = c(sample),
                names_from = target,
                values_from = c(droplets, positives),
                names_glue = "{target}_{.value}") %>%
    select(-c(HbA_droplets)) %>%
    dplyr::rename(vf_assay_droplets = HbS_droplets,
                  variant_positives = HbS_positives,
                  reference_positives = HbA_positives) %>%
    mutate(GE_level = GE_level)
  
  return(output)
}

# Bind the datasets together, then apply var_ref_calculations

lod_data_merged <- var_ref_calculations(rbind(
  merge_lod_wells(wells_3000_GE, 3000), 
  merge_lod_wells(wells_6000_GE, 6000),
  merge_lod_wells(wells_9000_GE, 9000),
  merge_lod_wells(wells_12000_GE, 12000))) %>%
  # Allows easier colour labeling
  mutate(sample = factor(sample, levels = 
                 c("SS 12%", "SS 10%", "SS 8%", "SS 6%", "SS 4%", "SS 2%",
                   "0%",
                   "AA 2%", "AA 4%", "AA 6%", "AA 8%", "AA 10%", "AA 12%")))

lod_plot_title <- expression(paste(italic("HBB"), " c.20A>T ddPCR limit of detection study"))

## Colour schemes

# In order of shade
# "#FFFFFF" = white
# "#CCCCCC" = grey 1
# "#9999CC" = grey 2
# "#999999" = grey 3
# "#666666" = grey 4
# "#333333" = grey 5
# "#000000" = black; 

black_and_white_scheme <- c(
  # SS 12% to 2%
  "#000000", "#333333", "#666666", "#999999", "#9999CC", "#CCCCCC",
  # 0%
  "#FFFFFF",
  # AA 2% to 12%
  "#CCCCCC", "#9999CC", "#999999", "#666666", "#333333","#000000")

red_and_blue_scheme <- c(
  # SS 12% to 2%
  "#660000", "#990000", "#CC0000", "#FF3333", "#FF6666", "#FF9999",
  # 0%
  "#FFFFFF",
  # AA 2% to 12%
  "#CCCCFF", "#9999FF", "#6666FF", "#3333CC", "#0000CC", "#000099")


# Limit of detection plot
lod_plot <- ggplot(lod_data_merged, 
                   aes(x = vf_assay_molecules, y = variant_percent)) +
  theme_bw() +
 theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.position = "right", 
    plot.title = element_text(size = 11),
    legend.title = element_blank(),
    legend.text = element_text(size = 9)) +
  geom_errorbar(aes(ymin = variant_percent_min, ymax = variant_percent_max),
                alpha = 0.2) +
  geom_errorbarh(aes(xmin = vf_assay_molecules_min, 
                     xmax = vf_assay_molecules_max),
                 alpha = 0.2) +
  scale_fill_manual(values = red_and_blue_scheme) +
  scale_shape_manual(values = c(24, 24, 24, 24, 24, 24, 21, 
                                25,25, 25, 25, 25, 25)) +
  geom_point(size = 1, aes(fill = sample, shape = sample)) +
  ylim(39, 61) +
  xlim(0, 15000) +
  geom_hline(yintercept = 50, linetype = "dashed") +
  labs(x = "Genome Equivalents (GE)",
       y = "Variant fraction (%)",
       title = lod_plot_title)

# 01/11/22: choose png as filetype instead of tiff (pngs have compressed file-sizes)
# A4 page is 21.6cm wide. 15 cm width for plots gives good results.

ggsave(plot = lod_plot, 
       filename = paste0("lod_plot_",
                         format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
       path = "plots/", 
       device='png', 
       dpi=1200,
       units = "cm",
       width = 15,
       height = 10)

#########################
# Plot themes
#########################

# Consistent theme and axes for plots
multiplot_theme <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), 
  legend.position = "none", 
  plot.title = element_text(size = 11),
  legend.title = element_blank(),
  legend.text = element_text(size = 9))

multiplot_y <- ylim(38, 62)

multiplot_x <- scale_x_continuous(limits = c(0,32000),
                                  breaks = c(0, 10000, 20000 ,30000))

## Z score lines

z3_line <- geom_hline(yintercept = het_gDNA_mean_vp+(3*het_gDNA_sd_vp),
                      linetype = "dashed", alpha = 0.5) 

zminus3_line <- geom_hline(yintercept = het_gDNA_mean_vp-(3*het_gDNA_sd_vp),
                           linetype = "dashed", alpha = 0.5)

z2_line <- geom_hline(yintercept = het_gDNA_mean_vp+(2*het_gDNA_sd_vp), 
                      linetype = "dashed", alpha = 0.5)

zminus2_line <- geom_hline(yintercept = het_gDNA_mean_vp-(2*het_gDNA_sd_vp),
                           linetype = "dashed", alpha = 0.5)

vertical_line <- geom_vline(xintercept = vf_assay_molecules_limit,
           linetype = "dashed", alpha = 0.5)

cfdna_fill <- scale_fill_manual(values=c("#FFFFFF", "#FF0000", "#999999"), 
                  guide = "none")

cfdna_alpha <- scale_alpha_manual(values = c(1, 1, 0.2), guide = "none")

cfdna_shape <- scale_shape_manual(values = c(24, 24, 21, 25, 25))

###################
# Plot cfDNA results
###################

samples_longer <- all_samples_unblinded %>%
  select(r_number, vf_assay_molecules, variant_percent, sprt_outcome, 
         mcmc_outcome, zscore_outcome,
         fetal_genotype) %>%
  pivot_longer(cols = c(sprt_outcome, mcmc_outcome, zscore_outcome),
               names_to = "analysis_type",
               values_to = "analysis_outcome") %>%
  mutate(fetal_genotype2 = case_when(
    fetal_genotype %in% c("hemizygous reference", "homozygous reference") ~ "hemi/homozygous reference",
    fetal_genotype %in% c("hemizygous variant", "homozygous variant") ~ "hemi/homozygous variant",
    fetal_genotype == "heterozygous" ~"heterozygous"),
    
    fetal_genotype2 = factor(fetal_genotype2, levels = c("hemi/homozygous variant",
                                                         "heterozygous",
                                                         "hemi/homozygous reference")),
    
    analysis_outcome = factor(analysis_outcome, levels = c("correct", "incorrect", 
                                                           "inconclusive")),
    analysis_type2 = case_when(
      analysis_type == "sprt_outcome" ~"SPRT analysis",
      # Change name based on Sarah's comment
      analysis_type == "mcmc_outcome" ~"Bayesian analysis",
      analysis_type == "zscore_outcome" ~"Z score analysis"),
    
    analysis_type2 = factor(analysis_type2, levels = c("SPRT analysis", 
                                                       "Bayesian analysis",
                                                       "Z score analysis")))

cfdna_plot_title <- paste0("ddPCR analysis for ", nrow(all_samples_unblinded), " cfDNA samples")

cfdna_plot <- ggplot(samples_longer, aes(x = vf_assay_molecules, 
                                         y = variant_percent)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.position = "bottom", 
    plot.title = element_text(size = 11),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9)) +
  ylim(38, 62) +
  scale_x_continuous(limits = c(0,38000),
                     breaks = c(0, 2000, 10000, 20000, 30000)) +
  scale_fill_manual(values=c("#FFFFFF", "#FF0000", "#999999"), 
                    guide = "none") +
  scale_alpha_manual(values = c(1, 1, 0.2), guide = "none") +
  scale_shape_manual(values = c(24, 21, 25)) +
  geom_point(size = 1, aes(fill = analysis_outcome,
                           alpha = analysis_outcome,
                           shape = fetal_genotype2),
             colour = "black") +
  labs(y = "Variant fraction (%)", x = "Genome equivalents (GE)",
       shape = "Fetal genotype (shape)",
       title = cfdna_plot_title) +
  geom_vline(aes(xintercept = 2000), linetype = "dashed") +
  facet_wrap(~analysis_type2, nrow = 3, ncol = 1) +
  
  geom_hline(data = subset(samples_longer, analysis_type2 == "Z score analysis"), 
             aes(yintercept = het_gDNA_mean_vp+(3*het_gDNA_sd_vp)), 
             linetype = "dashed") +
  
  geom_hline(data = subset(samples_longer, analysis_type2 == "Z score analysis"), 
             aes(yintercept = het_gDNA_mean_vp-(3*het_gDNA_sd_vp)),
             linetype = "dashed") +
  
  geom_hline(data = subset(samples_longer, analysis_type2 == "Z score analysis"), 
             aes(yintercept = het_gDNA_mean_vp+(2*het_gDNA_sd_vp)), 
             linetype = "dashed") +
  
  geom_hline(data = subset(samples_longer, analysis_type2 == "Z score analysis"),
             aes(yintercept = het_gDNA_mean_vp-(2*het_gDNA_sd_vp)), 
             linetype = "dashed") +
  
  geom_text(data = subset(samples_longer, analysis_type2 == "Z score analysis"),
            aes(x = 36000, y = 55), label = "z > 3") +
  
  geom_text(data = subset(samples_longer, analysis_type2 == "Z score analysis"),
            aes(x = 36000, y = 45), label = "z < -3") +
  
  geom_text(data = subset(samples_longer, analysis_type2 == "Z score analysis"),
            aes(x = 36000, y = 50), label = "2 > z > -2")

ggsave(plot = cfdna_plot, 
       filename = paste0("cfdna_plot_",
                         format(Sys.time(), "%Y%m%d_%H%M%S"),
                         ".png"),
       path = "plots/", 
       device='png', 
       dpi=600,
       units = "cm",
       width = 15,
       height = 18)

###################
# cfDNA analysis summary plot
###################

# Summary of how each analysis method performed.

analysis_summary <- rbind(
  all_samples_unblinded %>%
  count(zscore_outcome) %>%
  dplyr::rename(outcome = zscore_outcome) %>%
  mutate(analysis = "z score"),
  all_samples_unblinded %>%
    count(sprt_outcome) %>%
    dplyr::rename(outcome = sprt_outcome) %>%
    mutate(analysis = "sprt"),
  all_samples_unblinded %>%
    count(mcmc_outcome) %>%
    dplyr::rename(outcome = mcmc_outcome) %>%
    mutate(analysis = "mcmc")) %>%
  mutate(outcome = factor(outcome, 
                          levels = c("correct",
                                     "inconclusive",
                                     "incorrect")),
         analysis = factor(analysis,
                           levels = c("sprt", "mcmc", "z score")))

# Bar plot
analysis_summary_plot <- ggplot(analysis_summary, aes(x = outcome, y = n))+
  geom_col(fill = "white", colour = "black") +
  facet_wrap(~analysis) +
  ylim(0, 120) +
  labs(y = "Samples", x = "",
       title = "") +
  theme_bw() +
  geom_text(aes(x = outcome, y = n, label = n), vjust = -1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

###################
# ROC curve analysis for sickle cell disease
###################

# Plotting ROC curves can be achieved manually usually tidyverse functions,
# or using the pROC package.

# We need a binary classifier for the fetal genotype.

# It is best to restrict this analysis to sickle cell disease: this is the 
# largest cohort (88 samples) and trying to combine three different forms
# of inheritance with different dosage relationships into a binary 
# classification system doesn't work.

roc_binary_calls <- all_samples_unblinded %>%
  filter(vf_assay == "HBB c.20A>T") %>%
  # Convert invasive results to binary outcomes
  mutate(fetus_unbalanced = case_when(
      # As there are 3 possible genotypes, both homozygous variant and 
      # homozygous reference fetuses are coded as "affected" (TRUE), 
    # as they both involve an imbalance in cfDNA.
      fetal_genotype %in% c("homozygous variant", 
                            "homozygous reference")  ~"TRUE",
        fetal_genotype == "heterozygous" ~"FALSE"),
      # Convert "fetus_affected" column to Boolean vector
      fetus_unbalanced = as.logical(fetus_unbalanced),
      # SPRT likelihood score and z scores are already single predictors
      # We need to convert the MCMC calls to a single value
      mcmc_unbalanced_call = pmax(p_G1, p_G3),
      # Remove minus signs from z score
      zscore_unbalanced_call = abs(zscore))

# Create ROC objects for each analysis method
sprt_roc_object <- roc(
  # Response
  roc_binary_calls$fetus_unbalanced, 
  # Predictor
  roc_binary_calls$likelihood_ratio,
  direction="<",
  # Levels indicate "controls", then "cases"
  levels = c("FALSE", "TRUE"))

mcmc_roc_object <- roc(
  # Response
  roc_binary_calls$fetus_unbalanced, 
  # Predictor
  roc_binary_calls$mcmc_unbalanced_call,
  direction="<",
  levels = c("FALSE", "TRUE"))

zscore_roc_object <- roc(
  # Response
  roc_binary_calls$fetus_unbalanced, 
  # Predictor
  roc_binary_calls$zscore_unbalanced_call,
  direction="<",
  levels = c("FALSE", "TRUE"))

# SPRT plot
sprt_roc <- ggroc(sprt_roc_object, size = 1) +
  labs(y = "Sensitivity",
       x = "",
       title = "SPRT analysis") +
  theme_bw() +
  multiplot_theme +
  geom_abline(linetype = "dashed",
              intercept = 1) +
  annotate(geom = "text", x = 0.25, y = 0.3, 
           label = paste0("AUC = ", 
                    round(auc(roc_binary_calls$fetus_unbalanced, 
                              roc_binary_calls$likelihood_ratio),3)))

# MCMC plot
mcmc_roc <- ggroc(mcmc_roc_object, size = 1) +
  labs(y = "",
       x = "Specificity",
    title = "MCMC analysis") +
  theme_bw() +
  multiplot_theme +
  geom_abline(linetype = "dashed",
              intercept = 1)+
  annotate(geom = "text", x = 0.25, y = 0.3, 
           label = paste0("AUC = ", 
                          round(auc(roc_binary_calls$fetus_unbalanced, 
                                    roc_binary_calls$mcmc_unbalanced_call),3)))

# Z score plot
zscore_roc <- ggroc(zscore_roc_object, size = 1) +
  labs(x = "",
       y = "",
       title = "Z score analysis") +
  theme_bw() +
  multiplot_theme +
  geom_abline(linetype = "dashed",
              intercept = 1)+
  annotate(geom = "text", x = 0.25, y = 0.25, 
           label = paste0("AUC = ", 
                          round(auc(roc_binary_calls$fetus_unbalanced, 
                                    roc_binary_calls$zscore_unbalanced_call),3)))

# All ROC plots together
roc_plot <- ggpubr::ggarrange(sprt_roc, mcmc_roc, 
                              zscore_roc, 
                              ncol = 3, nrow = 1)

ggsave(plot = roc_plot, 
       filename = paste0("roc_plot_",
                         format(Sys.time(), "%Y%m%d_%H%M%S"),
                         ".tiff"),
       path = "plots/", device='tiff',
       units = "in",
       width = 10,
       height = 4)

###################
# Sensitivity and specificity
###################

sprt_scd <- sensitivity_metrics(df = all_samples_unblinded, 
                                prediction_binary = quo(sprt_prediction_binary), 
                                outcome = quo(sprt_outcome), 
                                cohort_input = "sickle cell disease",
                                cohort_name = "sickle cell disease") %>%
  dplyr::rename(sprt = analysis_method)

sprt_bespoke <- sensitivity_metrics(df = all_samples_unblinded, 
                                    prediction_binary = quo(sprt_prediction_binary), 
                                    outcome = quo(sprt_outcome), 
                                    cohort_input = "bespoke design",
                                    cohort_name = "bespoke design") %>%
  dplyr::rename(sprt = analysis_method)

sprt_all <- sensitivity_metrics(df = all_samples_unblinded, 
                                prediction_binary = quo(sprt_prediction_binary), 
                                outcome = quo(sprt_outcome), 
                                cohort_input = c("bespoke design", 
                                                 "sickle cell disease"),
                                cohort_name = "all") %>%
  dplyr::rename(sprt = analysis_method)

mcmc_scd <- sensitivity_metrics(df = all_samples_unblinded, 
                                prediction_binary = quo(mcmc_prediction_binary), 
                                outcome = quo(mcmc_outcome), 
                                cohort_input = "sickle cell disease",
                                cohort_name = "sickle cell disease") %>%
  dplyr::rename(mcmc = analysis_method)

mcmc_bespoke <- sensitivity_metrics(df = all_samples_unblinded, 
                                    prediction_binary = quo(mcmc_prediction_binary), 
                                    outcome = quo(mcmc_outcome), 
                                    cohort_input = "bespoke design",
                                    cohort_name = "bespoke design") %>%
  dplyr::rename(mcmc = analysis_method)

mcmc_all <- sensitivity_metrics(df = all_samples_unblinded, 
                                prediction_binary = quo(mcmc_prediction_binary), 
                                outcome = quo(mcmc_outcome), 
                                cohort_input = c("bespoke design", 
                                                 "sickle cell disease"),
                                cohort_name = "all") %>%
  dplyr::rename(mcmc = analysis_method)

zscore_scd <- sensitivity_metrics(df = all_samples_unblinded, 
                                  prediction_binary = quo(zscore_prediction_binary), 
                                  outcome = quo(zscore_outcome), 
                                  cohort_input = "sickle cell disease",
                                  cohort_name = "sickle cell disease") %>%
  dplyr::rename(zscore = analysis_method)

zscore_bespoke <- sensitivity_metrics(df = all_samples_unblinded, 
                                      prediction_binary = quo(zscore_prediction_binary), 
                                      outcome = quo(zscore_outcome), 
                                      cohort_input = "bespoke design",
                                      cohort_name = "bespoke design") %>%
  dplyr::rename(zscore = analysis_method)

zscore_all <- sensitivity_metrics(df = all_samples_unblinded, 
                                  prediction_binary = quo(zscore_prediction_binary), 
                                  outcome = quo(zscore_outcome), 
                                  cohort_input = c("bespoke design", 
                                                   "sickle cell disease"),
                                  cohort_name = "all") %>%
  dplyr::rename(zscore = analysis_method)


# Row  bind analysis methods together

sprt_metrics <- rbind(sprt_scd, sprt_bespoke, sprt_all)
mcmc_metrics <- rbind(mcmc_scd, mcmc_bespoke, mcmc_all)
zscore_metrics <- rbind(zscore_scd, zscore_bespoke, zscore_all)

analysis_metrics <- sprt_metrics %>%
  left_join(mcmc_metrics,
            by = c("group", "category")) %>%
  left_join(zscore_metrics,
            by = c("group", "category")) %>%
  # Tidy up names for paper
  dplyr::rename(
    `SPRT` = sprt,
    `MCMC` = mcmc,
    `Z score` = zscore) %>%
  mutate(`Result` = case_when(
    category == "true_positive" ~ "True positive",
    category == "true_negative" ~ "True negative",
    category == "false_positive" ~ "False positive",
    category == "false_negative" ~ "False negative",
    category == "sensitivity" ~ "Sensitivity (%)",
    category == "specificity" ~ "Specificity (%)",
    category == "inconclusive" ~ "Inconclusive"),
    `Cohort` = case_when(
      group == "sickle cell disease" ~ "Sickle cell disease",
      group == "bespoke design" ~ "Bespoke design",
      group == "all" ~"All samples")) %>%
  select(`Cohort`, `Result`, `SPRT`, `MCMC`, `Z score`)

spe_sen_table <- analysis_metrics %>%
  filter(Cohort == "All samples" &
           Result %in% c("Sensitivity (%)", "Specificity (%)")) %>%
  select(`Result`, `SPRT`, `MCMC`, `Z score`)

export_timestamp(spe_sen_table)

#########################
# SPRT analysis - gDNA
#########################

het_gdna_sprt <- het_gdna %>%
  # Remove replicates outside the DNA input range of cfDNA
  filter(vf_assay_molecules < cfdna_molecules_max) %>%
  # Assign an artificial fetal fraction of 4% for SPRT calculations
  mutate(fetal_fraction = 0.04,
         # Calculate the likelihood ratio
         likelihood_ratio = case_when(
           inheritance_chromosomal == "x_linked" ~ 
             calc_lr_x_linked(fetal_fraction, (major_allele_percent/100), 
                              vf_assay_molecules),
           inheritance_chromosomal == "autosomal" ~ 
             calc_lr_autosomal(fetal_fraction, (major_allele_percent/100), 
                               vf_assay_molecules)))

# Predict genotypes
het_gdna_sprt <- predict_sprt_genotypes(het_gdna_sprt, 8)

#########################
# MCMC analysis - gDNA
#########################

# Set an artificial fetal fraction of 4%
mcmc_fetal_fraction <- 0.04

gdna_data_mcmc <- het_gdna %>%
  # Remove replicates outside the DNA input range of cfDNA
  filter(vf_assay_molecules < cfdna_molecules_max) %>%
  mutate(
    var_ref_positives = variant_positives + reference_positives,
    var_ref_molecules = poisson_correct(vf_assay_droplets, var_ref_positives),
    
    # Create an artificial dataset for the fetal fraction.
    # Use the same number of total droplets as the variant fraction assay.
    ff_assay_droplets = vf_assay_droplets,
    
    #Calculate the number of fetal-specific and maternal molecules
    paternal_molecules = var_ref_molecules * mcmc_fetal_fraction,
    maternal_molecules = var_ref_molecules * (1-mcmc_fetal_fraction),
    
    # Calculate the expected number of positive partitions
    paternal_positives = reverse_poisson(paternal_molecules, 
                                         ff_assay_droplets),
    maternal_positives = reverse_poisson(maternal_molecules, 
                                         ff_assay_droplets)) %>%
  # Rename for pipeline
  dplyr::rename(n_K = vf_assay_droplets,
                K_M = variant_positives,
                K_N = reference_positives,
                Z_X = maternal_positives,
                Z_Y = paternal_positives,
                n_Z = ff_assay_droplets) %>%
  select(worksheet_well_sample, inheritance_chromosomal, inheritance_pattern,
         vf_assay, n_K, n_Z, K_N, K_M, n_Z, Z_X, Z_Y)

# 30/10/22: Analysing the 551 gDNA replicates takes ~1 hour 23 minutes
# 11/12/22: 1 hour 3 minutes.
gdna_mcmc_analysed <- run_mcmc(gdna_data_mcmc, 0.95)

# Export the file as a csv so you don't have to keep rerunning the analysis.
export_timestamp(gdna_mcmc_analysed)

###################
# SPRT and MCMC gDNA plots
###################
# Set fetal fraction and likelihood ratio thresholds
ff_for_graph <- 0.04
lr_for_graph <- 8

# SPRT function lines
het_upper_line <- geom_function(fun = "calc_het_upper_boundary",
              aes(x = vf_assay_molecules, y =
                    calc_het_upper_boundary(vf_assay_molecules, ff_for_graph, 
                                           lr_for_graph)),
              colour = "black",
              args = c(ff_for_graph, lr_for_graph),
              alpha = 0.5)
  
het_lower_line <- geom_function(fun = "calc_het_lower_boundary",
                aes(x = vf_assay_molecules, y =
                      calc_het_lower_boundary(vf_assay_molecules, ff_for_graph, 
                                             lr_for_graph)),
                colour = "black",
                args = c(ff_for_graph, lr_for_graph),
                alpha = 0.5)
  
hom_var_line <- geom_function(fun = "calc_hom_var_boundary",
                aes(x = vf_assay_molecules, y =
                      calc_hom_var_boundary(vf_assay_molecules, ff_for_graph, 
                                       lr_for_graph)),
                colour = "black",
                args = c(ff_for_graph, lr_for_graph),
                alpha = 0.5)
  
hom_ref_line <- geom_function(fun = "calc_hom_ref_boundary",
                aes(x = vf_assay_molecules, y =
                      calc_hom_ref_boundary(vf_assay_molecules, ff_for_graph, 
                                       lr_for_graph)),
                colour = "black",
                args = c(ff_for_graph, lr_for_graph),
                alpha = 0.5)

hemi_var_line <- geom_function(fun = "calc_hemi_var_boundary",
              aes(x = vf_assay_molecules, y =
                    calc_hemi_var_boundary(vf_assay_molecules, ff_for_graph, 
                                           lr_for_graph)),
              colour = "black",
              args = c(ff_for_graph, lr_for_graph),
              alpha = 0.5)
  
hemi_ref_line <- geom_function(fun = "calc_hemi_ref_boundary",
                aes(x = vf_assay_molecules, y =
                      calc_hemi_ref_boundary(vf_assay_molecules, ff_for_graph, 
                                             lr_for_graph)),
                colour = "black",
                args = c(ff_for_graph, lr_for_graph),
                alpha = 0.5)
  
het_region_label <- geom_text(aes(x = 28000, y = 50), size = 3,
                             label = "het")

hom_ref_region_label <- geom_text(aes(x = 28000, y = 47), size = 3, 
                          label = "hom ref")

hom_var_region_label <- geom_text(aes(x = 28000, y = 53), size = 3,
                          label = "hom var")

hemi_ref_region_label <- geom_text(aes(x = 28000, y = 47), size = 3,
                           label = "hemi ref")

hemi_var_region_label <- geom_text(aes(x = 28000, y = 53), size = 3, 
                           label = "hemi var")

gdna_plots_x <- xlim(0, 30000)
  
gdna_plots_y <- ylim(43, 57)

###########
# SPRT on gDNA for HBB c.20A>T
###########

# Factorise for plots
het_gdna_sprt_factored <- het_gdna_sprt %>%
  mutate(sprt_prediction = factor(sprt_prediction,
                                     levels = c("inconclusive",
                                                "heterozygous",
                                                "homozygous reference",
                                                "hemizygous reference",
                                                "homozygous variant",
                                                "hemizygous variant")),
         outcome = case_when(
           sprt_prediction == "inconclusive" ~"inconclusive",
           sprt_prediction %in% c("homozygous reference",
                                  "hemizygous reference",
                                  "homozygous variant",
                                  "hemizygous variant") ~"imbalance",
           sprt_prediction == "heterozygous" ~"balance"),
         outcome = factor(outcome, levels = c("balance",
                                              "inconclusive",
                                              "imbalance")))

gdna_plot_title <- expression(paste(italic("HBB"), 
                                   "c.20A>T gDNA: SPRT"))

plot_s5a <- ggplot(het_gdna_sprt_factored %>%
         filter(vf_assay == "HBB c.20A>T"), 
       aes(x = vf_assay_molecules, y = variant_percent)) +
  scale_fill_manual(values=c(
    # Balance - white
    "#FFFFFF",
    # Inconclusive - grey
    "#999999",
    # Imbalance - red
    "#FF0000")) +
  geom_point(size = 1,
             aes(fill = outcome), 
             pch=21,
             alpha = 0.8) +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.position = "bottom", 
    plot.title = element_text(size = 11),
    legend.title = element_blank(),
    legend.text = element_text(size = 9)) +
  het_upper_line +
  het_lower_line +
  hom_var_line +
  hom_ref_line +
  gdna_plots_x +
  gdna_plots_y +
  labs(x = "",
       y = "Variant fraction (%)",
       title = gdna_plot_title) +
  het_region_label +
  hom_ref_region_label +
  hom_var_region_label

###########
# SPRT on gDNA for autosomal assays
###########

plot_s5b <- ggplot(het_gdna_sprt_factored %>%
         filter(inheritance_chromosomal == "autosomal" &
                  vf_assay != "HBB c.20A>T"), 
       aes(x = vf_assay_molecules, y = variant_percent))+
  scale_fill_manual(values=c(
    # Balance - white
    "#FFFFFF",
    # Inconclusive - grey
    "#999999",
    # Imbalance - red
    "#FF0000")) +    
  geom_point(size = 1, 
             colour = "black", 
             aes(fill = outcome), 
                   pch=21,
                   alpha = 0.8) +
  theme_bw()+
  multiplot_theme +
  het_upper_line +
  het_lower_line +
  hom_var_line +
  hom_ref_line +
  gdna_plots_x +
  gdna_plots_y +
  labs(x = "",
       y = "Variant fraction (%)",
       title = "Autosomal variant gDNA: SPRT") +
  het_region_label +
  hom_ref_region_label +
  hom_var_region_label

###########
# SPRT on gDNA for X-linked assays
###########

plot_s5c <- ggplot(het_gdna_sprt_factored %>%
         filter(inheritance_chromosomal == "x_linked"), 
       aes(x = vf_assay_molecules, y = variant_percent))+
  scale_fill_manual(values=c(
    # Inconclusive - grey
    "#999999",
    # Imbalance - red
    "#FF0000")) +  
  geom_point(size = 1, 
             aes(fill = outcome),
             pch=21,
             alpha = 0.8) +
  theme_bw()+
  multiplot_theme +
  hemi_var_line +
  hemi_ref_line +
  gdna_plots_x +
  gdna_plots_y +
  labs(x = "Genome equivalents (GE)",
       y = "Variant fraction (%)",
       title = "X-linked variant gDNA: SPRT") +
  hemi_ref_region_label +
  hemi_var_region_label
 
###########
# MCMC on gDNA for HBB c.20A>T
###########

mcmc_gdna_for_plot <- gdna_mcmc_analysed %>%
  select(-c(inheritance_chromosomal, inheritance_pattern,
            vf_assay, fetal_fraction)) %>%
  left_join(het_gdna_sprt,
    by = "worksheet_well_sample") %>%
  # Factorise for plotting
  mutate(mcmc_prediction = factor(mcmc_prediction, levels = 
                             c("hemizygous variant",
                               "homozygous variant",
                               "heterozygous",
                               "homozygous reference",
                               "hemizygous reference",
                               "inconclusive")),
         outcome = case_when(
           mcmc_prediction == "inconclusive" ~"inconclusive",
           mcmc_prediction %in% c("homozygous reference",
                                  "hemizygous reference",
                                  "homozygous variant",
                                  "hemizygous variant") ~"imbalance",
           mcmc_prediction == "heterozygous" ~"balance"),
         outcome = factor(outcome, levels = c("balance",
                                              "inconclusive",
                                              "imbalance")))

gdna_plot_title <- expression(paste(italic("HBB"), 
                                   "c.20A>T gDNA: Bayesian"))

plot_s5d <- ggplot(mcmc_gdna_for_plot %>%
                   filter(vf_assay == "HBB c.20A>T"), 
                 aes(x = vf_assay_molecules, y = variant_percent)) +
  scale_fill_manual(values=c(
    # Balance - white
    "#FFFFFF",
    # Inconclusive - grey
    "#999999")) +    
  geom_point(size = 1,
             aes(fill = outcome), 
             pch=21,
             alpha = 0.8) +
  theme_bw()+
  multiplot_theme +
  gdna_plots_x +
  gdna_plots_y +
  labs(x = "",
       y = "",
       title = gdna_plot_title)

###########
# MCMC on gDNA for autosomal assays
###########

plot_s5e <- ggplot(mcmc_gdna_for_plot %>%
         filter(inheritance_chromosomal == "autosomal" &
                  vf_assay != "HBB c.20A>T"), 
       aes(x = vf_assay_molecules, y = variant_percent)) +
  scale_fill_manual(values=c(
    # Balance - white
    "#FFFFFF",
    # Inconclusive - grey
    "#999999")) +    
  geom_point(size = 1,
             aes(fill = outcome), 
             pch=21,
             alpha = 0.8) +
  theme_bw()+
  multiplot_theme +
  gdna_plots_x +
  gdna_plots_y +
  labs(x = "",
       y = "",
       title = "Autosomal variant gDNA: Bayesian")

###########
# MCMC on gDNA for X-linked assays
###########

plot_s5f <- ggplot(mcmc_gdna_for_plot %>%
         filter(inheritance_chromosomal == "x_linked"), 
       aes(x = vf_assay_molecules, y = variant_percent))+
  scale_fill_manual(values=c(
    # Inconclusive - grey
    "#999999",
    # Imbalance - red
    "#FF0000")) +    
  geom_point(size = 1, 
             aes(fill = outcome),
             pch=21,
             alpha = 0.8) +
  theme_bw()+
  multiplot_theme +
  gdna_plots_x +
  gdna_plots_y +
  labs(x = "Genome equivalents (GE)",
       y = "",
       title = "X-linked variant gDNA: Bayesian")

###########
# Arrange plots together
###########

gdna_results_plot <- ggpubr::ggarrange(plot_s5a, plot_s5d,
                                       plot_s5b, plot_s5e, 
                                       plot_s5c, plot_s5f,
                                       common.legend = TRUE,
                                       legend = "bottom",
                                  ncol = 2, nrow = 3, align = "v",
                                  labels = c("A", "D", "B",
                                             "E", "C", "F"))

ggsave(plot = gdna_results_plot,
       filename = paste0("gdna_results_plot_",
              format(Sys.time(), "%Y%m%d_%H%M%S"),
              ".png"),
       path = "plots/", device='png', dpi=600,
       units = "cm",
       width = 15,
       height = 18)

#########################
# Relative mutation dosage plots for each case
#########################

rmd_plots <-list()

# Plot an rmd plot for each case
for (i in sample_wells$cfdna_sample) {
  
  case <- sample_wells %>%
    filter(cfdna_sample == i)
  
  new_plot <- draw_rmd_plot(i, 
                            c(case[, 2], case[, 3]),
                            c(case[, 4], case[, 5]))
  rmd_plots <- list(rmd_plots, new_plot)
  rm(new_plot)
}

# Export RMD plots for every sample as a single pdf
ggexport(plotlist = rmd_plots, filename = paste0(
  "plots/rmd_plots_all_samples_", 
  format(Sys.time(), "%Y%m%d_%H%M%S"),
  ".pdf"), width=15, height=8, res=300)

#########################
# Why are is there a higher inconclusive rate in the bespoke cohort?
#########################

bespoke_inconclusive_plot <- all_samples_unblinded %>%
  mutate(combined_sample_id = paste0(r_number, "_", sample_id)) %>%
  filter(cohort == "bespoke design") %>%
  ggplot(aes(x = reorder(combined_sample_id, desc(vf_assay_molecules)), y = vf_assay_molecules, 
             colour = quality_filter))+
  scale_shape_manual(values = c(15, 17)) +
  scale_colour_manual(values = c("#CC0000", "#0000FF")) +
  geom_point(size = 3, aes(shape = inheritance_chromosomal)) +
  geom_hline(yintercept = 2000, linetype = "dashed") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0,20000),
                      breaks = c(0, 2000, 10000, 15000, 20000)) +
  labs(x = "", y = "Molecules detected in variant fraction assay",
       title = "Bespoke design ddPCR NIPD cohort")
  
ggsave(plot = bespoke_inconclusive_plot, 
       filename = paste0("bespoke_inconclusive_plot",
                         format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
       path = "plots/", 
       device='png', 
       dpi=1200,
       units = "cm",
       width = 20,
       height = 15)

#########################