################################################################################
## Analysis of Caswell et al 2020 cohort
## September 2021
## Joseph.Shaw@gosh.nhs.uk
################################################################################

#########################
# Load libraries and resources
#########################

## Load necessary packages
library(tidyverse)
library(readxl)
library(janitor)

## Set working directory
setwd("W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/ddPCR_R_Analysis")

# Source functions (this includes loading ddPCR data)
source("ddpcr_nipd/functions/ddPCR_nipd_functions.R")

# Load Caswell et al supplementary info

caswell_s1 <- read_excel(path = "archived_files/data/Paper_data_comparison/Caswell_et_al_supp_data.xlsx",
                         sheet = "supp_table_1",
                         skip = 3, 
                         n_max = 57) %>%
  janitor::clean_names() %>%
  mutate(sample_id = paste(family_number, pregnancy_number, weeks_gestation, 
                           sep = "_"))

caswell_s2 <- read_excel(path = "archived_files/data/Paper_data_comparison/Caswell_et_al_supp_data.xlsx",
                         sheet = "supp_table_2",
                         skip = 2) %>%
  janitor::clean_names() %>%
  mutate(sample_id = paste(family_number, pregnancy_number, weeks_gestation, 
                           sep = "_")) %>%
  select(sample_id, probabilty_for_homozygosity,
         predicted_genotype, actual_genotype)

#########################
# Analyse using SPRT
#########################

# Join s1 and s2 together, and rename for SPRT calculations
caswell_supp <- caswell_s1 %>%
  left_join(caswell_s2, by = "sample_id") %>%
  dplyr::rename(
    vf_assay_droplets = total_droplets_analysed_6,
    variant_positives = positive_alt,
    reference_positives = positive_ref,
    ff_assay_droplets = total_droplets_analysed_17,
    paternal_positives = positive_pat,
    maternal_positives = positive_mat)

lr_threshold <- 8

caswell_supp_sprt <- ff_calculations(var_ref_calculations(caswell_supp)) %>%
  mutate(likelihood_ratio = calc_lr_autosomal(fetal_fraction,
                                         major_allele_percent/100,
                                         vf_assay_molecules),
    sprt_genotype_prediction = case_when(
    likelihood_ratio > lr_threshold &
      major_allele == "variant allele"
    ~"hom",
    likelihood_ratio > lr_threshold &
      major_allele == "reference allele"
    ~"hom",
    likelihood_ratio < 1/lr_threshold
    ~"het",
    TRUE ~"no call")) %>%
  select(sample_id, predicted_genotype, sprt_genotype_prediction,
         actual_genotype)

#########################