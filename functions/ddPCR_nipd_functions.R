################################################################################
## Functions for ddPCR NIPD
## March 2021
## Joseph.Shaw@gosh.nhs.uk
################################################################################

#########################
# Load libraries
#########################

# Load packages
library(tidyverse)

# Load resources
controls <- readr::read_csv("resources/controls.csv")
ddpcr_target_panel <- readr::read_csv("resources/ddpcr_target_panel.csv") 

#########################
# ddPCR functions
#########################

# Function for Poisson correction. See Barrett et al 2012 
# supplementary information (PMID: 22451622)

Poisson_correct <- function(N, P) {
  num_molecules <- as.integer(-log((N-P)/N)*N)
  return(num_molecules)}

# Calculate the fetal fraction from paternal allele copies.

calc_ff <- function(Maternal_copies, Paternal_copies) {
  fetal_fraction <- (Paternal_copies*2) / (Paternal_copies + Maternal_copies)
  return(fetal_fraction)
}

# Calculates the 95% Poisson confidence intervals.
# Poisson intervals are lamda +/- 1.96 (sqrt(lamda/N)), 
# when lambda is copies per droplet and N is total droplets.

Poisson_max <- function(Copies_per_droplet, AcceptedDroplets) {
  Cpd_poisson_max = (Copies_per_droplet + 
                       1.96*(sqrt(Copies_per_droplet / AcceptedDroplets))) * AcceptedDroplets
  return(Cpd_poisson_max)
}

Poisson_min <- function(Copies_per_droplet, AcceptedDroplets) {
  Cpd_poisson_min = (Copies_per_droplet - 
                       1.96*(sqrt(Copies_per_droplet / AcceptedDroplets))) * AcceptedDroplets
  return(Cpd_poisson_min)
}

# Calculates the 95% Poisson fractional abundances based on the 95% max 
# and min values for one allele, but not the other allele.

Poisson_fraction_max <- function(copies_allele_max, copies_allele_other) {
  Allele_fraction_max = copies_allele_max / 
    (copies_allele_max + copies_allele_other)
  return(Allele_fraction_max)
}

Poisson_fraction_min <- function(copies_allele_min, copies_allele_other) {
  Allele_fraction_min = copies_allele_min / 
    (copies_allele_min + copies_allele_other)
  return(Allele_fraction_min)
}

#########################
# SPRT functions
#########################

# Calculates the likelihood ratio for a ddPCR test with X-linked inheritance.

calc_LR_X_linked <- function(fetal_fraction, overrep_fraction, total_copies) {
  q0 <- (1 - fetal_fraction) / (2 - fetal_fraction)
  q1 <- 1 / (2 - fetal_fraction)
  Delta <- (1- q1)/(1-q0)
  Gamma <- ((q1 * (1-q0))/ (q0*(1-q1)))
  LR <- exp((((overrep_fraction*log(Gamma)) + log(Delta))*total_copies))
  return(LR)
}

# Calculates the likelihood ratio for a ddPCR test when the variant is 
# on an autosome (recessive or dominant).

calc_LR_autosomal <- function(fetal_fraction, overrep_fraction, total_copies) {
  q0 = 0.5
  # I modified the q1 expression to make it easier to use. 
  # Fetal fraction should be in the right format
  # I.e. 0.05 not 5
  q1 <- 0.5+(fetal_fraction/2)
  Delta <- (1- q1)/(1-q0)
  Gamma <- ((q1 * (1-q0))/ (q0*(1-q1)))
  LR <- exp((((overrep_fraction*log(Gamma)) + log(Delta))*total_copies))
  return(LR)
}

# These functions calculate the SPRT thresholds assuming 
# a fetal fraction of 4% and a likelihood ratio of 8.

# Variables for each function
q0 <- 0.5
q1 <- 0.5+(0.04/2)
Delta <- (1- q1)/(1-q0)
Gamma <- ((q1 * (1-q0))/ (q0*(1-q1)))

calc_SS_boundary <- function(total_copies) {
  SS_boundary <- ((log(8)/total_copies) - log(Delta))/log(Gamma)
  # Convert to a percentage for output
  return(SS_boundary*100)
}

calc_AS_upper_boundary <- function(total_copies) {
  AS_upper_boundary <- ((log(1/8)/total_copies) - log(Delta))/log(Gamma)
  # Convert to a percentage for output
  return(AS_upper_boundary*100)
}

calc_AS_lower_boundary <- function(total_copies) {
  AS_upper_boundary <- ((log(1/8)/total_copies) - log(Delta))/log(Gamma)
  AS_lower_boundary <- 0.5-(AS_upper_boundary-0.5)
  # Convert to a percentage for output
  return(AS_lower_boundary*100)
}

calc_AA_boundary <- function(total_copies) {
  SS_boundary <- ((log(8)/total_copies) - log(Delta))/log(Gamma)
  AA_boundary <- 0.5-(SS_boundary-0.5)
  # Convert to a percentage for output
  return(AA_boundary*100)
}

#########################
# Grouped functions
#########################

## Naming convention
# For consistency, variables are named as: target_category_qualifier
# Target examples: variant, reference, paternal, maternal, fetal, difference, 
#                   major_allele, minor_allele, ff_assay, vf_assay, total
# Category examples: molecules, fraction, percent, positives
# Qualifier examples: max, min

# "maternal" and "paternal" refer to distinct maternally-inherited and 
# paternally-inherited alleles in the fetal fraction assay
# I.e. "maternal_positives" refers to the number of positive droplets for
# the maternally-inherited allele 

# Example: fetal_percent_max

# These functions calculate the number of molecules and fractional
# abundance for each target, including 95% Poisson confidence intervals.
# This is to prevent massive duplication of code when performing these
# calculations for gDNA, cfDNA and limit of detection datasets.

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
  # maternal_positives, paternal_positives, AcceptedDroplets_FetalFrac
  
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
# Read in ddPCR data 
#########################

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
# Wrangle cfDNA data into shape
#########################

# The data wrangling in this section is due to collating 
# data acquired over 3 years with various lab workflows.

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
              select(Target, Target_category), 
            by = "Target")

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
    dplyr::rename(variant_assay = Assay),
  # Join by
  by = "Target")

ref_table_ff <- left_join(
  # First table
  ddpcr_with_target %>%
    filter(Target_category == "ff_allele1"),
  # Second table
  ddpcr_target_panel %>%
    select(Assay, Target) %>%
    dplyr::rename(ff_assay = Assay),
  # Join by
  by = "Target")

cfdna_ddpcr_data <- pivotted_ddpcr %>%
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
  dplyr::rename(AcceptedDroplets_FetalFrac = AcceptedDroplets_ff_allele1,
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

# Get single well and merged well control data without NTCs
ddpcr_controls <- ddpcr_data %>%
  filter(Sample %in% controls$Sample & Sample != "NTC")

# Add on target category 
ddpcr_controls_with_target <- ddpcr_controls %>% 
  left_join(ddpcr_target_panel %>%
              select(Target, Target_category, Assay), by = "Target") %>%
  # Add on sample identity as "maternal", "paternal" or "generic"
  left_join(controls %>%
              select(Sample,	identity),
            by = "Sample")

# Use pivot_wider to get one row per well for 
# each well tested with a variant assay, and then calculate
# molecular counts using var_ref_calculations

parent_gDNA_var_ref <- var_ref_calculations(
  ddpcr_controls_with_target %>%
  filter(Target_category %in% c("variant", "reference")) %>%
  select(Worksheet_well, Sample, identity, Assay, 
         Target_category, AcceptedDroplets, Positives) %>%
  pivot_wider(id_cols = c(Worksheet_well, Sample, identity, Assay),
              names_from = Target_category,
              values_from = c(AcceptedDroplets, Positives)) %>%
  select(-AcceptedDroplets_reference) %>%
  dplyr::rename(AcceptedDroplets_Variant_assay = AcceptedDroplets_variant,
                variant_assay = Assay))

# Use pivot_wider to get one row per well for each well tested with 
# a fetal fraction assay, and then calculate
# molecular counts using ff_calculations

parent_gDNA_ff <- ff_calculations(ddpcr_controls_with_target %>%
  filter(Target_category %in% c("ff_allele2", "ff_allele1") &
           identity %in% c("maternal gDNA", "paternal gDNA")) %>%
  select(Worksheet_well, Sample, identity, Assay,
         Target_category, AcceptedDroplets, Positives) %>%
  pivot_wider(id_cols = c(Worksheet_well, Sample, identity, Assay),
              names_from = Target_category,
              values_from = c(AcceptedDroplets, Positives)) %>%
  dplyr::rename(AcceptedDroplets_FetalFrac = AcceptedDroplets_ff_allele1,
                ff_assay = Assay) %>%
  
  # Need to create the maternal_positives and paternal_positives columns.
  mutate(
    maternal_positives = case_when(
      identity == "paternal gDNA" ~pmin(Positives_ff_allele1, 
                                        Positives_ff_allele2),
      identity == "maternal gDNA" ~ pmax(Positives_ff_allele1, 
                                         Positives_ff_allele2)),
    
    paternal_positives = case_when(
      identity == "paternal gDNA" ~pmax(Positives_ff_allele1, 
                                        Positives_ff_allele2),
      identity == "maternal gDNA" ~ pmin(Positives_ff_allele1, 
                                         Positives_ff_allele2))) %>%
  select(-c(AcceptedDroplets_ff_allele2, Positives_ff_allele2, 
         Positives_ff_allele1)))
