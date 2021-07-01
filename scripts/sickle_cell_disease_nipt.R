#############################################################
## ddPCR for Non Invasive Prenatal Diagnosis (NIPD)
## July 2021
## Joseph.Shaw@gosh.nhs.uk
## This is an analysis script for the prediction of fetal 
## genotypes from cfDNA testing using ddPCR for sickle cell 
## disease, based on dosage experiments with heterozygous gDNA controls.
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

# This step selects every HbAS gDNA control which we have run, including
# single well data, and data merged by sample for each worksheet.
gDNA_scd_data <- var_ref_calculations(
  ddpcr_data %>%
    # Create a unique identifier for every well in the dataset
    mutate(worksheet_sample_well = paste0(Worksheet,"_",Sample,"_",Well)) %>%
   filter(Sample %in% controls$Sample & 
            # Remove any controls that aren't HbAS
            !Sample %in% c("NTC", "30139", "30130", 
                           "19RG-220G0191", "19RG-220G0193",
                           "21RG-120G0072") &
            Target %in% c("HbA", "HbS")) %>%
   pivot_wider(id_cols = c(worksheet_sample_well, Sample, Worksheet, Well),
               names_from = Target,
               values_from = 
                 c(AcceptedDroplets, Positives)) %>%
   select(-AcceptedDroplets_HbA) %>%
   dplyr::rename(
     AcceptedDroplets_Variant_assay =AcceptedDroplets_HbS,
     Positives_variant = Positives_HbA,
     Positives_reference = Positives_HbS)) %>%
  # Select relevant columns
  select(vf_assay_molecules, vf_assay_molecules_max, 
         vf_assay_molecules_min, variant_percent,
         variant_percent_min, variant_percent_max,
         worksheet_sample_well) %>%
  # New column to distinguish control gDNA
  mutate(sample_type = "gDNA") %>%
  dplyr::rename(Sample = worksheet_sample_well) %>%
  filter(vf_assay_molecules > 200)

# This function does not work yet.
# The aim is to classify samples as "affected" or "unaffected" depending
# on whether the variant fraction for the cfDNA sample is above or below
# the upper variant fraction limit for HbAS gDNA controls in the same 
# range of molecules tested.

dosage_predictor <- function(sample_variant_percent, limit, control_data) {
  
  # This filter step selects gDNA with equivalent numbers of molecules
  equivalent_gDNA <- control_data %>%
    dplyr::filter(vf_assay_molecules > limit)
  
  variant_percent_upper_limit <- max(equivalent_gDNA$variant_percent)
  
  prediction = case_when(
    
    sample_variant_percent > variant_percent_upper_limit
    ~ "affected", 
    
    sample_variant_percent < variant_percent_upper_limit
    ~ "unaffected")
  
  rm(equivalent_gDNA)
  return(prediction)
}

cfDNA_scd_data <- ff_calculations(
  var_ref_calculations(cfdna_ddpcr_data )) %>%
  # Select only sickle cell disease cases
  # And cases with fetal fraction over 4%
  filter(variant_assay == "HBB c.20A>T" & fetal_percent > 4) %>%
  # This is where the dosage_predictor function will go
  mutate(sample_type = "cfDNA") %>%
  select(vf_assay_molecules, vf_assay_molecules_max, 
         vf_assay_molecules_min, variant_percent,
         variant_percent_min, variant_percent_max, Sample,
         sample_type)

# Bind data frames together
gDNA_cfDNA_bind <- rbind(cfDNA_scd_data, gDNA_scd_data)

# Plot the gDNA controls and the cfDNA on the same plot.
# This gives a visual representation of how the samples will
# be classified as "affected" or "unaffected", based on lying above the
ggplot(gDNA_cfDNA_bind %>%
         arrange(vf_assay_molecules), aes(x = vf_assay_molecules, 
                                          y = variant_percent))+
  geom_point(aes(colour = sample_type))+
  geom_errorbar(aes(ymin = variant_percent_min, ymax = variant_percent_max),
                alpha = 0.2)+
  geom_errorbarh(aes(xmin = vf_assay_molecules_min, xmax = vf_assay_molecules_max),
                 alpha = 0.2)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  ylim(42, 58)+
  xlim(0, 42000)


