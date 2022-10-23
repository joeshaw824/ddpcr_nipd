################################################################################
## Load ddPCR data
## October 2021
## Joseph.Shaw@gosh.nhs.uk
################################################################################

#########################
# Set working directory
#########################

setwd("W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/ddPCR_R_Analysis/ddpcr_nipd")

#########################
# Load functions and packages
#########################

# Load functions (includes tidyverse)
source("functions/ddPCR_nipd_functions.R")

# Load janitor for cleaning header names
library(janitor)

# stringr for counting string patterns
library(stringr)

#########################
# Load resources
#########################

# Control gDNA information
controls <- readr::read_csv("resources/controls.csv")

# Target panel
ddpcr_target_panel <- readr::read_csv("resources/ddpcr_target_panel.csv") 

# Number and type of plasma extractions
plasma_extractions <- read.csv("resources/extraction_volumes.csv") %>%
  group_by(r_number) %>%
  summarise(plasma_volume_ml = (sum(tubes_removed))*2) 

plasma_replicates <- read.csv("resources/extraction_volumes.csv") %>%
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
  tmp_dat <- readr::read_csv(paste0(dataPath,dataFile), col_names = TRUE)
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
  
  # Rename sample to r_number to allow merge with RAPID Biobank data.
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
stopifnot(nrow(controls) == length(unique(controls$sample)))

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

#########################