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

library(epiR)

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
# Curve functions
#########################

# This step calculates the SPRT thresholds assuming a fetal fraction of 4% and
# a likelihood ratio of 8.
calc_SS_boundary <- function(total_copies) {
  q0 = 0.5
  q1 <- 0.5+(0.04/2)
  Delta <- (1- q1)/(1-q0)
  Gamma <- ((q1 * (1-q0))/ (q0*(1-q1)))
  SS_boundary <- ((log(8)/total_copies) - log(Delta))/log(Gamma)
  # Convert to a percentage for output
  return(SS_boundary*100)
}

calc_AS_upper_boundary <- function(total_copies) {
  q0 = 0.5
  q1 <- 0.5+(0.04/2)
  Delta <- (1- q1)/(1-q0)
  Gamma <- ((q1 * (1-q0))/ (q0*(1-q1)))
  AS_upper_boundary <- ((log(1/8)/total_copies) - log(Delta))/log(Gamma)
  # Convert to a percentage for output
  return(AS_upper_boundary*100)
}

calc_AS_lower_boundary <- function(total_copies) {
  q0 = 0.5
  q1 <- 0.5+(0.04/2)
  Delta <- (1- q1)/(1-q0)
  Gamma <- ((q1 * (1-q0))/ (q0*(1-q1)))
  AS_upper_boundary <- ((log(1/8)/total_copies) - log(Delta))/log(Gamma)
  AS_lower_boundary <- 0.5-(AS_upper_boundary-0.5)
  # Convert to a percentage for output
  return(AS_lower_boundary*100)
}

calc_AA_boundary <- function(total_copies) {
  q0 = 0.5
  q1 <- 0.5+(0.04/2)
  Delta <- (1- q1)/(1-q0)
  Gamma <- ((q1 * (1-q0))/ (q0*(1-q1)))
  SS_boundary <- ((log(8)/total_copies) - log(Delta))/log(Gamma)
  AA_boundary <- 0.5-(SS_boundary-0.5)
  # Convert to a percentage for output
  return(AA_boundary*100)
}

# This function calculates the boundary for the gDNA heterozygous controls.
# cfDNA samples with variant fractions above this will be predicted to be 
# affected.
control_boundary_SS <- function(input){
  output <- ((log(125)/input)+0.0409)/0.0805
  return(output*100)
}

control_boundary_AA <- function(input){
  output <- 50-(control_boundary_SS(input)-50)
  return(output)
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
# Select sickle cell disease samples
#########################

# Samples to exlude:
# 13262 - this sample had contamination
# 17004 - this sample was actually HbAC
# 20915 - this sample was from a twin pregnancy
samples_to_exclude <- c(13262, 17004, 20915)

cfDNA_scd_data <- ff_calculations(
  var_ref_calculations(cfdna_ddpcr_data)) %>%
  # Select only sickle cell disease cases
  # And cases with fetal fraction over 4%
  dplyr::rename(r_number = Sample) %>%
  filter(variant_assay == "HBB c.20A>T" &
           !r_number %in% samples_to_exclude)

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
            # Remove empty well 
            Positives != 0 &
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
  # New column to distinguish control gDNA
  mutate(sample_type = "gDNA",
        sprt_prediction = case_when(
          
         variant_percent > calc_SS_boundary(vf_assay_molecules) &
              major_allele == "variant allele"
            ~"HbSS",
         
         variant_percent < calc_AA_boundary(vf_assay_molecules) &
           major_allele == "reference allele"
         ~"HbAA",
         
         variant_percent > calc_AS_lower_boundary(vf_assay_molecules) &
          variant_percent < calc_AS_upper_boundary(vf_assay_molecules)
         ~"HbAS",
          TRUE ~"inconclusive")) %>%
  dplyr::rename(r_number = Sample)
  
#########################
# SPRT with gDNA controls
#########################

# This plot shows that the common form of the SPRT equation is 
# inappropriate for this dataset

ggplot(gDNA_scd_data, aes(x = vf_assay_molecules, y = variant_percent))+
  geom_point(size = 2, colour = "grey", pch=21)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none")+
  labs(x = "Total molecules (variant and reference)",
       y = "Variant fraction (%)") +
  xlim(200, 35000)+
  scale_y_continuous(breaks = c(45, 50, 55), 
                     limits = c(45, 55))+
  geom_function(fun = "calc_AS_upper_boundary",
                aes(x = vf_assay_molecules, y =
                      calc_AS_upper_boundary(vf_assay_molecules)),
                colour = "black")+
  geom_function(fun = "calc_AS_lower_boundary",
                aes(x = vf_assay_molecules, y =
                      calc_AS_lower_boundary(vf_assay_molecules)),
                colour = "black") +
  geom_function(fun = "calc_SS_boundary",
                aes(x = vf_assay_molecules, y =
                      calc_SS_boundary(vf_assay_molecules)),
                colour = "black")+
  geom_function(fun = "calc_AA_boundary",
                aes(x = vf_assay_molecules, y =
                      calc_AA_boundary(vf_assay_molecules)),
                colour = "black")+
  annotate(geom = "text", x = 30000, y = 50, 
                              label = "HbAS classification")+
  annotate(geom = "text", x = 30000, y = 47, 
           label = "HbAA classification") +
  annotate(geom = "text", x = 30000, y = 53, 
           label = "HbSS classification")

# This plot adds on more appropriate curves to show the range
# of variation in the gDNA dataset.

ggplot(gDNA_scd_data, aes(x = vf_assay_molecules, y = variant_percent))+
  geom_point(size = 2, colour = "grey", pch=21)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none")+
  labs(x = "Total molecules detected",
       y = "Variant fraction (%)") +
  xlim(200, 30000)+
  ylim(40, 60)+
  geom_function(fun = "control_boundary_SS",
                aes(x = vf_assay_molecules, y =
                      control_boundary_SS(vf_assay_molecules)),
                colour = "black")+
  geom_function(fun = "control_boundary_AA",
                aes(x = vf_assay_molecules, y =
                      control_boundary_AA(vf_assay_molecules)),
                colour = "black")

#########################
# Predict fetal genotypes using gDNA control data
#########################

# The aim is to classify samples as "affected" or "unaffected" depending
# on whether the variant fraction for the cfDNA sample is outside or within
# the variant fraction limits for HbAS gDNA controls in the same 
# range of molecules tested.

lr_threshold <- 8

cfDNA_scd_predictions <- cfDNA_scd_data %>%
  mutate(
    genotype_prediction = case_when(
      # Samples with fetal fractions below 4% are inconclusive
      fetal_percent < 4 ~"inconclusive",
      variant_percent > control_boundary_SS(vf_assay_molecules)
      ~"HbSS",
      variant_percent < control_boundary_AA(vf_assay_molecules)
      ~"HbAA",
      variant_percent < control_boundary_SS(vf_assay_molecules) &
      variant_percent > control_boundary_AA(vf_assay_molecules)
      ~"HbAS"),
    clinical_prediction = case_when(
      variant_percent > control_boundary_SS(vf_assay_molecules)
      ~"affected",
      variant_percent < control_boundary_SS(vf_assay_molecules)
      ~"unaffected"),
    # Calculate likelihood ratio
    likelihood_ratio = calc_LR_autosomal(fetal_fraction,
                                         major_allele_percent/100,
                                         vf_assay_molecules),
    sprt_prediction = case_when(
      likelihood_ratio > lr_threshold &
        major_allele == "variant allele"
      ~"HbSS",
      likelihood_ratio > lr_threshold &
        major_allele == "reference allele"
      ~"HbAA",
      likelihood_ratio < 1/lr_threshold
      ~"HbAS",
      TRUE ~"inconclusive"),
    
    sample_id = 
             paste0("HBB-", as.character(row.names(cfDNA_scd_data))))

# Bind the gDNA and cfDNA datasets together
gDNA_cfDNA_bind <- rbind(
  gDNA_scd_data %>%
    dplyr::rename(genotype_prediction = sample_type) %>%
    select(r_number, genotype_prediction, vf_assay_molecules, 
           vf_assay_molecules_max, 
           vf_assay_molecules_min, variant_percent, variant_percent_max, 
           variant_percent_min),
  cfDNA_scd_predictions %>%
    select(r_number, genotype_prediction, vf_assay_molecules, vf_assay_molecules_max, 
           vf_assay_molecules_min, variant_percent, variant_percent_max, 
           variant_percent_min)) %>%
  # Factorise for plotting
  mutate(genotype_prediction = factor(genotype_prediction, levels = c(
    "HbSS", "HbAS", "HbAA", "inconclusive", "gDNA")))

# Plot the gDNA controls and the cfDNA on the same plot.
# This gives a visual representation of how the samples are
# classified as "affected" or "unaffected", based on lying within or outside
# the gDNA control variation.

ggplot(gDNA_cfDNA_bind %>%
         filter(!genotype_prediction %in% c("gDNA",
                                            "inconclusive")), aes(x = vf_assay_molecules, 
                                            y = variant_percent))+
  # Add error bars only for the cfDNA data points
  # To stop the plot getting messy
  geom_errorbar(data = gDNA_cfDNA_bind %>%
                  filter(!genotype_prediction %in% c("gDNA",
                                                     "inconclusive")),
                aes(ymin = variant_percent_min, ymax = variant_percent_max),
                alpha = 0.2)+
  geom_errorbarh(data = gDNA_cfDNA_bind %>%
                   filter(!genotype_prediction %in% c("gDNA",
                                                    "inconclusive")),
                 aes(xmin = vf_assay_molecules_min, 
                     xmax = vf_assay_molecules_max),
                 alpha = 0.2)+
  # "#000000" = black; "#99CCFF" = light blue
  # "#0000FF" = dark blue; "#FFFFFF" = white
  #scale_fill_manual(values=c("#000000", "#99CCFF", "#0000FF", "#FFFFFF"))+
  #scale_alpha_manual(values=c(1, 1, 1, 0.2))+
  scale_fill_manual(values=c("#FFFFFF", "#FFFFFF", "#FFFFFF", "#000000"))+
  geom_point(size = 2, aes(fill = genotype_prediction),
             colour = "black", pch=21)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.title = element_blank(), 
    legend.position = "none")+
  labs(y = "Variant fraction (%)", x = "Total molecules detected",
       title = "cfDNA cohort")+
  #geom_function(fun = "control_boundary_SS",
                #aes(x = vf_assay_molecules, y =
                      #control_boundary_SS(vf_assay_molecules)),
                #colour = "black")+
  #geom_function(fun = "control_boundary_AA",
              #aes(x = vf_assay_molecules, y =
                      #control_boundary_AA(vf_assay_molecules)),
                #colour = "black")+
  ylim(40, 60)+
  xlim(200, 35000)+
  # Circle the incorrect predictions
  #geom_point(data=subset(gDNA_cfDNA_bind, r_number %in% 
                                 #c("10209", "17472",
                                   #"18836", "20763")),
             #pch=1,
             #size=5,
             #colour = "red",
             #alpha = 0.8)+
  #annotate(geom = "text", x = 30000, y = 50, 
           #label = "HbAS classification")+
  #annotate(geom = "text", x = 30000, y = 47, 
           #label = "HbAA classification") +
 # annotate(geom = "text", x = 30000, y = 53, 
           #label = "HbSS classification")+
  geom_hline(yintercept = 52.2, linetype = "dashed")+
  geom_hline(yintercept = 51.5, linetype = "dashed")+
  geom_hline(yintercept = 47.8, linetype = "dashed")+
  geom_hline(yintercept = 48.5, linetype = "dashed")+
  geom_vline(xintercept = 4000, linetype = "dashed")
  
# Plot of just the gDNA controls
ggplot(gDNA_cfDNA_bind %>%
         filter(genotype_prediction == "gDNA"), aes(x = vf_assay_molecules, 
                            y = variant_percent))+
  
  # "#000000" = black; "#99CCFF" = light blue
  # "#0000FF" = dark blue; "#FFFFFF" = white
  scale_fill_manual(values=c("#FFFFFF"))+
  geom_point(size = 2, aes(fill = genotype_prediction),
             alpha = 0.5,
             colour = "black", pch=21)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.title = element_blank(), 
    legend.position = "none")+
  labs(y = "Variant fraction (%)", x = "Total molecules detected",
       title = "Heterozygous gDNA controls")+
  #geom_function(fun = "control_boundary_SS",
                #aes(x = vf_assay_molecules, y =
                      #control_boundary_SS(vf_assay_molecules)),
              #colour = "black")+
  #geom_function(fun = "control_boundary_AA",
                #aes(x = vf_assay_molecules, y =
                      #control_boundary_AA(vf_assay_molecules)),
                #colour = "black")+
  ylim(40, 60)+
  xlim(200, 35000)+
  geom_hline(yintercept = 52.2, linetype = "dashed")+
  geom_hline(yintercept = 51.5, linetype = "dashed")+
  geom_hline(yintercept = 47.8, linetype = "dashed")+
  geom_hline(yintercept = 48.5, linetype = "dashed")+
  geom_vline(xintercept = 4000, linetype = "dashed")


#########################
# Z score
#########################

# Calculate the z score using the gDNA data

mean_vf <- mean(gDNA_scd_data$variant_percent)
sd_vf <- sd(gDNA_scd_data$variant_percent)

cfDNA_scd_z_score <- cfDNA_scd_data %>%
  mutate(zscore = (variant_percent - mean_vf)/sd_vf)

ggplot(cfDNA_scd_z_score, aes(x = vf_assay_molecules,
                              y = zscore))+
  geom_point()+
  theme_bw()+
  geom_hline(yintercept = 3, linetype = "dashed")+
  geom_hline(yintercept = -
               3, linetype = "dashed")

#########################
# Compare predictions against Biobank
#########################

cfDNA_scd_outcomes <- left_join(
  cfDNA_scd_predictions,
  RAPID_biobank %>%
    mutate(r_number = as.character(r_number)) %>%
    select(r_number, study_id, date_of_blood_sample, Gestation_total_weeks,
           gestation_character,
           mutation_genetic_info_fetus, Partner_sample_available),
  by = "r_number") %>%
  mutate(outcome = case_when(
    genotype_prediction == "inconclusive" ~"inconclusive",
    is.na(mutation_genetic_info_fetus) ~"awaiting result",
    genotype_prediction == mutation_genetic_info_fetus
    ~"correct",
    TRUE ~"incorrect"),
    sprt_outcome = case_when(
      sprt_prediction == "inconclusive" ~"inconclusive",
      is.na(mutation_genetic_info_fetus) ~"awaiting result",
      sprt_prediction == mutation_genetic_info_fetus
      ~"correct",
      TRUE ~"incorrect"))

# This method calls 69/72 genotypes correct (95.8%)
count(cfDNA_scd_outcomes, genotype_prediction, mutation_genetic_info_fetus)
# This method calls 73/73 clinical outcomes correct
# 8/73 samples had fetal fractions below 4% (11% inconclusive rate)

# 10209, 18836 and 20763 are incorrect

count(cfDNA_scd_outcomes, genotype_prediction)

ggplot(cfDNA_scd_outcomes, aes(x = mutation_genetic_info_fetus,
                                       y = variant_percent))+
  geom_jitter()

# False negative rate is 1/61, 
# False positive rate is 2/61
# True positive rate is 12/12

# This plot shows the balanced vs imbalanced samples
ggplot(cfDNA_scd_outcomes %>%
         mutate(deviation = abs(variant_percent-50)+50), 
       aes(x = vf_assay_molecules, 
           y = deviation))+
  # "#000000" = black; "#99CCFF" = light blue
  # "#0000FF" = dark blue; "#FFFFFF" = white
  scale_fill_manual(values=c("#000000", "#99CCFF", "#0000FF"))+
  geom_point(size = 2, aes(fill = mutation_genetic_info_fetus),
             colour = "black", pch=21)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.title = element_blank(), 
    legend.position = "bottom")+
  labs(y = "Variant fraction (%)", x = "Total molecules detected")+
  geom_function(fun = "control_boundary_SS",
                aes(x = vf_assay_molecules, y =
                      control_boundary_SS(vf_assay_molecules)),
                colour = "black")+
  geom_function(fun = "control_boundary_AA",
                aes(x = vf_assay_molecules, y =
                      control_boundary_AA(vf_assay_molecules)),
                colour = "black")+
  ylim(50, 60)+
  xlim(200, 35000)

#########################
# Results table
#########################

scd_cohort_table <- cfDNA_scd_outcomes %>%
  select(sample_id, r_number, study_id, date_of_blood_sample, 
         gestation_character, Partner_sample_available, variant_assay, 
         AcceptedDroplets_Variant_assay, Positives_variant,
         Positives_reference, ff_assay, AcceptedDroplets_FetalFrac, 
         maternal_positives, paternal_positives, variant_molecules, 
         reference_molecules, maternal_molecules, paternal_molecules,
         variant_percent, fetal_percent, genotype_prediction, 
         clinical_prediction, 
         mutation_genetic_info_fetus, outcome) %>%
  dplyr::rename(
    gestation = gestation_character, 
    "invasive genotype" = mutation_genetic_info_fetus)

view(scd_cohort_table)

#########################
# Limit of detection
#########################

LOD_data <- read_csv("data/20-1557_LOD.csv", col_names = TRUE)

LOD_data_longer <- var_ref_calculations(LOD_data %>%
  mutate(unique_identifier = paste(Input_molecules, Sample)) %>%
  pivot_wider(id_cols = c(unique_identifier, Sample, Input_molecules, Mass_molecules, fetal_fraction),
              names_from = Target,
              values_from = c(AcceptedDroplets, Positives)) %>%
  select(-c(AcceptedDroplets_HbA)) %>%
  dplyr::rename(AcceptedDroplets_Variant_assay = AcceptedDroplets_HbS,
                Positives_variant = Positives_HbS,
                Positives_reference = Positives_HbA,
                sample_type = Sample,
                r_number = unique_identifier)) %>%
  # Factorise for plot
  mutate(sample_type = factor(sample_type, levels = c(
    "0%",
    "AA 2%", "SS 2%",
    "AA 4%", "SS 4%",
    "AA 6%", "SS 6%",
    "AA 8%", "SS 8%",
    "AA 10%", "SS 10%",
    "AA 12%", "SS 12%")))

gDNA_cfDNA_lod_bind <- rbind(gDNA_scd_data %>%
                        select(r_number, sample_type, vf_assay_molecules, 
                               vf_assay_molecules_max, vf_assay_molecules_min, 
                               variant_percent, variant_percent_max, 
                               variant_percent_min),
                       LOD_data_longer %>%
                         select(r_number, sample_type, vf_assay_molecules, 
                                vf_assay_molecules_max, vf_assay_molecules_min, 
                                variant_percent, variant_percent_max, 
                                variant_percent_min),
                       cfDNA_scd_predictions %>%
                         dplyr:: rename(sample_type = genotype_prediction) %>%
                         select(r_number, sample_type, vf_assay_molecules, 
                                vf_assay_molecules_max, vf_assay_molecules_min, 
                                variant_percent, variant_percent_max, 
                                variant_percent_min)) %>%
  mutate(sample_type = factor(sample_type, levels = c(
    "0%",
    "AA 2%", "SS 2%",
    "AA 4%", "SS 4%",
    "AA 6%", "SS 6%",
    "AA 8%", "SS 8%",
    "AA 10%", "SS 10%",
    "AA 12%", "SS 12%", 
    "gDNA",
    "HbSS", "HbAA", "HbAS")))


ggplot(gDNA_cfDNA_lod_bind %>%
         filter(!sample_type %in% c("HbSS", "HbAA", "HbAS")), 
       aes(x = vf_assay_molecules, y = variant_percent))+
  #geom_errorbar(aes(ymin = variant_percent_min, ymax = variant_percent_max),
                #alpha = 0.2)+
  #geom_errorbarh(aes(xmin = vf_assay_molecules_min, 
                     #xmax = vf_assay_molecules_max))+
  # In order of shade
  # "#FFFFFF" = white
  # "#CCCCCC" = grey 1
  # "#9999CC" = grey 2
  # "#999999" = grey 3
  # "#666666" = grey 4
  # "#333333" = grey 5
  # "#000000" = black; 
  scale_fill_manual(values = c("#FFFFFF", 
                               "#CCCCCC", "#CCCCCC",
                               "#9999CC", "#9999CC",
                               "#999999", "#999999",
                               "#666666", "#666666",
                               "#333333", "#333333",
                               "#000000", "#000000",
                               "#FFFFFF"))+
  geom_point(size = 3, aes(fill = sample_type), pch=21, alpha = 0.5)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom")+
  labs(x = "Total molecules detected",
       y = "Variant fraction (%)") +
  xlim(200, 15000)+
  ylim(40, 60)+
  #geom_function(fun = "control_boundary_SS",
                #aes(x = vf_assay_molecules, y =
                      #control_boundary_SS(vf_assay_molecules)),
                #colour = "black")+
  #geom_function(fun = "control_boundary_AA",
                #aes(x = vf_assay_molecules, y =
                    #control_boundary_AA(vf_assay_molecules)),
                #colour = "black")+
  geom_hline(yintercept = 50, linetype = "dashed")+
  geom_hline(yintercept = 52, linetype = "dashed")+
  geom_hline(yintercept = 52.5, linetype = "dashed")+
  geom_hline(yintercept = 48, linetype = "dashed")+
  geom_hline(yintercept = 47.5, linetype = "dashed")+
  geom_vline(xintercept = 3000, linetype = "dashed")

nrow(cfDNA_scd_data %>%
       filter(fetal_percent > 6))


# Median is 50.09% 
median(gDNA_scd_data$variant_percent)

cfDNA_scd_predictions %>%
  filter(genotype_prediction == "HbAA") %>%
  select(r_number, fetal_percent, variant_percent, vf_assay_molecules) %>%
  arrange(vf_assay_molecules)

nrow(cfDNA_scd_data %>%
  filter(vf_assay_molecules > 3000 & fetal_percent > 4))

# Use the epiR package to calculate sensitivity
count(cfDNA_scd_outcomes, genotype_prediction, mutation_genetic_info_fetus)

# Balanced vs unbalanced
# True positives, false positives, false negatives, true negatives
scd_data <- as.table(matrix(c(24, 3, 1, 44), nrow = 2, byrow = TRUE))
scd_metrics <- epi.tests(scd_data, conf.level = 0.95)

# Affected vs unaffected
# True positives, false positives, false negatives, true negatives
scd_data <- as.table(matrix(c(12, 0, 0, 60), nrow = 2, byrow = TRUE))
scd_metrics <- epi.tests(scd_data, conf.level = 0.95)


cfDNA_scd_outcomes %>%
  filter(outcome == "incorrect") %>%
  select(r_number, sample_id)
