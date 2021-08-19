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
library(janitor)
library(stringi)

# Load resources
controls <- readr::read_csv("resources/controls.csv")
ddpcr_target_panel <- readr::read_csv("resources/ddpcr_target_panel.csv") 

#########################
# ddPCR functions
#########################

# Function for Poisson correction. See Barrett et al 2012 
# supplementary information (PMID: 22451622)

poisson_correct <- function(N, P) {
  num_molecules <- as.integer(-log((N-P)/N)*N)
  return(num_molecules)}

# Calculate the fetal fraction from paternal allele copies.

calc_ff <- function(maternal_copies, paternal_copies) {
  fetal_fraction <- (paternal_copies*2) / (paternal_copies + maternal_copies)
  return(fetal_fraction)
}

# Calculates the 95% Poisson confidence intervals.
# Poisson intervals are lamda +/- 1.96 (sqrt(lamda/N)), 
# when lambda is copies per droplet (cpd) and N is total droplets.

poisson_max <- function(copies_per_droplet, droplets) {
  copies_poisson_max = (copies_per_droplet + 
                       1.96*(sqrt(copies_per_droplet / droplets))) * droplets
  return(copies_poisson_max)
}

poisson_min <- function(copies_per_droplet, droplets) {
  copies_poisson_min = (copies_per_droplet - 
                       1.96*(sqrt(copies_per_droplet / droplets))) * droplets
  return(copies_poisson_min)
}

# Calculates the 95% Poisson fractional abundances based on the 95% max 
# and min values for one allele, but not the other allele.

poisson_fraction_max <- function(copies_allele_max, copies_allele_other) {
  allele_fraction_max = copies_allele_max / 
    (copies_allele_max + copies_allele_other)
  return(allele_fraction_max)
}

poisson_fraction_min <- function(copies_allele_min, copies_allele_other) {
  allele_fraction_min = copies_allele_min / 
    (copies_allele_min + copies_allele_other)
  return(allele_fraction_min)
}

#########################
# SPRT functions
#########################

# Calculates the likelihood ratio (lr) for a ddPCR test with 
# X-linked inheritance.

calc_lr_x_linked <- function(fetal_fraction, overrep_fraction, total_copies) {
  q0 <- (1 - fetal_fraction) / (2 - fetal_fraction)
  q1 <- 1 / (2 - fetal_fraction)
  delta <- (1- q1)/(1-q0)
  gamma <- ((q1 * (1-q0))/ (q0*(1-q1)))
  lr <- exp((((overrep_fraction*log(gamma)) + log(delta))*total_copies))
  return(lr)
}

# Calculates the likelihood ratio for a ddPCR test when the variant is 
# on an autosome (recessive or dominant).

calc_lr_autosomal <- function(fetal_fraction, overrep_fraction, total_copies) {
  q0 = 0.5
  # I modified the q1 expression to make it easier to use. 
  # Fetal fraction should be in the right format
  # I.e. 0.05 not 5
  q1 <- 0.5+(fetal_fraction/2)
  delta <- (1- q1)/(1-q0)
  gamma <- ((q1 * (1-q0))/ (q0*(1-q1)))
  lr <- exp((((overrep_fraction*log(gamma)) + log(delta))*total_copies))
  return(lr)
}

# These functions calculate the SPRT thresholds assuming 
# a fetal fraction of 4% and a likelihood ratio of 8.

# Variables for each function
q0 <- 0.5
q1 <- 0.5+(0.04/2)
delta <- (1- q1)/(1-q0)
gamma <- ((q1 * (1-q0))/ (q0*(1-q1)))

calc_SS_boundary <- function(total_copies) {
  SS_boundary <- ((log(8)/total_copies) - log(delta))/log(gamma)
  # Convert to a percentage for output
  return(SS_boundary*100)
}

calc_AS_upper_boundary <- function(total_copies) {
  AS_upper_boundary <- ((log(1/8)/total_copies) - log(delta))/log(gamma)
  # Convert to a percentage for output
  return(AS_upper_boundary*100)
}

calc_AS_lower_boundary <- function(total_copies) {
  AS_upper_boundary <- ((log(1/8)/total_copies) - log(delta))/log(gamma)
  AS_lower_boundary <- 0.5-(AS_upper_boundary-0.5)
  # Convert to a percentage for output
  return(AS_lower_boundary*100)
}

calc_AA_boundary <- function(total_copies) {
  SS_boundary <- ((log(8)/total_copies) - log(delta))/log(gamma)
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

# Example: fetal_percent_max

# These functions calculate the number of molecules and fractional
# abundance for each target, including 95% Poisson confidence intervals.
# This is to prevent massive duplication of code when performing these
# calculations for gDNA, cfDNA and limit of detection datasets.

var_ref_calculations <- function(data_input) {
  
  # The data_input needs these variables for this function to work:
  # vf_assay_droplets, variant_positives, reference_positives
  
  stopifnot(c("vf_assay_droplets", "variant_positives", 
              "reference_positives")
            %in% colnames(data_input))
  
  data_output <- data_input %>%
    
    mutate(
      # Perform Poisson correction to determine the total number of molecules 
      # detected for each target.
      variant_molecules = poisson_correct(vf_assay_droplets,
                                          variant_positives), 
      reference_molecules = poisson_correct(vf_assay_droplets,
                                            reference_positives),   
      
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
      variant_molecules_max = poisson_max((variant_molecules/ 
                                             vf_assay_droplets), 
                                          vf_assay_droplets),
      
      variant_molecules_min = poisson_min((variant_molecules/
                                             vf_assay_droplets),
                                          vf_assay_droplets),
      
      reference_molecules_max = poisson_max((reference_molecules/
                                               vf_assay_droplets),
                                            vf_assay_droplets),
      
      reference_molecules_min = poisson_min((reference_molecules/
                                               vf_assay_droplets), 
                                            vf_assay_droplets),
      
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
      variant_percent_max = (poisson_fraction_max(
        variant_molecules_max, reference_molecules))*100,
      
      variant_percent_min = (poisson_fraction_min(
        variant_molecules_min, reference_molecules))*100,
      
      reference_percent_max = (poisson_fraction_max(
        reference_molecules_max, variant_molecules))*100,
      
      reference_percent_min = (poisson_fraction_min(
        reference_molecules_min, variant_molecules))*100,
      
      major_allele_percent_max = case_when(
        major_allele == "variant allele" ~variant_percent_max,
        major_allele == "reference allele" ~reference_percent_max),
      
      major_allele_percent_min = case_when(
        major_allele == "variant allele" ~variant_percent_min,
        major_allele == "reference allele" ~reference_percent_min),
      
      # Calculate the 95% confidence intervals of the numbers of molecules 
      # detected by each assay
      vf_assay_molecules_max = poisson_max((
        vf_assay_molecules/vf_assay_droplets), 
        vf_assay_droplets),
      
      vf_assay_molecules_min = poisson_min((
        vf_assay_molecules/vf_assay_droplets),
        vf_assay_droplets))
  
  return(data_output)
}

ff_calculations <- function(data_input) {
  
  # The data_input needs these variables for this function to work:
  # maternal_positives, paternal_positives, ff_assay_droplets
  
  stopifnot(c("maternal_positives", "paternal_positives",
              "ff_assay_droplets") %in% colnames(data_input))
  
  data_output <- data_input %>%
    mutate(
      
      # Calculate the same variables for the fetal fraction assay
      
      maternal_molecules = poisson_correct(
        ff_assay_droplets,maternal_positives),
      
      paternal_molecules = poisson_correct(
        ff_assay_droplets,paternal_positives),
      
      # Calculate the fetal fraction
      fetal_fraction = calc_ff(maternal_molecules, paternal_molecules),
      
      fetal_percent = fetal_fraction*100,
      
      ff_assay_molecules = maternal_molecules + paternal_molecules,
      
      paternal_molecules_max = poisson_max((
        paternal_molecules / ff_assay_droplets), 
        ff_assay_droplets),
      
      paternal_molecules_min = poisson_min((
        paternal_molecules / ff_assay_droplets), 
        ff_assay_droplets),
      
      maternal_molecules_max = poisson_max((
        maternal_molecules / ff_assay_droplets), 
        ff_assay_droplets),
      
      maternal_molecules_min = poisson_min((
        maternal_molecules / ff_assay_droplets), 
        ff_assay_droplets),
      
      ff_assay_molecules_max = poisson_max((
        ff_assay_molecules / ff_assay_droplets), 
        ff_assay_droplets),
      
      ff_assay_molecules_min = poisson_min((
        ff_assay_molecules / ff_assay_droplets), 
        ff_assay_droplets),
      
      # When calculating for the fetal fraction, the paternal fraction 
      # must be multiplied by 2.
      fetal_percent_max = 200* (poisson_fraction_max(
        paternal_molecules_max, maternal_molecules)),
      
      fetal_percent_min = 200* (poisson_fraction_min(
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
  mutate(num_wells = str_count(merged_wells, ",")+1)

# Reshape the data frame to sum all values by Target.
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

# Use pivot_wider to get one row per sample
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
  # This removes several samples which didn't have both assays performed
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
  # fetal fraction calculation
  mutate(maternal_positives = pmax(ff_allele1_positives, ff_allele2_positives),
         paternal_positives = pmin(ff_allele1_positives, ff_allele2_positives))

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

parent_gDNA_ff <- ff_calculations(ddpcr_controls_with_target %>%
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
  mutate(
    maternal_positives = case_when(
      identity == "paternal gDNA" ~pmin(ff_allele1_positives, 
                                        ff_allele2_positives),
      identity == "maternal gDNA" ~ pmax(ff_allele1_positives, 
                                         ff_allele2_positives)),
    
    paternal_positives = case_when(
      identity == "paternal gDNA" ~pmax(ff_allele1_positives, 
                                        ff_allele2_positives),
      identity == "maternal gDNA" ~ pmin(ff_allele1_positives, 
                                         ff_allele2_positives))) %>%
  select(-c(ff_allele2_droplets, ff_allele2_positives, 
         ff_allele1_positives)))

#########################
# RMD plot function
#########################

# This is the function for plotting relative mutation dosage (RMD) 
# results for ddPCR, including parental gDNA controls. The amount of 
# information on this plot can be modified to suit user preference.

# Materal sample only example
# draw_rmd_plot("30065", "21-1863.csv_M07", "21-1862.csv_M10")

# Maternal and paternal samples example
# draw_rmd_plot("20238", 
              # c("21-1413.csv_M02", "21-1413.csv_M08"), 
              # c("21-1412.csv_M08",	"21-1412.csv_M10"))

draw_rmd_plot <- function(cfdna_sample, parent_vf_wells, parent_ff_wells) {
  
  # Get cfDNA data
  cfDNA_rmd <- ddpcr_sprt_unblinded %>%
    filter(r_number %in% cfdna_sample) %>%
    dplyr::rename(sample = r_number) %>%
    mutate(identity = "cfDNA") %>%
    select(sample, identity, sprt_prediction, fetal_percent,
           variant_percent, 
           vf_assay, ff_assay, mutation_genetic_info_fetus,
           variant_molecules, reference_molecules,
           variant_molecules_max, variant_molecules_min, 
           reference_molecules_max, reference_molecules_min, 
           maternal_molecules, paternal_molecules,
           paternal_molecules_max, paternal_molecules_min,
           maternal_molecules_max, maternal_molecules_min)
  
  # Get parent data and cfDNA data for a single NIPT case and pivot into
  # the format for plotting
  parents_rmd <- left_join(
    parent_gDNA_var_ref %>%
      # First table
      filter(worksheet_well_sample %in% parent_vf_wells) %>%
      select(sample, identity, vf_assay, variant_molecules, 
             reference_molecules, variant_molecules_max, variant_molecules_min, 
             reference_molecules_max, reference_molecules_min),
    # Second table
    parent_gDNA_ff %>%
      filter(worksheet_well_sample %in% parent_ff_wells) %>%
      select(sample, ff_assay, maternal_molecules, paternal_molecules,
             paternal_molecules_max, paternal_molecules_min,
             maternal_molecules_max, maternal_molecules_min),
    by = "sample") %>%
    select(sample, identity, vf_assay, ff_assay,
           variant_molecules, reference_molecules,
           variant_molecules_max, variant_molecules_min, 
           reference_molecules_max, reference_molecules_min, 
           maternal_molecules, paternal_molecules,
           paternal_molecules_max, paternal_molecules_min,
           maternal_molecules_max, maternal_molecules_min)
  
  # Bind with cfDNA data
  case_rmd <- rbind(parents_rmd,
                    cfDNA_rmd %>%
                      select(sample, identity, vf_assay, 
                             ff_assay, variant_molecules, 
                             reference_molecules, variant_molecules_max, 
                             variant_molecules_min, reference_molecules_max, 
                             reference_molecules_min, maternal_molecules, 
                             paternal_molecules, paternal_molecules_max, 
                             paternal_molecules_min, maternal_molecules_max, 
                             maternal_molecules_min))
  
  # The fetal fraction and variant assay for all samples in a case 
  # should be the same
  # na.omit used because some cases do not have both assays tested for 
  # parental samples
  stopifnot(length(unique(na.omit(case_rmd$ff_assay)))==1)
  stopifnot(length(unique(na.omit(case_rmd$vf_assay)))==1)
  
  # This section can probably be simplified with a clever pivot
  # to get the max and min values in separate columns
  
  case_data_variant <- case_rmd %>%
    select(sample, identity, variant_molecules,
           variant_molecules_max, variant_molecules_min) %>%
    mutate(target_type = "Variant") %>%
    rename(molecules = variant_molecules,
           molecules_max = variant_molecules_max,
           molecules_min = variant_molecules_min)
  
  case_data_reference <- case_rmd %>%
    select(sample, identity, reference_molecules,
           reference_molecules_max, reference_molecules_min) %>%
    mutate(target_type = "Reference") %>%
    rename(molecules = reference_molecules,
           molecules_max = reference_molecules_max,
           molecules_min = reference_molecules_min)
  
  case_data_maternal <- case_rmd %>%
    select(sample, identity, maternal_molecules,
           maternal_molecules_max, maternal_molecules_min) %>%
    mutate(target_type = "Shared allele") %>%
    rename(molecules = maternal_molecules,
           molecules_max = maternal_molecules_max,
           molecules_min = maternal_molecules_min)
  
  case_data_paternal <- case_rmd %>%
    select(sample, identity, paternal_molecules,
           paternal_molecules_max, paternal_molecules_min) %>%
    mutate(target_type = "Fetal-specific allele") %>%
    rename(molecules = paternal_molecules,
           molecules_max = paternal_molecules_max,
           molecules_min = paternal_molecules_min)
  
  case_data_long <- rbind(case_data_variant,
                          case_data_reference,
                          case_data_maternal,
                          case_data_paternal)%>%
    
    # Control factor order for plot
    mutate(identity = factor(identity, levels = c("maternal gDNA",
                                                  "paternal gDNA",
                                                  "cfDNA")),
           target_type = factor(target_type, levels = c("Shared allele",
                                                        "Fetal-specific allele",
                                                        "Reference",
                                                        "Variant")))
  
  rmd_plot <- ggplot(case_data_long, 
                     aes(x = identity, y = molecules, fill = target_type)) +
    geom_col(position = position_dodge(width = 0.9), 
             colour="black", alpha = 0.6)+
    geom_errorbar(aes(ymin = molecules_min, ymax = molecules_max, 
                      width = 0.3), position = position_dodge(width = 0.9))+
    scale_fill_manual(values = c("#99FFFF", "#FFCC99", "#3366FF", "#FF0000"))+
    theme_bw()+
    theme(axis.text=element_text(size=18), axis.title = element_text(size=18),
          legend.position = "bottom", legend.title = element_blank(),
          legend.text = element_text(size= 14),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    labs(x = "", y = "Molecules", 
         title = paste("cfDNA:", cfdna_sample),
         subtitle = paste("Variant fraction assay:", cfDNA_rmd[1,6], "  ",
                          "Fetal fraction assay:", cfDNA_rmd[1,7],
                          "  ", "SPRT prediction: ",
                          cfDNA_rmd[1,3],
                          "  ", "Fetal fraction: ",
                          round(cfDNA_rmd[1,4], 1),"%  ",
                          "Variant fraction: ",
                          round(cfDNA_rmd[1,5], 1),"%  ",
                          "Invasive result: ",
                          cfDNA_rmd[1,8]))+
    geom_text(aes(x = identity, y = molecules_max, label = molecules), 
              position = position_dodge(width = 0.9), vjust = -1)
  
  return(rmd_plot)
}

#########################
# Sequence functions
#########################

reverse_complement <- function(input_sequence){
  
  rev_comp <- stri_reverse(chartr("ATGC","TACG",input_sequence))
  
  return(rev_comp)
  
}
