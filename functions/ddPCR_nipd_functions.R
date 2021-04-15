#############################################################
## Functions for ddPCR NIPD
## March 2021
## Joseph.Shaw@gosh.nhs.uk
#############################################################

#############################################################
# Load libraries
#############################################################

# Load packages
library(tidyverse)
library(readxl)
library(lubridate)

#############################################################
# ddPCR functions
#############################################################

# Function for Poisson correction. See Barrett et al 2012 supplementary information (PMID: 22451622)

Poisson_correct <- function(N, P) {
  num_molecules <- as.integer(-log((N-P)/N)*N)
  return(num_molecules)}

# Calculate the fetal fraction from paternal allele copies.

calc_ff <- function(Maternal_copies, Paternal_copies) {
  fetal_fraction <- (Paternal_copies*2) / (Paternal_copies + Maternal_copies)
  return(fetal_fraction)
}

# Calculates the likelihood ratio for a ddPCR test with X-linked inheritance.

calc_LR_X_linked <- function(fetal_fraction, overrep_fraction, total_copies) {
  q0 <- (1 - fetal_fraction) / (2 - fetal_fraction)
  q1 <- 1 / (2 - fetal_fraction)
  Delta <- (1- q1)/(1-q0)
  Gamma <- ((q1 * (1-q0))/ (q0*(1-q1)))
  LR <- exp((((overrep_fraction*log(Gamma)) + log(Delta))*total_copies))
  return(LR)
}

# Calculates the likelihood ratio for a ddPCR test when the variant is on an autosome (recessive or dominant).

calc_LR_autosomal <- function(fetal_fraction, overrep_fraction, total_copies) {
  q0 = 0.5
  # I modified the q1 expression to make it easier to use. Fetal fraction should be in the right format
  # I.e. 0.05 not 5
  q1 <- 0.5+(fetal_fraction/2)
  Delta <- (1- q1)/(1-q0)
  Gamma <- ((q1 * (1-q0))/ (q0*(1-q1)))
  LR <- exp((((overrep_fraction*log(Gamma)) + log(Delta))*total_copies))
  return(LR)
}

# Calculates the 95% Poisson confidence intervals.
# Poisson intervals are lamda +/- 1.96 (sqrt(lamda/N)), when lambda is copies per droplet and N is total droplets.

Poisson_max <- function(Copies_per_droplet, AcceptedDroplets) {
  Cpd_poisson_max = (Copies_per_droplet + 1.96*(sqrt(Copies_per_droplet / AcceptedDroplets))) * AcceptedDroplets
  return(Cpd_poisson_max)
}

Poisson_min <- function(Copies_per_droplet, AcceptedDroplets) {
  Cpd_poisson_min = (Copies_per_droplet - 1.96*(sqrt(Copies_per_droplet / AcceptedDroplets))) * AcceptedDroplets
  return(Cpd_poisson_min)
}

# Calculates the 95% Poisson fractional abundances based on the 95% max and min values for 
# one allele, but not the other allele.

Poisson_fraction_max <- function(copies_allele_max, copies_allele_other) {
  Allele_fraction_max = copies_allele_max / (copies_allele_max + copies_allele_other)
  return(Allele_fraction_max)
}

Poisson_fraction_min <- function(copies_allele_min, copies_allele_other) {
  Allele_fraction_min = copies_allele_min / (copies_allele_min + copies_allele_other)
  return(Allele_fraction_min)
}

#############################################################
# ddPCR analysis: composite functions
#############################################################

# This function prepares the data for the SPRT calculation by calculating the numbers of 
# molecules of each species and calculating the different fractional abundances.
calc_molecules <- function(worksheet) {
  
  worksheet_analysed <- worksheet %>%
    
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
    mutate(Over_represented_fraction_percent = (pmax(Reference_fraction, Variant_fraction))*100)
  
  return(worksheet_analysed)
}

# This function adds on the the 95% Poisson confidence intervals to the input.
calc_conf_intervals <- function(worksheet) {
  
  worksheet_analysed <- worksheet %>%
    
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
  
  return(worksheet_analysed)
}

# This function calculates the likelihood ratio for SPRT and classifies samples with autosomal and x-linked
# inheritance patterns according to a user-defined likelihood ratio threshold.
calc_SPRT <- function(worksheet, LR_threshold) {
  
  worksheet_analysed <- worksheet %>%
    
    # Perform SPRT and return the likelihood ratio
    mutate(Likelihood_ratio = case_when(
      Inheritance == "x_linked" ~ calc_LR_X_linked(Fetal_fraction, Over_represented_fraction, Molecules_variant_assay),
      Inheritance == "autosomal" ~ calc_LR_autosomal(Fetal_fraction, Over_represented_fraction, Molecules_variant_assay))) %>%
    
    # Classify based on likelihood ratio threshold supplied
    mutate(SPRT_prediction = case_when(
      Inheritance == "autosomal" & Likelihood_ratio > LR_threshold & Over_represented_fraction == Reference_fraction ~ "homozygous reference",
      Inheritance == "autosomal" & Likelihood_ratio > LR_threshold & Over_represented_fraction == Variant_fraction ~ "homozygous variant",
      Inheritance == "autosomal" & Likelihood_ratio < (1/LR_threshold) ~ "heterozygous",
      Inheritance == "autosomal" & Likelihood_ratio < LR_threshold & Likelihood_ratio > (1/LR_threshold) ~ "no call",
      Inheritance == "x_linked" & Likelihood_ratio > LR_threshold & Over_represented_fraction == Reference_fraction ~ "hemizygous reference",
      Inheritance == "x_linked" & Likelihood_ratio > LR_threshold & Over_represented_fraction == Variant_fraction ~ "hemizygous variant",
      Inheritance == "x_linked" & Likelihood_ratio < LR_threshold ~ "no call"))
  
  return(worksheet_analysed)
}

#############################################################
# ggplot2 plot theme
#############################################################

# Theme for plots
ddPCR_plot_theme <- theme_bw()+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18), 
        plot.title = element_text(size=20),
        legend.position = c(0.85, 0.85), legend.background = element_rect(fill="white"), 
        legend.title = element_text(size= 18), legend.text = element_text(size= 15))

#Add a horizontal line to show 50%
fifty_percent_line <- geom_hline(yintercept=50, linetype="dashed", size = 1)


#############################################################
# Load RAPID Biobank
#############################################################

biobank_filepath <- "I:/Genetics/RAPID/Biobank Library/AA current biobank/"

# Use filename pattern recognition with a wildcard for the date to select the current biobank.
biobank_current <- list.files(path= biobank_filepath, pattern="^RAPID_Biobank_database_.*.xlsx$")

# Loads the most recent version of the RAPID biobank.
RAPID_biobank <- read_excel(paste0(biobank_filepath,biobank_current),
                            sheet = "database",
                            # Have to specify the data types for each column as not all cells have data
                            # and read_excel cannot guess correctly.
                            col_types = c(
                              #r_number
                              "numeric",
                              #study_id
                              "text",
                              #site
                              "text",
                              #date_of_blood_sample
                              "date",
                              #time_of_blood sample
                              "date",
                              #maternal_DOB
                              "date",
                              #expected_delivery_date
                              "date",
                              #gestation_weeks
                              "numeric",
                              #gestation_days
                              "numeric",
                              #date_of_first_spin
                              "date",
                              #time_first_spin
                              "date",
                              #first_interval
                              "text",
                              #date_of_2nd_spin
                              "date",
                              #time_of_2nd_spin
                              "date",
                              #second_ Interval
                              "text",
                              #vacutainer
                              "text",
                              #original_plasma_vol
                              "numeric",
                              #current_volume
                              "numeric",
                              #removed_for_plasma
                              "text",
                              #date_removed
                              "date",
                              #blood_Pellets
                              "numeric",
                              #maternal_episode_number
                              "text",
                              #maternal_episode_removed_for
                              "text",
                              #maternal_episode_date_removed
                              "date",
                              #paternal_episode_number
                              "text",
                              #paternal_episode_removed_for
                              "text",
                              #paternal_episode_date_removed
                              "date",
                              #fetal_episode_number
                              "text",
                              #fetal_episode_removed for_excluding_PAGE
                              "text",
                              #fetal_episode_date_removed
                              "date",
                              #indication
                              "text",
                              #risk
                              "text",
                              #confirmed_diagnosis
                              "text",
                              #reference
                              "text",
                              #mutation_genetic_info_fetus
                              "text",
                              #maternal_mutation
                              "text",
                              #paternal_mutation
                              "text",
                              #confirmed_fetal_sex
                              "text",
                              #ethnicity
                              "text",
                              #maternal_weight_kg
                              "numeric",
                              #maternal_height_cm
                              "numeric",
                              #diabetic
                              "text",
                              #smoking
                              "text",
                              #before_after_invasive
                              "text",
                              #partner_date_of_birth
                              "date",
                              #additional_comments
                              "text",
                              #Recruited_to_EACH
                              "text",
                              #study_name
                              "skip",
                              #Recruited_to_PAGE
                              "skip",
                              #Recruited_to_BOOSTB4
                              "skip",
                              #Recruited_to_Qiagen
                              "text",
                              #EU_Consent
                              "text",
                              #Long_Term_follow_up
                              "text",
                              #DNA_storage_at_Sanger
                              "text",
                              #Checked_By_GOSH
                              "text",
                              #Fetal_Sex_Result_by_Lab_qPCR_dPCR
                              "text",
                              #Fetal_sex_double_checked_locally
                              "text",
                              #action_plasma
                              "text",
                              #tubes_plasma_original
                              "numeric",
                              #tubes_plasma_current
                              "numeric",
                              #Plasma_location_tray_number
                              "numeric",
                              #Plasma_location_Y
                              "numeric",
                              #Plasma_location_X
                              "numeric",
                              #action_blood
                              "text",
                              #tubes_blood_original
                              "numeric",
                              #tubes_blood_current
                              "numeric",
                              #Blood_location_tray_number
                              "numeric",
                              #Blood_location_box
                              "numeric",
                              #Reorganisation_notes
                              "text")) %>%
  # Convert gestation into numeric form
  mutate(gestation_days_as_weeks = gestation_days/7) %>%
  mutate(Gestation_total_weeks = gestation_weeks + gestation_days_as_weeks) %>%
  mutate(Partner_sample_available = ifelse(paste0("Partner",study_id) %in% study_id, "Yes", "No")) %>%
  mutate(year = year(date_of_blood_sample))
