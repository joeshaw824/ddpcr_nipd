################################################################################
## ddPCR variant fraction assays
## joseph.shaw3@nhs.net
## Collation of all variant fraction assays used in the ddPCR project.
################################################################################

#########################
# Load resources
#########################

library(tidyverse)
library(readxl)
library(stringr)
library(stringi)

# Working directory
setwd("W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/ddPCR_R_Analysis/ddpcr_nipd")

# Load functions for "reverse_complement"
source("functions/ddPCR_nipd_functions.R")

#########################
# Load sequences
#########################

# Read in amplicons
amplicons <- read_excel("resources/ddpcr_assay_primer_probe_sequences.xlsx",
                        sheet = "amplicons") %>%
  mutate(
    reference_amplicon_length = case_when(
    reference_amplicon == "Not provided" ~"Not provided",
    TRUE ~as.character(nchar(reference_amplicon))),
    
    variant_amplicon_length = case_when(
      variant_amplicon == "Not provided" ~"Not provided",
      TRUE ~as.character(nchar(variant_amplicon)))) %>%
  # Remove ZFXY assay
  filter(assay_name != "ZFXY")

primer_probes <- read_excel("resources/ddpcr_assay_primer_probe_sequences.xlsx",
                            sheet = "sequences") %>%
  # Remove "+" from locked nucleic acid probes
  mutate(variant_probe_no_lnas = gsub("+","",variant_probe, fixed = TRUE),
         reference_probe_no_lnas = gsub("+","",reference_probe, fixed = TRUE),
         
         # Create reverse complement primers
         reverse_primer_rc = reverse_complement(reverse_primer),
         forward_primer_rc = reverse_complement(forward_primer),
         annealing_temperature = round(annealing_temperature,0)) %>%
  # Remove ZFXY assay
  filter(assay_name != "ZFXY")

#########################
# Collate information and sense check
#########################

amplicons_primers_probes <- full_join(
  amplicons,
  primer_probes,
  by = "assay_name") %>%
  # Checks
  mutate(
    var_in_var_amp = str_detect(variant_amplicon, variant_probe_no_lnas),
    ref_in_ref_amp = str_detect(reference_amplicon, reference_probe_no_lnas),
    var_in_ref_amp = str_detect(variant_amplicon, reference_probe_no_lnas),
    ref_in_var_amp = str_detect(reference_amplicon, variant_probe_no_lnas),
    # Check forward primer
    forward_in_ref_amp = str_detect(reference_amplicon, forward_primer),
    forward_rc_in_ref_amp = str_detect(reference_amplicon, forward_primer_rc),
    forward_in_var_amp = str_detect(variant_amplicon, forward_primer),
    forward_rc_in_var_amp = str_detect(variant_amplicon, forward_primer_rc),
    # Check reverse primer
    reverse_in_ref_amp = str_detect(reference_amplicon, reverse_primer),
    reverse_rc_in_ref_amp = str_detect(reference_amplicon, reverse_primer_rc),
    reverse_in_var_amp = str_detect(variant_amplicon, reverse_primer),
    reverse_rc_in_var_amp = str_detect(variant_amplicon, reverse_primer_rc),
    
    # Sense check
    sense_check = (case_when(
      
      # Option 1
      var_in_var_amp == "TRUE" &
        ref_in_ref_amp ==  "TRUE" &
        var_in_ref_amp == "FALSE" &
        ref_in_var_amp == "FALSE" &
        forward_in_ref_amp == "TRUE" &
        forward_rc_in_ref_amp == "FALSE" &
        forward_in_var_amp == "TRUE" &
        forward_rc_in_var_amp == "FALSE" &
        reverse_in_ref_amp == "FALSE" &
        reverse_rc_in_ref_amp == "TRUE" &
        reverse_in_var_amp == "FALSE" &
        reverse_rc_in_var_amp == "TRUE"
      ~ "pass",
      
      # Option 2
      var_in_var_amp == "TRUE" &
        ref_in_ref_amp ==  "TRUE" &
        var_in_ref_amp == "FALSE" &
        ref_in_var_amp == "FALSE" &
        forward_in_ref_amp == "FALSE" &
        forward_rc_in_ref_amp == "TRUE" &
        forward_in_var_amp == "FALSE" &
        forward_rc_in_var_amp == "TRUE" &
        reverse_in_ref_amp == "TRUE" &
        reverse_rc_in_ref_amp == "FALSE" &
        reverse_in_var_amp == "TRUE" &
        reverse_rc_in_var_amp == "FALSE"
      ~ "pass",
      
      TRUE ~"fail")))

#########################
# Format for paper
#########################

ddpcr_assay_sequences <- amplicons_primers_probes %>%
  select(assay_name, forward_primer, reverse_primer, 
         reference_probe, variant_probe, fluorophores, 
         modifications, manufacturer, annealing_temperature, 
         reference_amplicon, variant_amplicon,
         reference_amplicon_length, variant_amplicon_length, 
         context_sequence) %>%
  mutate(manufacturer = factor(manufacturer, 
                               levels = c("IDT PrimeTime qPCR",
                                          "BioRad ddPCR assay",
                                          "IDT Eclipse MGB",
                                          "IDT Locked Nucleic Acid",
                                          "ThermoFisher TaqMan"))) %>%
  arrange(manufacturer)

export_timestamp(ddpcr_assay_sequences)

#########################