################################################################################
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
################################################################################

#########################
# Load libraries and resources
#########################

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

# Source functions (this includes loading ddPCR data)
source("functions/ddPCR_nipd_functions.R")
source("functions/RAPID_biobank.R")

#########################
# cfDNA SPRT analysis
#########################

# Set likelihood ratio threshold
LR_threshold <- 250

ddpcr_sprt_analysed <- ff_calculations(
  var_ref_calculations(cfdna_ddpcr_data)) %>% 
  
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

#########################
# cfDNA MCMC analysis
#########################

# Prepare ddPCR data for MCMC

# n_K	= number of droplets tested for variant assay
# K_M	= number of droplets positive for variant (mutant) allele
# K_N	= number of droplets positive for normal (reference) allele
# n_Z	= number of droplets tested for fetal fraction assay
# Z_X	= number of droplets positive for maternal homozygous allele
# Z_Y	= number of droplets positive for paternal allele

ddpcr_data_mcmc <- cfdna_ddpcr_data %>%
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

#########################
# Collation of results
#########################

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
    select(r_number, study_id, gestation_weeks, gestation_days, 
           Gestation_total_weeks, gestation_character, 
           date_of_blood_sample, vacutainer, mutation_genetic_info_fetus, 
           Partner_sample_available, original_plasma_vol, 
           tubes_plasma_current, report_acquired),
  by = "r_number")

#########################
# Sickle cell disease analysis for paper
#########################

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


# New algorithm idea

variant_percent_limit <- 51.1

new_attempt <- scd_ddpcr %>%
  mutate(new_algorithm = case_when(
    
    vf_assay_molecules < 5999 &
      fetal_percent > 4 &
      variant_percent < variant_percent_limit
    ~"unaffected",
    
    vf_assay_molecules < 5999 & 
      fetal_percent > 4 &
      variant_percent > variant_percent_limit
    ~"affected",
    
    vf_assay_molecules > 6000 &
      variant_percent < variant_percent_limit
    ~"unaffected",
    
    vf_assay_molecules > 6000 &
      variant_percent > variant_percent_limit      
    ~"affected",
    
    TRUE ~"inconclusive"),
    
    outcome = case_when(
      mutation_genetic_info_fetus == "HbAS" |
      mutation_genetic_info_fetus == "HbAA" 
       ~"unaffected",
      mutation_genetic_info_fetus == "HbSS"
      ~"affected"),
    outcome_check = ifelse(outcome == new_algorithm, "TRUE", "FALSE"))
  
ggplot(new_attempt, aes(x = fetal_percent, y = variant_percent,
                        colour = new_algorithm))+
         geom_point()

new_attempt_check <- new_attempt %>%
  filter(new_algorithm != "inconclusive") %>%
  select(r_number, new_algorithm, outcome, outcome_check)

#########################
# RMD plot function
#########################

# This is the function for plotting relative mutation dosage (rmd) 
# results for ddPCR, including parental gDNA controls. The amount of 
# information on this plot can be modified to suit user preference.

# Input parental gDNA wells as vectors

colnames(ddpcr_sprt_analysed)

plot_rmd_graph <- function(cfdna_sample, parent_vf_wells, parent_ff_wells) {
  
  # Get cfDNA data
  cfDNA_single_sample <- ddpcr_sprt_analysed %>%
    filter(r_number == cfdna_sample) %>%
    dplyr::rename(Sample = r_number) %>%
    mutate(identity = "cfDNA") %>%
    select(Sample, identity, SPRT_prediction, variant_assay, ff_assay,
           variant_molecules, reference_molecules,
           variant_molecules_max, variant_molecules_min, 
           reference_molecules_max, reference_molecules_min, 
           maternal_molecules, paternal_molecules,
           paternal_molecules_max, paternal_molecules_min,
           maternal_molecules_max, maternal_molecules_min)
  
  # Get parent data and cfDNA data for a single NIPT case and pivot into
  # the format for plotting
  case_data <- left_join(
    parent_gDNA_var_ref %>%
    # First table
      filter(Worksheet_well %in% parent_vf_wells) %>%
      select(Sample, identity, variant_assay, variant_molecules, 
             reference_molecules, variant_molecules_max, variant_molecules_min, 
             reference_molecules_max, reference_molecules_min,
             variant_assay),
    # Second table
      parent_gDNA_ff %>%
      filter(Worksheet_well %in% parent_ff_wells) %>%
      select(Sample, ff_assay, maternal_molecules, paternal_molecules,
             paternal_molecules_max, paternal_molecules_min,
             maternal_molecules_max, maternal_molecules_min),
      by = "Sample") %>%
    select(Sample, identity, variant_assay, ff_assay,
         variant_molecules, reference_molecules,
         variant_molecules_max, variant_molecules_min, 
         reference_molecules_max, reference_molecules_min, 
         maternal_molecules, paternal_molecules,
         paternal_molecules_max, paternal_molecules_min,
         maternal_molecules_max, maternal_molecules_min) %>%
    # Bind with cfDNA data
    rbind(cfDNA_single_sample %>%
            select(-SPRT_prediction))
  
  # The fetal fraction and variant assay for all samples in a case 
  # should be the same
  stopifnot(length(unique(case_data$ff_assay))==1)
  stopifnot(length(unique(case_data$variant_assay))==1)
  
  # This section can probably be simplified with a clever pivot
  # to get the max and min values in separate columns
  
  case_data_variant <- case_data %>%
    select(Sample, identity, variant_molecules,
           variant_molecules_max, variant_molecules_min) %>%
    mutate(Target_type = "Variant") %>%
    rename(molecules = variant_molecules,
           molecules_max = variant_molecules_max,
           molecules_min = variant_molecules_min)
  
  case_data_reference <- case_data %>%
    select(Sample, identity, reference_molecules,
           reference_molecules_max, reference_molecules_min) %>%
    mutate(Target_type = "Reference") %>%
    rename(molecules = reference_molecules,
           molecules_max = reference_molecules_max,
           molecules_min = reference_molecules_min)
  
  case_data_maternal <- case_data %>%
    select(Sample, identity, maternal_molecules,
           maternal_molecules_max, maternal_molecules_min) %>%
    mutate(Target_type = "Shared allele") %>%
    rename(molecules = maternal_molecules,
           molecules_max = maternal_molecules_max,
           molecules_min = maternal_molecules_min)
  
  case_data_paternal <- case_data %>%
    select(Sample, identity, paternal_molecules,
           paternal_molecules_max, paternal_molecules_min) %>%
    mutate(Target_type = "Fetal-specific allele") %>%
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
           Target_type = factor(Target_type, levels = c("Shared allele",
                                              "Fetal-specific allele",
                                              "Reference",
                                              "Variant")))
  
  rmd_plot <- ggplot(case_data_long, 
                     aes(x = identity, y = molecules, fill = Target_type))+
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
         subtitle = paste("Variant assay:", cfDNA_single_sample[1,4], "  ",
                          "Fetal fraction assay:", cfDNA_single_sample[1,5],
                          "  ", "SPRT prediction: ", 
                          cfDNA_single_sample[1,3]))+
    geom_text(aes(x = identity, y = molecules_max, label = molecules), 
              position = position_dodge(width = 0.9), vjust = -1)

return(list(rmd_plot, case_data))
}

#plot_rmd_graph("30065", "21-1863.csv_M07", "21-1862.csv_M10")
