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
# ggpubr for exporting graphs
library(ggpubr)

# bayesplot and posterior can be added for MCMC diagnostic modelling
#library(bayesplot)
#library(posterior)

## Set working directory
setwd("W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/ddPCR_R_Analysis/ddpcr_nipd")

# Load resources
sample_wells <- read.csv("resources/sample_wells.csv",
                         colClasses = c("character", "character", 
                                        "character", "character",
                                        "character"))

# Source functions (this includes loading ddPCR data)
source("functions/ddPCR_nipd_functions.R")
source("functions/RAPID_biobank.R")

#########################
# cfDNA SPRT analysis
#########################

# Set likelihood ratio threshold
lr_threshold <- 250

ddpcr_sprt_analysed <- ff_calculations(
  var_ref_calculations(cfdna_ddpcr_data)) %>% 
  
  # Rename sample to r_number to allow merge with RAPID Biobank data.
  rename(r_number = sample) %>%
  
  mutate(
    total_molecules = vf_assay_molecules + ff_assay_molecules,

  # Perform SPRT and return the likelihood ratio
  likelihood_ratio = case_when(
     inheritance_chromosomal == "x_linked" ~ 
       calc_lr_x_linked(fetal_fraction, (major_allele_percent/100), 
                        vf_assay_molecules),
     inheritance_chromosomal == "autosomal" ~ 
       calc_lr_autosomal(fetal_fraction, (major_allele_percent/100), 
                         vf_assay_molecules)),
   
  # Classify based on likelihood ratio threshold supplied
  # Fetal genotype predictions are named consistently as 
  # "inconclusive", "heterozygous", "homozygous/hemizygous reference/variant"
  sprt_prediction = case_when(
      
      inheritance_chromosomal == "autosomal" &
      likelihood_ratio > lr_threshold &
      major_allele == "reference allele" 
          ~ "homozygous reference",
      
      inheritance_chromosomal == "autosomal" &
      likelihood_ratio > lr_threshold &
      major_allele == "variant allele" 
          ~ "homozygous variant",
      
      inheritance_chromosomal == "autosomal" &
      likelihood_ratio < (1/lr_threshold)
          ~ "heterozygous",
      
      inheritance_chromosomal == "autosomal" &
      likelihood_ratio < lr_threshold &
      likelihood_ratio > (1/lr_threshold) 
          ~ "inconclusive",
      
      inheritance_chromosomal == "x_linked" &
      likelihood_ratio > lr_threshold &
      major_allele == "reference allele" 
          ~ "hemizygous reference",
      
      inheritance_chromosomal == "x_linked" &
      likelihood_ratio > lr_threshold &
      major_allele == "variant allele" 
          ~ "hemizygous variant",
      
      inheritance_chromosomal == "x_linked" &
      likelihood_ratio < lr_threshold 
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
  inheritance_chromosomal == "autosomal" &
  inheritance_pattern == "dominant" ~map(
    data, ~ dominant_model$sample(data = .,
                                  init = initialise_chains_dominant,
                                  step_size = 0.2,
                                  parallel_chains = parallel::detectCores())),
   
  inheritance_chromosomal == "autosomal" & 
  inheritance_pattern == "recessive" ~ map(
    data, ~ recessive_model$sample(data = .,
                                  init = initialise_chains_recessive,
                                  step_size = 0.2,
                                  parallel_chains = parallel::detectCores())),
   
  inheritance_chromosomal == "x_linked" ~map(
    data, ~ x_linked_model$sample(data = .,
                                 init = initialise_chains_xlinked,
                                 step_size = 0.2,
                                 parallel_chains = parallel::detectCores()))),
 
 results = case_when(
    inheritance_chromosomal == "autosomal" &
    inheritance_pattern == "dominant" ~map(
      fit,  ~ setNames(pivot_wider(.$summary(c("pG", "rho"), "mean"),
                                  names_from = "variable",
                                  values_from = "mean"),
                       c("p_G1", "p_G2", "rho_est"))),
   
   inheritance_chromosomal == "autosomal" & 
    inheritance_pattern == "recessive" ~map(
      fit,  ~ setNames(pivot_wider(.$summary(c("pG", "rho"), "mean"),
                                 names_from = "variable",
                                 values_from = "mean"),
                     c("p_G1", "p_G2", "p_G3", "rho_est"))),
   
   inheritance_chromosomal == "x_linked" ~map(
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
    inheritance_chromosomal == "autosomal" & 
    inheritance_pattern == "dominant" & 
    p_G1 > mcmc_threshold ~"heterozygous",
    
    inheritance_chromosomal == "autosomal" &
    inheritance_pattern == "dominant" & 
    p_G2 > mcmc_threshold ~"homozygous reference",
    
    inheritance_chromosomal == "autosomal" &
    inheritance_pattern == "dominant" & 
    p_G1 < mcmc_threshold & 
    p_G2 < mcmc_threshold ~"inconclusive",
    
    # Recessive predictions
    inheritance_chromosomal == "autosomal" &
    inheritance_pattern == "recessive" & 
    p_G1 > mcmc_threshold ~"homozygous reference",
    
    inheritance_chromosomal == "autosomal" &
    inheritance_pattern == "recessive" &
    p_G2 > mcmc_threshold ~"heterozygous",
    # ok
    
    inheritance_chromosomal == "autosomal" &
    inheritance_pattern == "recessive" &
    p_G3 > mcmc_threshold ~"homozygous variant",
    
    inheritance_chromosomal == "autosomal" &
    inheritance_pattern == "recessive" & 
    p_G1 < mcmc_threshold & 
    p_G2 < mcmc_threshold ~"inconclusive",
    
    # X-linked predictions
    inheritance_chromosomal == "x_linked" &
    p_G0 > mcmc_threshold ~"hemizygous reference",
    
    inheritance_chromosomal == "x_linked" &
    p_G1 > mcmc_threshold ~"hemizygous variant",
    
    inheritance_chromosomal == "x_linked" &
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

write.csv(ddpcr_nipd_unblinded, "analysis_outputs/ddpcr_nipd_unblinded_20210811.csv",
          row.names = FALSE)

# If you don't want to run the MCMC pipeline every time, 
# just compare using SPRT

ddpcr_sprt_unblinded <- left_join(
  ddpcr_sprt_analysed,
  RAPID_biobank %>%
    # Change r_number to a character to match ddpcr_analysed
    mutate(r_number = as.character(r_number)) %>%
    filter(r_number %in% ddpcr_sprt_analysed$r_number) %>%
    select(r_number, study_id, gestation_weeks, gestation_days, 
           Gestation_total_weeks, gestation_character, 
           date_of_blood_sample, vacutainer, mutation_genetic_info_fetus, 
           Partner_sample_available, original_plasma_vol, 
           tubes_plasma_current, report_acquired),
  by = "r_number")

#########################
# Plot bespoke cohort cases
#########################

bespoke_cases <- ddpcr_sprt_analysed %>%
  filter(vf_assay != "HBB c.20A>T")

bespoke_wells <- sample_wells %>%
  filter(cfdna_sample %in% bespoke_cases$r_number)

bespoke_plots <-list()

# Plot an rmd plot for each case
for (i in bespoke_wells$cfdna_sample) {
  
  case <- bespoke_wells %>%
    filter(cfdna_sample == i)
  
  new_plot <- draw_rmd_plot(i, 
                      c(case[, 2], case[, 3]),
                      c(case[, 4], case[, 5]))
  bespoke_plots <- list(bespoke_plots, new_plot)
  rm(new_plot)
}

# Export bespoke cohort plots a single pdf
ggexport(plotlist = bespoke_plots, filename = "plots/bespoke_cohort.pdf",
         width=15, height=8, res=300)

#########################
# Plot sickle cell disease cohort cases
#########################

secondary_cohort <- c("14182", "19868", "20238", "20611", 
                      "20874", "30063", "30068", "30113", "30142", 
                      "30206", "30228", "30230", "30078", "30065", 
                      "13402", "20939", "30215", "30203",
                      "20911", "30236", "30112")

scd_cases_secondary <- ddpcr_sprt_analysed %>%
  filter(vf_assay == "HBB c.20A>T" & 
           r_number %in% secondary_cohort)

scd_wells_secondary <- sample_wells %>%
  filter(cfdna_sample %in% scd_cases_secondary$r_number)

scd_plots <-list()

# Plot an rmd plot for each case
for (i in scd_wells_secondary$cfdna_sample) {
  
  case <- scd_wells_secondary %>%
    filter(cfdna_sample == i)
  
  new_plot <- draw_rmd_plot(i, 
                            c(case[, 2], case[, 3]),
                            c(case[, 4], case[, 5]))
  scd_plots <- list(scd_plots, new_plot)
  rm(new_plot)
}

# Export sickle cell plots as a single pdf
ggexport(plotlist = scd_plots, 
         filename = "plots/sickle_cell_secondary_cohort.pdf",
         width=15, height=8, res=300)

#########################
# Fetal fraction versus variant fraction
#########################

ggplot(ddpcr_sprt_analysed %>%
         filter(!r_number %in% samples_to_exclude &
                  inheritance_chromosomal == "autosomal"), aes(x = fetal_percent,
                                y = major_allele_percent,
                                colour = sprt_prediction)) +
  geom_errorbar(aes(ymin = major_allele_percent_min,
                    ymax = major_allele_percent_max),
                alpha = 0.2) +
  geom_errorbarh(aes(xmin = fetal_percent_min,
                     xmax = fetal_percent_max),
                 alpha = 0.2) +
  geom_point(pch=21, fill = "white", size = 3) +
  theme_bw() +
  theme(
    panel.grid.major =  element_blank(),
    panel.grid.minor = element_blank()) +
  labs(y = "Over-represented allele (%)",
       x = "Fetal fraction (%)",
       title = "ddPCR for maternally-inherited autosomal variants") +
  geom_abline(slope = 0.5, intercept = 50, linetype = "dashed") +
  geom_hline(yintercept = 50, linetype = "dashed") +
  xlim(0, 21) +
  ylim(50, 60)

#########################
# Worksheet 21-2418
#########################

# IDS, MAGED2 and ADA assays

cfDNA_samples <- c("19611", "20980", "14142")

parental_samples <- c("21RG-027G0070", "21RG-027G0010", "21RG-047G0089",
                      "21RG-047G0093")

parents_and_cfDNA <- rbind(parent_gDNA_var_ref %>%
                             filter(sample %in% parental_samples &
                                      worksheet_well_sample %in% 
                                      grep("21-2418.csv", 
                                           parent_gDNA_var_ref$worksheet_well_sample,
                                           value = TRUE)) %>%
  select(sample, vf_assay, vf_assay_molecules, vf_assay_molecules_max, 
          vf_assay_molecules_min, variant_percent, variant_percent_max, 
         variant_percent_min) %>%
  mutate(sample_type = "gDNA"),
  
  ddpcr_sprt_analysed %>%
    filter(r_number %in% cfDNA_samples) %>%
    select(r_number, vf_assay, vf_assay_molecules, vf_assay_molecules_max, 
           vf_assay_molecules_min, variant_percent, variant_percent_max, 
           variant_percent_min) %>%
    dplyr::rename(sample = r_number) %>%
    mutate(sample_type = "cfDNA"))


cfDNA_z_scores <- left_join(parents_and_cfDNA %>%
  filter(sample_type == "cfDNA"),
  parents_and_cfDNA %>%
    filter(sample_type == "gDNA")%>%
    dplyr::group_by(vf_assay) %>%
    dplyr::summarise(mean_vp = mean(variant_percent),
                     stand_dev_vp = sd(variant_percent)),
  by = "vf_assay") %>%
  mutate(z_score = (variant_percent - mean_vp) / stand_dev_vp)

plot_title <- expression(paste(italic("IDS"), " c.182_189del ddPCR"))

ggplot(parents_and_cfDNA %>%
         filter(vf_assay == "IDS c.182_189del"), aes(
           x = vf_assay_molecules,
          y = variant_percent)) +
  geom_errorbar(aes(ymin = variant_percent_min, 
                    ymax = variant_percent_max), alpha = 0.2) +
  geom_errorbarh(aes(xmin = vf_assay_molecules_min, 
                    xmax = vf_assay_molecules_max), alpha = 0.2) +
  geom_point(size = 3, aes(shape = sample_type)) +
  scale_shape_manual(values = c(19, 1)) +
  ylim(40, 60) +
  xlim(0, 10000) +
  theme_bw() +
  theme(legend.title = element_blank(),
        panel.grid = element_blank()) +
  geom_hline(yintercept = 50, linetype = "dashed") +
  theme(legend.position = "bottom") +
  labs(x = "Molecules detected",
       y = "Variant percent (%)",
       title = plot_title)

#########################
# IDS incorrect result
#########################

# Sample was at 0.788ng/ul
# Worksheet 21-0588. 7 reps for IDS assay, 4 reps for the ZFXY assay.

IDS_only <- ddpcr_sprt_analysed %>%
  filter(r_number == "19611") %>%
  mutate(IDS_expected = (0.788*1000/3.3)*(5*7),
         ZFXY_expected = (0.788*1000/3.3)*(5*4),
         amplifiability_IDS = (vf_assay_molecules/IDS_expected)*100,
         amplifiability_ZFXY = (ff_assay_molecules/ZFXY_expected)*100) %>%
  select(vf_assay_molecules, IDS_expected, ff_assay_molecules,
         ZFXY_expected, amplifiability_IDS, amplifiability_ZFXY)
