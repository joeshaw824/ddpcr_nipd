################################################################################
## ddPCR for Non Invasive Prenatal Diagnosis (NIPD)
## November 2021
## Joseph.Shaw@gosh.nhs.uk
## This is an analysis script for the prediction of fetal genotypes from 
## cfDNA testing using ddPCR for maternally-inherited
## variants. This includes modules for analysis of ddPCR data
## using the sequential probability ratio test (SPRT) (Lo et al,
## 2008; PMID:  17664418), MonteCarlo Markov Chain 
## (MCMC) analysis (Caswell et al, 2020; PMID:  32533152) and z score
## analysis (Chiu et al, 2008; PMID:  19073917).
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
source("W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/RAPID_project_biobank/scripts/RAPID_biobank.R")

RAPID_biobank <- load_biobank()

#########################
# cfDNA SPRT analysis
#########################

# Set likelihood ratio threshold
lr_threshold <- 8

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
# (12/10/2021: this takes approximately 14.5 minutes)
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
# Collation of SPRT and MCMC results for all samples
#########################

ddpcr_sprt_mcmc <- left_join(
  ddpcr_sprt_analysed,
  ddpcr_mcmc_analysed %>%
    select(r_number, p_G0, p_G1, p_G2, p_G3, mcmc_prediction),
  by = "r_number") %>%
  left_join(RAPID_biobank %>%
              # Change r_number to a character to match ddpcr_sprt_analysed
              mutate(r_number = as.character(r_number)) %>%
              filter(r_number %in% ddpcr_sprt_analysed$r_number) %>%
              select(r_number, study_id, gestation_weeks, gestation_days, 
                     gestation_total_weeks, gestation_character, 
                     date_of_blood_sample, vacutainer, mutation_genetic_info_fetus, 
                     partner_sample_available, original_plasma_vol, 
                     tubes_plasma_current, report_acquired),
            by = "r_number") %>%
  
  # Change the nomenclature of the outcomes for sickle cell disease
  mutate(fetal_genotype = case_when(
    mutation_genetic_info_fetus == "HbSS" ~"homozygous variant",
    mutation_genetic_info_fetus == "HbAS" ~"heterozygous",
    mutation_genetic_info_fetus == "HbAA" ~"homozygous reference",
    TRUE ~mutation_genetic_info_fetus),
    
    # Compare predictions to invasive testing
    outcome_sprt = case_when(
      sprt_prediction == "inconclusive" ~"inconclusive",
      sprt_prediction == fetal_genotype
      ~"correct",
      TRUE ~"incorrect"),
    
    outcome_mcmc = case_when(
      mcmc_prediction == "inconclusive" ~"inconclusive",
      mcmc_prediction == fetal_genotype
      ~"correct",
      TRUE ~"incorrect"))

write.csv(ddpcr_sprt_mcmc, 
          file = (paste0("analysis_outputs/total_cohort_sprt_mcmc",
                         format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")),
          row.names = FALSE)

#########################
# Select bespoke case parental data for z score analysis
#########################

bespoke_gdna <- parent_gDNA_var_ref %>%
  filter(sample %in% controls$sample & 
           # Remove sickle cell disease samples
           vf_assay != "HBB c.20A>T") %>%
  mutate(sample_type = "het gDNA") %>%
  filter(variant_percent > 40 &
           vf_assay_molecules < 25000)

# vp is "variant percent"
bespoke_gdna_mean_vp <- mean(bespoke_gdna$variant_percent)            

# Find standard deviation of variant percent
bespoke_gdna_sd_vp <- sd(bespoke_gdna$variant_percent)

# Coefficient of variation
bespoke_gdna_cv_vp <- (bespoke_gdna_sd_vp/bespoke_gdna_mean_vp) * 100

#########################
# cfDNA z score analysis
#########################

z_score_imbalance_threshold <- 3
z_score_balance_threshold <- 2

bespoke_z_score <- ddpcr_sprt_mcmc %>%
  filter(vf_assay != "HBB c.20A>T"  &
         # Remove 30116 as we don't have an outcome (yet)
         r_number != "30116") %>%
  mutate(z_score = (variant_percent-bespoke_gdna_mean_vp) / bespoke_gdna_sd_vp,
         
         z_score_prediction = case_when(
           
           # Exclude cases with low cffDNA
           fetal_percent < 4 ~"low fetal fraction",
           
           # Predict fetal genotypes for X-linked cases
           inheritance_chromosomal == "x_linked" &
             z_score > z_score_imbalance_threshold ~ "hemizygous variant",
           
           inheritance_chromosomal == "x_linked" &
             z_score < -z_score_imbalance_threshold ~ "hemizygous reference",
           
           inheritance_chromosomal == "x_linked" &
             z_score > -z_score_imbalance_threshold &
             z_score < z_score_imbalance_threshold ~ "inconclusive",
           
           # Predict fetal genotypes for autosomal cases
           inheritance_chromosomal == "autosomal" &
             z_score > z_score_imbalance_threshold ~ "homozygous variant",
           
           inheritance_chromosomal == "autosomal" &
             z_score < -z_score_imbalance_threshold ~ "homozygous reference",
           
           inheritance_chromosomal == "autosomal" &
             z_score < z_score_balance_threshold &
             z_score > -z_score_balance_threshold ~ "heterozygous",
           
           inheritance_chromosomal == "autosomal" &
             z_score < z_score_imbalance_threshold &
             z_score > z_score_balance_threshold ~ "inconclusive",
           
           inheritance_chromosomal == "autosomal" &
             z_score > -z_score_imbalance_threshold &
             z_score < -z_score_balance_threshold ~ "inconclusive"),
         
         # Compare z score results to invasive testing
         
         outcome_zscore = case_when(
           z_score_prediction == "inconclusive" ~"inconclusive",
           z_score_prediction == "low fetal fraction" ~"inconclusive",
           z_score_prediction == mutation_genetic_info_fetus
           ~"correct",
           TRUE ~"incorrect"),
         
         outcome_zscore = factor(outcome_zscore, levels = c(
           "correct", "incorrect", "inconclusive")),
         
         mutation_genetic_info_fetus = 
           paste0(mutation_genetic_info_fetus, " fetus"),
         
         mutation_genetic_info_fetus = factor(mutation_genetic_info_fetus,
                                    levels = c("hemizygous variant fetus",
                                               "homozygous variant fetus",
                                               "heterozygous fetus",
                                               "homozygous reference fetus",
                                               "hemizygous reference fetus")))

#########################
# Individual plots for bespoke cases
#########################

# "Bespoke" is any assay other than the HBB c.20A>T assay for sickle cell
# disease.

bespoke_cases <- ddpcr_sprt_analysed %>%
  # Apply quality control filters to select only cases with 
  # sufficient total cfDNA (over 3000 molecules detected)
  filter(vf_assay != "HBB c.20A>T" &
           vf_assay_molecules >= 3000)

bespoke_wells <- sample_wells %>%
  filter(cfdna_sample %in% bespoke_cases$r_number) %>%
  arrange(cfdna_sample)

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
ggexport(plotlist = bespoke_plots, filename = "plots/bespoke_cohort_qc.pdf",
         width=15, height=8, res=300)

#########################
# cfDNA z score plots
#########################

# Consistent theme and axes
multiplot_theme <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.title = element_text(size = 11),
  axis.text.x = element_text(size = 9),
  legend.title = element_blank(),
  legend.text = element_text(size = 9),
  legend.position = "bottom")

multiplot_y <- ylim(37, 63)

multiplot_x <- scale_x_continuous(limits = c(0,20000),
                                  breaks = c(0, 10000, 20000))

## Lines

z3_line <- geom_hline(yintercept = bespoke_gdna_mean_vp+(3*bespoke_gdna_sd_vp),
                      linetype = "dashed", alpha = 0.5) 

zminus3_line <- geom_hline(yintercept = bespoke_gdna_mean_vp-(3*bespoke_gdna_sd_vp),
                           linetype = "dashed", alpha = 0.5) 

z2_line <- geom_hline(yintercept = bespoke_gdna_mean_vp+(2*bespoke_gdna_sd_vp), 
                      linetype = "dashed", alpha = 0.5)

zminus2_line <- geom_hline(yintercept = bespoke_gdna_mean_vp-(2*bespoke_gdna_sd_vp),
                           linetype = "dashed", alpha = 0.5)

# Plot 1: heterozygous gDNA

plot1 <- ggplot(bespoke_gdna, aes(x = vf_assay_molecules, y = variant_percent))+
  scale_shape_manual(values = c(21)) +
  z3_line +
  zminus3_line +
  z2_line +
  zminus2_line +
  geom_point(size = 2, colour = "black", fill = "white", 
             aes(shape = sample_type),
             alpha = 0.8) +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none")+
  labs(x = "",
       y = "Variant fraction (%)") +
  theme_bw() +
  multiplot_theme +
  multiplot_y +
  multiplot_x +
  labs(x = "", y = "Variant fraction (%)",
       title = "Heterozygous parental gDNA") +
  guides(shape=guide_legend(ncol=1))

# Plot 2: autosomal cases cfDNA

plot2 <- ggplot(bespoke_z_score %>%
         filter(inheritance_chromosomal == "autosomal"),
       aes(x = vf_assay_molecules, 
               y = variant_percent)) +
  z3_line +
  zminus3_line +
  z2_line +
  zminus2_line +
  
  # No incorrect results, so black option removed
  scale_fill_manual(values=c("#FFFFFF", "#999999"), guide = "none") +
  scale_alpha_manual(values = c(1, 0.2), guide = "none") +
  scale_shape_manual(values = c(24, 21, 25)) +
  geom_point(size = 2, aes(fill = outcome_zscore,
                           alpha = outcome_zscore,
                           shape = mutation_genetic_info_fetus),
             colour = "black") +
  theme_bw() +
  multiplot_theme +
  multiplot_y +
  multiplot_x +
  labs(x = "Genome equivalents (GE)", y = "",
       title = "cfDNA: autosomal variant cases") +
  guides(shape=guide_legend(ncol=1))

# Plot 3: x-linked cases cfDNA

plot3 <- ggplot(bespoke_z_score  %>%
         filter(inheritance_chromosomal == "x_linked"),
       aes(x = vf_assay_molecules, 
           y = variant_percent)) +
  z3_line +
  zminus3_line +
  z2_line +
  zminus2_line +
  scale_fill_manual(values=c("#FFFFFF", "#000000", "#999999"), guide = "none") +
  scale_alpha_manual(values = c(1, 1, 0.2), guide = "none") +
  # No "heterozygous" genotype option
  scale_shape_manual(values = c(24, 25)) +
  geom_point(size = 2, aes(fill = outcome_zscore,
                           alpha = outcome_zscore,
                           shape = mutation_genetic_info_fetus),
             colour = "black") +
  theme_bw() +
  multiplot_theme +
  multiplot_y +
  multiplot_x +
  labs(x = "", y = "",
       title = "cfDNA: X-linked variant cases") +
  theme(legend.position = "bottom") +
  guides(shape=guide_legend(ncol=1))

triptych <- ggpubr::ggarrange(plot1, plot2, plot3,
                                ncol = 3, nrow = 1, align = "h")

ggsave(plot = triptych, 
       filename = "bespoke_cohort_triptych.tiff",
       path = "plots/", device='tiff', dpi=600,
       units = "in",
       width = 10,
       height = 5)

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

#########################



