################################################################################
## ddPCR NIPT Analysis for All Samples
## October 2021
## Joseph.Shaw@gosh.nhs.uk
## This analysis script performs 3 different analyses on the entire cohort
## of samples tested using ddPCR: the sequential probability ratio test (SPRT)
## (Lo et al, 2008; PMID:  17664418), MonteCarlo Markov Chain 
## (MCMC) analysis (Caswell et al, 2020; PMID:  32533152) and z score
## analysis (Chiu et al, 2008; PMID:  19073917).
## The MCMC analysis section  is compiled from scripts and
## models supplied by Tristan Snowsill (Exeter).
################################################################################

#########################
# Load packages, functions and data
#########################

# Source functions
source("functions/ddPCR_nipd_functions.R")

# Load ddPCR data
source("scripts/load_ddpcr_data.R")

# epiR for sensitivity calculations
library(epiR)

#ggpubr for plot collation
library(ggpubr)

# cmdstanr required for running stan models
library(cmdstanr)

# Source RAPID biobank
source("W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/RAPID_project_biobank/scripts/RAPID_biobank.R")

# Load ddPCR SNP panel (from Camunas Soler et al 2018)
ddpcr_snp_panel <- read.csv("W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/ddPCR SNP Panel/Camunas_Soler_panel_47_GnomAD_frequencies.csv")

gene_info <- read.csv("resources/vf_assay_gene_information.csv")

#########################
# SPRT analysis
#########################

# Set likelihood ratio threshold
lr_threshold <- 8

all_samples_sprt <- cfdna_ddpcr_data_molecules %>%
  mutate(
    
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
# MCMC analysis
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

all_samples_mcmc <- ddpcr_with_fits %>%
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
# Z score analysis
#########################
###################
# Heterozygous gDNA cohort
###################

# Upper limit of number of cfDNA molecules per case
cfdna_molecules_max <- max(cfdna_ddpcr_data_molecules$vf_assay_molecules)

# Lower limit for analysis 

vf_assay_molecules_limit <- 2000

het_controls <- controls %>%
  filter(identity %in% c("maternal gDNA", "paternal gDNA") &
           zygosity == "heterozygous")

# Select heterozygous gDNA data
het_gdna <- parent_gDNA_var_ref %>%
  filter(sample %in% het_controls$sample &
           # Remove empty well 
           reference_positives != 0) %>%
  
  # New column to distinguish control gDNA
  mutate(sample_type = "gDNA") %>%
  dplyr::rename(r_number = sample)

het_gdna_range <- het_gdna %>%
  # Filter the gDNA dataset to a similar range to the cfDNA dataset
  filter(vf_assay_molecules < (cfdna_molecules_max +1000) &
           vf_assay_molecules > vf_assay_molecules_limit)

###################
# Heterozygous gDNA variation
###################

# vp is variant percent
het_gDNA_mean_vp <- mean(het_gdna_range$variant_percent)            

# sd is standard deviation

het_gDNA_sd_vp <- sd(het_gdna_range$variant_percent)

# cv is coefficient of variation
het_gDNA_cv_vp <- (het_gDNA_sd_vp/het_gDNA_mean_vp) * 100

###################
# cfDNA z score analysis
###################

z_score_imbalance_threshold <- 3
z_score_balance_threshold <- 2

all_samples_zscore <- cfdna_ddpcr_data_molecules %>%
  mutate(z_score = (variant_percent - het_gDNA_mean_vp) / het_gDNA_sd_vp,
         
         z_score_prediction = case_when(
           
           # Exclude cases with low cffDNA or low cfDNA
           fetal_percent < 4 ~"inconclusive",
           vf_assay_molecules < vf_assay_molecules_limit ~ "inconclusive",
           
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
             z_score < -z_score_balance_threshold ~ "inconclusive"))

#########################
# Collate analyses 
#########################

# Samples to exclude
# 13262 - this sample had contamination
# 17004 - father was HbAC
# 20915 - twin pregnancy
samples_to_exclude <- c("13262", "20915", "17004")

all_samples_blinded <- left_join(
  all_samples_sprt,
  all_samples_mcmc %>%
    select(r_number, p_G0, p_G1, p_G2, p_G3, mcmc_prediction),
  by = "r_number") %>%
  left_join(
    all_samples_zscore %>%
      select(r_number, z_score, z_score_prediction),
    by = "r_number") %>%
  filter(!r_number %in% samples_to_exclude)

#########################
# Compare predictions against Biobank
#########################

all_samples_unblinded <- all_samples_blinded %>%
  left_join(RAPID_biobank %>%
            # Change r_number to a character
            mutate(r_number = as.character(r_number)) %>%
            select(r_number, study_id, site, gestation_weeks, gestation_days, 
                   gestation_total_weeks, gestation_character,
                   hours_to_first_spin, days_to_storage, vacutainer, 
                   mutation_genetic_info_fetus, 
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
      TRUE ~"incorrect"),
    
    outcome_zscore = case_when(
      z_score_prediction %in% c("inconclusive", "low fetal fraction")
      ~"inconclusive",
      z_score_prediction == fetal_genotype
      ~"correct",
      TRUE ~"incorrect")) %>%
  select(-mutation_genetic_info_fetus) %>%
  
  mutate(fetal_genotype = paste0(fetal_genotype, " fetus"),
         
         # Factorise for plotting
         fetal_genotype = factor(fetal_genotype,
                                     levels = c("hemizygous variant fetus",
                                                "homozygous variant fetus",
                                                "heterozygous fetus",
                                                "homozygous reference fetus",
                                                "hemizygous reference fetus")),
         outcome_zscore = factor(outcome_zscore,
                                 levels = c("correct", "incorrect", 
                                            "inconclusive")),
         outcome_mcmc = factor(outcome_mcmc,
                               levels = c("correct", "incorrect", 
                                          "inconclusive")),
         outcome_sprt = factor(outcome_sprt,
                               levels = c("correct", "incorrect", 
                                          "inconclusive")))

#########################
# Summary of all results
#########################

all_samples_arranged <- all_samples_unblinded %>%
  mutate(cohort = ifelse(vf_assay == "HBB c.20A>T", "sickle cell disease",
                       "bespoke design"),
         cohort = factor(cohort, levels = c("sickle cell disease",
                                          "bespoke design"))) %>%
  arrange(cohort, inheritance_chromosomal, inheritance_pattern, study_id)
  
families <- data.frame(
  study_id = unique(all_samples_arranged$study_id))

families <- mutate(families, 
                  family_number = rownames(families))

supplementary_table <- families %>%
  full_join(all_samples_arranged,
            by = "study_id") %>%
  mutate(sample_id = 
         paste0("cfDNA-", as.character(row.names(all_samples_arranged))),
         
         # Specify the method used to determine the fetal fraction
         ff_determination = case_when(
           inheritance_chromosomal == "x_linked" ~ "ZFXY",
           inheritance_chromosomal == "autosomal" & ff_assay %in%
             ddpcr_snp_panel$dbSNP ~ "ddPCR SNP panel - workflow 2",
           TRUE ~"NGS SNP panel - workflow 1")) %>%
  
  # Reorder columns
  select(
    # Sample identifiers
    sample_id, r_number, family_number, study_id, 
    # Sampling information
    hours_to_first_spin, days_to_storage,
    vacutainer, gestation_character,
    diagnostic_sampling, 
    # Extraction information
    plasma_volume_ml, extraction_replicates, 
    # Variant information
    cohort, inheritance_chromosomal, inheritance_pattern, vf_assay,
    ff_determination, ff_assay, vf_assay_num_wells, 
    vf_assay_droplets, reference_positives,
    variant_positives, ff_assay_num_wells, ff_assay_droplets, 
    maternal_positives, paternal_positives,
    # Molecules
    variant_molecules, reference_molecules, vf_assay_molecules,
    variant_percent,
    maternal_molecules, paternal_molecules, fetal_percent,
    ff_assay_molecules, total_molecules, totalGE_ml_plasma,
    # Analysis
    likelihood_ratio, sprt_prediction, p_G0, p_G1, p_G2,
    p_G3, mcmc_prediction, z_score, z_score_prediction,
    # Outcome
    fetal_genotype, outcome_sprt, outcome_mcmc, outcome_zscore)

# Export the results
write.csv(supplementary_table, 
          file = (paste0("analysis_outputs/all_samples_ddpcr ",
               format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")),
          row.names = FALSE)

###################
# Gene information table
###################

xl_count_table <- all_samples_arranged %>%
  filter(inheritance_chromosomal == "x_linked") %>%
  count(vf_assay)

ad_count_table <- all_samples_arranged %>%
  filter(inheritance_chromosomal == "autosomal" &
           inheritance_pattern == "dominant") %>%
  count(vf_assay)

ar_count_table <- all_samples_arranged %>%
  filter(inheritance_chromosomal == "autosomal" &
           inheritance_pattern == "recessive") %>%
  count(vf_assay)

vf_assay_count_table <- rbind(xl_count_table, 
                              ad_count_table,
                              ar_count_table) %>%
  dplyr::rename(samples = n) %>%
  
  # Bind to gene information
  left_join(gene_info,
            by = "vf_assay") %>%
  
  # Add inheritance information
  left_join(
    ddpcr_target_panel %>%
      select(assay, inheritance_chromosomal, inheritance_pattern) %>%
      #Remove duplicate rows
      distinct(.keep_all = TRUE) %>%
      dplyr::rename(vf_assay = assay),
    by = "vf_assay") %>%
  
  # Add on abbreviation of inheritance patterns
  mutate(inheritance_abbreviation = case_when(
    inheritance_chromosomal == "autosomal" &
      inheritance_pattern == "recessive" ~"AR",
    inheritance_chromosomal == "autosomal" &
      inheritance_pattern == "dominant" ~"AD",
    inheritance_chromosomal == "x_linked" &
      inheritance_pattern == "dominant" ~"XLD",
    inheritance_chromosomal == "x_linked" &
      inheritance_pattern == "recessive" ~"XLR"),
    
    # Factorise for ordering
    inheritance_abbreviation = factor(inheritance_abbreviation,
                                      levels = c("AR", "AD", "XLR", "XLD"))) %>%
  
  arrange(inheritance_abbreviation, vf_assay) %>%
  
  select(inheritance_abbreviation, gene, transcript, variant_dna, 
         variant_protein, condition,
         samples) %>%
  
  dplyr::rename(Inheritance = inheritance_abbreviation,
                Gene = gene,
                Transcript = transcript,
                DNA = variant_dna,
                Protein = variant_protein,
                Condition = condition,
                Samples = samples)

write.csv(vf_assay_count_table, 
          "analysis_outputs/vf_assay_count_table.csv",
          row.names = FALSE)

###################
# Sensitivity and specificity 4x4 tables
###################

sensitivity_scores <- function(df, prediction, outcome) {
  
  unbalanced_genotypes <- c("homozygous variant",
                            "hemizygous variant",
                            "homozygous reference",
                            "hemizygous reference")
  
  true_positives <- nrow(df %>%
    filter(!!prediction %in% unbalanced_genotypes &
             !!outcome == "correct"))
  
  true_negatives <- nrow(df %>%
    filter(!!prediction == "heterozygous" &
             !!outcome == "correct"))
  
  false_positives <- nrow(df %>%
    filter(!!prediction %in% unbalanced_genotypes &
             !!outcome == "incorrect"))
  
  false_negatives <- nrow(df %>%
    filter(!!prediction == "heterozygous" &
             !!outcome == "incorrect"))
  
  inconclusives <- nrow(df %>%
    filter(!!prediction == "inconclusive"))
  
  output <- data.frame(
    "analysis" = c(as_label(prediction)),
    "true_positives" = c(true_positives),
    "true_negatives" = c(true_negatives),
    "false_positives" = c(false_positives),
    "false_negatives" = c(false_negatives),
    "inconclusive" = c(inconclusives))
  
  return(output)
}

sensitivity_table <- rbind(
  # SPRT results
  sensitivity_scores(df = supplementary_table,
                   prediction = quo(sprt_prediction), 
                   outcome = quo(outcome_sprt)),
  # MCMC results
  sensitivity_scores(df = supplementary_table,
                     prediction = quo(mcmc_prediction), 
                     outcome = quo(outcome_mcmc)),
  # Z score results
  sensitivity_scores(df = supplementary_table,
                     prediction = quo(z_score_prediction), 
                     outcome = quo(outcome_zscore))) %>%
  
  # Calculate sensitivity and specificity
  mutate(
    sensitivity = round((true_positives / (true_positives + false_negatives))*100, 1),
    specificity = round((true_negatives / (false_positives + true_negatives))*100, 1))


test_table <- supplementary_table %>%
  mutate(sprt_binary = case_when(
    # Autosomal dominant inheritance
    inheritance_chromosomal == "autosomal" &
      inheritance_pattern == "dominant" &
      sprt_prediction == "heterozygous" ~"positive",
    inheritance_chromosomal == "autosomal" &
      inheritance_pattern == "dominant" &
      sprt_prediction == "homozygous reference" ~"negative",
    
    # Autosomal recessive inheritance
    inheritance_chromosomal == "autosomal" &
      inheritance_pattern == "recessive" &
      sprt_prediction %in% c("homozygous variant", 
                             "homozygous reference")  ~"positive",
    inheritance_chromosomal == "autosomal" &
      inheritance_pattern == "recessive" &
      sprt_prediction == "heterozygous" ~"negative",
    
    # X linked inheritance
    inheritance_chromosomal == "x_linked" &
      sprt_prediction == "hemizygous variant"  ~"positive",
    sprt_prediction == "hemizygous reference"  ~"negative",
    TRUE ~"inconclusive")) %>%
  select(r_number, cohort, inheritance_chromosomal, 
         inheritance_pattern, sprt_prediction, sprt_binary,
         outcome_sprt)

count(test_table, sprt_binary, outcome_sprt) %>%
  mutate(cohort = "all") %>%
  pivot_wider(id_cols = c(cohort),
              names_from = c(sprt_binary, outcome_sprt),
              values_from = n) %>%
  dplyr::rename(true_negative = negative_correct,
                false_negative = negative_incorrect,
                true_positive = positive_correct,
                false_positive = positive_incorrect,
                inconclusive = inconclusive_inconclusive) %>%
  select(cohort, true_positive, true_negative, false_positive, 
         false_negative, inconclusive) %>%
  mutate(
    sensitivity = round((true_positive / (true_positive + false_negative))*100, 1),
    specificity = round((true_negative / (false_positive + true_negative))*100, 1)) %>%
  pivot_longer(cols = c(-cohort),
               names_to = "category",
                 values_to = "sprt")

###################
# Text for paper
###################

paste0("Overall, we analysed ", nrow(all_samples_unblinded), " samples from ", 
       nrow(families),
       " families, for ", length(unique(all_samples_unblinded$vf_assay)), 
       " pathogenic variants.")

paste0("Samples were collected from between ", 
       min(all_samples_unblinded$gestation_weeks),
       " and ", max(all_samples_unblinded$gestation_weeks),
       " weeks gestation (median: ", 
       round(median(all_samples_unblinded$gestation_total_weeks), 1),
       ") with fetal fractions ranging from ",
       round(min(all_samples_unblinded$fetal_percent), 1), "-", 
       round(max(all_samples_unblinded$fetal_percent), 1),
       "% (median: ", round(median(all_samples_unblinded$fetal_percent), 1),
       "%).")


unique(all_samples_unblinded$fetal_genotype)


paste0("When compared to the results of invasive testing, SPRT and MCMC analysis generated ",
       nrow(all_samples_unblinded %>%
              filter(vf_assay == "HBB c.20A>T" &
                       outcome_sprt == "correct")),
       " and ",
       nrow(all_samples_unblinded %>%
              filter(vf_assay == "HBB c.20A>T" &
                       outcome_mcmc == "correct")),
       " correct fetal genotype predictions, respectively, including ",
       nrow(all_samples_unblinded %>%
              filter(vf_assay == "HBB c.20A>T" &
                       fetal_genotype == "homozygous variant fetus",
                       outcome_sprt == "correct")),
       " and ",
       nrow(all_samples_unblinded %>%
              filter(vf_assay == "HBB c.20A>T" &
                       fetal_genotype == "homozygous variant fetus",
                     outcome_mcmc == "correct")),
       " fetuses affected with sickle cell disease. However, each method also generated ",
       nrow(all_samples_unblinded %>%
              filter(vf_assay == "HBB c.20A>T" &
                       outcome_sprt == "incorrect")),
       " and ",
       nrow(all_samples_unblinded %>%
              filter(vf_assay == "HBB c.20A>T" &
                       outcome_mcmc == "incorrect")),
       " incorrect fetal gentoype predictions.")


nrow(all_samples_unblinded %>%
       filter(vf_assay != "HBB c.20A>T" &
                outcome_sprt == "correct"))

nrow(all_samples_unblinded %>%
       filter(vf_assay != "HBB c.20A>T" &
                outcome_mcmc == "correct"))

nrow(all_samples_unblinded %>%
       filter(vf_assay != "HBB c.20A>T" &
                outcome_sprt == "incorrect"))

nrow(all_samples_unblinded %>%
       filter(vf_assay != "HBB c.20A>T" &
                outcome_mcmc == "incorrect"))

all_samples_unblinded %>%
  filter(outcome_mcmc == "incorrect" |
           outcome_sprt == "incorrect") %>%
  select(r_number, z_score_prediction, outcome_zscore)

#########################
# Plot results
#########################

# Consistent theme and axes for plots
multiplot_theme <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), 
  legend.position = "none", 
  plot.title = element_text(size = 11),
  legend.title = element_blank(),
  legend.text = element_text(size = 9))

multiplot_y <- ylim(38, 62)

multiplot_x <- scale_x_continuous(limits = c(0,32000),
                                  breaks = c(0, 10000, 20000 ,30000))

## Z score lines

z3_line <- geom_hline(yintercept = het_gDNA_mean_vp+(3*het_gDNA_sd_vp),
                      linetype = "dashed", alpha = 0.5) 

zminus3_line <- geom_hline(yintercept = het_gDNA_mean_vp-(3*het_gDNA_sd_vp),
                           linetype = "dashed", alpha = 0.5)

z2_line <- geom_hline(yintercept = het_gDNA_mean_vp+(2*het_gDNA_sd_vp), 
                      linetype = "dashed", alpha = 0.5)

zminus2_line <- geom_hline(yintercept = het_gDNA_mean_vp-(2*het_gDNA_sd_vp),
                           linetype = "dashed", alpha = 0.5)

vertical_line <- geom_vline(xintercept = vf_assay_molecules_limit,
           linetype = "dashed", alpha = 0.5)

cfdna_fill <- scale_fill_manual(values=c("#FFFFFF", "#000000", "#999999"), 
                  guide = "none")

cfdna_alpha <- scale_alpha_manual(values = c(1, 1, 0.2), guide = "none")

cfdna_shape <- scale_shape_manual(values = c(24, 24, 21, 25, 25))

###################
# Plot 1: heterozygous gDNA controls
###################

plot_1 <- ggplot(het_gdna %>%
         mutate(sample_type = "het gDNA"), 
                    aes(x = vf_assay_molecules, 
                        y = variant_percent)) +
  z3_line +
  zminus3_line +
  z2_line +
  zminus2_line +
  vertical_line +
  multiplot_y +
  multiplot_x +
  
  scale_shape_manual(values = c(21)) +
  geom_point(size = 2, aes(shape = sample_type), fill= "white",
             colour = "black") +
  theme_bw() +
  multiplot_theme +
  labs(y = "Variant fraction (%)", x = "",
       title = "ddPCR for 82 heterozygous gDNA controls")


ggplot(het_gdna %>%
         filter(vf_assay == "HBB c.20A>T") %>%
                   mutate(sample_type = "het gDNA"), 
                 aes(x = vf_assay_molecules, 
                     y = variant_percent)) +
  geom_hline(yintercept = 50,
             linetype = "dashed", alpha = 0.5) +
  ylim(43, 57)+
  scale_x_continuous(limits = c(0,30000),
                     breaks = c(0, 10000, 20000 ,30000)) +
  scale_shape_manual(values = c(21)) +
  geom_point(size = 3, aes(shape = sample_type), fill= "white",
             colour = "black", alpha =0.7) +
  theme_bw() +
  multiplot_theme +
  labs(y = "Variant fraction (%)", x = "DNA input (molecules)",
       title = "ddPCR for 42 HBB c.20A>T gDNA controls")

###################
# Plot 2: z score analysis
###################

plot_2 <- ggplot(all_samples_unblinded, aes(x = vf_assay_molecules, 
             y = variant_percent)) +
  theme_bw() +
  multiplot_theme + 
  z3_line +
  zminus3_line +
  z2_line +
  zminus2_line +
  vertical_line +
  multiplot_y +
  multiplot_x +
  cfdna_fill +
  cfdna_alpha +
  cfdna_shape +
  geom_point(size = 2, aes(fill = outcome_zscore,
                           alpha = outcome_zscore,
                           shape = fetal_genotype),
             colour = "black") +
  
  labs(y = "Variant fraction (%)", x = "Genome equivalents (GE)", 
       title = "ddPCR for 124 cfDNA samples with z score classification")

###################
# Plot 3: SPRT analysis
###################

plot_3 <- ggplot(all_samples_unblinded, aes(x = vf_assay_molecules, 
                                  y = variant_percent)) +
  theme_bw() +
  multiplot_theme + 
  #z3_line +
  #zminus3_line +
  #z2_line +
  #zminus2_line +
  #vertical_line +
  multiplot_y +
  multiplot_x +
  cfdna_fill +
  cfdna_alpha +
  cfdna_shape +
  geom_point(size = 2, aes(fill = outcome_sprt,
                           alpha = outcome_sprt,
                           shape = fetal_genotype),
             colour = "black") +
  
  labs(y = "Variant fraction (%)", x = "Genome equivalents (GE)", 
       title = "ddPCR for 124 cfDNA samples with SPRT classification")

###################
# Plot 4: MCMC analysis
###################

plot_4 <- ggplot(all_samples_unblinded, aes(x = vf_assay_molecules, 
                                  y = variant_percent)) +
  theme_bw() +
  multiplot_theme + 
  #z3_line +
  #zminus3_line +
  #z2_line +
  #zminus2_line +
  #vertical_line +
  multiplot_y +
  multiplot_x +
  cfdna_fill +
  cfdna_alpha +
  cfdna_shape +
  geom_point(size = 2, aes(fill = outcome_mcmc,
                           alpha = outcome_mcmc,
                           shape = fetal_genotype),
             colour = "black") +
  
  labs(y = "Variant fraction (%)", x = "Genome equivalents (GE)", 
       title = "ddPCR for 124 cfDNA samples with MCMC classification")

###################
# All plots together
###################

ddpcr_cohort <- ggpubr::ggarrange(plot_1, plot_2, 
                                plot_3, plot_4,
                                ncol = 2, nrow = 2, align = "v",
                                labels = c("A", "B", "C", "D"))

ggsave(plot = ddpcr_cohort, 
       filename = "ddpcr_cohort.tiff",
       path = "plots/", device='tiff', dpi=600,
       units = "in",
       width = 12.5,
       height = 7)

###################
# Analysis summary plot
###################

# Summary of how each analysis method performed.

analysis_summary <- rbind(
  all_samples_unblinded %>%
  count(outcome_zscore) %>%
  dplyr::rename(outcome = outcome_zscore) %>%
  mutate(analysis = "z score"),
  all_samples_unblinded %>%
    count(outcome_sprt) %>%
    dplyr::rename(outcome = outcome_sprt) %>%
    mutate(analysis = "sprt"),
  all_samples_unblinded %>%
    count(outcome_mcmc) %>%
    dplyr::rename(outcome = outcome_mcmc) %>%
    mutate(analysis = "mcmc")) %>%
  mutate(outcome = factor(outcome, 
                          levels = c("correct",
                                     "inconclusive",
                                     "incorrect")),
         analysis = factor(analysis,
                           levels = c("sprt", "mcmc", "z score")))

# Bar plot
ggplot(analysis_summary, aes(x = outcome, y = n))+
  geom_col(fill = "white", colour = "black") +
  facet_wrap(~analysis) +
  ylim(0, 120) +
  labs(y = "Samples", x = "",
       title = "") +
  theme_bw() +
  geom_text(aes(x = outcome, y = n, label = n), vjust = -1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

###################
# ROC curve analysis
###################

roc_binary_calls <- all_samples_unblinded %>%
  # Convert invasive results to binary outcomes
  mutate(unbalanced = case_when(
    fetal_genotype %in% c("homozygous reference fetus", 
                          "hemizygous variant fetus",
                          "homozygous variant fetus",
                          "hemizygous reference fetus") ~"TRUE",
    fetal_genotype == "heterozygous fetus" ~"FALSE"),
    
    # Convert "unbalanced" column to Boolean vector
    unbalanced = as.logical(unbalanced),
    # Convert the MCMC calls to a binary outcome
    # Annoyingly, pG2 has different meanings for recessive and dominant
    # cases.
    mcmc_unbalanced_call = case_when(
      inheritance_chromosomal == "autosomal" &
        inheritance_pattern == "recessive" ~pmax(p_G1, p_G3),
      inheritance_chromosomal == "autosomal" &
        inheritance_pattern == "dominant" ~p_G2,
      inheritance_chromosomal == "x_linked" ~pmax(p_G0, p_G1)),
    # Remove minus signs from z score
    zscore_unbalanced_call = abs(z_score))

sprt_roc <- roc_binary_calls %>%
  arrange(desc(likelihood_ratio)) %>%
  # tpr is "true positive rate" and fpr is "false positive rate"
  mutate(sprt_tpr = cumsum(unbalanced)/sum(unbalanced)) %>%
  mutate(sprt_fpr = cumsum(!unbalanced)/sum(!unbalanced)) %>%
  select(r_number, likelihood_ratio, sprt_tpr, sprt_fpr, unbalanced)

mcmc_roc <- roc_binary_calls %>%
  arrange(desc(mcmc_unbalanced_call)) %>%
  mutate(mcmc_tpr = cumsum(unbalanced)/sum(unbalanced)) %>%
  mutate(mcmc_fpr = cumsum(!unbalanced)/sum(!unbalanced)) %>%
  select(r_number, mcmc_unbalanced_call, mcmc_tpr, mcmc_fpr)

zscore_roc <- roc_binary_calls %>%
  arrange(desc(zscore_unbalanced_call)) %>%
  mutate(zscore_tpr = cumsum(unbalanced)/sum(unbalanced)) %>%
  mutate(zscore_fpr = cumsum(!unbalanced)/sum(!unbalanced)) %>%
  select(r_number, zscore_tpr, zscore_fpr)

sprt_roc_plot <- ggplot(sprt_roc, aes(x = sprt_fpr, y = sprt_tpr))+
  geom_line(size = 2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(linetype = "dashed")+
  ylim(0, 1)+
  xlim(0,1)+
  labs(x = "",
       y = "True positive rate", 
       title = "SPRT analysis")

mcmc_roc_plot <- ggplot(mcmc_roc, aes(x = mcmc_fpr, y = mcmc_tpr))+
  geom_line(size = 2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(linetype = "dashed")+
  ylim(0, 1)+
  xlim(0,1)+
  labs(x = "False positive rate",
       y = "", 
       title = "MCMC analysis")

zscore_roc_plot <- ggplot(zscore_roc, aes(x = zscore_fpr, y = zscore_tpr))+
  geom_line(size = 2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(linetype = "dashed")+
  ylim(0, 1)+
  xlim(0,1)+
  labs(x = "",
       y = "", 
       title = "Z score analysis")

# All plots together
roc_plot <- ggpubr::ggarrange(sprt_roc_plot, mcmc_roc_plot, 
                  zscore_roc_plot, 
                  ncol = 3, nrow = 1)

ggsave(plot = roc_plot, 
       filename = "roc_plot.tiff",
       path = "plots/", device='tiff')


unique(all_samples_unblinded$site)

###################
# SPRT on gDNA controls
###################

# This plot shows that the common form of the SPRT equation is 
# inappropriate for this dataset

ff_for_graph <- 0.04
lr_for_graph <- 8

het_dna_sprt <- ggplot(het_gdna, aes(x = vf_assay_molecules, y = variant_percent))+
  geom_point(size = 2, colour = "black", fill = "white", pch=21,
             alpha = 0.8) +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none")+
  labs(x = "Genome equivalents (GE)",
       y = "Variant fraction (%)") +
  geom_function(fun = "calc_AS_upper_boundary",
                aes(x = vf_assay_molecules, y =
                      calc_AS_upper_boundary(vf_assay_molecules, ff_for_graph, 
                                             lr_for_graph)),
                colour = "black",
                args = c(ff_for_graph, lr_for_graph)) +
  
  geom_function(fun = "calc_AS_lower_boundary",
                aes(x = vf_assay_molecules, y =
                      calc_AS_lower_boundary(vf_assay_molecules, ff_for_graph, 
                                             lr_for_graph)),
                colour = "black",
                args = c(ff_for_graph, lr_for_graph)) +
  
  geom_function(fun = "calc_SS_boundary",
                aes(x = vf_assay_molecules, y =
                      calc_SS_boundary(vf_assay_molecules, ff_for_graph, 
                                       lr_for_graph)),
                colour = "black",
                args = c(ff_for_graph, lr_for_graph))+
  
  geom_function(fun = "calc_AA_boundary",
                aes(x = vf_assay_molecules, y =
                      calc_AA_boundary(vf_assay_molecules, ff_for_graph, 
                                       lr_for_graph)),
                colour = "black",
                args = c(ff_for_graph, lr_for_graph)) +
  scale_y_continuous(limits = c(43,57)) + 
  annotate(geom = "text", x = 28000, y = 50, label = "heterozygous") +
  annotate(geom = "text", x = 28000, y = 47, label = "homozygous reference") +
  annotate(geom = "text", x = 28000, y = 53, label = "homozygous variant") +
  xlim(0, 30000)

ggsave(plot = het_dna_sprt, 
       filename = "het_dna_sprt.tiff",
       path = "plots/", device='tiff', dpi=600,
       units = "in",
       width = 12.5,
       height = 7)

###################