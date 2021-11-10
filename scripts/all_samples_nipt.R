################################################################################
## ddPCR NIPT Analysis for All Samples
## October 2021
## Joseph.Shaw@gosh.nhs.uk
## This analysis script performs 3 different analyses on the entire cohort
## of samples tested using ddPCR: the sequential probability ratio test (SPRT)
## (Lo et al, 2008; PMID:  17664418), MonteCarlo Markov Chain 
## (MCMC) analysis (Caswell et al, 2020; PMID:  32533152) and z score
## analysis.
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

# pROC for plotting ROC curves
library(pROC)

# Source RAPID biobank
source("W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/RAPID_project_biobank/scripts/RAPID_biobank.R")

# Load ddPCR SNP panel (from Camunas Soler et al 2018)
ddpcr_snp_panel <- read.csv("W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/ddPCR SNP Panel/Camunas_Soler_panel_47_GnomAD_frequencies.csv")

gene_info <- read.csv("resources/vf_assay_gene_information.csv")

# Run the script to collate ddPCR assay primer and probe sequence information
source("scripts/ddpcr_vf_assays.R")

#########################
# Heterozygous gDNA cohort
#########################

# Upper limit of number of cfDNA molecules per case
cfdna_molecules_max <- max(cfdna_ddpcr_data_molecules$vf_assay_molecules)

# Lower limit for analysis (gDNA variation is high below this)
vf_assay_molecules_limit <- 2000

# Remove duplicates in ddPCR target panel
ddpcr_targets <- ddpcr_target_panel[!duplicated(ddpcr_target_panel$assay),]

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
  dplyr::rename(r_number = sample) %>%
  # Add inheritance patterns
  left_join(ddpcr_targets %>%
              dplyr::rename(vf_assay = assay) %>%
              select(vf_assay, inheritance_chromosomal, inheritance_pattern),
            by = "vf_assay")

het_gdna_range <- het_gdna %>%
  # Filter the gDNA dataset to a similar range to the cfDNA dataset
  filter(vf_assay_molecules < (cfdna_molecules_max +1000) &
           vf_assay_molecules > vf_assay_molecules_limit)

#########################
# SPRT analysis - cfDNA
#########################

all_samples_sprt <- cfdna_ddpcr_data_molecules %>%
  mutate(
    
    # Perform SPRT and return the likelihood ratio
    likelihood_ratio = case_when(
      inheritance_chromosomal == "x_linked" ~ 
        calc_lr_x_linked(fetal_fraction, (major_allele_percent/100), 
                         vf_assay_molecules),
      inheritance_chromosomal == "autosomal" ~ 
        calc_lr_autosomal(fetal_fraction, (major_allele_percent/100), 
                          vf_assay_molecules)))

# Predict genotypes with a likelihood ratio threshold of 8
all_samples_sprt <- predict_sprt_genotypes(all_samples_sprt, 8)

# Convert SPRT predictions into a binary format, for assistance with
# sensitivity calculations later on
all_samples_sprt <- binary_predictions(df = all_samples_sprt, 
                              prediction = quo(sprt_prediction)) %>%
  dplyr::rename(sprt_binary = binary_call)

#########################
# SPRT analysis - gDNA
#########################

het_gdna_sprt <- het_gdna %>%
  # Remove replicates outside the DNA input range of cfDNA
  filter(vf_assay_molecules < 30000) %>%
  # Assign an artificial fetal fraction of 4% for SPRT calculations
  mutate(fetal_fraction = 0.04,
         # Calculate the likelihood ratio
         likelihood_ratio = case_when(
           inheritance_chromosomal == "x_linked" ~ 
             calc_lr_x_linked(fetal_fraction, (major_allele_percent/100), 
                              vf_assay_molecules),
           inheritance_chromosomal == "autosomal" ~ 
             calc_lr_autosomal(fetal_fraction, (major_allele_percent/100), 
                               vf_assay_molecules)))

# Predict genotypes
het_gdna_sprt <- predict_sprt_genotypes(het_gdna_sprt, 8)

#########################
# MCMC analysis - cfDNA
#########################

# Prepare ddPCR data for MCMC

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

# 07/11/2021 - 15 minutes for 127 cfDNA samples
all_samples_mcmc <- run_mcmc(ddpcr_data_mcmc, 0.95)

all_samples_mcmc <- binary_predictions(
  df = all_samples_mcmc,
  prediction = quo(mcmc_prediction)) %>%
  dplyr::rename(mcmc_binary = binary_call)

#########################
# MCMC analysis - gDNA
#########################

# Set an artificial fetal fraction of 4%
mcmc_fetal_fraction <- 0.04

gdna_data_mcmc <- het_gdna %>%
  # Remove replicates outside the DNA input range of cfDNA
  filter(vf_assay_molecules < 30000) %>%
  mutate(
    var_ref_positives = variant_positives + reference_positives,
    var_ref_molecules = poisson_correct(vf_assay_droplets, var_ref_positives),
    
    # Create an artificial dataset for the fetal fraction.
    # Use the same number of total droplets as the variant fraction assay.
    ff_assay_droplets = vf_assay_droplets,
    
    #Calculate the number of fetal-specific and maternal molecules
    paternal_molecules = var_ref_molecules * mcmc_fetal_fraction,
    maternal_molecules = var_ref_molecules * (1-mcmc_fetal_fraction),
    
    # Calculate the expected number of positive partitions
    paternal_positives = reverse_poisson(paternal_molecules, 
                                         ff_assay_droplets),
    maternal_positives = reverse_poisson(maternal_molecules, 
                                         ff_assay_droplets)) %>%
  # Rename for pipeline
  dplyr::rename(n_K = vf_assay_droplets,
                K_M = variant_positives,
                K_N = reference_positives,
                Z_X = maternal_positives,
                Z_Y = paternal_positives,
                n_Z = ff_assay_droplets) %>%
  select(worksheet_well_sample, inheritance_chromosomal, inheritance_pattern,
         vf_assay, n_K, n_Z, K_N, K_M, n_Z, Z_X, Z_Y)

# 527 replicates of gDNA takes ~1 hour 4 mins to run
gdna_mcmc_analysed <- run_mcmc(gdna_data_mcmc, 0.95)

# Export the file as a csv so you don't have to keep rerunning the analysis.

write.csv(gdna_mcmc_analysed, paste0("analysis_outputs/gdna_mcmc_analysed ",
                 format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
          row.names = FALSE)

#########################
# Heterozygous gDNA variation
#########################

# vp is variant percent
het_gDNA_mean_vp <- mean(het_gdna_range$variant_percent)            

# sd is standard deviation

het_gDNA_sd_vp <- sd(het_gdna_range$variant_percent)

# cv is coefficient of variation
het_gDNA_cv_vp <- (het_gDNA_sd_vp/het_gDNA_mean_vp) * 100

#########################
# Z score analysis - cfDNA
#########################

zscore_imbalance_threshold <- 3
zscore_balance_threshold <- 2

all_samples_zscore <- cfdna_ddpcr_data_molecules %>%
  mutate(zscore = (variant_percent - het_gDNA_mean_vp) / het_gDNA_sd_vp,
         
         zscore_prediction = case_when(
           
           # Exclude cases with low cffDNA or low cfDNA
           fetal_percent < 4 ~"inconclusive",
           vf_assay_molecules < vf_assay_molecules_limit ~"inconclusive",
           
           # Predict fetal genotypes for X-linked cases
           inheritance_chromosomal == "x_linked" &
             zscore > zscore_imbalance_threshold ~"hemizygous variant",
           
           inheritance_chromosomal == "x_linked" &
             zscore < -zscore_imbalance_threshold ~"hemizygous reference",
           
           inheritance_chromosomal == "x_linked" &
             zscore > -zscore_imbalance_threshold &
             zscore < zscore_imbalance_threshold ~"inconclusive",
           
           # Predict fetal genotypes for autosomal cases
           inheritance_chromosomal == "autosomal" &
             zscore > zscore_imbalance_threshold ~"homozygous variant",
           
           inheritance_chromosomal == "autosomal" &
             zscore < -zscore_imbalance_threshold ~"homozygous reference",
           
           inheritance_chromosomal == "autosomal" &
             zscore < zscore_balance_threshold &
             zscore > -zscore_balance_threshold ~"heterozygous",
           
           inheritance_chromosomal == "autosomal" &
             zscore < zscore_imbalance_threshold &
             zscore > zscore_balance_threshold ~"inconclusive",
           
           inheritance_chromosomal == "autosomal" &
             zscore > -zscore_imbalance_threshold &
             zscore < -zscore_balance_threshold ~ "inconclusive"))

all_samples_zscore <- binary_predictions(
  df = all_samples_zscore,
  prediction = quo(zscore_prediction)) %>%
  dplyr::rename(zscore_binary = binary_call)

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
    select(r_number, p_G0, p_G1, p_G2, p_G3, mcmc_prediction, mcmc_binary),
  by = "r_number") %>%
  left_join(
    all_samples_zscore %>%
      select(r_number, zscore, zscore_prediction, zscore_binary),
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
    sprt_outcome = case_when(
      sprt_prediction == "inconclusive" ~"inconclusive",
      sprt_prediction == fetal_genotype ~"correct",
      TRUE ~"incorrect"),
    
    mcmc_outcome = case_when(
      mcmc_prediction == "inconclusive" ~"inconclusive",
      mcmc_prediction == fetal_genotype ~"correct",
      TRUE ~"incorrect"),
    
    zscore_outcome = case_when(
      zscore_prediction == "inconclusive" ~"inconclusive",
      zscore_prediction == fetal_genotype ~"correct",
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
         zscore_outcome = factor(zscore_outcome,
                                 levels = c("correct", "incorrect", 
                                            "inconclusive")),
         mcmc_outcome = factor(mcmc_outcome,
                               levels = c("correct", "incorrect", 
                                          "inconclusive")),
         sprt_outcome = factor(sprt_outcome,
                               levels = c("correct", "incorrect", 
                                          "inconclusive")))

#########################
# Summary of all results: Supplementary Data
#########################

all_samples_arranged <- all_samples_unblinded %>%
  mutate(cohort = ifelse(vf_assay == "HBB c.20A>T", "sickle cell disease",
                       "bespoke design"),
         cohort = factor(cohort, levels = c("sickle cell disease",
                                          "bespoke design"))) %>%
  arrange(cohort, inheritance_chromosomal, inheritance_pattern, study_id)
  
families <- data.frame(
  study_id = unique(all_samples_arranged$study_id))

families <- families %>%
  mutate(family_number = rownames(families))

supplementary_table <- families %>%
  full_join(all_samples_arranged,
            by = "study_id") %>%
  # Add on condition name
  left_join(gene_info %>%
              select(vf_assay, condition),
            by = "vf_assay") %>%
  mutate(sample_id = 
         paste0("cfDNA-", as.character(row.names(all_samples_arranged))),
         
         # Specify the method used to determine the fetal fraction
         ff_determination = case_when(
           inheritance_chromosomal == "x_linked" ~"ZFXY",
           inheritance_chromosomal == "autosomal" & ff_assay %in%
             ddpcr_snp_panel$dbSNP ~ "ddPCR SNP panel - workflow 2",
           TRUE ~"NGS SNP panel - workflow 1")) %>%
  
  # Reorder columns
  select(
    # Sample identifiers (r_number and study_id removed)
    sample_id, family_number, 
    # Sampling information
    hours_to_first_spin, days_to_storage,
    vacutainer, gestation_character,
    diagnostic_sampling, 
    # Extraction information
    plasma_volume_ml, extraction_replicates, 
    # Variant information
    cohort, inheritance_chromosomal, inheritance_pattern, condition,
    vf_assay, ff_determination, ff_assay, vf_assay_num_wells, 
    vf_assay_droplets, reference_positives,
    variant_positives, ff_assay_num_wells, ff_assay_droplets, 
    maternal_positives, paternal_positives,
    # Molecules
    variant_molecules, reference_molecules, vf_assay_molecules,
    variant_percent,
    maternal_molecules, paternal_molecules, fetal_percent,
    ff_assay_molecules, total_molecules, totalGE_ml_plasma,
    # Poisson confidence limits
    fetal_percent_min,
    fetal_percent_max,
    variant_percent_min,
    variant_percent_max,
    vf_assay_molecules_min,
    vf_assay_molecules_max,
    # Analysis
    likelihood_ratio, sprt_prediction, sprt_binary, 
    p_G0, p_G1, p_G2, p_G3, mcmc_prediction, mcmc_binary, 
    zscore, zscore_prediction, zscore_binary,
    # Outcome
    fetal_genotype, sprt_outcome, mcmc_outcome, zscore_outcome)

# Export the results
write.csv(supplementary_table, 
          file = (paste0("analysis_outputs/Supplementary_data_",
               format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")),
          row.names = FALSE)

###################
# Gene information table: Table 1
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
  
  # Rename without underscores
  dplyr::rename(Inheritance = inheritance_abbreviation,
                Gene = gene,
                Transcript = transcript,
                DNA = variant_dna,
                Protein = variant_protein,
                Condition = condition,
                Samples = samples)

write.csv(vf_assay_count_table, 
          file = (paste0("analysis_outputs/Table 1 ",
                         format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")),
          row.names = FALSE)

#########################
# HBB c.20A>T limit of detection study
#########################

lod_data <- read_csv("data/ddpcr_data/20-1557.csv", col_names = TRUE) %>%
  janitor::clean_names() %>%
  dplyr::rename(droplets = accepted_droplets)

# There are 4 replicates per spike-in. Merging them together gives
# outputs of 3000, 6000, 9000 and 12000 genome equivalents.

# Wells for 3000 GE dataset
wells_3000_GE <- c("A02", "A03", "E03", "A04", "E04", "A05", "E05", 
                   "A06", "E06", "A07", "E07", "A08", "E08")

wells_6000_GE <- c(wells_3000_GE, "B02", "B03", "F03", "B04", "F04", 
                   "B05", "F05", "B06", "F06", "B07", "F07", "B08", "F08")

wells_9000_GE <- c(wells_6000_GE, "C02", "C03", "G03", "C04", "G04", 
                   "C05", "G05", "C06", "G06", "C07", "G07", "C08", "G08")

wells_12000_GE <- c(wells_9000_GE, "D02", "D03", "H03", "D04", "H04", 
                    "D05", "H05", "D06", "H06", "D07", "H07", "D08", "H08")

# Merge values together with a function

merge_lod_wells <- function(input_wells, GE_level){
  
  output <- lod_data %>% 
    select(well, sample, target, positives, droplets) %>%
    filter(well %in% input_wells) %>%
    group_by(sample, target) %>% 
    summarise(positives = sum(positives),
              droplets = sum(droplets),
              .groups="drop") %>%
    pivot_wider(id_cols = c(sample),
                names_from = target,
                values_from = c(droplets, positives),
                names_glue = "{target}_{.value}") %>%
    select(-c(HbA_droplets)) %>%
    dplyr::rename(vf_assay_droplets = HbS_droplets,
                  variant_positives = HbS_positives,
                  reference_positives = HbA_positives) %>%
    mutate(GE_level = GE_level)
  
  return(output)
}

# Bind the datasets together, then apply var_ref_calculations

lod_data_merged <- var_ref_calculations(rbind(
  merge_lod_wells(wells_3000_GE, 3000), 
  merge_lod_wells(wells_6000_GE, 6000),
  merge_lod_wells(wells_9000_GE, 9000),
  merge_lod_wells(wells_12000_GE, 12000))) %>%
  # Allows easier colour labeling
  mutate(sample = factor(sample, levels = 
                 c("SS 12%", "SS 10%", "SS 8%", "SS 6%", "SS 4%", "SS 2%",
                   "0%",
                   "AA 2%", "AA 4%", "AA 6%", "AA 8%", "AA 10%", "AA 12%")))

# Supplemental Figure 4
lod_plot <- ggplot(lod_data_merged, 
                   aes(x = vf_assay_molecules, y = variant_percent)) +
  theme_bw() +
 theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.position = "right", 
    plot.title = element_text(size = 11),
    legend.title = element_blank(),
    legend.text = element_text(size = 9)) +
  geom_errorbar(aes(ymin = variant_percent_min, ymax = variant_percent_max),
                alpha = 0.2) +
  geom_errorbarh(aes(xmin = vf_assay_molecules_min, 
                     xmax = vf_assay_molecules_max),
                 alpha = 0.2) +
  # In order of shade
  # "#FFFFFF" = white
  # "#CCCCCC" = grey 1
  # "#9999CC" = grey 2
  # "#999999" = grey 3
  # "#666666" = grey 4
  # "#333333" = grey 5
  # "#000000" = black; 
  scale_fill_manual(values = c(
    # SS 12% to 2%
    "#000000", "#333333", "#666666", "#999999", "#9999CC", "#CCCCCC",
    # 0%
    "#FFFFFF",
    # 2% to 12%
    "#CCCCCC", "#9999CC", "#999999", "#666666", "#333333","#000000")) +
  scale_shape_manual(values = c(24, 24, 24, 24, 24, 24, 21, 
                                25,25, 25, 25, 25, 25)) +
  geom_point(size = 2, aes(fill = sample, shape = sample)) +
  ylim(39, 61) +
  xlim(0, 15000) +
  geom_hline(yintercept = 50, linetype = "dashed") +
  labs(x = "Genome Equivalents (GE)",
       y = "Variant fraction (%)",
       title = "")

ggsave(plot = lod_plot, 
       filename = paste0("lod_plot_",
                         format(Sys.time(), "%Y%m%d_%H%M%S"), ".tiff"),
       path = "plots/", 
       device='tiff', 
       dpi=300,
       units = "in",
       width = 8,
       height = 6)

###################
# Sensitivity and specificity table - Table 2
###################

sprt_scd <- sensitivity_metrics(df = supplementary_table, 
                    prediction_binary = quo(sprt_binary), 
                    outcome = quo(sprt_outcome), 
                    cohort_input = "sickle cell disease",
                    cohort_name = "sickle cell disease") %>%
  dplyr::rename(sprt = analysis_method)

sprt_bespoke <- sensitivity_metrics(df = supplementary_table, 
                                prediction_binary = quo(sprt_binary), 
                                outcome = quo(sprt_outcome), 
                                cohort_input = "bespoke design",
                                cohort_name = "bespoke design") %>%
  dplyr::rename(sprt = analysis_method)

sprt_all <- sensitivity_metrics(df = supplementary_table, 
                                    prediction_binary = quo(sprt_binary), 
                                    outcome = quo(sprt_outcome), 
                                    cohort_input = c("bespoke design", 
                                                     "sickle cell disease"),
                                    cohort_name = "all") %>%
  dplyr::rename(sprt = analysis_method)

mcmc_scd <- sensitivity_metrics(df = supplementary_table, 
                                prediction_binary = quo(mcmc_binary), 
                                outcome = quo(mcmc_outcome), 
                                cohort_input = "sickle cell disease",
                                cohort_name = "sickle cell disease") %>%
  dplyr::rename(mcmc = analysis_method)

mcmc_bespoke <- sensitivity_metrics(df = supplementary_table, 
                                    prediction_binary = quo(mcmc_binary), 
                                    outcome = quo(mcmc_outcome), 
                                    cohort_input = "bespoke design",
                                    cohort_name = "bespoke design") %>%
  dplyr::rename(mcmc = analysis_method)

mcmc_all <- sensitivity_metrics(df = supplementary_table, 
                                prediction_binary = quo(mcmc_binary), 
                                outcome = quo(mcmc_outcome), 
                                cohort_input = c("bespoke design", 
                                                 "sickle cell disease"),
                                cohort_name = "all") %>%
  dplyr::rename(mcmc = analysis_method)

zscore_scd <- sensitivity_metrics(df = supplementary_table, 
                                prediction_binary = quo(zscore_binary), 
                                outcome = quo(zscore_outcome), 
                                cohort_input = "sickle cell disease",
                                cohort_name = "sickle cell disease") %>%
  dplyr::rename(zscore = analysis_method)

zscore_bespoke <- sensitivity_metrics(df = supplementary_table, 
                                    prediction_binary = quo(zscore_binary), 
                                    outcome = quo(zscore_outcome), 
                                    cohort_input = "bespoke design",
                                    cohort_name = "bespoke design") %>%
  dplyr::rename(zscore = analysis_method)

zscore_all <- sensitivity_metrics(df = supplementary_table, 
                                prediction_binary = quo(zscore_binary), 
                                outcome = quo(zscore_outcome), 
                                cohort_input = c("bespoke design", 
                                                 "sickle cell disease"),
                                cohort_name = "all") %>%
  dplyr::rename(zscore = analysis_method)


# Row  bind analysis methods together

sprt_metrics <- rbind(sprt_scd, sprt_bespoke, sprt_all)
mcmc_metrics <- rbind(mcmc_scd, mcmc_bespoke, mcmc_all)
zscore_metrics <- rbind(zscore_scd, zscore_bespoke, zscore_all)

analysis_metrics <- sprt_metrics %>%
  left_join(mcmc_metrics,
            by = c("group", "category")) %>%
  left_join(zscore_metrics,
            by = c("group", "category")) %>%
  # Tidy up names for paper
  dplyr::rename(
    `SPRT` = sprt,
    `MCMC` = mcmc,
    `Z score` = zscore) %>%
  mutate(`Result` = case_when(
            category == "true_positive" ~ "True positive",
            category == "true_negative" ~ "True negative",
            category == "false_positive" ~ "False positive",
            category == "false_negative" ~ "False negative",
            category == "sensitivity" ~ "Sensitivity (%)",
            category == "specificity" ~ "Specificity (%)",
            category == "inconclusive" ~ "Inconclusive"),
         `Cohort` = case_when(
           group == "sickle cell disease" ~ "Sickle cell disease",
           group == "bespoke design" ~ "Bespoke design",
           group == "all" ~"All samples")) %>%
  select(`Cohort`, `Result`, `SPRT`, `MCMC`, `Z score`)

write.csv(x = analysis_metrics, 
          file = paste0("analysis_outputs/Table_2_",
                  format(Sys.time(), "%Y%m%d_%H%M%S"),
                  ".csv"),
          row.names = FALSE)

#########################
# Plot cfDNA results - Figure 2
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
# Plot 2A: SPRT analysis
###################

plot_2a <- ggplot(all_samples_unblinded, aes(x = vf_assay_molecules, 
                                            y = variant_percent)) +
  theme_bw() +
  multiplot_theme + 
  multiplot_y +
  multiplot_x +
  cfdna_fill +
  cfdna_alpha +
  cfdna_shape +
  geom_point(size = 2, aes(fill = sprt_outcome,
                           alpha = sprt_outcome,
                           shape = fetal_genotype),
             colour = "black") +
  
  labs(y = "Variant fraction (%)", x = "", 
       title = "ddPCR for 124 cfDNA samples with SPRT classification")

###################
# Plot 2B: MCMC analysis
###################

plot_2b <- ggplot(all_samples_unblinded, aes(x = vf_assay_molecules, 
                                            y = variant_percent)) +
  theme_bw() +
  multiplot_theme + 
  multiplot_y +
  multiplot_x +
  cfdna_fill +
  cfdna_alpha +
  cfdna_shape +
  geom_point(size = 2, aes(fill = mcmc_outcome,
                           alpha = mcmc_outcome,
                           shape = fetal_genotype),
             colour = "black") +
  
  labs(y = "", x = "", 
       title = "ddPCR for 124 cfDNA samples with MCMC classification")

###################
# Plot 2C: heterozygous gDNA controls
###################

plot_2c <- ggplot(het_gdna %>%
         mutate(sample_type = "het gDNA"), 
                    aes(x = vf_assay_molecules, 
                        y = variant_percent)) +
  z3_line +
  zminus3_line +
  z2_line +
  zminus2_line +
  vertical_line +
  multiplot_y +
  scale_x_continuous(limits = c(0,32000),
                     breaks = c(0, 2000, 10000, 20000 ,30000)) +
  scale_shape_manual(values = c(21)) +
  geom_point(size = 2, aes(shape = sample_type), fill= "white",
             colour = "black") +
  theme_bw() +
  multiplot_theme +
  labs(y = "Variant fraction (%)", x = "Genome equivalents (GE)",
       title = "ddPCR for 82 heterozygous gDNA controls")

###################
# Plot 2D: z score analysis
###################

plot_2d <- ggplot(all_samples_unblinded, aes(x = vf_assay_molecules, 
             y = variant_percent)) +
  theme_bw() +
  multiplot_theme + 
  z3_line +
  zminus3_line +
  z2_line +
  zminus2_line +
  vertical_line +
  multiplot_y +
  scale_x_continuous(limits = c(0,32000),
                     breaks = c(0, 2000, 10000, 20000 ,30000)) +
  cfdna_fill +
  cfdna_alpha +
  cfdna_shape +
  geom_point(size = 2, aes(fill = zscore_outcome,
                           alpha = zscore_outcome,
                           shape = fetal_genotype),
             colour = "black") +
  
  labs(y = "", x = "Genome equivalents (GE)", 
       title = "ddPCR for 124 cfDNA samples with z score classification")

###################
# All plots together
###################

ddpcr_cohort <- ggpubr::ggarrange(plot_2a, plot_2b,
                                  plot_2c, plot_2d,
                                  ncol = 2, nrow = 2, align = "v",
                                  labels = c("A", "B", "C", "D"))

ggsave(plot = ddpcr_cohort, 
       filename = paste0("figure_2_",
                     format(Sys.time(), "%Y%m%d_%H%M%S"),
                     ".tiff"),
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
  count(zscore_outcome) %>%
  dplyr::rename(outcome = zscore_outcome) %>%
  mutate(analysis = "z score"),
  all_samples_unblinded %>%
    count(sprt_outcome) %>%
    dplyr::rename(outcome = sprt_outcome) %>%
    mutate(analysis = "sprt"),
  all_samples_unblinded %>%
    count(mcmc_outcome) %>%
    dplyr::rename(outcome = mcmc_outcome) %>%
    mutate(analysis = "mcmc")) %>%
  mutate(outcome = factor(outcome, 
                          levels = c("correct",
                                     "inconclusive",
                                     "incorrect")),
         analysis = factor(analysis,
                           levels = c("sprt", "mcmc", "z score")))

# Bar plot
analysis_summary_plot <- ggplot(analysis_summary, aes(x = outcome, y = n))+
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
# ROC curve analysis - Supplemental Figure 6
###################

# Plotting ROC curves can be achieved manually usually tidyverse functions,
# or using the pROC package.

# We need a binary classifier for the fetal genotype.

roc_binary_calls <- all_samples_unblinded %>%
  # Convert invasive results to binary outcomes based on whether the 
  # fetus is affected or not
  mutate(fetus_affected = case_when(
      # Autosomal dominant inheritance
      inheritance_chromosomal == "autosomal" &
        inheritance_pattern == "dominant" &
        fetal_genotype == "heterozygous fetus" ~"TRUE",
      inheritance_chromosomal == "autosomal" &
        inheritance_pattern == "dominant" &
        fetal_genotype == "homozygous reference fetus" ~"FALSE",
      
      # Autosomal recessive inheritance
      # As there are 3 possible genotypes, both homozygous variant and 
      # homozygous reference fetuses are coded as "affected", as they both
      # involve an imbalance in cfDNA.
      inheritance_chromosomal == "autosomal" &
        inheritance_pattern == "recessive" &
        fetal_genotype %in% c("homozygous variant fetus", 
                            "homozygous reference fetus")  ~"TRUE",
      inheritance_chromosomal == "autosomal" &
        inheritance_pattern == "recessive" &
        fetal_genotype == "heterozygous fetus" ~"FALSE",
      
      # X linked inheritance
      inheritance_chromosomal == "x_linked" &
        fetal_genotype == "hemizygous variant fetus"  ~"TRUE",
      inheritance_chromosomal == "x_linked" &
        fetal_genotype == "hemizygous reference fetus"  ~"FALSE"),
      
      # Convert "fetus_affected" column to Boolean vector
      fetus_affected = as.logical(fetus_affected),
      # SPRT likelihood score and z scores are already single predictors
      # We need to convert the MCMC calls to a single value
      # Annoyingly, pG2 has different meanings for recessive and dominant
      # cases.
      mcmc_unbalanced_call = case_when(
        inheritance_chromosomal == "autosomal" &
          inheritance_pattern == "recessive" ~pmax(p_G1, p_G3),
        inheritance_chromosomal == "autosomal" &
          inheritance_pattern == "dominant" ~p_G1,
        inheritance_chromosomal == "x_linked" ~p_G1),
      # Remove minus signs from z score
      zscore_unbalanced_call = abs(zscore))

# 12585 has an "infinite" likelihood ratio, which messes up the ROC
# calculations, so need to replace with a finite number.

roc_binary_calls[roc_binary_calls$r_number == "12585", "likelihood_ratio"] <- 
  # Make it one larger than the largest finite likelihood ratio in the 
  # dataset
  (1.064861e+60 + 1)

# Create ROC objects for each analysis method
sprt_roc_object <- roc(
  # Response
  roc_binary_calls$fetus_affected, 
  # Predictor
  roc_binary_calls$likelihood_ratio)

mcmc_roc_object <- roc(
  # Response
  roc_binary_calls$fetus_affected, 
  # Predictor
  roc_binary_calls$mcmc_unbalanced_call)

zscore_roc_object <- roc(
  # Response
  roc_binary_calls$fetus_affected, 
  # Predictor
  roc_binary_calls$zscore_unbalanced_call)

# SPRT plot
sprt_roc <- ggroc(sprt_roc_object, size = 2) +
  labs(x = "", 
       y = "Sensitivity", 
       title = "SPRT analysis") +
  theme_bw() +
  multiplot_theme +
  geom_abline(linetype = "dashed",
              intercept = 1)+
  annotate(geom = "text", x = 0.25, y = 0.25, 
           label = paste0("AUC = ", 
                          round(auc(roc_binary_calls$fetus_affected, 
                                    roc_binary_calls$likelihood_ratio),3)))
# MCMC plot
mcmc_roc <- ggroc(mcmc_roc_object, size = 2) +
  labs(x = "Specificity", 
       y = "", 
       title = "MCMC analysis") +
  theme_bw() +
  multiplot_theme +
  geom_abline(linetype = "dashed",
              intercept = 1)+
  annotate(geom = "text", x = 0.25, y = 0.25, 
           label = paste0("AUC = ", 
                          round(auc(roc_binary_calls$fetus_affected, 
                                    roc_binary_calls$mcmc_unbalanced_call),3)))

# Z score plot
zscore_roc <- ggroc(zscore_roc_object, size = 2) +
  labs(x = "", 
       y = "", 
       title = "Z score analysis") +
  theme_bw() +
  multiplot_theme +
  geom_abline(linetype = "dashed",
              intercept = 1)+
  annotate(geom = "text", x = 0.25, y = 0.25, 
           label = paste0("AUC = ", 
                          round(auc(roc_binary_calls$fetus_affected, 
                                    roc_binary_calls$zscore_unbalanced_call),3)))

# All ROC plots together
roc_plot <- ggpubr::ggarrange(sprt_roc, mcmc_roc, 
                              zscore_roc, 
                              ncol = 3, nrow = 1)

ggsave(plot = roc_plot, 
       filename = paste0("roc_plot_",
                         format(Sys.time(), "%Y%m%d_%H%M%S"),
                         ".tiff"),
       path = "plots/", device='tiff',
       units = "in",
       width = 10,
       height = 7)

###################
# SPRT and MCMC gDNA plots - Supplemental Figure 5
###################

# Set fetal fraction and likelihood ratio thresholds
ff_for_graph <- 0.04
lr_for_graph <- 8

# SPRT function lines
het_upper_line <- geom_function(fun = "calc_het_upper_boundary",
              aes(x = vf_assay_molecules, y =
                    calc_het_upper_boundary(vf_assay_molecules, ff_for_graph, 
                                           lr_for_graph)),
              colour = "black",
              args = c(ff_for_graph, lr_for_graph),
              alpha = 0.5)
  
het_lower_line <- geom_function(fun = "calc_het_lower_boundary",
                aes(x = vf_assay_molecules, y =
                      calc_het_lower_boundary(vf_assay_molecules, ff_for_graph, 
                                             lr_for_graph)),
                colour = "black",
                args = c(ff_for_graph, lr_for_graph),
                alpha = 0.5)
  
hom_var_line <- geom_function(fun = "calc_hom_var_boundary",
                aes(x = vf_assay_molecules, y =
                      calc_hom_var_boundary(vf_assay_molecules, ff_for_graph, 
                                       lr_for_graph)),
                colour = "black",
                args = c(ff_for_graph, lr_for_graph),
                alpha = 0.5)
  
hom_ref_line <- geom_function(fun = "calc_hom_ref_boundary",
                aes(x = vf_assay_molecules, y =
                      calc_hom_ref_boundary(vf_assay_molecules, ff_for_graph, 
                                       lr_for_graph)),
                colour = "black",
                args = c(ff_for_graph, lr_for_graph),
                alpha = 0.5)

hemi_var_line <- geom_function(fun = "calc_hemi_var_boundary",
              aes(x = vf_assay_molecules, y =
                    calc_hemi_var_boundary(vf_assay_molecules, ff_for_graph, 
                                           lr_for_graph)),
              colour = "black",
              args = c(ff_for_graph, lr_for_graph),
              alpha = 0.5)
  
hemi_ref_line <- geom_function(fun = "calc_hemi_ref_boundary",
                aes(x = vf_assay_molecules, y =
                      calc_hemi_ref_boundary(vf_assay_molecules, ff_for_graph, 
                                             lr_for_graph)),
                colour = "black",
                args = c(ff_for_graph, lr_for_graph),
                alpha = 0.5)
  
het_region_label <- annotate(geom = "text", x = 28000, y = 50, 
                             label = "heterozygous")

hom_ref_region_label <- annotate(geom = "text", x = 28000, y = 47, 
                          label = "homozygous reference")

hom_var_region_label <- annotate(geom = "text", x = 28000, y = 53, 
                          label = "homozygous variant")

hemi_ref_region_label <- annotate(geom = "text", x = 28000, y = 47, 
                           label = "hemizygous reference")

hemi_var_region_label <- annotate(geom = "text", x = 28000, y = 53, 
                           label = "hemizygous variant")

gdna_plots_x <- xlim(0, 30000)
  
gdna_plots_y <- ylim(43, 57)

###########
# Plot S5A: SPRT on gDNA for HBB c.20A>T
###########

gdna_plot_title <- expression(paste(italic("HBB"), 
                                    "c.20A>T gDNA: SPRT results"))
plot_s5a <- ggplot(het_gdna_sprt %>%
         filter(vf_assay == "HBB c.20A>T"), 
       aes(x = vf_assay_molecules, y = variant_percent)) +
  # Four colours required
  scale_fill_manual(values=c("#000000", "#FFFFFF", "#000000", "#999999"), 
                    guide = "none") +
  geom_point(size = 2,
             aes(fill = sprt_prediction), 
             pch=21,
             alpha = 0.8) +
  theme_bw()+
  multiplot_theme +
  het_upper_line +
  het_lower_line +
  hom_var_line +
  hom_ref_line +
  gdna_plots_x +
  gdna_plots_y +
  labs(x = "",
       y = "Variant fraction (%)",
       title = gdna_plot_title)
  
###########
# Plot S5B: SPRT on gDNA for autosomal assays
###########

plot_s5b <- ggplot(het_gdna_sprt %>%
         filter(inheritance_chromosomal == "autosomal" &
                  vf_assay != "HBB c.20A>T"), 
       aes(x = vf_assay_molecules, y = variant_percent))+
  # Three colours required
  scale_fill_manual(values=c("#000000", "#FFFFFF", "#999999"), 
                    guide = "none") +      
  geom_point(size = 2, 
             colour = "black", 
             aes(fill = sprt_prediction), 
                   pch=21,
                   alpha = 0.8) +
  theme_bw()+
  multiplot_theme +
  het_upper_line +
  het_lower_line +
  hom_var_line +
  hom_ref_line +
  gdna_plots_x +
  gdna_plots_y +
  labs(x = "",
       y = "Variant fraction (%)",
       title = "Autosomal variant gDNA: SPRT results")

###########
# Plot S5C: SPRT on gDNA for X-linked assays
###########

plot_s5c <- ggplot(het_gdna_sprt %>%
         filter(inheritance_chromosomal == "x_linked"), 
       aes(x = vf_assay_molecules, y = variant_percent))+
  # Three colours required
  scale_fill_manual(values=c("#000000", "#000000", "#999999"), 
                    guide = "none") +
  geom_point(size = 2, 
             aes(fill = sprt_prediction),
             pch=21,
             alpha = 0.8) +
  theme_bw()+
  multiplot_theme +
  hemi_var_line +
  hemi_ref_line +
  gdna_plots_x +
  gdna_plots_y +
  labs(x = "Genome equivalents (GE)",
       y = "Variant fraction (%)",
       title = "X-linked variant gDNA: SPRT results")
 
###########
# Plot S5D: MCMC on gDNA for HBB c.20A>T
###########

mcmc_gdna_for_plot <- gdna_mcmc_analysed %>%
  select(-c(inheritance_chromosomal, inheritance_pattern,
            vf_assay, fetal_fraction)) %>%
  left_join(het_gdna_sprt,
    by = "worksheet_well_sample") %>%
  # Factorise for plotting
  mutate(mcmc_prediction = factor(mcmc_prediction, levels = 
                             c("hemizygous variant",
                               "homozygous variant",
                               "heterozygous",
                               "homozygous reference",
                               "hemizygous reference",
                               "inconclusive")))

gdna_plot_title <- expression(paste(italic("HBB"), 
                                   "c.20A>T gDNA: MCMC results"))

plot_s5d <- ggplot(mcmc_gdna_for_plot %>%
                   filter(vf_assay == "HBB c.20A>T"), 
                 aes(x = vf_assay_molecules, y = variant_percent)) +
  # Two colours required
  scale_fill_manual(values=c("#FFFFFF", "#999999"), 
                    guide = "none") +
  geom_point(size = 2,
             aes(fill = mcmc_prediction), 
             pch=21,
             alpha = 0.8) +
  theme_bw()+
  multiplot_theme +
  gdna_plots_x +
  gdna_plots_y +
  labs(x = "",
       y = "",
       title = gdna_plot_title)

###########
# Plot S5E: MCMC on gDNA for autosomal assays
###########

plot_s5e <- ggplot(mcmc_gdna_for_plot %>%
         filter(inheritance_chromosomal == "autosomal" &
                  vf_assay != "HBB c.20A>T"), 
       aes(x = vf_assay_molecules, y = variant_percent)) +
  # Two colours required
  scale_fill_manual(values=c("#FFFFFF", "#999999"), 
                    guide = "none") +
  geom_point(size = 2,
             aes(fill = mcmc_prediction), 
             pch=21,
             alpha = 0.8) +
  theme_bw()+
  multiplot_theme +
  gdna_plots_x +
  gdna_plots_y +
  labs(x = "",
       y = "",
       title = "Autosomal variant gDNA: MCMC results")

###########
# Plot S5F: MCMC on gDNA for X-linked assays
###########

plot_s5f <- ggplot(mcmc_gdna_for_plot %>%
         filter(inheritance_chromosomal == "x_linked"), 
       aes(x = vf_assay_molecules, y = variant_percent))+
  # Three colours required
  scale_fill_manual(values=c("#000000", "#000000", "#999999"), 
                    guide = "none") +
  geom_point(size = 2, 
             aes(fill = mcmc_prediction),
             pch=21,
             alpha = 0.8) +
  theme_bw()+
  multiplot_theme +
  gdna_plots_x +
  gdna_plots_y +
  labs(x = "Genome equivalents (GE)",
       y = "",
       title = "X-linked variant gDNA: MCMC results")

###################
# Arrange plots together
###################

gdna_results_plot <- ggpubr::ggarrange(plot_s5a, plot_s5d,
                                       plot_s5b, plot_s5e, 
                                       plot_s5c, plot_s5f,
                                  ncol = 2, nrow = 3, align = "v",
                                  labels = c("A", "D", "B",
                                             "E", "C", "F"))
ggsave(plot = gdna_results_plot,
       filename = paste0("gdna_results_plot_",
              format(Sys.time(), "%Y%m%d_%H%M%S"),
              ".tiff"),
       path = "plots/", device='tiff', dpi=600,
       units = "in",
       width = 8,
       height = 10)

###################