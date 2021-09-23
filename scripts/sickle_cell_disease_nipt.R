################################################################################
## Accuracy of Non-Invasive Fetal Genotyping for Sickle Cell Disease 
## with Droplet Digital PCR
## September 2021
## Joseph.Shaw@gosh.nhs.uk
## This is an analysis script for the prediction of fetal 
## genotypes from cfDNA testing using ddPCR for sickle cell 
## disease, based on dosage experiments with heterozygous gDNA controls.
################################################################################

#########################
# Load libraries and resources
#########################

## Load necessary packages
library(tidyverse)
library(epiR)
library(ggpubr)

## Set working directory
setwd("W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/ddPCR_R_Analysis/ddpcr_nipd")

# Source functions (this includes loading ddPCR data)
source("functions/ddPCR_nipd_functions.R")
source("W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/RAPID_project_biobank/scripts/RAPID_biobank.R")

#########################
# Select sickle cell disease cfDNA ddPCR data
#########################

# Sample to exclude:
# 13262 - this sample had contamination
samples_to_exclude <- c(13262)

# The "secondary cohort" refers to the most recent phase of the 
# ddPCR sickle cell disease project, when all samples were extracted 
# using a 6ml protocol.

secondary_cohort <- c("14182", "19868", "20238", "20611", 
                      "20874", "30063", "30068", "30113", "30142", 
                      "30206", "30228", "30230", "30078", "30065", 
                      "13402", "20939", "30215", "30203",
                      "20911", "30236", "30112", "30251", "30257",
                      "30216", "30232")

cfDNA_scd_data <- ff_calculations(
  var_ref_calculations(cfdna_ddpcr_data)) %>%
  dplyr::rename(r_number = sample) %>%
  filter(vf_assay == "HBB c.20A>T" &
           !r_number %in% samples_to_exclude) %>%
  mutate(cohort = ifelse(r_number %in% secondary_cohort,
                         "secondary", "primary"),
         sample_type = "cfDNA")

#########################
# Prepare HbAS gDNA data
#########################

# This step selects every HbAS gDNA control which we have run, including
# single well data, and data merged by sample for each worksheet.

gDNA_scd_data <- parent_gDNA_var_ref %>%
   filter(sample %in% controls$sample & 
            # Remove any controls that aren't HbAS gDNA
            !sample %in% c("NTC", "30139", "30130", 
                           "19RG-220G0191", "19RG-220G0193",
                           "21RG-120G0072") &
            # Remove empty well 
            reference_positives != 0 &
            vf_assay == "HBB c.20A>T") %>%
  
  # New column to distinguish control gDNA
  mutate(sample_type = "gDNA",
        sprt_prediction = case_when(
          
         variant_percent > calc_SS_boundary(vf_assay_molecules,
                                            0.04, 8) &
              major_allele == "variant allele"
            ~"HbSS",
         
         variant_percent < calc_AA_boundary(vf_assay_molecules,
                                            0.04, 8) &
           major_allele == "reference allele"
         ~"HbAA",
         
         variant_percent > calc_AS_lower_boundary(vf_assay_molecules,
                                                  0.04, 8) &
          variant_percent < calc_AS_upper_boundary(vf_assay_molecules,
                                                   0.04, 8)
         ~"HbAS",
          TRUE ~"inconclusive")) %>%
  dplyr::rename(r_number = sample)

#########################
# Supplemental Figure 4: SPRT with HbAS gDNA controls
#########################

# This plot shows that the common form of the SPRT equation is 
# inappropriate for this dataset

ff_for_graph <- 0.04
lr_for_graph <- 8

ggplot(gDNA_scd_data, aes(x = vf_assay_molecules, y = variant_percent))+
  
  # scale_fill_manual(values=c("#FFFFFF", "#FFFFFF", "#FFFFFF", "#999999")) +
  # scale_alpha_manual(values = c(1, 1, 1, 0.4)) +
  geom_point(size = 2, colour = "black", fill = "white", pch=21,
             alpha = 0.8) +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none")+
  labs(x = "Genome equivalents (GE)",
       y = "Variant fraction (%)"
       #, title = "Heterozygous gDNA samples"
       ) +
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
  scale_y_continuous(breaks = c(45, 50, 55), 
                     limits = c(45,55)) + 
  annotate(geom = "text", x = 45000, y = 50, label = "HbAS") +
  annotate(geom = "text", x = 45000, y = 47, label = "HbAA") +
  annotate(geom = "text", x = 45000, y = 53, label = "HbSS") +
  
  # Add arrows to show the incorrect predictions
  geom_segment(aes(x = 5500, y = 51.9, 
                   xend = 4500, yend = 51.9),
               arrow = arrow(length = unit(0.2, "cm"))) +
  
  geom_segment(aes(x = 4000, y = 53, 
                   xend = 3000, yend = 53),
               arrow = arrow(length = unit(0.2, "cm"))) +
  
  geom_segment(aes(x = 4000, y = 46.5, 
                   xend = 3000, yend = 46.5),
               arrow = arrow(length = unit(0.2, "cm"))) +
  
  geom_segment(aes(x = 7000, y = 47.7, 
                   xend = 6000, yend = 47.7),
               arrow = arrow(length = unit(0.2, "cm")))

#########################
# Variation in HbAS gDNA controls
#########################

vf_assay_molecules_limit <- 4000

# This section shows the variation in HbAS gDNA controls over 
# 4000 GE, in preparation for z score analysis in the next section

gDNA_scd_data_4000 <- gDNA_scd_data %>%
  filter(vf_assay_molecules > vf_assay_molecules_limit)

gDNA_scd_data_sub4000 <- gDNA_scd_data %>%
  filter(vf_assay_molecules < vf_assay_molecules_limit)

# Get max and min values for the paper writeup
max(gDNA_scd_data_4000$variant_percent)
min(gDNA_scd_data_4000$variant_percent)
max(gDNA_scd_data_sub4000$variant_percent)
min(gDNA_scd_data_sub4000$variant_percent)

# vp is "variant percent"
gDNA_mean_vp <- mean(gDNA_scd_data_4000$variant_percent)            

# Find the standard deviation, by using the "sd" function in stats package. 
# The standard deviation is the square root of the variance.
# The variance is the average of the squared differences.

gDNA_stand_dev_vp <- sd(gDNA_scd_data_4000$variant_percent)

# Coefficient of variation
gDNA_cv_vp <- (gDNA_stand_dev_vp/gDNA_mean_vp) * 100

# Z score thresholds
gDNA_mean_vp+(3*gDNA_stand_dev_vp)
gDNA_mean_vp-(3*gDNA_stand_dev_vp)
gDNA_mean_vp+(2*gDNA_stand_dev_vp)
gDNA_mean_vp-(2*gDNA_stand_dev_vp)

#########################
# Limit of detection study
#########################

lod_data <- read_csv("data/20-1557.csv", col_names = TRUE) %>%
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
  # Allows easier colour labelling
  mutate(sample = factor(sample, levels = 
                c("SS 12%", "SS 10%", "SS 8%", "SS 6%", "SS 4%", "SS 2%",
                  "0%",
                  "AA 2%", "AA 4%", "AA 6%", "AA 8%", "AA 10%", "AA 12%")),
         z_score = (variant_percent - gDNA_mean_vp) / gDNA_stand_dev_vp)

#########################
# Predict fetal genotypes
#########################

# The aim is to classify samples using z score analysis, in the same
# way as NIPT for trisomy 21.

lr_threshold <- 8
z_score_imbalance_threshold <- 3
z_score_balance_threshold <- 2

cfDNA_scd_predictions <- cfDNA_scd_data %>%
  mutate(
    z_score = (variant_percent - gDNA_mean_vp) / gDNA_stand_dev_vp,
    z_score_major = (major_allele_percent - gDNA_mean_vp) / gDNA_stand_dev_vp,
    
    z_score_genotype_prediction = case_when(
      # Samples with fetal fractions below 4% are inconclusive
      
      z_score > z_score_imbalance_threshold &
      fetal_percent > 4 &
      vf_assay_molecules > vf_assay_molecules_limit
      ~"HbSS",
      
      z_score < -z_score_imbalance_threshold &
      fetal_percent > 4 &
      vf_assay_molecules > vf_assay_molecules_limit
      ~"HbAA",
      
      z_score < z_score_balance_threshold &
      z_score > -z_score_balance_threshold &
      fetal_percent > 4 &
      vf_assay_molecules > vf_assay_molecules_limit
      ~"HbAS",
      
      fetal_percent < 4 ~"insufficient data",
      vf_assay_molecules < vf_assay_molecules_limit ~"insufficient data",
      
      TRUE ~"inconclusive"),
    
    z_score_clinical_prediction = case_when(
      
      z_score > z_score_imbalance_threshold &
        fetal_percent > 4 &
        vf_assay_molecules > vf_assay_molecules_limit
      ~"affected",
      z_score < z_score_balance_threshold &
        fetal_percent > 4 &
        vf_assay_molecules > vf_assay_molecules_limit
      ~ "unaffected",
      TRUE ~"inconclusive"),
    
    # Factorise for plot
    z_score_genotype_prediction = factor(z_score_genotype_prediction, levels = 
                                   c("HbSS", "HbAS", "HbAA",
                                     "inconclusive", "insufficient data")),
    
    # Calculate likelihood ratio
    likelihood_ratio = calc_lr_autosomal(fetal_fraction,
                                         major_allele_percent/100,
                                         vf_assay_molecules),
    
    # Calculate the upper and lower SPRT boundaries to check that my 
    # rearrangement of the SPRT likelihood ratio calculation doesn't contain
    # any errors
    upper_sprt_boundary = calc_SS_boundary(vf_assay_molecules,
                     fetal_fraction, 8),
    
    lower_sprt_boundary = calc_AS_upper_boundary(vf_assay_molecules,
                           fetal_fraction, 8),
    
    sprt_genotype_prediction = case_when(
      likelihood_ratio > lr_threshold &
        major_allele == "variant allele"
      ~"HbSS",
      likelihood_ratio > lr_threshold &
        major_allele == "reference allele"
      ~"HbAA",
      likelihood_ratio < 1/lr_threshold
      ~"HbAS",
      TRUE ~"inconclusive"))

#########################
# Extraction volumes and invasive sampling
#########################

plasma_extractions <- read.csv("resources/extraction_volumes.csv") %>%
  group_by(r_number) %>%
  summarise(plasma_volume_ml = (sum(tubes_removed))*2) 

plasma_replicates <- read.csv("resources/extraction_volumes.csv") %>%
  group_by(r_number) %>%
  summarise(extraction_replicates = n())

invasive_sampling <- read.csv("resources/confirmation_testing.csv")

#########################
# Compare predictions against Biobank
#########################

cfDNA_scd_outcomes <- left_join(
  cfDNA_scd_predictions,
  RAPID_biobank %>%
    mutate(r_number = as.character(r_number)) %>%
    select(r_number, study_id, site, 
           gestation_total_weeks,
           gestation_character, vacutainer,
           mutation_genetic_info_fetus, partner_sample_available,
           report_acquired, hours_to_first_spin,
           days_to_storage, sampling_date_time),
  by = "r_number") %>%
  mutate(outcome_zscore = case_when(
    z_score_genotype_prediction == "inconclusive" ~"inconclusive",
    z_score_genotype_prediction == "insufficient data" ~"insufficient data",
    is.na(mutation_genetic_info_fetus) ~"awaiting result",
    z_score_genotype_prediction == mutation_genetic_info_fetus
    ~"correct",
    TRUE ~"incorrect"),
    outcome_sprt = case_when(
      sprt_genotype_prediction == "inconclusive" ~"inconclusive",
      is.na(mutation_genetic_info_fetus) ~"awaiting result",
      sprt_genotype_prediction == mutation_genetic_info_fetus
      ~"correct",
      TRUE ~"incorrect")) %>%
  left_join(plasma_extractions %>%
              mutate(r_number = as.character(r_number)), by = "r_number") %>%
  left_join(plasma_replicates %>%
              mutate(r_number = as.character(r_number)), by = "r_number") %>%
  left_join(invasive_sampling %>%
              mutate(r_number = as.character(r_number)), by = "r_number") %>%

  # Calculate the number of cfDNA molecules per ml plasma
  
  mutate(total_molecules = vf_assay_molecules + ff_assay_molecules,
         
         GE_ml_plasma = total_molecules / plasma_volume_ml,
         
         report_acquired = ifelse(is.na(report_acquired), "No", 
                                  report_acquired))

# Numbers for paper
min(cfDNA_scd_outcomes$gestation_total_weeks)
max(cfDNA_scd_outcomes$gestation_total_weeks)
median(cfDNA_scd_outcomes$gestation_total_weeks)
min(cfDNA_scd_outcomes$fetal_percent)
max(cfDNA_scd_outcomes$fetal_percent)
median(cfDNA_scd_outcomes$fetal_percent)

#########################
# Results table
#########################

scd_cohort_table <- cfDNA_scd_outcomes %>%
  arrange(sampling_date_time) %>%
  mutate(sample_id = 
           paste0("HBB-", as.character(row.names(cfDNA_scd_outcomes)))) %>%
  select(sample_id, r_number, study_id, site, 
         sampling_date_time, hours_to_first_spin,
         days_to_storage,
         vacutainer, gestation_character, cohort, plasma_volume_ml, 
         extraction_replicates, 
         partner_sample_available, vf_assay, vf_assay_num_wells,  
         vf_assay_droplets, variant_positives,
         reference_positives, ff_assay, ff_assay_num_wells, ff_assay_droplets, 
         maternal_positives, paternal_positives, variant_molecules, 
         reference_molecules, vf_assay_molecules,
         maternal_molecules, paternal_molecules, total_molecules, GE_ml_plasma,
         variant_percent, fetal_percent, 
         major_allele_percent, upper_sprt_boundary,
         lower_sprt_boundary, likelihood_ratio,
         sprt_genotype_prediction, z_score, z_score_genotype_prediction, 
         diagnostic_sampling,
         mutation_genetic_info_fetus,
         report_acquired,
         outcome_zscore, outcome_sprt) %>%
  # Rename columns for ease of reading
  dplyr::rename(
    gestation = gestation_character, 
    "invasive_genotype" = mutation_genetic_info_fetus,
    "HBB_GE" = vf_assay_molecules,
    "variant_fraction_assay" = vf_assay,
    "fetal_fraction_assay" = ff_assay,
    "fetal_fraction" = fetal_percent,
    "variant_fraction" = variant_percent,
    "variant_allele_positives" = variant_positives,
    "reference_allele_positives" = reference_positives,
    "maternal_allele_positives" = maternal_positives,
    "paternal_allele_positives" = paternal_positives,
    "variant_allele_molecules" = variant_molecules,
    "reference_allele_molecules" = reference_molecules,
    "maternal_allele_molecules" = maternal_molecules,
    "paternal_allele_molecules" = paternal_molecules)


# Change the table to exclude HbAC and twin pregnancy
scd_cohort_table[scd_cohort_table$r_number == 17004, "outcome_zscore"] <- 
  "excluded: HbAC sample"
scd_cohort_table[scd_cohort_table$r_number == 17004, "outcome_sprt"] <- 
  "excluded: HbAC sample"
scd_cohort_table[scd_cohort_table$r_number == 20915, "outcome_zscore"] <- 
  "excluded: twin pregnancy"
scd_cohort_table[scd_cohort_table$r_number == 20915, "outcome_sprt"] <- 
  "excluded: twin pregnancy"

# Export table with time stamp
write.csv(scd_cohort_table, 
          file = (paste0("analysis_outputs/Supplementary_data ",
                         format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")),
          row.names = FALSE)

#########################
# Sensitivity and specificity
#########################

# Z score: balanced vs unbalanced
count(cfDNA_scd_outcomes %>%
        filter(vf_assay_molecules > 4000 &
                 fetal_percent > 4), z_score_genotype_prediction, 
      mutation_genetic_info_fetus)

# True positives (8+10), false positives (1), false negatives (1),
# true negatives (37)
scd_zscore_data <- as.table(matrix(c(18, 1, 1, 37), nrow = 2, byrow = TRUE))
scd_zscore_metrics <- epi.tests(scd_zscore_data, conf.level = 0.95)

# SPRT: balanced vs unbalanced
count(cfDNA_scd_outcomes, sprt_genotype_prediction, 
      mutation_genetic_info_fetus)

# True positives (15+12), false positives (3), false negatives (1),
# true negatives (50)
scd_sprt_data <- as.table(matrix(c(27, 3, 1, 50), nrow = 2, byrow = TRUE))
scd_sprt_metrics <- epi.tests(scd_sprt_data, conf.level = 0.95)

#########################
# Figure 2
#########################

# Consistent theme and axes for plots
multiplot_theme <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), 
  legend.position = "right", 
  plot.title = element_text(size = 11),
  legend.title = element_blank(),
  legend.text = element_text(size = 9))

multiplot_y <- ylim(39, 61)

multiplot_x <- scale_x_continuous(limits = c(0,30500),
                     breaks = c(0, 4000, 10000, 20000 ,30000))

## Lines
vertical_line <- geom_vline(xintercept = 4000, 
                            linetype = "dashed", alpha = 0.5)

z3_line <- geom_hline(yintercept = gDNA_mean_vp+(3*gDNA_stand_dev_vp),
                        linetype = "dashed", alpha = 0.5) 

zminus3_line <- geom_hline(yintercept = gDNA_mean_vp-(3*gDNA_stand_dev_vp),
                           linetype = "dashed", alpha = 0.5)

z2_line <- geom_hline(yintercept = gDNA_mean_vp+(2*gDNA_stand_dev_vp), 
           linetype = "dashed", alpha = 0.5)

zminus2_line <- geom_hline(yintercept = gDNA_mean_vp-(2*gDNA_stand_dev_vp),
             linetype = "dashed", alpha = 0.5)

# Plot 1: HbAS gDNA controls
gdna_plot_title <- expression(paste("ddPCR for 41 ", italic("HBB"), 
                  "c.20A>T gDNA controls"))

gdna_plot <- ggplot(gDNA_scd_data %>%
                      mutate(sample_type = "HbAS gDNA"), 
                    aes(x = vf_assay_molecules, 
                        y = variant_percent)) +
  vertical_line +
  z3_line +
  zminus3_line +
  z2_line +
  zminus2_line +
  multiplot_y +
  multiplot_x +
  
  scale_shape_manual(values = c(21)) +
  geom_point(size = 2, aes(shape = sample_type), fill= "white",
             colour = "black") +
  theme_bw() +
  multiplot_theme +
  labs(y = "Variant fraction (%)", x = "",
       title = gdna_plot_title)


# Plot 2: limit of detection gDNA mixtures
lod_plot <- ggplot(lod_data_merged, 
                   aes(x = vf_assay_molecules, y = variant_percent)) +
  theme_bw() +
  multiplot_theme + 
  vertical_line +
  z3_line +
  zminus3_line +
  z2_line +
  zminus2_line +
  multiplot_y +
  multiplot_x +
  
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
  
  labs(x = "",
       y = "Variant fraction (%)",
       title = "ddPCR for gDNA mixtures")

# Plot 3: cfDNA samples with z score analysis

cfdna_z_score_plot <- cfDNA_scd_outcomes %>%
  filter(!r_number %in% c("20915", "17004")) %>%
           mutate(mutation_genetic_info_fetus = 
                    paste0(mutation_genetic_info_fetus, " fetus"),
                  mutation_genetic_info_fetus = factor(mutation_genetic_info_fetus,
                                 levels = c("HbSS fetus","HbAS fetus",
                                            "HbAA fetus")),
                  outcome_zscore = factor(outcome_zscore, levels = c(
                    "correct", "incorrect", "insufficient data",
                    "inconclusive"))
                  ) %>%
                    ggplot(aes(x = vf_assay_molecules, 
                     y = variant_percent)) +
  
  theme_bw() +
  multiplot_theme + 
  vertical_line +
  z3_line +
  zminus3_line +
  z2_line +
  zminus2_line +
  multiplot_y +
  multiplot_x +
  
  scale_fill_manual(values=c("#FFFFFF", "#000000", "#999999",
                             "#999999"), guide = "none") +
  scale_alpha_manual(values = c(1, 1, 0.2, 0.2), guide = "none") +
  scale_shape_manual(values = c(24, 21, 25)) +
  geom_point(size = 2, aes(fill = outcome_zscore,
                           alpha = outcome_zscore,
                           shape = mutation_genetic_info_fetus),
             colour = "black") +
  
  labs(y = "Variant fraction (%)", x = "Genome equivalents (GE)", 
       title = "ddPCR for 88 cfDNA samples with z score classification")

# Plot 4: cfDNA samples with SPRT analysis
cfdna_sprt_plot <- cfDNA_scd_outcomes %>%
  filter(!r_number %in% c("20915", "17004")) %>%
  mutate(mutation_genetic_info_fetus = 
           paste0(mutation_genetic_info_fetus, " fetus"),
         mutation_genetic_info_fetus = factor(mutation_genetic_info_fetus,
                                              levels = c("HbSS fetus","HbAS fetus",
                                                         "HbAA fetus")),
         outcome_sprt = factor(outcome_sprt, levels = 
                                 c("correct", "incorrect", 
                                   "inconclusive"))) %>%
  ggplot(aes(x = vf_assay_molecules, 
             y = variant_percent)) +
  
  theme_bw() +
  multiplot_theme + 
  vertical_line +
  z3_line +
  zminus3_line +
  z2_line +
  zminus2_line +
  multiplot_y +
  multiplot_x +
  
  scale_fill_manual(values=c("#FFFFFF", "#000000", "#999999"), 
                    guide = "none") +
  scale_alpha_manual(values = c(1, 1, 0.2), guide = "none") +
  scale_shape_manual(values = c(24, 21, 25)) +
  geom_point(size = 2, aes(fill = outcome_sprt,
                           alpha = outcome_sprt,
                           shape = mutation_genetic_info_fetus),
             colour = "black") +
  
  labs(y = "Variant fraction (%)", x = "Genome equivalents (GE)", 
       title = "ddPCR for 88 cfDNA samples with SPRT classification") 

# Display 4 plots together

multi_plot <- ggpubr::ggarrange(gdna_plot, lod_plot, 
                                cfdna_z_score_plot,
                                cfdna_sprt_plot,
                                labels = c("2A", "2B", "2C", "2D"),
                                ncol = 2, nrow = 2, align = "v")

ggsave(plot = multi_plot, 
       filename = "figure_2.tiff",
       path = "plots/", device='tiff', dpi=600,
       units = "in",
       width = 12.5,
       height = 7)