################################################################################
## ddPCR for Non Invasive Prenatal Testing of Sickle Cell Disease
## July 2021
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

## Set working directory
setwd("W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/ddPCR_R_Analysis/ddpcr_nipd")

# Source functions (this includes loading ddPCR data)
source("functions/ddPCR_nipd_functions.R")
source("functions/RAPID_biobank.R")

#########################
# Select sickle cell disease cfDNA ddPCR data
#########################

# Samples to exclude:
# 13262 - this sample had contamination
# 17004 - this sample was actually HbAC
# 20915 - this sample was from a twin pregnancy
samples_to_exclude <- c(13262, 17004, 20915)

# The "secondary cohort" refers to the most recent phase of the 
# ddPCR sickle cell disease project, when all samples were extracted 
# using a 6ml protocol.

secondary_cohort <- c("14182", "19868", "20238", "20611", 
                      "20874", "30063", "30068", "30113", "30142", 
                      "30206", "30228", "30230", "30078", "30065", 
                      "13402", "20939", "30215", "30203",
                      "20911", "30236", "30112")

cfDNA_scd_data <- ff_calculations(
  var_ref_calculations(cfdna_ddpcr_data)) %>%
  dplyr::rename(r_number = sample) %>%
  filter(vf_assay == "HBB c.20A>T" &
           !r_number %in% samples_to_exclude) %>%
  mutate(extraction_volume = ifelse(r_number %in% secondary_cohort,
                         "6ml", "2 or 4ml"),
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
  dplyr::rename(r_number = sample)

#########################
# SPRT with HbAS gDNA controls
#########################

# This plot shows that the common form of the SPRT equation is 
# inappropriate for this dataset

ggplot(gDNA_scd_data, aes(x = vf_assay_molecules, y = variant_percent))+
  geom_point(size = 2, pch=21, alpha =0.6)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none")+
  labs(x = "Genome equivalents (GE)",
       y = "Variant fraction (%)",
       title = "Heterozygous gDNA samples") +
  xlim(0, 20000)+
  scale_y_continuous(breaks = c(40, 50, 60), 
                     limits = c(40,60))+
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
  annotate(geom = "text", x = 20000, y = 50, label = "HbAS")+
  annotate(geom = "text", x = 20000, y = 47, label = "HbAA") +
  annotate(geom = "text", x = 20000, y = 53, label = "HbSS")

#########################
# Variation in HbAS gDNA controls
#########################

vf_assay_molecules_limit <- 4000

# This section shows the variation in HbAS gDNA controls over 
# 4000 GE, in preparation for z score analysis in the next section

gDNA_scd_data_4000 <- gDNA_scd_data %>%
  filter(vf_assay_molecules > 4000)

# vp is "variant percent"
gDNA_mean_vp <- mean(gDNA_scd_data_4000$variant_percent)            

# Find the standard deviation, by using the "sd" function in stats package. 
# The stanard deviation is the square root of the variance.
# The variance is the average of the squared differences.

gDNA_stand_dev_vp <- sd(gDNA_scd_data_4000$variant_percent)

# Coefficient of variation
gDNA_cv_vp <- (gDNA_stand_dev_vp/gDNA_mean_vp) * 100

# Plot controls with normal distribution curve
gDNA_scd_data %>%
  mutate(variant_percent = round(variant_percent/0.5)*0.5) %>%
  filter(vf_assay_molecules > 4000) %>%
  ggplot(aes(x = variant_percent, y = )) +
  geom_bar(aes(y = (..count..)/sum(..count..))) +
  xlim(45, 55) +
  stat_function(fun = dnorm, n = 171, args = list(mean = 49.9, sd = 0.57),
                linetype = "dashed") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  labs(y = "Percentage of controls", x = "Variant percent (%)",
       title = "Distrubution of HbAS gDNA controls")

# gDNA controls with z score limits
ggplot(gDNA_scd_data, aes(x = vf_assay_molecules, 
                               y = variant_percent)) +
  geom_point(pch = 21, size = 3) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.position = "bottom",
    legend.title = element_blank()) +
  labs(y = "Variant fraction (%)", x = "Genome equivalents (GE)",
       title = "Heterozygous gDNA samples") +
  ylim(40, 60) +
  xlim(0, 21000) +
  geom_vline(xintercept = vf_assay_molecules_limit, linetype = "dashed") +
  geom_hline(yintercept = gDNA_mean_vp+(2*gDNA_stand_dev_vp), 
             linetype = "dashed") +
  geom_hline(yintercept = gDNA_mean_vp-(2*gDNA_stand_dev_vp),
             linetype = "dashed") +
  geom_hline(yintercept = gDNA_mean_vp+(3*gDNA_stand_dev_vp),
             linetype = "dashed") +
  geom_hline(yintercept = gDNA_mean_vp-(3*gDNA_stand_dev_vp),
             linetype = "dashed") +
  annotate(geom = "text", x = 21000, y = 50, 
           label = "HbAS") +
  annotate(geom = "text", x = 21000, y = 47, 
           label = "HbAA") +
  annotate(geom = "text", x = 21000, y = 53, 
           label = "HbSS")


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

# Plot the LOD data with Z scores and Z score thresholds
ggplot(lod_data_merged, 
       aes(x = vf_assay_molecules, y = variant_percent)) +
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
  geom_point(size = 3, aes(fill = sample), pch=21) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right", 
    legend.title = element_blank()) + 
  labs(x = "Genome equivalents (GE)",
       y = "Variant fraction (%)",
       title = "Limit of detection experiment Z score") +
  geom_vline(xintercept = vf_assay_molecules_limit, linetype = "dashed") +
  geom_hline(yintercept = gDNA_mean_vp+(2*gDNA_stand_dev_vp), 
             linetype = "dashed") +
  geom_hline(yintercept = gDNA_mean_vp-(2*gDNA_stand_dev_vp),
             linetype = "dashed") +
  geom_hline(yintercept = gDNA_mean_vp+(3*gDNA_stand_dev_vp),
             linetype = "dashed") +
  geom_hline(yintercept = gDNA_mean_vp-(3*gDNA_stand_dev_vp),
             linetype = "dashed") +
  annotate(geom = "text", x = 21000, y = 50, 
           label = "HbAS") +
  annotate(geom = "text", x = 21000, y = 47, 
           label = "HbAA") +
  annotate(geom = "text", x = 21000, y = 53, 
           label = "HbSS") +
  xlim(0, 21000) +
  ylim(38, 62)

#########################
# Predict fetal genotypes using Z score analysis
#########################

# The aim is to classify samples using z score analysis, in the same
# way as NIPT for trisomy 21.

lr_threshold <- 8

colnames(cfDNA_scd_predictions)

cfDNA_scd_predictions <- cfDNA_scd_data %>%
  mutate(
    z_score = (variant_percent - gDNA_mean_vp) / gDNA_stand_dev_vp,
    z_score_major = (major_allele_percent - gDNA_mean_vp) / gDNA_stand_dev_vp,
    
    z_score_genotype_prediction = case_when(
      # Samples with fetal fractions below 4% are inconclusive
      
      z_score > 3 &
      fetal_percent > 4 &
      vf_assay_molecules > vf_assay_molecules_limit
      ~"HbSS",
      
      z_score < -3 &
      fetal_percent > 4 &
      vf_assay_molecules > vf_assay_molecules_limit
      ~"HbAA",
      
      z_score < 2 &
      z_score > -2 &
      fetal_percent > 4 &
      vf_assay_molecules > vf_assay_molecules_limit
      ~"HbAS",
      
      TRUE ~"inconclusive"),
    
    z_score_clinical_prediction = case_when(
      
      z_score > 3 &
        fetal_percent > 4 &
        vf_assay_molecules > vf_assay_molecules_limit
      ~"affected",
      z_score < 2 &
        fetal_percent > 4 &
        vf_assay_molecules > vf_assay_molecules_limit
      ~ "unaffected",
      TRUE ~"inconclusive"),
    
    # Factorise for plot
    z_score_genotype_prediction = factor(z_score_genotype_prediction, levels = 
                                   c("HbSS", "HbAS", "HbAA",
                                     "inconclusive")),
    
    # Calculate likelihood ratio
    likelihood_ratio = calc_lr_autosomal(fetal_fraction,
                                         major_allele_percent/100,
                                         vf_assay_molecules),
    sprt_genotype_prediction = case_when(
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

#########################
# Compare predictions against Biobank
#########################

cfDNA_scd_outcomes <- left_join(
  cfDNA_scd_predictions,
  RAPID_biobank %>%
    mutate(r_number = as.character(r_number)) %>%
    select(r_number, study_id, site, maternal_DOB, 
           original_plasma_vol,
           date_of_blood_sample, Gestation_total_weeks,
           gestation_character, vacutainer,
           mutation_genetic_info_fetus, Partner_sample_available,
           report_acquired),
  by = "r_number") %>%
  mutate(outcome = case_when(
    z_score_genotype_prediction == "inconclusive" ~"inconclusive",
    is.na(mutation_genetic_info_fetus) ~"awaiting result",
    z_score_genotype_prediction == mutation_genetic_info_fetus
    ~"correct",
    TRUE ~"incorrect"))

#########################
# Plot cfDNA results
#########################

ggplot(cfDNA_scd_outcomes, aes(x = vf_assay_molecules, 
                                  y = variant_percent)) +
  geom_errorbar(aes(ymin = variant_percent_min, ymax = variant_percent_max),
                alpha = 0.2) +
  geom_errorbarh(aes(xmin = vf_assay_molecules_min, 
                     xmax = vf_assay_molecules_max),
                 alpha = 0.2) +
  # "#000000" = black; "#99CCFF" = light blue
  # "#0000FF" = dark blue; "#999999" = grey
  scale_fill_manual(values=c("#000000", "#99CCFF", "#0000FF", "#999999")) +
  scale_alpha_manual(values = c(1, 1, 1, 0.4)) +
  geom_point(size = 2, aes(fill = z_score_genotype_prediction,
                           alpha = z_score_genotype_prediction),
             colour = "black", pch=21) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.position = "bottom",
    legend.title = element_blank()) +
  labs(y = "Variant fraction (%)", x = "Genome equivalents (GE)",
       title = "cfDNA fetal genotype predictions") +
  ylim(40, 60) +
  xlim(0, 31000) +
  geom_vline(xintercept = vf_assay_molecules_limit, linetype = "dashed") +
  geom_hline(yintercept = gDNA_mean_vp+(2*gDNA_stand_dev_vp), 
                          linetype = "dashed") +
  geom_hline(yintercept = gDNA_mean_vp-(2*gDNA_stand_dev_vp),
                          linetype = "dashed") +
  geom_hline(yintercept = gDNA_mean_vp+(3*gDNA_stand_dev_vp),
                          linetype = "dashed") +
  geom_hline(yintercept = gDNA_mean_vp-(3*gDNA_stand_dev_vp),
                          linetype = "dashed") +
  
  # Circle the incorrect predictions
  geom_point(data=subset(cfDNA_scd_outcomes, 
                         r_number %in% c("20763")),
             pch=1, size=5, colour = "red", alpha = 0.8)

#########################
# Results table
#########################

scd_cohort_table <- cfDNA_scd_outcomes %>%
  select(sample_id, r_number, study_id, site, 
         date_of_blood_sample, maternal_DOB, original_plasma_vol,
         vacutainer, gestation_character, extraction_volume, 
         Partner_sample_available, vf_assay, 
         vf_assay_droplets, variant_positives,
         reference_positives, ff_assay, ff_assay_droplets, 
         maternal_positives, paternal_positives, variant_molecules, 
         reference_molecules, vf_assay_molecules,
         maternal_molecules, paternal_molecules,
         variant_percent, fetal_percent, z_score,
         z_score_genotype_prediction, z_score_clinical_prediction,
         mutation_genetic_info_fetus, report_acquired, outcome) %>%
  dplyr::rename(
    gestation = gestation_character, 
    "invasive genotype" = mutation_genetic_info_fetus,
    "HBB GE measured" = vf_assay_molecules)

write.csv(scd_cohort_table, "analysis_outputs/Supplementary Table 1.csv",
          row.names = FALSE)

#########################
# Sensitivity and specificity
#########################

# Use the epiR package to calculate sensitivity
count(cfDNA_scd_outcomes, z_score_genotype_prediction, 
      mutation_genetic_info_fetus)

# Balanced vs unbalanced
# True positives (9+10), false positives (1), false negatives (0),
# true negatives (32)
scd_data <- as.table(matrix(c(19, 1, 0, 32), nrow = 2, byrow = TRUE))
scd_metrics <- epi.tests(scd_data, conf.level = 0.95)

#########################
# Z score plotting
#########################

gDNA_cfDNA <- rbind(gDNA_scd_data_4000 %>%
  mutate(z_score = (variant_percent - gDNA_mean_vp) / gDNA_stand_dev_vp) %>%
  select(worksheet_well_sample, z_score, sample_type) %>%
  dplyr::rename(sample = worksheet_well_sample),
  
  cfDNA_scd_outcomes %>%
    filter(vf_assay_molecules > 4000 &
             fetal_percent > 4) %>%
    select(r_number, z_score_major, mutation_genetic_info_fetus) %>%
    filter(!is.na(mutation_genetic_info_fetus)) %>%
    dplyr::rename(sample_type = mutation_genetic_info_fetus,
                  sample = r_number,
                  z_score = z_score_major)) %>%
    
    mutate(sample_type = factor(sample_type, levels = c("gDNA", "HbAS",
                                                        "HbAA", "HbSS"))) %>%
    arrange(sample_type) %>%
    mutate(sample = factor(sample))

ggplot(gDNA_cfDNA %>%
         filter(sample_type != "gDNA"), aes(x = sample_type, y = z_score))+
  geom_jitter(size = 3, aes(shape = sample_type)) + 
  scale_shape_manual(values = c(1, 0, 3)) +
  theme_bw() +
  labs(y = "Z statistic", x = "") +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.text.x = element_blank(), 
    legend.title = element_blank()) +
  ylim(-5, 20)

#########################

