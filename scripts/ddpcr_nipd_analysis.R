#############################################################
## ddPCR for Non Invasive Prenatal Diagnosis (NIPD)
## January 2021
## Joseph.Shaw@gosh.nhs.uk
## Analysis script for prediction of fetal genotypes from 
## cfDNA testing using ddPCR for maternally-inherited
## variants.
#############################################################

#############################################################
# Load libraries and resources
#############################################################

## Load necessary packages
library(tidyverse)

## Set working directory
setwd("W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/ddPCR_R_Analysis/ddpcr_nipd")

# Load resources
controls <- readr::read_csv("resources/controls.csv")
ddpcr_target_panel <- readr::read_csv("resources/ddpcr_target_panel.csv") 

# Source functions
source("functions/ddPCR_nipd_functions.R")

#############################################################
# Read in cfDNA data and wrangle into shape
#############################################################

## Load in all the csv files exported from QuantaSoft.

dataPath <- "data/ddpcr_data/"

ddpcr_files <- list.files(dataPath)

#Empty data frame
ddpcr_data <- data.frame()

# Read and collate each worksheet csv
for (dataFile in ddpcr_files){
  tmp_dat <- read_csv(paste0(dataPath,dataFile), col_names = TRUE)
  # Add a worksheet identifier to make things easier later on.
  tmp_dat_ws <- tmp_dat %>%
    mutate(Worksheet = dataFile) %>%
    mutate(Worksheet_well = paste(Worksheet, Well, sep = "_"))
  ddpcr_data <-rbind(ddpcr_data, tmp_dat_ws)
  rm(tmp_dat_ws)
  rm(tmp_dat)
}

# Remove single wells and controls
ddpcr_data_merged_samples <- ddpcr_data %>%
  # Remove single wells
  drop_na("MergedWells") %>%
  # Remove controls
  filter(!(Sample %in% controls$Controls))

# Reshape the data frame to sum all values by Target
ddpcr_data_reshape <- ddpcr_data_merged_samples %>% 
  group_by(Sample, Target) %>% 
  summarise(Positives = sum(Positives),
            AcceptedDroplets = sum(AcceptedDroplets),
            .groups="drop")

ddpcr_with_target <- ddpcr_data_reshape %>% 
  left_join(ddpcr_target_panel %>%
              select(Target, Target_category), by = "Target")

# Use pivot_wider to get one row per sample
pivotted_ddpcr <- ddpcr_with_target %>% 
  pivot_wider(id_cols = Sample,
              names_from = Target_category,
              values_from = c(AcceptedDroplets, Positives))

# Add on assay and inheritance pattern to the table.
ref_table_var <- left_join(
  # First table
  ddpcr_with_target %>%
    filter(Target_category == "variant"),
  # Second table
  ddpcr_target_panel %>%
    select(Assay, Inheritance, Target) %>%
    # Rename the assay column to avoid clash with fetal fraction assay
    rename(variant_assay = Assay),
  # Join by
  by = "Target")
  
ref_table_ff <- left_join(
  # First table
  ddpcr_with_target %>%
    filter(Target_category == "ff_allele1"),
  # Second table
  ddpcr_target_panel %>%
    select(Assay, Target) %>%
    rename(ff_assay = Assay),
  # Join by
  by = "Target")

ddpcr_data_tbl <- pivotted_ddpcr %>%
  left_join(ref_table_var %>%
              select(Sample, Inheritance, variant_assay), by = "Sample",
            .groups="drop") %>% 
  left_join(ref_table_ff %>%
              select(Sample, ff_assay), by = "Sample",
            .groups="drop") %>%
  
  # Remove any samples which haven't had both assays performed.
  filter(!is.na(Positives_variant) & !is.na(Positives_ff_allele1)) %>%
  
  # Remove duplicate columns and rename to be compatible with functions
  select(-c(AcceptedDroplets_ff_allele2, AcceptedDroplets_reference)) %>%
  rename(AcceptedDroplets_FetalFrac = AcceptedDroplets_ff_allele1) %>%
  rename(AcceptedDroplets_Variant_assay = AcceptedDroplets_variant) %>%
  
  # Determine the maternal and paternally-inherited alleles for the fetal fraction calculation
  mutate(Positives_maternal = pmax(Positives_ff_allele1, Positives_ff_allele2)) %>%
  mutate(Positives_paternal = pmin(Positives_ff_allele1, Positives_ff_allele2))

#############################################################
# Read in gDNA data and wrangle
#############################################################

# Get single well controls only without NTC
ddpcr_controls <- ddpcr_data %>%
  filter(Sample %in% controls$Controls & Sample != "NTC" & !is.na(CopiesPer20uLWell))

ddpcr_controls_with_target <- ddpcr_controls %>% 
  left_join(ddpcr_target_panel %>%
              select(Target, Target_category), by = "Target")

# No need to pivot to get one row per sample because we want the control
# data as individual wells.
# Use pivot_wider to get one row per well for each well tested with 
# a variant assay
pivotted_controls_var <- ddpcr_controls_with_target %>%
  filter(Target_category %in% c("variant", "reference")) %>%
  pivot_wider(id_cols = c(Worksheet_well, Sample),
              names_from = Target_category,
              values_from = c(AcceptedDroplets, Positives)) %>%
  select(-AcceptedDroplets_reference)

# Use pivot_wider to get one row per well for each well tested with 
# a fetal fraction assay
pivotted_controls_ff <- ddpcr_controls_with_target %>%
  filter(Target_category %in% c("ff_allele2", "ff_allele1")) %>%
  pivot_wider(id_cols = c(Worksheet_well, Sample),
              names_from = Target_category,
              values_from = c(AcceptedDroplets, Positives)) %>%
  select(-AcceptedDroplets_ff_allele2)

# Add on assay and inheritance pattern to the variant control table.
control_table_var <- left_join(
  # First table
  ddpcr_controls_with_target %>%
    filter(Target_category == "variant") %>%
    select(Worksheet_well, Sample, Target),
  # Second table
  ddpcr_target_panel %>%
    select(Assay, Inheritance, Target) %>%
    # Rename the assay column to avoid clash with fetal fraction assay
    rename(variant_assay = Assay),
  # Join by
  by = "Target")

# Add on assay to the fetal fraction control table.
control_table_ff <- left_join(
  # First table
  ddpcr_controls_with_target %>%
    filter(Target_category == "ff_allele1") %>%
    select(Worksheet_well, Sample, Target),
  # Second table
  ddpcr_target_panel %>%
    select(Assay, Target) %>%
    rename(ff_assay = Assay),
  # Join by
  by = "Target")

# Get the control information for the control wells tested with a variant assay
ddpcr_control_tbl_var <- pivotted_controls_var %>%
  left_join(control_table_var %>%
              select(Worksheet_well, Inheritance, variant_assay), by = "Worksheet_well",
            .groups="drop") %>%
  # Remove duplicate columns and rename to be compatible with functions
  rename(AcceptedDroplets_Variant_assay = AcceptedDroplets_variant) %>%
  mutate(Molecules_variant = Poisson_correct(AcceptedDroplets_Variant_assay,Positives_variant)) %>%
  mutate(Molecules_reference = Poisson_correct(AcceptedDroplets_Variant_assay,Positives_reference))

# Get the control information for the control wells tested with a fetal fraction assay
ddpcr_control_tbl_ff <- pivotted_controls_ff %>%
  left_join(control_table_ff %>%
              select(Worksheet_well, ff_assay), by = "Worksheet_well",
            .groups="drop") %>%
  # Remove duplicate columns and rename to be compatible with functions
  rename(AcceptedDroplets_FetalFrac = AcceptedDroplets_ff_allele1) %>%
  mutate(Positives_maternal = pmax(Positives_ff_allele1, Positives_ff_allele2)) %>%
  mutate(Positives_paternal = pmin(Positives_ff_allele1, Positives_ff_allele2)) %>%
  mutate(Molecules_maternal = Poisson_correct(AcceptedDroplets_FetalFrac,Positives_maternal)) %>%
  mutate(Molecules_paternal = Poisson_correct(AcceptedDroplets_FetalFrac,Positives_paternal))

#############################################################
# Sickle cell disease cfDNA analysis
#############################################################

# Ammend the fetal genotype prediction for samples which failed QC steps
samples_failing_qc <- c(13262, 20763, 20810)

# Analyse all the samples using Poisson correction and SPRT.
sickle_cell_analysed <- calc_SPRT(calc_conf_intervals(calc_molecules(ddpcr_data_tbl %>%
                                                                     filter(variant_assay == "HBB c.20A>T"))), 250) %>%
  # Add in a column for the overall prediction based on the data.
  mutate(overall_prediction = ifelse(Sample %in% samples_failing_qc, "no call", SPRT_prediction)) %>%
  
  # Add in a "call" column to allow comparison of inconclusive and conclusive results
  mutate(call = ifelse(overall_prediction == "no call", "no call", "call")) %>%
  
  # Remove sample 13262 as Natalie thinks this should not go in the paper.
  filter(Sample != "13262") %>%
  
  # Rename sample to r_number to allow merge with RAPID Biobank data in next step.
  # Have to specify dplyr for rename.
  dplyr::rename(r_number = Sample) 

# Ammend the dataframe for the sample from a twin pregnancy (20915).

sickle_cell_analysed$SPRT_prediction[sickle_cell_analysed$r_number == 20915] <- "twin pregnancy"
sickle_cell_analysed$overall_prediction[sickle_cell_analysed$r_number == 20915] <- "twin pregnancy"

# Add on a generic identifier for the samples.
Identifier <- paste0("HBB", "-", rownames(sickle_cell_analysed))

sickle_cell_blinded <- cbind(Identifier, sickle_cell_analysed)

# Compare results to those in the RAPID biobank.
sickle_cell_cohort <- RAPID_biobank %>%
  filter(r_number %in% sickle_cell_blinded$r_number) %>%
  select(r_number, study_id, gestation_weeks, gestation_days, date_of_blood_sample, vacutainer, mutation_genetic_info_fetus, 
         Partner_sample_available)

# Need to convert the biobank mutation information into consistent strings

HbAA_genotypes <- c("HbAA", "HbAA. Normal PCR/karyotype.", "HbAA_HBB c.20A hom, 46,XX", 
                    "HbAA_HBB c.20A hom, 46,XY", "HbAA_HBB c.20A hom")

HbAS_genotypes <- c("HbAS", "HbAS_HBB c.20A>T het", "HbAS_HBB c.20A>T het, 46,XY", 
                    "HbAS_HBB c.20A>T het, 46,XX", "HbAS, normal PCR", "HbAS_HBB c.20A>T het, mos 46,XY,?inv(12)(p11.2q13)[5]/46,XY[40] on CVS",
                    "HbAS, 46,XY", "HbAS, PCR normal", "HbAS_HBB c.20A>T het, PCR normal",
                    "47,XXY (also affected with Klinefelter's), HbAS_HBB c.20A>T het")

HbSS_genotypes <- c("HbSS_HBB c.20A>T hom, 46,XX", "HbSS_HBB c.20A>T hom", "HbSS", "Affected")

HbAC_genotypes <- c("HbAC")

# Merge the two dataframes together by r_number.

sickle_cell_unblinded <- merge(sickle_cell_blinded, sickle_cell_cohort, by = "r_number") %>%
  
  # Add a consistent genotype column
  mutate(invasive_result = case_when(
    mutation_genetic_info_fetus %in% HbAA_genotypes ~"HbAA",
    mutation_genetic_info_fetus %in% HbAS_genotypes ~"HbAS",
    mutation_genetic_info_fetus %in% HbSS_genotypes ~"HbSS",
    mutation_genetic_info_fetus %in% HbAC_genotypes ~"HbAC",
    TRUE ~mutation_genetic_info_fetus)) %>%
  
  # Reorganise the columns to make it easier to look at.
  select(Identifier, r_number, study_id, gestation_weeks, gestation_days, date_of_blood_sample, 
         Partner_sample_available, vacutainer, AcceptedDroplets_Variant_assay, Positives_variant,
         Positives_reference, ff_assay, AcceptedDroplets_FetalFrac, Positives_maternal, Positives_paternal, 
         Molecules_variant, Molecules_reference, Molecules_maternal, Molecules_paternal, Molecules_variant_assay,
         Fetal_fraction_max_percent, Fetal_fraction_percent, Fetal_fraction_min_percent,
         Variant_fraction_max_percent, Variant_fraction_percent, Variant_fraction_min_percent,
         Likelihood_ratio, SPRT_prediction, overall_prediction, mutation_genetic_info_fetus, invasive_result)
  
# Export the inconclusives for the paper table
inconclusives_only <- sickle_cell_unblinded %>%
  filter(overall_prediction == "no call") %>%
  select(Identifier, r_number, gestation_weeks, gestation_days, Molecules_variant_assay, 
         Variant_fraction_percent, Fetal_fraction_percent, 
         Likelihood_ratio, SPRT_prediction, overall_prediction, invasive_result)

# Rename the columns to be more helpful
colnames(inconclusives_only) <- c("Sample", "Research number", "Gestation_weeks", "Gestation_days",	"DNA molecules at HBB c.20",	"Variant fraction (%)",	"Fetal fraction (%)",
                                  "Likelihood ratio",	"SPRT prediction",	"Overall prediction",	"Invasive result")

# Are the fetal fractions significantly different in inconclusive and conclusive samples?
call_fetal_fractions <- sickle_cell_blinded %>%
  filter(call == "call")

no_call_fetal_fractions <- sickle_cell_blinded %>%
  filter(call == "no call")

t.test(call_fetal_fractions$Fetal_fraction_percent, no_call_fetal_fractions$Fetal_fraction_percent,
       paired = FALSE)

#############################################################
# Sickle cell gDNA analysis
#############################################################

# Worksheets containing 20ng 19RG-220G0190 per replicate
worksheets_20ng <- c("20-4046.csv", "20-3824.csv", "20-3756_FAM_HEX.csv", "20-4321.csv", 
                     "20-4317.csv", "20-4227.csv")

het_DNA_ids <- read.csv("resources/het_DNA_ids.csv")

mat_het_gDNA <- ddpcr_data %>%
  filter(Sample == "19RG-220G0190") %>%
  filter(is.na(MergedWells)) %>%
  filter(Worksheet %in% worksheets_20ng) %>%
  filter(Target == "HbS") %>%
  select(Worksheet_well, Positives, FractionalAbundance, 
         PoissonFractionalAbundanceMax, PoissonFractionalAbundanceMin)

pat_cfDNA <- ddpcr_data %>%
  filter(Sample %in% c(30139, 30130)) %>%
  filter(!is.na(MergedWells)) %>%
  filter(Target == "HbS") %>%
  select(Worksheet_well, Positives, FractionalAbundance, 
         PoissonFractionalAbundanceMax, PoissonFractionalAbundanceMin)

het_DNA <- rbind(mat_het_gDNA, pat_cfDNA)

het_gDNA_arranged <- left_join(het_DNA, het_DNA_ids, by = "Worksheet_well") %>%
  select(Control_id, FractionalAbundance, PoissonFractionalAbundanceMax, PoissonFractionalAbundanceMin)

colnames(het_gDNA_arranged) <- c("Identifier", "Variant_fraction_percent", "Variant_fraction_max_percent", 
                                 "Variant_fraction_min_percent")

mat_cfDNA_het <- sickle_cell_unblinded %>%
  filter(invasive_result == "HbAS") %>%
  arrange(overall_prediction) %>%
  select(Identifier, Variant_fraction_percent, Variant_fraction_max_percent, Variant_fraction_min_percent)

all_het_DNA <- rbind(het_gDNA_arranged, mat_cfDNA_het)

het_DNA_plot <- ggplot(het_gDNA_arranged, aes(x = Identifier, y = Variant_fraction_percent))+
  geom_point(size = 4)+
  geom_pointrange(aes(ymin = Variant_fraction_min_percent, ymax = Variant_fraction_max_percent))+
  geom_hline(aes(yintercept=50), linetype="dashed", size = 1)+
  geom_hline(aes(yintercept=51.2), linetype="dashed", colour = "grey", size = 1)+
  geom_hline(aes(yintercept=48.9), linetype="dashed", colour = "grey", size = 1)+
  ylim(40, 60)+
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        plot.title = element_text(size=20),
        legend.position = c(0.85, 0.85), legend.background = element_rect(fill="white"), 
        legend.title = element_text(size= 18), legend.text = element_text(size= 15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "", y = "Variant fractional abundance (%)")+
  theme(axis.text.x = element_text(angle = 90))

#############################################################
# Sickle cell disease limit of detection study
#############################################################

LOD_data <- read_csv("data/20-1557_LOD.csv", col_names = TRUE)

LOD_data_longer <- LOD_data %>%
  mutate(unique_identifier = paste(Input_molecules, Sample)) %>%
  pivot_wider(id_cols = c(unique_identifier, Sample, Input_molecules, Mass_molecules, fetal_fraction),
              names_from = Target,
              values_from = c(AcceptedDroplets, Positives, FractionalAbundance, PoissonFractionalAbundanceMax, 
                              PoissonFractionalAbundanceMin)) %>%
  select(-c(AcceptedDroplets_HbS, FractionalAbundance_HbA, PoissonFractionalAbundanceMax_HbA, PoissonFractionalAbundanceMin_HbA)) %>%
  rename(AcceptedDroplets_Variant_assay = AcceptedDroplets_HbA) %>%
  mutate(HbS_molecules = Poisson_correct(AcceptedDroplets_Variant_assay, Positives_HbS)) %>%
  mutate(HbA_molecules = Poisson_correct(AcceptedDroplets_Variant_assay, Positives_HbA)) %>%
  mutate(total_DNA_molecules = HbS_molecules + HbA_molecules) %>%
  mutate(Reference_fraction = HbA_molecules / total_DNA_molecules) %>%
  mutate(Variant_fraction = HbS_molecules / total_DNA_molecules) %>%
  mutate(Over_represented_fraction = pmax(Reference_fraction, Variant_fraction)) %>%
  mutate(Likelihood_ratio = calc_LR_autosomal(fetal_fraction, Over_represented_fraction, total_DNA_molecules)) %>%
  mutate(SPRT_prediction = case_when(
    Likelihood_ratio > 250 & Over_represented_fraction == Reference_fraction ~ "homozygous reference",
    Likelihood_ratio > 250 & Over_represented_fraction == Variant_fraction ~ "homozygous variant",
    Likelihood_ratio < (1/250) ~ "heterozygous",
    Likelihood_ratio < 250 & Likelihood_ratio > (1/250) ~ "no call"))

LOD_data_even_longer <- LOD_data_longer %>%
  select(Sample, Mass_molecules, Input_molecules, HbS_molecules, HbA_molecules, AcceptedDroplets_Variant_assay) %>%
  pivot_longer(cols = c(HbS_molecules, HbA_molecules), names_to = "Target", values_to = "molecules") %>%
  # Order the factors
  mutate(Sample = factor(Sample, levels = c("AA 12%", "AA 10%", "AA 8%", "AA 6%", "AA 4%","AA 2%",
                                            "0%", "SS 2%", "SS 4%", "SS 6%", 
                                            "SS 8%", "SS 10%", "SS 12%"))) %>%
  mutate(Mass_molecules = factor(Mass_molecules, levels = c("9.9ng (~3,000 molecules)", "19.8ng (~6,000 molecules)",
                                                            "29.7ng (~9,000 molecules)", "39.6ng (~12,000 molecules)"))) %>%
  
  # Add in the 95% Poisson confidence intervals
  
  # Calculate copies per droplet (cpd) for each allele
  mutate(Cpd = molecules / AcceptedDroplets_Variant_assay) %>%
  
  # Calculate the 95% confidence intervals of the numbers of molecules for each allele
  mutate(Molecules_max = Poisson_max(Cpd, AcceptedDroplets_Variant_assay)) %>%
  mutate(Molecules_min = Poisson_min(Cpd, AcceptedDroplets_Variant_assay)) %>%
  mutate(Percentage = case_when(
    Sample %in% c("SS 2%", "AA 2%") ~"2%",
    Sample %in% c("SS 4%", "AA 4%") ~"4%",
    Sample %in% c("SS 6%", "AA 6%") ~"6%",
    Sample %in% c("SS 8%", "AA 8%") ~"8%",
    Sample %in% c("SS 10%", "AA 10%") ~"10%",
    Sample %in% c("SS 12%", "AA 12%") ~"12%",
    Sample %in% c("0") ~"0%"))

# Plot the results for presentation
ggplot(LOD_data_even_longer %>%
         filter(Input_molecules %in% c(12000)), aes(x = Sample, y = molecules, fill = Target))+
  geom_col(width = 0.5, position = position_dodge(width =0.5), colour = "black", alpha = 0.6)+
  geom_errorbar(aes(ymin = Molecules_min, ymax = Molecules_max), width = .2, position=position_dodge(width=0.5))+
  scale_fill_manual(values=c("#3366FF", "#FF0000"), labels= c("Reference", "Variant"))+
  theme_bw()+
  theme(axis.text=element_text(size=10), axis.title = element_text(size=14), 
        plot.title = element_text(size=20), legend.position = "bottom", legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "", y = "Molecules detected by ddPCR")+
  ylim(0, 8000)+
  scale_x_discrete(labels=c("12%","10%","8%","6%","4%","2%","0%",
                            "2%","4%","6%","8%","10%","12%"))

# A different plot - not as nice to look at.
ggplot(LOD_data_longer, aes(x = total_DNA_molecules, y = FractionalAbundance_HbS, colour = Sample))+
  geom_point(size = 4)+
  geom_pointrange(aes(ymin = PoissonFractionalAbundanceMin_HbS, ymax = PoissonFractionalAbundanceMax_HbS))+
  geom_hline(yintercept=50, linetype="dashed", size = 1)

#############################################################
# Bespoke cohort analysis
#############################################################

bespoke_cohort <- ddpcr_data_tbl %>%
  filter(variant_assay != "HBB c.20A>T") %>%
  dplyr::rename(r_number = Sample)

bespoke_cohort$r_number <- as.numeric(bespoke_cohort$r_number)

AD_samples <- c(12990, 13519, 14116, 14491, 14522, 19261, 19711)
AR_samples <- c(17531, 19102)
XLR_samples <- c(10280, 11928, 12585, 13625, 13965, 14247, 
                 14917, 16319, 16468, 16881, 18164, 18385, 18891,
                 20817, 19611, 20980, 17667)
XLD_samples <- c(12945)

# Can remove the sample with a primer binding SNP if need be
# filter(r_number != 19397) %>%

# Perform Poisson correction and SPRT analysis
bespoke_cohort_analysed <- calc_SPRT(calc_conf_intervals(calc_molecules(bespoke_cohort)), 8) %>%

  mutate(Inheritance_pattern = case_when(
    r_number %in% AD_samples ~"AD",
    r_number %in% AR_samples ~"AR",
    r_number %in% XLR_samples ~"XLR",
    r_number %in% XLD_samples ~"XLD"))

#Arrange into a nice order
bespoke_cohort_analysed <- arrange(bespoke_cohort_analysed, Inheritance_pattern)

# Add on a generic identifier for the samples.
Identifier <- paste0("cfDNA", "-", rownames(bespoke_cohort_analysed))

bespoke_cohort_blinded <- cbind(Identifier, bespoke_cohort_analysed)

bespoke_cohort_unblinded <- left_join(bespoke_cohort_blinded,
                                      RAPID_biobank %>%
                                        select(r_number, study_id, gestation_weeks, gestation_days, 
                                               Gestation_total_weeks, date_of_blood_sample, 
                                               vacutainer, mutation_genetic_info_fetus),
                                      by = "r_number") %>%
  mutate(Call = ifelse(SPRT_prediction == "no call", "no call", "call"))


#############################################################
# cfDNA genomic equivalents per ml plasma (GE/ml)
#############################################################

# cfDNA concentration in plasma depends on the extraction volume,
# number of replicates and input volume, all of which was changed
# for different samples.

# 12585: ZFXY tested on 20-4420, Test9 and 19-3581, F8 tested on
# 20-4515 and Test7
# 16319: ZFXY tested on 19-3581 and Test9
# 16881: ZFXY tested on 20-4420, Test9 and 19-3581; F8 tested on
# 20-4515 and Test8
# 18164: ZFXY tested on Test9 and 20-3156
# 18385: ZFXY tested on 20-4420 and Test27; ATP7A tested on 
# 20-4468 and Test26
# 18891: ZFXY tested on 20-4126 and Test27; FOXP3 tested on
# 20-4126 and Test22

# Specify the sample_worksheet combinations to exclude
rows_to_exclude <- c("12585_Test9.csv", "12585_Test7.csv", "12585_19-3581.csv",
                     "16319_Test9.csv", "16881_Test9.csv", "16881_19-3581.csv",
                     "16881_Test8.csv", "18385_Test26.csv",
                     "18164_Test9.csv", "18385_Test27.csv", "18891_Test27.csv",
                     "18891_Test22.csv")

# This tibble can then replace the ddpcr_data_reshape tibble from earlier.
# This should get all the bespoke ddPCR cohort with one worksheet each for the
# variant and fetal fraction assays, so that the GE/ml can be calculated.

for_GE_calc <- ddpcr_data_merged_samples %>%
  filter(Sample %in% bespoke_cohort_analysed$r_number) %>%
  # Change concentration to a numeric
  mutate(Concentration = as.numeric(Concentration))%>%
  mutate(PoissonConfMax = as.numeric(PoissonConfMax))%>%
  mutate(PoissonCNVMin = as.numeric(PoissonCNVMin))%>%
  mutate(replicates = str_count(MergedWells, ",")+1) %>%
  mutate(Sample_worksheet = paste(Sample, Worksheet, sep = "_"))%>%
  # Filter out the repeat samples
  filter(!Sample_worksheet %in% rows_to_exclude) %>%
  arrange(Sample)

# Which samples/worksheets had odd elution or extraction volumes?
# 13625 on Test23 (variant assay) would have had a 2ml extraction
# 16319 on Test8 (variant assay) would have had a 2ml extraction
# 16468 on Test24 (variant assay) would have had a 2ml extraction
# 18891 on 20-4126 had 9ul input into each replicate for the variant assay
# 20817 would have had a 2ml extraction on both assays

# Add on the target category
for_GE_calc_with_target <- for_GE_calc %>% 
  left_join(ddpcr_target_panel %>%
              select(Target, Target_category), by = "Target") %>%
  select(Well, Sample, Target, Target_category, Concentration, PoissonConfMax, PoissonConfMin,
         Positives, AcceptedDroplets, Sample_worksheet, replicates) %>%
  arrange(Sample) %>%
  pivot_wider(id_cols = c(Sample),
              names_from = Target_category,
              values_from = c(Concentration, PoissonConfMax, PoissonConfMin, Positives,
                              AcceptedDroplets, replicates)) %>%
  mutate(ff_assay_replicate_input_ul = ifelse(Sample == 18891, 10, 5)) %>%
  mutate(variant_assay_replicate_input_ul = ifelse(Sample == 18891, 9, 5)) %>%
  mutate(ff_extraction_input_ml = ifelse(Sample %in% c(13625, 16319, 16468, 20817), 2, 4)) %>%
  mutate(variant_extraction_input_ml = ifelse(Sample %in% c(13625, 16319, 16468, 20817), 2, 4)) %>%
  
  # Determine which targets are the fetals signal and the shared maternal allele
  mutate(Positives_maternal = pmax(Positives_ff_allele1, Positives_ff_allele2)) %>%
  mutate(Positives_paternal = pmin(Positives_ff_allele1, Positives_ff_allele2)) %>%
  mutate(Concentration_maternal = pmax(Concentration_ff_allele1, Concentration_ff_allele2)) %>%
  mutate(Concentration_paternal = pmin(Concentration_ff_allele1, Concentration_ff_allele2)) %>%
  mutate(PoissonConfMax_maternal = pmax(PoissonConfMax_ff_allele1, PoissonConfMax_ff_allele2)) %>%
  mutate(PoissonConfMax_maternal = pmax(PoissonConfMax_ff_allele1, PoissonConfMax_ff_allele2)) %>%
  mutate(PoissonConfMin_maternal = pmax(PoissonConfMin_ff_allele1, PoissonConfMin_ff_allele2)) %>%
  mutate(PoissonConfMax_paternal = pmin(PoissonConfMax_ff_allele1, PoissonConfMax_ff_allele2)) %>%
  mutate(PoissonConfMin_paternal = pmin(PoissonConfMin_ff_allele1, PoissonConfMin_ff_allele2)) %>%
  
  # Get rid of duplicated columns 
  select(-c("AcceptedDroplets_ff_allele2", "AcceptedDroplets_reference", "replicates_ff_allele2",
            "replicates_reference")) %>%
  
  # Calculate the total cfDNA detected by each assay
  mutate(Positives_ff_assay = Positives_maternal + Positives_paternal) %>%
  mutate(Positives_variant_assay = Positives_reference + Positives_variant) %>%
  
  # Calculate concentration in copies per ul for each assay. Each droplet is 0.00085ul (0.85nl)
  mutate(Concentration_ff_assay = -log((AcceptedDroplets_ff_allele1 - Positives_ff_assay)/ 
                                         AcceptedDroplets_ff_allele1) / 0.00085) %>%
  mutate(Concentration_variant_assay = -log((AcceptedDroplets_variant - Positives_variant_assay)/ 
                                         AcceptedDroplets_variant) / 0.00085) %>%
  
  mutate(Concentration_ff_assay_max = Concentration_ff_assay + 1.96*(sqrt(Concentration_ff_assay / AcceptedDroplets_ff_allele1*0.00085))) %>%
  mutate(Concentration_ff_assay_min = Concentration_ff_assay - 1.96*(sqrt(Concentration_ff_assay / AcceptedDroplets_ff_allele1*0.00085))) %>%
  mutate(Concentration_variant_assay_max = Concentration_variant_assay + 1.96*(sqrt(Concentration_variant_assay / AcceptedDroplets_variant*0.00085))) %>%
  mutate(Concentration_variant_assay_min = Concentration_variant_assay - 1.96*(sqrt(Concentration_variant_assay / AcceptedDroplets_variant*0.00085))) %>%
  
  # Calculate the copies per ul of the cfDNA extraction elution
  mutate(cfDNA_elution_concentration_ff_assay = Concentration_ff_assay * (replicates_ff_allele1*22) /
           (replicates_ff_allele1*ff_assay_replicate_input_ul)) %>%
  mutate(cfDNA_GE_ml_plasma_ff_assay = (cfDNA_elution_concentration_ff_assay * 60) / ff_extraction_input_ml) %>%
  
  mutate(cfDNA_elution_concentration_ff_assay_max = Concentration_ff_assay_max * (replicates_ff_allele1*22) /
           (replicates_ff_allele1*ff_assay_replicate_input_ul)) %>%
  mutate(cfDNA_GE_ml_plasma_ff_assay_max = (cfDNA_elution_concentration_ff_assay_max * 60) / ff_extraction_input_ml) %>%
  
  mutate(cfDNA_elution_concentration_ff_assay_min = Concentration_ff_assay_min * (replicates_ff_allele1*22) /
           (replicates_ff_allele1*ff_assay_replicate_input_ul)) %>%
  mutate(cfDNA_GE_ml_plasma_ff_assay_min = (cfDNA_elution_concentration_ff_assay_min * 60) / ff_extraction_input_ml) %>%
  
  mutate(cfDNA_elution_concentration_variant_assay = Concentration_variant_assay * (replicates_variant*22) /
           (replicates_variant*variant_assay_replicate_input_ul)) %>%
  mutate(cfDNA_GE_ml_plasma_variant_assay = (cfDNA_elution_concentration_variant_assay * 60) / variant_extraction_input_ml) %>%
  
  mutate(cfDNA_elution_concentration_variant_assay_max = Concentration_variant_assay_max * (replicates_variant*22) /
           (replicates_variant*variant_assay_replicate_input_ul)) %>%
  mutate(cfDNA_GE_ml_plasma_variant_assay_max = (cfDNA_elution_concentration_variant_assay_max * 60) / variant_extraction_input_ml) %>%
  
  mutate(cfDNA_elution_concentration_variant_assay_min = Concentration_variant_assay_min * (replicates_variant*22) /
           (replicates_variant*variant_assay_replicate_input_ul)) %>%
  mutate(cfDNA_GE_ml_plasma_variant_assay_min = (cfDNA_elution_concentration_variant_assay_min * 60) / variant_extraction_input_ml) %>%
  
  mutate(cffDNA_elution_concentration_ff_assay = (Concentration_paternal * (replicates_ff_allele1*22)) /
           (replicates_ff_allele1*ff_assay_replicate_input_ul)) %>%
  
  mutate(cffDNA_elution_concentration_ff_assay_max = PoissonConfMax_paternal * (replicates_ff_allele1*22) / 
           (replicates_ff_allele1*ff_assay_replicate_input_ul)) %>%
  
  mutate(cffDNA_elution_concentration_ff_assay_min = PoissonConfMin_paternal * (replicates_ff_allele1*22) / 
           (replicates_ff_allele1*ff_assay_replicate_input_ul)) %>%
  
  mutate(cffDNA_GE_ml_plasma_ff_assay = (cffDNA_elution_concentration_ff_assay * 60) / ff_extraction_input_ml) %>%
  
  mutate(cffDNA_GE_ml_plasma_ff_assay_max = (cffDNA_elution_concentration_ff_assay_max * 60) / ff_extraction_input_ml) %>%
  
  mutate(cffDNA_GE_ml_plasma_ff_assay_min = (cffDNA_elution_concentration_ff_assay_min * 60) / ff_extraction_input_ml)
  
cfDNA_conc_gestations <- left_join(
  bespoke_cohort_unblinded, 
  for_GE_calc_with_target %>%
    mutate(r_number = as.numeric(Sample)),
  by = "r_number")

ggplot(cfDNA_conc_gestations, aes(x = cfDNA_GE_ml_plasma_variant_assay, y = cfDNA_GE_ml_plasma_ff_assay, 
                                  shape = Call))+
  geom_point(size = 5, alpha = 0.6)+
  scale_shape_manual(values = c(19, 1))+
  geom_errorbarh(aes(xmin = cfDNA_GE_ml_plasma_variant_assay_min, xmax = cfDNA_GE_ml_plasma_variant_assay_max)) +
  geom_errorbar(aes(ymin = cfDNA_GE_ml_plasma_ff_assay_min, ymax = cfDNA_GE_ml_plasma_ff_assay_max))+
  ylim(0, 7000)+
  xlim(0, 7000)+
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none")+
  geom_abline(linetype = "dashed", alpha = 0.5)+
  labs(x = "cfDNA concentration from variant/reference ddPCR (GE/ml plasma)",
       y = "cfDNA concentration from fetal fraction ddPCR (GE/ml plasma)")


ggplot(cfDNA_conc_gestations, aes(x = Gestation_total_weeks, y = cffDNA_GE_ml_plasma_ff_assay,
                                  shape = Call))+
  geom_point(size = 3)+
  scale_shape_manual(values = c(19, 1))+
  geom_errorbar(aes(ymin = cffDNA_GE_ml_plasma_ff_assay_min, ymax = cffDNA_GE_ml_plasma_ff_assay_max))+
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none")+
  xlim(0, 40)+
  labs(x = "Gestational Age (weeks)", y = "Fetal specific allele (GE/ml plasma)")

# Perform a t test to see if the fetal fraction and variant fraction assays amplify cfDNA equally
t.test(cfDNA_conc_gestations$cfDNA_GE_ml_plasma_variant_assay, cfDNA_conc_gestations$cfDNA_GE_ml_plasma_ff_assay,
       paired = TRUE)

#############################################################
# Sickle cell disease twin case (20915)
#############################################################

twin_sample <- sickle_cell_blinded %>%
  filter(r_number %in% c(20915, 13751))

colnames(sickle_cell_blinded)

sample_20915 <- sickle_cell_blinded %>%
  filter(r_number == 20915)

sample_13751 <- sickle_cell_blinded %>%
  filter(r_number == 13751)

ff_error_20915 <- sample_20915$Fetal_fraction_max_percent - sample_20915$Fetal_fraction_min_percent
ff_20915 <- sample_20915$Fetal_fraction_percent

ff_error_13751 <- sample_13751$Fetal_fraction_max_percent - sample_13751$Fetal_fraction_min_percent
ff_13751 <- sample_13751$Fetal_fraction_percent

# Plot the results of singleton vs twin pregnancy
ggplot(twin_sample, aes(x = r_number, y = Variant_fraction_percent))+
  geom_point()+
  geom_errorbar(aes(ymin = Variant_fraction_min_percent, 
                    ymax = Variant_fraction_max_percent), width =0.1)+
  theme_bw()+
  
  ## Twin sample classification areas
  
  # Region for classification as both heterozygous
  geom_rect(aes(xmin = 1.5, xmax = Inf, ymin = 50+(ff_error_20915/2), ymax = 50-(ff_error_20915/2)),
            alpha = 0.1,
            fill = "red")+
  # Region for classification as one heterozygous, one homozygous variant
  geom_rect(aes(xmin = 1.5, xmax = Inf, ymin = 50+(ff_20915/4) -(ff_error_20915/2), 
            ymax = 50+(ff_20915/4) +(ff_error_20915/2)),
            alpha = 0.1,
            fill = "red")+
  # Region for classification as both homozygous variant
  geom_rect(aes(xmin = 1.5, xmax = Inf, ymin = 50+(ff_20915/2) -(ff_error_20915/2), 
                ymax = 50+(ff_20915/2) +(ff_error_20915/2)),
            alpha = 0.1,
            fill = "red")+
  # Region for classification as one heterozygous, one homozygous reference
  geom_rect(aes(xmin = 1.5, xmax = Inf, ymin = 50-(ff_20915/4) -(ff_error_20915/2), 
                ymax = 50-(ff_20915/4) +(ff_error_20915/2)),
            alpha = 0.1,
            fill = "red")+
  # Region for classification as both homozygous reference
  geom_rect(aes(xmin = 1.5, xmax = Inf, ymin = 50-(ff_20915/2) -(ff_error_20915/2), 
                ymax = 50-(ff_20915/2) +(ff_error_20915/2)),
            alpha = 0.1,
            fill = "red")+
  
  ## Singleton sample classification areas
  ## Twin sample classification areas
  
  # Region for classification as both heterozygous
  geom_rect(aes(xmin = -Inf, xmax = 1.5, ymin = 50+(ff_error_13751/2), ymax = 50-(ff_error_13751/2)),
            alpha = 0.1,
            fill = "red")+
  # Region for classification as homozygous variant
  geom_rect(aes(xmin = -Inf, xmax = 1.5, ymin = 50+(ff_13751/2) -(ff_error_13751/2), 
                ymax = 50+(ff_13751/2) +(ff_error_13751/2)),
            alpha = 0.1,
            fill = "red")+
  # Region for classification as homozygous reference
  geom_rect(aes(xmin = -Inf, xmax = 1.5, ymin = 50-(ff_13751/2) -(ff_error_13751/2), 
                ymax = 50-(ff_13751/2) +(ff_error_13751/2)),
            alpha = 0.1,
            fill = "red")+
  labs(x = "", y = "Variant fraction (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(aes(xintercept=1.5), linetype="dashed", size = 1)+
  
  # Add on labels
  annotate(geom="text", x=1.65, y=58, label="Twin 1: HbSS
Twin 2: HbSS")+
  annotate(geom="text", x=1.65, y=54, label="Twin 1: HbAS
Twin 2: HbSS")+
  annotate(geom="text", x=1.7, y=50, label="Twin 1: HbAS or HbSS
Twin 2: HbAS or HbAA")+
  annotate(geom="text", x=1.65, y=46, label="Twin 1: HbAS
Twin 2: HbAA")+
  annotate(geom="text", x=1.65, y=42, label="Twin 1: HbAA
Twin 2: HbAA") +
  annotate(geom="text", x=0.55, y=46, label="Fetus HbAA")+
  annotate(geom="text", x=0.55, y=50, label="Fetus HbAS")+
  annotate(geom="text", x=0.55, y=54, label="Fetus HbSS")+
  annotate(geom="text", x=2.1, y=60, label="Sample 20915: dichorionic diamniotic twin 
pregnancy sample with fetal fraction of 16.1%", size = 4, fontface = "bold")+
  annotate(geom="text", x=1, y=60, label="Sample 13751: singleton pregnancy 
sample with fetal fraction of 9.0%", size = 4, fontface = "bold")

# Set the likelihood ratio
LR <- 250

sample_20915_twin_SPRT <- sample_20915 %>%
  mutate(q0 = 0.5) %>%
  # Prediction for one homozygous fetus, one heterozygous
  mutate(q1 = 0.5+(Fetal_fraction/4)) %>%
  # Prediction for both fetuses homozygous
  mutate(q2 = 0.5+(Fetal_fraction/2)) %>%
  mutate(Delta_1 = (1- q1)/(1-q0)) %>%
  mutate(Delta_2 = (1- q2)/(1-q1)) %>%
  mutate(Gamma_1 = ((q1 * (1-q0))/ (q0*(1-q1)))) %>%
  mutate(Gamma_2 = ((q2 * (1-q1))/ (q1*(1-q2)))) %>%
  # Likelihood ratio for one homozygous
  mutate(LR_1 = exp((((Over_represented_fraction*log(Gamma_1)) + log(Delta_1))*Molecules_variant_assay))) %>%
  # Likelihood ratio for both homozygous
  mutate(LR_2 = exp((((Over_represented_fraction*log(Gamma_2)) + log(Delta_2))*Molecules_variant_assay))) %>%
  # Calculate the threshold variant fractions for each classification
  mutate(threshold_AS_SS_lower = (((log(LR)/Molecules_variant_assay) - log(Delta_1)) / log(Gamma_1))*100) %>%
  mutate(threshold_SS_SS = (((log(LR)/Molecules_variant_assay) - log(Delta_2)) / log(Gamma_2))*100) %>%
  mutate(threshold_AS_AS_upper = (((log(1/LR)/Molecules_variant_assay) - log(Delta_1)) / log(Gamma_1))*100) %>%
  mutate(threshold_AS_SS_upper = (((log(1/LR)/Molecules_variant_assay) - log(Delta_2)) / log(Gamma_2))*100) %>%
  mutate(threshold_AA_AA = 50-(threshold_SS_SS-50)) %>%
  mutate(threshold_AS_AA_upper = 50- (threshold_AS_SS_lower-50)) %>%
  mutate(threshold_AS_AS_lower = 50-(threshold_AS_AS_upper-50)) %>%
  mutate(threshold_AS_AA_lower = 50-(threshold_AS_SS_upper-50))


ggplot(sample_20915_twin_SPRT, aes(x = Identifier, y = Over_represented_fraction_percent))+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = Over_represented_fraction_min_percent, 
                    ymax = Over_represented_fraction_max_percent), width =0)+
  # Region for classification as both homozygous variant
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = threshold_SS_SS, ymax = Inf),
            alpha = 0.1,
            fill = "grey")+
  # Region for classification as one fetus homozygous variant
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = threshold_AS_SS_lower, 
                ymax = threshold_AS_SS_upper),
            alpha = 0.1,
            fill = "grey")+
  # Region for classification heterozygous
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = threshold_AS_AS_lower, 
                ymax = threshold_AS_AS_upper),
            alpha = 0.1,
            fill = "grey")+
  # Region for classification as one fetus homozygous reference
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = threshold_AS_AA_lower, 
                ymax = threshold_AS_AA_upper),
            alpha = 0.1,
            fill = "grey")+
  # Region for classification as both homozygous reference
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = threshold_AA_AA),
            alpha = 0.1,
            fill = "grey")+
  ylim(42, 58)+
  labs(x = "", y = "")+
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18), 
        plot.title = element_text(size=20),
        legend.position = "bottom", legend.background = element_rect(fill="white"), 
        legend.title = element_blank(), legend.text = element_text(size= 15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  # Add on labels
  annotate(geom="text", x=0.65, y=58, label="HbSS:HbSS", size = 4)+
  annotate(geom="text", x=0.65, y=54, label="HbAS:HbSS", size = 4)+
  annotate(geom="text", x=0.65, y=50, label="HbAS:HbAS
or
HbAA:HbSS", size = 4)+
  annotate(geom="text", x=0.65, y=46, label="HbAS:HbAA", size = 4)+
  annotate(geom="text", x=0.65, y=42, label="HbAA:HbAA",size = 4)


# Plot for presentation
ggplot(sample_20915_twin_SPRT, aes(x = Identifier, y = Over_represented_fraction_percent))+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = Over_represented_fraction_min_percent, 
                    ymax = Over_represented_fraction_max_percent), width =0)+
  # Region for classification as both homozygous variant
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = threshold_SS_SS, ymax = Inf),
            alpha = 0.2,
            fill = "blue")+
  # Region for classification as one fetus homozygous variant
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = threshold_AS_SS_lower, 
                ymax = threshold_AS_SS_upper),
            alpha = 0.2,
            fill = "blue")+
  # Region for classification heterozygous
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = threshold_AS_AS_lower, 
                ymax = threshold_AS_AS_upper),
            alpha = 0.2,
            fill = "blue")+
  # Region for classification as one fetus homozygous reference
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = threshold_AS_AA_lower, 
                ymax = threshold_AS_AA_upper),
            alpha = 0.2,
            fill = "blue")+
  # Region for classification as both homozygous reference
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = threshold_AA_AA),
            alpha = 0.2,
            fill = "blue")+
  ylim(42, 58)+
  labs(x = "", y = "Variant fraction (%)")+
  theme_bw()+
  theme(axis.text.y=element_text(size=15), axis.text.x=element_blank(),
        axis.title = element_text(size=18), 
        plot.title = element_text(size=20),
        legend.position = "bottom", legend.background = element_rect(fill="white"), 
        legend.title = element_blank(), legend.text = element_text(size= 15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


## Empty plot
ggplot(sample_20915_twin_SPRT %>%
         filter(Fetal_fraction_percent > 100), aes(x = Identifier, y = Over_represented_fraction_percent))+
  geom_point(size = 3)+
  ylim(42, 58)+
  labs(x = "", y = "Variant fraction (%)")+
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18), 
        plot.title = element_text(size=20),
        legend.position = "bottom", legend.background = element_rect(fill="white"), 
        legend.title = element_blank(), legend.text = element_text(size= 15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#############################################################
# Plotting individual sample graphs (function)
#############################################################

## This function plots a graph of the results for a single case,
# with maternal gDNA individual replicates and cfDNA merged results. 

plot_sample_graph <- function(R_number, control_number){
  
  sample_variant <- bespoke_cohort_analysed %>%
    filter(r_number == R_number) %>%
    select(r_number, Molecules_variant, variant_assay,
           Molecules_reference, AcceptedDroplets_Variant_assay) %>%
    pivot_longer(
      cols = !c(r_number, variant_assay, AcceptedDroplets_Variant_assay),
      names_to = "Target",
      values_to = "count"
    ) %>%
    rename(assay = variant_assay, AcceptedDroplets = AcceptedDroplets_Variant_assay)
  
  sample_ff <- bespoke_cohort_analysed %>%
    filter(r_number == R_number) %>%
    select(r_number, ff_assay, AcceptedDroplets_FetalFrac, Molecules_maternal,
           Molecules_paternal) %>%
    pivot_longer(
      cols = !c(r_number, ff_assay, AcceptedDroplets_FetalFrac),
      names_to = "Target",
      values_to = "count"
    ) %>%
    rename(assay = ff_assay, AcceptedDroplets = AcceptedDroplets_FetalFrac)
  
  control_variant <- ddpcr_control_tbl_var %>%
    filter(Sample == control_number) %>%
    select(Worksheet_well, Molecules_variant, variant_assay, AcceptedDroplets_Variant_assay,
           Molecules_reference) %>%
    pivot_longer(
      cols = !c(Worksheet_well, variant_assay, AcceptedDroplets_Variant_assay),
      names_to = "Target",
      values_to = "count"
    ) %>%
    rename(r_number = Worksheet_well, assay = variant_assay,
           AcceptedDroplets = AcceptedDroplets_Variant_assay)
  
  control_ff <- ddpcr_control_tbl_ff %>%
    filter(Sample == control_number) %>%
    select(Worksheet_well, ff_assay, AcceptedDroplets_FetalFrac, Molecules_maternal,
           Molecules_paternal) %>%
    pivot_longer(
      cols = !c(Worksheet_well, ff_assay, AcceptedDroplets_FetalFrac),
      names_to = "Target",
      values_to = "count") %>%
    rename(r_number = Worksheet_well, assay = ff_assay, 
           AcceptedDroplets = AcceptedDroplets_FetalFrac)
  
  # Bind the tables in the order they will appear in the plot
  sample_control <- rbind(control_ff, control_variant, sample_ff, sample_variant) %>%
    mutate(id = paste(r_number, assay)) %>%
    # This trick fixes the table so that the plot is ordered correctly
    mutate(id = factor(id, levels = unique(id))) %>%
    mutate(Cpd = count/AcceptedDroplets) %>%
    mutate(count_max = Poisson_max(Cpd, AcceptedDroplets)) %>%
    mutate(count_min = Poisson_min(Cpd, AcceptedDroplets)) %>%
    mutate(sample_type = ifelse(r_number %in% ddpcr_data_tbl$Sample, "cfDNA merged", "mat gDNA"))
  
  # To have useful x labels for the plot, I had to do this hack.
  x_axis_labels <- c(rep("mat gDNA 
single", (nrow(sample_control)-4)/2), rep("cfDNA 
merged", 2))
  
  sample_plot <- ggplot(sample_control, aes(x = id, y = count, fill = Target))+
    geom_col(position = position_dodge(width = 0.9), colour="black")+
    geom_errorbar(aes(ymin = count_min, ymax = count_max, width = 0.3), position = position_dodge(width = 0.9))+
    theme_bw()+
    theme(axis.text=element_text(size=15), axis.title = element_text(size=18),
          axis.text.x = element_text(angle = 45, vjust = 0.5))+
    labs(x = "", y = "Molecules detected", title = paste("Sample:", sample_variant$r_number, ";" ,"Variant assay:", 
                                                         sample_variant$assay, ";",
                                                         "Fetal fraction assay:", sample_ff$assay))+
    theme(axis.text.x = element_blank())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    
    # Label the x axis with the sample type to make things simpler
    scale_x_discrete(labels=x_axis_labels) +
    theme(axis.text.x=element_text(vjust=.5, size = 9))
  
  return(sample_plot)
}

# VHL
plot_sample_graph(14522, "20RG-148G0075")
# TCOF
plot_sample_graph(14116, "20RG-148G0076")
# FGFR3
plot_sample_graph(13519, "20RG-307G0060")
# NF1
plot_sample_graph(14491, "20RG-148G0081")
# ALPL
plot_sample_graph(19261, "20RG-307G0064")
# FGFR2
plot_sample_graph(19711, "20RG-148G0080")
# ADAR
plot_sample_graph(19102, "20RG-336G0066")
# ABCD1
plot_sample_graph(19397, "18G10773")
# MAGED2
plot_sample_graph(20980, "21RG-027G0010")
# IDS
plot_sample_graph(19611, "21RG-027G0070")
# OTC
plot_sample_graph(17667, "21RG-027G0004")

#############################################################
# Plotting individual sample graphs for paper
#############################################################

R_number <- 14491
control_number <- "20RG-148G0081"

sample_variant <- bespoke_cohort_analysed %>%
  filter(r_number == R_number) %>%
  # RENAME
  select(r_number, Molecules_variant, variant_assay,
         Molecules_reference, AcceptedDroplets_Variant_assay) %>%
  pivot_longer(
    cols = !c(r_number, variant_assay, AcceptedDroplets_Variant_assay),
    names_to = "Target",
    values_to = "count"
  ) %>%
  rename(assay = variant_assay, AcceptedDroplets = AcceptedDroplets_Variant_assay)

sample_ff <- bespoke_cohort_analysed %>%
  filter(r_number == R_number) %>%
  # RENAME
  select(r_number, ff_assay, AcceptedDroplets_FetalFrac, Molecules_maternal,
         Molecules_paternal) %>%
  pivot_longer(
    cols = !c(r_number, ff_assay, AcceptedDroplets_FetalFrac),
    names_to = "Target",
    values_to = "count"
  ) %>%
  rename(assay = ff_assay, AcceptedDroplets = AcceptedDroplets_FetalFrac)

control_variant <- ddpcr_control_tbl_var %>%
  filter(Sample == control_number) %>%
  # RENAME
  select(Worksheet_well, Molecules_variant, variant_assay, AcceptedDroplets_Variant_assay,
         Molecules_reference) %>%
  pivot_longer(
    cols = !c(Worksheet_well, variant_assay, AcceptedDroplets_Variant_assay),
    names_to = "Target",
    values_to = "count"
  ) %>%
  rename(r_number = Worksheet_well, assay = variant_assay,
         AcceptedDroplets = AcceptedDroplets_Variant_assay)

control_ff<- ddpcr_control_tbl_ff %>%
  filter(Sample == control_number) %>%
  select(Worksheet_well, ff_assay, AcceptedDroplets_FetalFrac, Molecules_maternal,
         Molecules_paternal) %>%
  pivot_longer(
    cols = !c(Worksheet_well, ff_assay, AcceptedDroplets_FetalFrac),
    names_to = "Target",
    values_to = "count") %>%
  rename(r_number = Worksheet_well, assay = ff_assay, 
         AcceptedDroplets = AcceptedDroplets_FetalFrac)

# Bind the tables in the order they will appear in the plot
sample_control <- rbind(control_ff, sample_ff, control_variant, sample_variant) %>%
  mutate(id = paste(r_number, assay)) %>%
  # This trick fixes the table so that the plot is ordered correctly
  mutate(id = factor(id, levels = unique(id))) %>%
  mutate(Cpd = count/AcceptedDroplets) %>%
  mutate(count_max = Poisson_max(Cpd, AcceptedDroplets)) %>%
  mutate(count_min = Poisson_min(Cpd, AcceptedDroplets)) %>%
  mutate(sample_type = ifelse(r_number %in% ddpcr_data_tbl$Sample, "cfDNA merged", "mat gDNA")) %>%
  mutate(new_target = case_when(
    Target == "Molecules_maternal" ~"Shared allele",
    Target == "Molecules_paternal" ~"Fetal allele",
    Target == "Molecules_variant" ~"Variant allele",
    Target == "Molecules_reference" ~"Reference allele"))

# Add factor levels to order the graph correctly.
sample_control$new_target <- factor(sample_control$new_target, levels = c("Shared allele", "Fetal allele",
                                                                          "Reference allele", "Variant allele"))

# To have useful x labels for the plot, I had to do this hack.
x_axis_labels <- c("mat gDNA 
rep 1", "mat gDNA 
rep 2", "mat gDNA 
rep 3", "cfDNA 
merged", "mat gDNA 
rep 1", "mat gDNA 
rep 2", "mat gDNA 
rep 3", "cfDNA 
merged")

sample_plot <- ggplot(sample_control, aes(x = id, y = count, fill = new_target))+
  geom_col(position = position_dodge(width = 0.9), colour="black", alpha = 0.6)+
  
  scale_fill_manual(values = c("#CCCCCC", "#00FF00", "#3366FF", "#FF0000"))+
  geom_errorbar(aes(ymin = count_min, ymax = count_max, width = 0.3), position = position_dodge(width = 0.9))+
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18),
        legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "", y = "Molecules", title = paste("Variant assay:", 
                                              sample_variant$assay, ";",
                                              "Fetal fraction assay:", sample_ff$assay))+
  # Label the x axis with the sample type to make things simpler
  scale_x_discrete(labels=x_axis_labels) +
  theme(axis.text.x=element_text(vjust=.5, size = 12))


# FGFR3
plot_sample_graph(13519, "20RG-307G0060")
# NF1
plot_sample_graph(14491, "20RG-148G0081")
# ALPL
plot_sample_graph(19261, "20RG-307G0064")
# FGFR2
plot_sample_graph(19711, "20RG-148G0080")




#############################################################
# Plotting graphs of sample cohorts
#############################################################

# Plot the results of fetal fraction versus variant fraction, including inconclusive results.
ggplot(sickle_cell_unblinded, aes(x = Fetal_fraction_percent, y = Variant_fraction_percent))+
  geom_point(aes(shape = overall_prediction, alpha = overall_prediction), size = 4)+
  # Change the shape of the shapes displayed
  scale_shape_manual(values = c(15, 17, 18, 19, 19))+
  scale_colour_manual(values = c("#FFFFFF"))+
  geom_errorbarh(aes(xmin = Fetal_fraction_min_percent, xmax = Fetal_fraction_max_percent,
                     alpha = overall_prediction)) +
  geom_errorbar(aes(ymin = Variant_fraction_min_percent, ymax = Variant_fraction_max_percent,
                    alpha = overall_prediction))+
  scale_alpha_manual(values = c(0.8,0.8,0.8,0.2, 0.8))+
  fifty_percent_line +
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18), 
        plot.title = element_text(size=20),
        legend.position = "bottom", legend.background = element_rect(fill="white"), 
        legend.title = element_blank(), legend.text = element_text(size= 15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Fetal fraction (%)", y = "Variant fraction (%)")+
  ylim(42.5, 57.5)+
  xlim(0, 22)+
  geom_text(data=subset(sickle_cell_blinded, Identifier %in% c("HBB-57", "HBB-62", "HBB-64")),
            aes(label=Identifier, vjust = 1.5, hjust = 1.2))

# Plot the results without calls or legend
ggplot(sickle_cell_unblinded, aes(x = Fetal_fraction_percent, y = Variant_fraction_percent))+
  geom_point(size = 8, alpha = 0.5)+
  # Change the shape of the shapes displayed
  geom_errorbarh(aes(xmin = Fetal_fraction_min_percent, xmax = Fetal_fraction_max_percent)) +
  geom_errorbar(aes(ymin = Variant_fraction_min_percent, ymax = Variant_fraction_max_percent))+
  fifty_percent_line +
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18), 
        plot.title = element_text(size=20),
        legend.position = "none", legend.background = element_rect(fill="white"), 
        legend.title = element_blank(), legend.text = element_text(size= 15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Fetal fraction (%)", y = "Variant fraction (%)")+
  ylim(42.5, 57.5)+
  xlim(0, 22)


# Empty plot without results
ggplot(sickle_cell_unblinded %>%
         filter(Fetal_fraction_percent > 100), aes(x = Fetal_fraction_percent, y = Variant_fraction_percent))+
  geom_point(size = 8, alpha = 0.5)+
  # Change the shape of the shapes displayed
  geom_errorbarh(aes(xmin = Fetal_fraction_min_percent, xmax = Fetal_fraction_max_percent)) +
  geom_errorbar(aes(ymin = Variant_fraction_min_percent, ymax = Variant_fraction_max_percent))+
  fifty_percent_line +
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18), 
        plot.title = element_text(size=20),
        legend.position = "none", legend.background = element_rect(fill="white"), 
        legend.title = element_blank(), legend.text = element_text(size= 15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Fetal fraction (%)", y = "Variant fraction (%)")+
  ylim(42.5, 57.5)+
  xlim(0, 22)


# Plot the results without calls or legend or confidence intervals
ggplot(sickle_cell_unblinded, aes(x = Fetal_fraction_percent, y = Variant_fraction_percent))+
  geom_point(size = 8, alpha = 0.5)+
  # Change the shape of the shapes displayed
  fifty_percent_line +
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18), 
        plot.title = element_text(size=20),
        legend.position = "none", legend.background = element_rect(fill="white"), 
        legend.title = element_blank(), legend.text = element_text(size= 15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Fetal fraction (%)", y = "Variant fraction (%)")+
  ylim(42.5, 57.5)+
  xlim(0, 22)


# Plot the SCD results in colour without error bars
ggplot(sickle_cell_unblinded, aes(x = Fetal_fraction_percent, y = Variant_fraction_percent))+
# Add colours for the classifications. Also put them in the right order for the legend: homozygous
# variant at the top
# 0000FF is dark blue
# 999999 is grey
# 000000 is black
# 99CCFF is light blue
# Factor order is alphabetical: het, hom normal, hom variant, no call
  scale_fill_manual(values=c("#000000", "#99CCFF", "#999999", "#0000FF","#FF0000"),
                  breaks=c("homozygous variant", "heterozygous", "no call", "homozygous reference", "twin pregnancy"))+
  geom_point(size = 8, aes(fill = overall_prediction), colour="black", pch=21, alpha = 0.8)+
  scale_colour_manual(values = c("#FFFFFF"))+
  fifty_percent_line +
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18), 
        plot.title = element_text(size=20),
        legend.background = element_rect(fill="white"),
        legend.position = c(0.85, 0.85), 
        legend.title = element_blank(), legend.text = element_text(size= 15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Fetal fraction (%)", y = "Variant fraction (%)")+
  ylim(42.5, 57.5)+
  xlim(0, 22)


cfDNA_for_graph <- sickle_cell_unblinded %>%
  filter(overall_prediction != "no call") %>%
  arrange(overall_prediction, Identifier) %>%
  mutate(Identifier=factor(Identifier, levels=Identifier)) %>%
  select(Identifier, Variant_fraction_percent, Variant_fraction_max_percent, Variant_fraction_min_percent)

all_DNA_for_graph <- rbind(het_gDNA_arranged, cfDNA_for_graph) %>%
  mutate(Identifier=factor(Identifier, levels=Identifier))

ggplot(sickle_cell_unblinded %>%
         filter(overall_prediction != "no call") %>%
         arrange(overall_prediction, Identifier) %>%
         mutate(Identifier=factor(Identifier, levels=Identifier)), aes(x = Identifier, y = Variant_fraction_percent))+
  geom_point(stat="identity")+
  geom_errorbar(aes(ymin = Variant_fraction_min_percent, ymax = Variant_fraction_max_percent))+
  fifty_percent_line +
  ddPCR_plot_theme+  
  labs(x = "", y = "Variant fraction (%)")

variant_fraction_plot <- ggplot(all_DNA_for_graph, aes(x = Identifier, y = Variant_fraction_percent))+
  geom_point(stat="identity")+
  geom_errorbar(aes(ymin = Variant_fraction_min_percent, ymax = Variant_fraction_max_percent))+
  geom_hline(aes(yintercept=50), linetype="dashed", size = 1)+
  geom_vline(aes(xintercept=23.5), linetype="dashed", colour = "grey", size = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 15), axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Samples", y = "Variant fraction (%)")+
  ylim(42, 58)

variant_fraction_plot + annotate(geom="text", x=10, y=57, label="HbAS gDNA and paternal cfDNA") +
  annotate(geom="text", x=35, y=57, label="cfDNA from pregnant HbAS carriers")

# Plot the results of all samples with HbAS fetuses by fetal and variant fraction.
ggplot(sickle_cell_unblinded %>%
         filter(invasive_result == "HbAS"), aes(x = Fetal_fraction_percent, y = Variant_fraction_percent,
                                                colour = overall_prediction))+
  geom_point()+
  geom_errorbar(aes(ymin = Variant_fraction_min_percent, ymax = Variant_fraction_max_percent))+
  geom_errorbarh(aes(xmin = Fetal_fraction_min_percent, xmax = Fetal_fraction_max_percent))+
  fifty_percent_line +
  ddPCR_plot_theme+  
  labs(x = "Total DNA molecules", y = "Variant fraction (%)")+
  ylim(47, 53)

# Plot the results of all samples with HbAS fetuses by variant fraction and total number of molecules
ggplot(sickle_cell_unblinded %>%
         filter(invasive_result == "HbAS"), aes(x = Molecules_variant_assay, y = Variant_fraction_percent,
                                                colour = overall_prediction))+
  geom_point()+
  geom_errorbar(aes(ymin = Variant_fraction_min_percent, ymax = Variant_fraction_max_percent))+
  fifty_percent_line +
  ddPCR_plot_theme+  
  labs(x = "Total DNA molecules", y = "Variant fraction (%)")+
  ylim(47, 53)

ggplot(bespoke_cohort_analysed, aes(x = Fetal_fraction_percent, y = Variant_fraction_percent))+
  geom_point(size = 4)+
  geom_errorbarh(aes(xmin = Fetal_fraction_min_percent, xmax = Fetal_fraction_max_percent)) +
  geom_errorbar(aes(ymin = Variant_fraction_min_percent, ymax = Variant_fraction_max_percent))+
  fifty_percent_line +
  ddPCR_plot_theme+
  labs(x = "Fetal fraction (%)", y = "Variant fraction (%)")

# Plot graph of variant imbalance vs fetal fraction.
autosomal_cohort <- ddpcr_data_tbl %>%
  filter(Inheritance == "autosomal")

autosomal_cohort_analysed <- calc_SPRT(calc_conf_intervals(calc_molecules(autosomal_cohort)), 8)

autosomal_cohort_imbalance <- autosomal_cohort_analysed %>%
  filter(SPRT_prediction != "heterozygous")

ggplot(autosomal_cohort_imbalance %>%
         filter(Sample != 20915), aes(x = Fetal_fraction_percent, y = Over_represented_fraction_percent))+
  geom_point(aes(shape = variant_assay))+
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18), 
        plot.title = element_text(size=20),
        legend.position = "bottom", legend.background = element_rect(fill="white"), 
        legend.title = element_text(size= 18), legend.text = element_text(size= 15))+
  geom_errorbarh(aes(xmin = Fetal_fraction_min_percent, xmax = Fetal_fraction_max_percent)) +
  geom_errorbar(aes(ymin = Over_represented_fraction_min_percent, ymax = Over_represented_fraction_max_percent))+
  ylim(48, 60)+
  xlim(0, 20)+
  geom_abline(intercept = 50, linetype="dashed", colour = "grey", size = 1, slope = 0.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Fetal fraction %", y = "Over represented fraction %",
       title = "Over-represented fraction versus fetal fraction
       in ddPCR cases with homozygous fetal genotypes")+
  geom_smooth(se = FALSE, method = "lm")

#############################################################
# Export csvs with time stamps
#############################################################

current_time <- Sys.time()

write.csv(sickle_cell_unblinded, 
          file = paste0("analysis_outputs/sickle_cell_cohort_analysed_unblinded", 
                        format(current_time, "%Y%m%d_%H%M%S"), ".csv"),
          row.names = FALSE)

write.csv(bespoke_cohort_unblinded, 
          file = paste0("analysis_outputs/bespoke_cohort_unblinded", 
                        format(current_time, "%Y%m%d_%H%M%S"), ".csv"),
          row.names = FALSE)


