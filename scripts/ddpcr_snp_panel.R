#############################################################
## ddPCR SNP panel analysis
## May 2021
## joseph.shaw3@nhs.net
## These scripts deal with the ddPCR SNP panel used for 
## determining the fetal fraction of cfDNA.
## This panel was originally described by Camunas-Soler et
## al (2018) (PMID: 29097507)
#############################################################

#############################################################
# Load packages and panel GnomAD data
#############################################################

library(tidyverse)
library(janitor)

## Set working directory
setwd("W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/ddPCR_R_Analysis/ddpcr_nipd")

# Load resources
controls <- readr::read_csv("resources/controls.csv")
ddpcr_target_panel <- readr::read_csv("resources/ddpcr_target_panel.csv") 

snp_panel <- read.csv("W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/ddPCR SNP Panel/Camunas_Soler_panel_47_GnomAD_frequencies.csv")

#############################################################
# Reshape and calculate probabilities for each SNP
#############################################################

snp_calc <- snp_panel %>% 
  select(dbSNP, African, Ashkenazi.Jewish, East.Asian, European..Finnish.,
         European..Non.Finnish., Latino, Other, GOSH_ID_ds) %>%
  pivot_longer(
    cols = !c(dbSNP, GOSH_ID_ds),
    names_to = "population",
    values_to = "MAF1") %>%
  # Remove SNPs not on GOSH panel
  filter(!is.na(GOSH_ID_ds)) %>%
  mutate(MAF2 = 1-MAF1,
         # A is major allele, B is minor allele
         frequency_A = pmax(MAF1, MAF2),
         frequency_B = pmin(MAF1, MAF2),
         # Calculate the probability of each genotype according to 
         # Hardy-Weinberg equilibrium
         p_AA = frequency_A^2,
         p_BB = frequency_B^2,
         p_AB = 2*(frequency_A*frequency_B),
         # Check it all adds up to 1.
         check = p_AA + p_BB + p_AB,
         # Calculate the probability of each parental genotype combination. 
         # Mother first, then father.
         p_AA_AA = p_AA*p_AA,
         p_AA_AB = p_AA*p_AB,
         p_AA_BB = p_AA*p_BB,
         p_AB_AA = p_AB*p_AA,
         p_AB_AB = p_AB*p_AB,
         p_AB_BB = p_AB*p_BB,
         p_BB_AA = p_BB*p_AA,
         p_BB_AB = p_BB*p_AB,
         p_BB_BB = p_BB*p_BB,
         # Check it all adds up to 1.
         check2 = p_AA_AA + p_AA_AB + p_AA_BB + p_AB_AA + 
           p_AB_AB + p_AB_BB + p_BB_AA + p_BB_AB + p_BB_BB,
         
         # Calculate the probability of an uninformative result with each 
         # parental genotype combination.
         # Probability is 1 for a type 2 SNP.
         p_AA_AA_uninf = p_AA_AA * 1,
         
         # If dad is heterozygous there's a 50% chance the fetus won't inherit 
         # the paternal-specific allele.
         # Note this step differs depending on whether you are dealing 
         # with prenatal testing (fetus inherits one paternal allele)
         # or transplant monitoring (both donor alleles will be present 
         # in the sample).
         p_AA_AB_uninf = p_AA_AB * 0.5,
         # Probability is 0 for a type 1 SNP,
         p_AA_BB_uninf = 0,
         # Any locus where the mother is heterozygous will not be informative
         p_AB_AA_uninf = p_AB_AA * 1,
         p_AB_AB_uninf = p_AB_AB * 1,
         p_AB_BB_uninf = p_AB_BB * 1,
         p_BB_AA_uninf = 0,
         p_BB_AB_uninf = p_BB_AB * 0.5,
         p_BB_BB_uninf = p_BB_BB * 1,
         
         # Calculate the total chance of the SNP being uninformative
         SNP_uninf_total = p_AA_AA_uninf + p_AA_AB_uninf + 
           p_AA_BB_uninf + p_AB_AA_uninf +
           p_AB_AB_uninf + p_AB_BB_uninf + 
           p_BB_AA_uninf + p_BB_AB_uninf + p_BB_BB_uninf,
         # Calculate the total chance of the SNP being a type 1 SNP.
         SNP_type1 = p_AA_BB + p_BB_AA,
         SNP_not_type1 = 1-SNP_type1)

#############################################################
# Uninformativeness of 24 SNP panel
#############################################################

# Check the probability of the 24 SNP panel being uninformative when 
# using pre-amplified cfDNA (type 1 and type 3 SNPs) 
snp_24_type1and3 <- snp_calc %>%
  # Select only the first 24 SNPs
  filter(GOSH_ID_ds < 25) %>%
  select(dbSNP, population, SNP_uninf_total) %>%
  # Reshape by pivot_wider
  pivot_wider(id_cols = population,
              names_from = dbSNP,
              values_from = SNP_uninf_total) %>%
  mutate(chance_uninformative = rs565522 * rs2737654 * rs2576241 * 
           rs1160680 * rs1914748 * rs12694624 * rs1399629 *rs2276702 * 
           rs6781236 * rs7653090 * rs357485 * rs9290003 * rs17017347 * 
           rs10027026 * rs4975819 * rs6899022 * rs1185246 * rs6877199 * 
           rs13218440 * rs6924733 * rs2535290 * rs172275 * rs4644087 * 
           rs12690832,
         uninformative_1_in  = 1/chance_uninformative) %>%
  select(population, chance_uninformative, uninformative_1_in)

# Check the probability of the 24 SNP panel being uninformative when 
# using parental gDNA (type 1 SNPs only) 
snp_24_type1_only <- snp_calc %>%
  # Select only the first 24 SNPs
  filter(GOSH_ID_ds < 25) %>%
  select(dbSNP, population, SNP_not_type1) %>%
  # Reshape by pivot_wider
  pivot_wider(id_cols = population,
              names_from = dbSNP,
              values_from = SNP_not_type1) %>%
  mutate(chance_uninformative = rs565522 * rs2737654 * rs2576241 * 
           rs1160680 * rs1914748 * rs12694624 * rs1399629 *rs2276702 * 
           rs6781236 * rs7653090 * rs357485 * rs9290003 * rs17017347 * 
           rs10027026 * rs4975819 * rs6899022 * rs1185246 * rs6877199 * 
           rs13218440 * rs6924733 * rs2535290 * rs172275 * rs4644087 * 
           rs12690832,
         uninformative_1_in  = 1/chance_uninformative) %>%
  select(population, chance_uninformative, uninformative_1_in)

#############################################################
# Uninformativeness of 40 SNP panel
#############################################################

# Check the probability of the 40 SNP panel being uninformative when 
# using pre-amplified cfDNA (type 1 and type 3 SNPs) 
snp_40_type1and3 <- snp_calc %>%
  select(dbSNP, population, SNP_uninf_total) %>%
  # Reshape by pivot_wider
  pivot_wider(id_cols = population,
              names_from = dbSNP,
              values_from = SNP_uninf_total) %>%
  mutate(chance_uninformative = rs565522 * rs2737654 * rs2576241 * 
           rs1160680 * rs1914748 * 
           rs12694624 * rs1399629 *rs2276702 * rs6781236 * rs7653090 * 
           rs357485 * rs9290003 * rs17017347 * rs10027026 * rs4975819 * 
           rs6899022 * rs1185246 * rs6877199 * rs13218440 * rs6924733 * 
           rs2535290 * rs172275 * rs4644087 * rs12690832 * rs7827391 * 
           rs2319150 * rs1410059 * rs10821808 * rs2370764 * rs1498553 * 
           rs2256111 * rs7325978 * rs3742560 * rs12148532 * rs249290 * 
           rs1544724 * rs3760269 * rs7233004 * rs4801945 * rs271981,
         uninformative_1_in  = 1/chance_uninformative) %>%
  select(population, chance_uninformative, uninformative_1_in)

# Check the probability of the 40 SNP panel being uninformative when using using parental gDNA 
# (type 1 SNPs only)
snp_40_type1_only <- snp_calc %>%
  select(dbSNP, population, SNP_not_type1) %>%
  # Reshape by pivot_wider
  pivot_wider(id_cols = population,
              names_from = dbSNP,
              values_from = SNP_not_type1) %>%
  mutate(chance_uninformative = rs565522 * rs2737654 * rs2576241 * 
           rs1160680 * rs1914748 * 
           rs12694624 * rs1399629 *rs2276702 * rs6781236 * rs7653090 * 
           rs357485 * rs9290003 * rs17017347 * rs10027026 * rs4975819 * 
           rs6899022 * rs1185246 * rs6877199 * rs13218440 * rs6924733 * 
           rs2535290 * rs172275 * rs4644087 * rs12690832 * rs7827391 * 
           rs2319150 * rs1410059 * rs10821808 * rs2370764 * rs1498553 * 
           rs2256111 * rs7325978 * rs3742560 * rs12148532 * rs249290 * 
           rs1544724 * rs3760269 * rs7233004 * rs4801945 * rs271981,
         uninformative_1_in  = 1/chance_uninformative) %>%
  select(population, chance_uninformative, uninformative_1_in)

#############################################################
# GnomAD SNP frequencies plot
#############################################################

# Plot of GnomAD SNP frequencies
ggplot(snp_calc, aes(x = GOSH_ID_ds, y = frequency_B))+
  geom_point(size = 5, alpha = 0.5, aes(colour =population))+
  ylim(0, 0.5)+
  labs(x = "SNPs", y = "GnomAD minor allele frequency", 
       title = "GnomAD minor allele frequencies of Camunas-Soler et al (2018) ddPCR SNP panel")+
  theme_bw()

#############################################################
# Collating SNP genotypes from cfDNA data
#############################################################

snp_panel <- read_csv("W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/ddPCR SNP Panel/Camunas_Soler_panel_47_GnomAD_frequencies.csv")

dataPath <- "data/ddPCR_SNP_genotyping/"

ddpcr_files <- list.files(dataPath)

snp_panel_24 <- c(unique(snp_panel %>%
                           # Selct only the first 24 SNPs
                           filter(GOSH_ID_ds < 25) %>%
                           select(dbSNP)))
#Empty data frame
SNP_data <- data.frame()

# Read and collate each worksheet csv
for (dataFile in ddpcr_files){
  tmp_dat <- read_csv(paste0(dataPath,dataFile), col_names = TRUE,
                      show_col_types = FALSE) %>%
    janitor::clean_names()
  SNP_data <-rbind(SNP_data, tmp_dat)
  rm(tmp_dat)
}

#############################################################
# Function for calculating fetal fraction for cfDNA samples
#############################################################

digital_snp <- function(cf_sample){
  sample_snp_table <- SNP_data %>%
    filter(sample == cf_sample) %>%
    select(sample, target, target_type, concentration, fractional_abundance) %>%
    
    left_join(ddpcr_target_panel %>%
              select(target, assay), by = "target") %>%
    # Convert the concentration column to a numeric and convert 
    # "no call" to zero
    mutate(copies_per_ul = as.integer(ifelse(concentration == "No Call", 
                                             0, concentration))) %>%
    # Determine fetal fraction (assuming it is between 1 and 20%)
    dplyr::rename(fraction_a = fractional_abundance) %>%
    mutate(fraction_b = 100 - fraction_a,
           minor_fraction = pmin(fraction_b, fraction_a),
           fetal_fraction = case_when(minor_fraction > 1 & 
                                        minor_fraction < 20 
                                      ~minor_fraction)) %>%
  # Filter to get one row per well
  filter(target_type == "Ch1Unknown") %>%
  select(c(assay, fetal_fraction))
  
  mean_ff <- mean(sample_snp_table$fetal_fraction)
  
  return(sample_snp_table)
}

#########################
# Select pre-amplification results
#########################

# This step can be used to select pre-amplified data for comparison
# with ddPCR on the original cfDNA sample. 

pre_amplification_fractions <- SNP_data %>%
    dplyr::rename(r_number = sample,
                  fraction_a = fractional_abundance,
                  fractionmax_a = poisson_fractional_abundance_max,
                  fractionmin_a = poisson_fractional_abundance_min) %>%
    
    mutate(fraction_b = 100 - fraction_a,
           fractionmax_b = 100 - fractionmin_a,
           fractionmin_b = 100 - fractionmax_a,
           
           minor_fraction = pmin(fraction_a, fraction_b),
           minor_fractionmax = pmin(fractionmax_a, fractionmax_b),
           minor_fractionmin = pmin(fractionmin_a, fractionmin_b)) %>%
    
    left_join(ddpcr_target_panel %>%
                select(target, assay),
              by = "target") %>%
  
  select(c(r_number, target, assay, minor_fraction, 
             minor_fractionmax, minor_fractionmin)) %>%
  filter(r_number != "NTC") %>%
  # Rename the columns to allow graph plotting
  dplyr::rename(preamp_minor_fraction = minor_fraction,
                preamp_minor_fractionmax = minor_fractionmax,
                preamp_minor_fractionmin = minor_fractionmin)

#############################################################