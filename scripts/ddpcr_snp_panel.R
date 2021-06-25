#############################################################
## ddPCR SNP panel analysis
## May 2021
## Joseph.Shaw@gosh.nhs.uk 
## These scripts deal with the ddPCR SNP panel used for 
## determining the fetal fraction of cfDNA.
## This panel was originally described by Camunas-Soler et
## al (2018) (PMID: 29097507)
#############################################################

#############################################################
# Load packages and panel GnomAD data
#############################################################

library(tidyverse)

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
  mutate(MAF2 = 1-MAF1) %>%
  # A is major allele, B is minor allele
  mutate(frequency_A = pmax(MAF1, MAF2)) %>%
  mutate(frequency_B = pmin(MAF1, MAF2)) %>%
  # Calculate the probability of each genotype according to Hardy-Weinberg equilibrium
  mutate(p_AA = frequency_A^2) %>%
  mutate(p_BB = frequency_B^2) %>%
  mutate(p_AB = 2*(frequency_A*frequency_B)) %>%
  # Check it all adds up to 1.
  mutate(check = p_AA + p_BB + p_AB) %>%
  # Calculate the probability of each parental genotype combination. Mother first, then father.
  mutate(p_AA_AA = p_AA*p_AA) %>%
  mutate(p_AA_AB = p_AA*p_AB) %>%
  mutate(p_AA_BB = p_AA*p_BB) %>%
  mutate(p_AB_AA = p_AB*p_AA) %>%
  mutate(p_AB_AB = p_AB*p_AB) %>%
  mutate(p_AB_BB = p_AB*p_BB) %>%
  mutate(p_BB_AA = p_BB*p_AA) %>%
  mutate(p_BB_AB = p_BB*p_AB) %>%
  mutate(p_BB_BB = p_BB*p_BB) %>%
  # Check it all adds up to 1.
  mutate(check2 = p_AA_AA + p_AA_AB + p_AA_BB + p_AB_AA + p_AB_AB + p_AB_BB + p_BB_AA + p_BB_AB + p_BB_BB) %>%
  
  # Calculate the probability of an uninformative result with each parental genotype combination.
  
  # Probability is 1 for a type 2 SNP.
  mutate(p_AA_AA_uninf = p_AA_AA * 1) %>%
  # If dad is heterozygous there's a 50% chance the fetus won't inherit the paternal-specific allele.
  # Note this step differs depending on whether you are dealing with prenatal testing (fetus inherits one paternal allele)
  # or transplant monitoring (both donor alleles will be present in the sample).
  mutate(p_AA_AB_uninf = p_AA_AB * 0.5) %>%
  # Probability is 0 for a type 1 SNP
  mutate(p_AA_BB_uninf = 0) %>%
  # Any locus where the mother is heterozygous will not be informative
  mutate(p_AB_AA_uninf = p_AB_AA * 1) %>%
  mutate(p_AB_AB_uninf = p_AB_AB * 1) %>%
  mutate(p_AB_BB_uninf = p_AB_BB * 1) %>%
  mutate(p_BB_AA_uninf = 0) %>%
  mutate(p_BB_AB_uninf = p_BB_AB * 0.5) %>%
  mutate(p_BB_BB_uninf = p_BB_BB * 1) %>%
  
  # Calculate the total chance of the SNP being uninformative
  mutate(SNP_uninf_total = p_AA_AA_uninf + p_AA_AB_uninf + p_AA_BB_uninf + p_AB_AA_uninf +
           p_AB_AB_uninf + p_AB_BB_uninf + p_BB_AA_uninf + p_BB_AB_uninf + p_BB_BB_uninf) %>%
  
  # Calculate the total chance of the SNP being a type 1 SNP.
  mutate(SNP_type1 = p_AA_BB + p_BB_AA) %>%
  mutate(SNP_not_type1 = 1-SNP_type1)

#############################################################
# Uninformativeness of 24 SNP panel
#############################################################

# Check the probability of the 24 SNP panel being uninformative when using pre-amplified cfDNA 
# (type 1 and type 3 SNPs) 
snp_24_type1and3 <- snp_calc %>%
  # Selct only the first 24 SNPs
  filter(GOSH_ID_ds < 25) %>%
  select(dbSNP, population, SNP_uninf_total) %>%
  # Reshape by pivot_wider
  pivot_wider(id_cols = population,
              names_from = dbSNP,
              values_from = SNP_uninf_total) %>%
  mutate(chance_uninformative = rs565522 * rs2737654 * rs2576241 * rs1160680 * rs1914748 * 
           rs12694624 * rs1399629 *rs2276702 * rs6781236 * rs7653090 * 
           rs357485 * rs9290003 * rs17017347 * rs10027026 * rs4975819 * 
           rs6899022 * rs1185246 * rs6877199 * rs13218440 * rs6924733 * 
           rs2535290 * rs172275 * rs4644087 * rs12690832) %>%
  mutate(uninformative_1_in  = 1/chance_uninformative) %>%
  select(population, chance_uninformative, uninformative_1_in)

(median(snp_24_type1and3$chance_uninformative))*100

# Check the probability of the 24 SNP panel being uninformative when using parental gDNA 
# (type 1 SNPs only) 
snp_24_type1_only <- snp_calc %>%
  # Selct only the first 24 SNPs
  filter(GOSH_ID_ds < 25) %>%
  select(dbSNP, population, SNP_not_type1) %>%
  # Reshape by pivot_wider
  pivot_wider(id_cols = population,
              names_from = dbSNP,
              values_from = SNP_not_type1) %>%
  mutate(chance_uninformative = rs565522 * rs2737654 * rs2576241 * rs1160680 * rs1914748 * 
           rs12694624 * rs1399629 *rs2276702 * rs6781236 * rs7653090 * 
           rs357485 * rs9290003 * rs17017347 * rs10027026 * rs4975819 * 
           rs6899022 * rs1185246 * rs6877199 * rs13218440 * rs6924733 * 
           rs2535290 * rs172275 * rs4644087 * rs12690832) %>%
  mutate(uninformative_1_in  = 1/chance_uninformative) %>%
  select(population, chance_uninformative, uninformative_1_in)

mean(snp_24_type1_only$uninformative_1_in)

#############################################################
# Uninformativeness of 40 SNP panel
#############################################################

# Check the probability of the 40 SNP panel being uninformative when using pre-amplified cfDNA 
# (type 1 and type 3 SNPs) 
snp_40_type1and3 <- snp_calc %>%
  select(dbSNP, population, SNP_uninf_total) %>%
  # Reshape by pivot_wider
  pivot_wider(id_cols = population,
              names_from = dbSNP,
              values_from = SNP_uninf_total) %>%
  mutate(chance_uninformative = rs565522 * rs2737654 * rs2576241 * rs1160680 * rs1914748 * 
           rs12694624 * rs1399629 *rs2276702 * rs6781236 * rs7653090 * 
           rs357485 * rs9290003 * rs17017347 * rs10027026 * rs4975819 * 
           rs6899022 * rs1185246 * rs6877199 * rs13218440 * rs6924733 * 
           rs2535290 * rs172275 * rs4644087 * rs12690832 * rs7827391 * 
           rs2319150 * rs1410059 * rs10821808 * rs2370764 * rs1498553 * 
           rs2256111 * rs7325978 * rs3742560 * rs12148532 * rs249290 * 
           rs1544724 * rs3760269 * rs7233004 * rs4801945 * rs271981) %>%
  mutate(uninformative_1_in  = 1/chance_uninformative) %>%
  select(population, chance_uninformative, uninformative_1_in)

mean(snp_40_type1and3$uninformative_1_in)

# Check the probability of the 40 SNP panel being uninformative when using using parental gDNA 
# (type 1 SNPs only)
snp_40_type1_only <- snp_calc %>%
  select(dbSNP, population, SNP_not_type1) %>%
  # Reshape by pivot_wider
  pivot_wider(id_cols = population,
              names_from = dbSNP,
              values_from = SNP_not_type1) %>%
  mutate(chance_uninformative = rs565522 * rs2737654 * rs2576241 * rs1160680 * rs1914748 * 
           rs12694624 * rs1399629 *rs2276702 * rs6781236 * rs7653090 * 
           rs357485 * rs9290003 * rs17017347 * rs10027026 * rs4975819 * 
           rs6899022 * rs1185246 * rs6877199 * rs13218440 * rs6924733 * 
           rs2535290 * rs172275 * rs4644087 * rs12690832 * rs7827391 * 
           rs2319150 * rs1410059 * rs10821808 * rs2370764 * rs1498553 * 
           rs2256111 * rs7325978 * rs3742560 * rs12148532 * rs249290 * 
           rs1544724 * rs3760269 * rs7233004 * rs4801945 * rs271981) %>%
  mutate(uninformative_1_in  = 1/chance_uninformative) %>%
  select(population, chance_uninformative, uninformative_1_in)

mean(snp_40_type1_only$uninformative_1_in)

#############################################################
# Plot graph
#############################################################

ggplot(snp_calc, aes(x = GOSH_ID_ds, y = frequency_B))+
  geom_point(size = 5, alpha = 0.5, aes(colour =population))+
  ylim(0, 0.5)+
  labs(x = "SNPs", y = "GnomAD minor allele frequency", title = "GnomAD minor allele frequencies of Camunas-Soler et al (2018) ddPCR SNP panel")+
  theme_bw()



SNP_data_table <- SNP_data %>%
  select(Sample, Target, TargetType, Concentration) %>%
  filter(Sample != "NTC") %>%
  left_join(ddpcr_target_panel %>%
              select(Target, Assay), by = "Target") %>%
  # Convert the concentration column to a numeric and convert "no call" to zero
  mutate(Copies_per_ul = as.integer(ifelse(Concentration == "No Call", 0, Concentration))) %>%
  # Add fluorophores
  mutate(Label = case_when(
    TargetType == "Ch1Unknown" ~ "FAM",
    TargetType == "Ch2Unknown" ~ "VIC")) %>%
  select(!c(Target, Concentration, TargetType)) %>%
  pivot_wider(
    id_cols = Assay,
    names_from = c(Sample, Label),
    values_from = Copies_per_ul)

view(SNP_data_plotting)

SNP_data_plotting <- SNP_data %>%
  filter(Sample != "NTC") %>%
  select(Sample, Target, TargetType, Concentration) %>%
  left_join(ddpcr_target_panel %>%
              select(Target, Assay), by = "Target") %>%
  # Convert the concentration column to a numeric and convert "no call" to zero
  mutate(Copies_per_ul = as.integer(ifelse(Concentration == "No Call", 0, Concentration))) %>%
  # Add fluorophores
  mutate(Label = case_when(
    TargetType == "Ch1Unknown" ~ "FAM",
    TargetType == "Ch2Unknown" ~ "VIC")) %>%
  select(!c(Target, Concentration, TargetType)) %>%
  pivot_wider(id_cols = c(Sample,Assay),
              names_from = Label,
              values_from = Copies_per_ul) %>%
  mutate(genotype = case_when(
    FAM > 50 & VIC < 50 ~"hom FAM",
    FAM > 50 & VIC > 50 ~"het",
    FAM < 50 & VIC > 50 ~"hom VIC",)) %>%
  filter(!is.na(genotype)) %>%
  mutate(
    fluorophore_VIC = case_when(
      genotype != "hom FAM" ~"Yes",
      genotype == "hom FAM" ~"No"),
    fluorophore_FAM= case_when(
      genotype != "hom VIC" ~"Yes",
      genotype == "hom VIC" ~"No"))


sample_30065_SNP <- SNP_data_plotting %>%
  filter(Sample %in% c("30065", "21RG-126G0126")) %>%
  select(-c(fluorophore_VIC, fluorophore_FAM)) %>%
  pivot_wider(
    id_cols = Assay,
    names_from = Sample,
    values_from = c(FAM, VIC, genotype)) %>%
  select(Assay, FAM_30065, VIC_30065, `FAM_21RG-126G0126`,
         `VIC_21RG-126G0126`, genotype_30065, `genotype_21RG-126G0126`)

write.csv(sample_30065_SNP, "analysis_outputs/sample_30065_SNP.csv", row.names = FALSE)



?pivot_wider()



# Plots Fluidigm-like graphs of different genotype clusters
ggplot(SNP_data_plotting %>%
         filter(!is.na(genotype)) %>%
         filter(!Sample %in% c("21RG-070G0033", "20611", "30113") ), aes(x = FAM, y = VIC, colour = genotype))+
  geom_point()+
  facet_wrap(~Sample)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "FAM target copies/ul", y = "VIC target copies/ul",
       title = "SNP Genotyping for Parental Samples")

ggplot(SNP_data_plotting %>%, aes(x = genotype, y =  Assay, colour = genotype))+
  geom_point(size = 3, alpha = 0.7)+
  theme_bw()+
  facet_wrap(~Sample)+
  labs(x = "")

# Individual plots
ggplot(SNP_data_plotting %>%
        filter(Sample == "21RG-083G0126"), aes(x = genotype, y = Assay, colour = genotype))+
  geom_point(size = 3, alpha = 0.7)+
  theme_bw()


write.csv(SNP_data_plotting, "analysis_outputs/SNP_data_plotting.csv", row.names = FALSE)

SNP_data_plotting %>%
  filter(Assay == "rs9290003" & genotype == "hom VIC")

#############################################################
# Function for calculating fetal fraction for cfDNA samples
#############################################################

digital_snp <- function(cf_sample){
  sample_snp_table <- SNP_data %>%
    filter(Sample == cf_sample) %>%
    select(Sample, Target, TargetType, Concentration, FractionalAbundance) %>%
    
    left_join(ddpcr_target_panel %>%
              select(Target, Assay), by = "Target") %>%
  # Convert the concentration column to a numeric and convert "no call" to zero
  mutate(Copies_per_ul = as.integer(ifelse(Concentration == "No Call", 0, Concentration))) %>%
  
  # Determine fetal fraction (assuming it is between 1 and 20%)
  dplyr::rename(Fraction_a = FractionalAbundance) %>%
  mutate(Fraction_b = 100 - Fraction_a) %>%
  mutate(minor_fraction = pmin(Fraction_b, Fraction_a)) %>%
  mutate(fetal_fraction = case_when(
    minor_fraction > 1 & minor_fraction < 20 ~minor_fraction)) %>%
  
  # Filter to get one row per well
  filter(TargetType == "Ch1Unknown") %>%
  
  select(c(Assay, fetal_fraction))

  mean_ff <- mean(sample_snp_table$fetal_fraction)
  
  return(sample_snp_table)
}

new_samples <- c("21RG-112G0065", "21RG-112G0027", "21RG-112G0029", "21RG-112G0098",
                 "21RG-112G0102", "21RG-112G0034")

snp_table_new <- SNP_data %>%
  filter(Sample %in% new_samples) %>%
  select(Sample, Target, TargetType, Concentration, FractionalAbundance) %>%
  left_join(ddpcr_target_panel %>%
              select(Target, Assay), by = "Target") %>%
  # Convert the concentration column to a numeric and convert "no call" to zero
  mutate(Copies_per_ul = as.integer(ifelse(Concentration == "No Call", 0, Concentration))) %>%
  
  # Determine fetal fraction (assuming it is between 1 and 20%)
  dplyr::rename(Fraction_a = FractionalAbundance) %>%
  mutate(Fraction_b = 100 - Fraction_a) %>%
  mutate(minor_fraction = pmin(Fraction_b, Fraction_a)) %>%
  mutate(fetal_fraction = case_when(
    minor_fraction > 1 & minor_fraction < 20 ~minor_fraction)) %>%
  
  # Filter to get one row per well
  filter(TargetType == "Ch1Unknown") %>%
  
  filter(!is.na(fetal_fraction)) %>%
  
  select(c(Sample, Assay, fetal_fraction)) %>%
  arrange(Assay)

view(snp_table_new)

# Find the median fetal fraction for each sample.
median_ff <- snp_table_new %>%
  select(Sample, fetal_fraction) %>%
  group_by(Sample) %>% 
  summarise_all(median)

ggplot(snp_table_new, aes(x = Assay, y = fetal_fraction))+
  geom_point()+
  facet_wrap(~Sample)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  ylim(0, 10)



     