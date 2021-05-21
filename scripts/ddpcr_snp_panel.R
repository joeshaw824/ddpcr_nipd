#############################################################
## ddPCR SNP panel analysis
## May 2021
## Joseph.Shaw@gosh.nhs.uk 
## This script is for calculating the probability of the 
## ddPCR SNP panel being uninformative for fetal fraction.
## This panel was originally described by Camunas-Soler et
## al (2018) (PMID: 29097507)
#############################################################

#############################################################
# Load packages and panel GnomAD data
#############################################################

library(tidyverse)

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
  
  # Probability is 0 for a type 2 SNP.
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

#############################################################
# Plot graph
#############################################################

ggplot(snp_calc, aes(x = GOSH_ID_ds, y = frequency_B))+
  geom_point(size = 5, alpha = 0.5, aes(colour =population))+
  ylim(0, 0.5)+
  labs(x = "SNPs", y = "GnomAD minor allele frequency", title = "GnomAD minor allele frequencies of Camunas-Soler et al (2018) ddPCR SNP panel")+
  theme_bw()

