################################################################################
## Tables for Paper Redraft
################################################################################

# The original analysis was run on CGEN-D0W15R2 with:
# R v4.1.0
# Cmdstanr v0.4.0
# Cmdstan v2.26.1

# For table redrafts I am going to use the original results generated on 24/11/2021.

library(tidyverse)
library(readxl)

setwd("W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/ddPCR_R_Analysis/ddpcr_nipd")

supplementary_table <- read_csv("analysis_outputs/supplementary_table20211124_142542.csv")

######################################
# Sickle cell results table
######################################

sickle_cell_results <- supplementary_table %>%
  filter(vf_assay == "HBB c.20A>T")

###################
# SPRT
###################
sprt_hbss_correct <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous variant" 
                                         & sickle_cell_results$sprt_outcome == "correct", "sample_id"])
                          
sprt_hbss_inconclusive <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous variant" 
                                              & sickle_cell_results$sprt_outcome == "inconclusive", "sample_id"])                          
                          
sprt_hbss_incorrect <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous variant" 
                                                   & sickle_cell_results$sprt_outcome == "incorrect", "sample_id"])                          

sprt_hbas_correct <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "heterozygous" 
                                              & sickle_cell_results$sprt_outcome == "correct", "sample_id"])

sprt_hbas_inconclusive <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "heterozygous" 
                                                   & sickle_cell_results$sprt_outcome == "inconclusive", "sample_id"])                          

sprt_hbas_incorrect <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "heterozygous" 
                                                & sickle_cell_results$sprt_outcome == "incorrect", "sample_id"])  

sprt_hbaa_correct <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous reference" 
                                              & sickle_cell_results$sprt_outcome == "correct", "sample_id"])

sprt_hbaa_inconclusive <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous reference" 
                                                   & sickle_cell_results$sprt_outcome == "inconclusive", "sample_id"])                          

sprt_hbaa_incorrect <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous reference" 
                                                & sickle_cell_results$sprt_outcome == "incorrect", "sample_id"])  

###################
# MCMC
###################
mcmc_hbss_correct <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous variant" 
                                              & sickle_cell_results$mcmc_outcome == "correct", "sample_id"])

mcmc_hbss_inconclusive <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous variant" 
                                              & sickle_cell_results$mcmc_outcome == "inconclusive", "sample_id"])

mcmc_hbss_incorrect <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous variant" 
                                                   & sickle_cell_results$mcmc_outcome == "incorrect", "sample_id"])

mcmc_hbas_correct <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "heterozygous" 
                                              & sickle_cell_results$mcmc_outcome == "correct", "sample_id"])

mcmc_hbas_inconclusive <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "heterozygous" 
                                                   & sickle_cell_results$mcmc_outcome == "inconclusive", "sample_id"])

mcmc_hbas_incorrect <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "heterozygous" 
                                                & sickle_cell_results$mcmc_outcome == "incorrect", "sample_id"])

mcmc_hbaa_correct <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous reference" 
                                              & sickle_cell_results$mcmc_outcome == "correct", "sample_id"])

mcmc_hbaa_inconclusive <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous reference" 
                                                   & sickle_cell_results$mcmc_outcome == "inconclusive", "sample_id"])

mcmc_hbaa_incorrect <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous reference" 
                                                & sickle_cell_results$mcmc_outcome == "incorrect", "sample_id"])

###################
# Z score
###################

zscore_hbss_correct <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous variant" 
                                              & sickle_cell_results$zscore_outcome == "correct", "sample_id"])

zscore_hbss_inconclusive <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous variant" 
                                                & sickle_cell_results$zscore_outcome == "inconclusive", "sample_id"])

zscore_hbss_incorrect <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous variant" 
                                                     & sickle_cell_results$zscore_outcome == "incorrect", "sample_id"])

zscore_hbas_correct <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "heterozygous" 
                                                & sickle_cell_results$zscore_outcome == "correct", "sample_id"])

zscore_hbas_inconclusive <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "heterozygous" 
                                                     & sickle_cell_results$zscore_outcome == "inconclusive", "sample_id"])

zscore_hbas_incorrect <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "heterozygous" 
                                                  & sickle_cell_results$zscore_outcome == "incorrect", "sample_id"])

zscore_hbaa_correct <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous reference" 
                                                & sickle_cell_results$zscore_outcome == "correct", "sample_id"])

zscore_hbaa_inconclusive <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous reference" 
                                                     & sickle_cell_results$zscore_outcome == "inconclusive", "sample_id"])

zscore_hbaa_incorrect <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous reference" 
                                                  & sickle_cell_results$zscore_outcome == "incorrect", "sample_id"])

###################
# New table
###################

sickle_table_new <- data.frame(
  Outcome = c("Correct", "Inconclusive", "Incorrect"),
  HbSS = c(sprt_hbss_correct, sprt_hbss_inconclusive, sprt_hbss_incorrect),
  HbAS = c(sprt_hbas_correct, sprt_hbas_inconclusive, sprt_hbas_incorrect),
  HbAA = c(sprt_hbaa_correct, sprt_hbaa_inconclusive, sprt_hbaa_incorrect),
  HbSS = c(mcmc_hbss_correct, mcmc_hbss_inconclusive, mcmc_hbss_incorrect),
  HbAS = c(mcmc_hbas_correct, mcmc_hbas_inconclusive, mcmc_hbas_incorrect),
  HbAA = c(mcmc_hbaa_correct, mcmc_hbaa_inconclusive, mcmc_hbaa_incorrect),
  HbSS = c(zscore_hbss_correct, zscore_hbss_inconclusive, zscore_hbss_incorrect),
  HbAS = c(zscore_hbas_correct, zscore_hbas_inconclusive, zscore_hbas_incorrect),
  HbAA = c(zscore_hbaa_correct, zscore_hbaa_inconclusive, zscore_hbaa_incorrect))

######################################
# Bespoke results table
######################################

bespoke_sprt <- supplementary_table %>%
  filter(vf_assay != "HBB c.20A>T") %>%
  group_by(sprt_outcome) %>%
  summarise(sprt = n()) %>%
  dplyr::rename(outcome = sprt_outcome)

bespoke_mcmc <- supplementary_table %>%
  filter(vf_assay != "HBB c.20A>T") %>%
  group_by(mcmc_outcome) %>%
  summarise(mcmc = n()) %>%
  select(-mcmc_outcome)

bespoke_zscore <- supplementary_table %>%
  filter(vf_assay != "HBB c.20A>T") %>%
  group_by(zscore_outcome) %>%
  summarise("Z score" = n()) %>%
  select(-zscore_outcome)

bespoke_table_new <- cbind(bespoke_sprt, bespoke_mcmc, bespoke_zscore)

######################################
# Streamlined supplementary table
######################################

# Restructure table to make it simpler.

supp_2 <- supplementary_table %>%
  select(-c("variant_percent_max", "variant_percent_min", "vf_assay_molecules_max", 
            "vf_assay_molecules_min", "fetal_percent_max", "fetal_percent_min", 
            "ff_assay_molecules_max", "ff_assay_molecules_min",
            "vf_assay_molecules", "ff_assay_molecules")) %>%
  dplyr::rename(variant_fraction = variant_percent) %>%
  dplyr::rename(fetal_fraction = fetal_percent)

################################################################################