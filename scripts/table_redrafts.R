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

write.csv(bespoke_table_new, 
          file = (paste0("analysis_outputs/bespoke_table_new",
                         format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")),
          row.names = FALSE)

######################################
# Streamlined supplementary table
######################################

# Restructure table to make it simpler following Lizzie's comment.

supp_2 <- supplementary_table %>%
  select(-c("variant_percent_max", "variant_percent_min", "vf_assay_molecules_max", 
            "vf_assay_molecules_min", "fetal_percent_max", "fetal_percent_min", 
            "ff_assay_molecules_max", "ff_assay_molecules_min",
            "vf_assay_molecules", "ff_assay_molecules")) %>%
  dplyr::rename(`variant_fraction_%` = variant_percent,
                `fetal_fraction_%` = fetal_percent,
                family = family_number,
                partner_sample_available = partner_sample,
                variant_fraction_assay = vf_assay,
                fetal_fraction_determination = ff_determination,
                fetal_fraction_assay =  ff_assay,
                sprt_likelihood_ratio = likelihood_ratio,
                variant_positive_droplets = variant_positives,
                reference_positive_droplets = reference_positives,
                maternal_marker_positive_droplets = maternal_positives,
                paternal_marker_positive_droplets = paternal_positives,
                maternal_marker_molecules = maternal_molecules,
                paternal_marker_molecules = paternal_molecules) %>%
  select(sample_id, family, cohort, inheritance, condition, gestation, partner_sample_available,
         variant_fraction_assay, fetal_fraction_determination, fetal_fraction_assay, `variant_fraction_%`, `fetal_fraction_%`,
         sprt_prediction, mcmc_prediction, zscore_prediction, fetal_genotype,
         sprt_outcome, mcmc_outcome, zscore_outcome, 
         # sampling information
         vacutainer, hours_to_first_spin, days_to_storage, diagnostic_sampling,
         plasma_volume_ml, extraction_replicates,
         # ddpcr information,
         vf_assay_num_wells, ff_assay_num_wells, variant_positive_droplets, reference_positive_droplets,
         vf_assay_droplets, maternal_marker_positive_droplets, paternal_marker_positive_droplets, variant_molecules,
         reference_molecules, maternal_marker_molecules, paternal_marker_molecules, 
         totalGE_ml_plasma,
         # analysis information
         sprt_likelihood_ratio, p_G0, p_G1, p_G2, p_G3, zscore)

draft_filepath <- "W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/Digital PCR Paper Drafts/ddPCR cohort paper/ddPCR cohort paper v5"

write.csv(supp_2, 
          file = (paste0(draft_filepath, "supplementary_table_rearranged_",
                         format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")),
          row.names = FALSE)

######################################
# Adding transcript IDs to bespoke results table
######################################

bespoke_results <- read.csv("analysis_outputs/bespoke_results_table20211124_142605.csv")

vf_assay_info <- read.csv("resources/vf_assay_gene_information.csv")

new_bespoke_results_table <- bespoke_results %>%
  left_join(vf_assay_info %>%
              dplyr::rename(DNA = variant_dna) %>%
              select(transcript, DNA), by = "DNA") %>%
  select(Inheritance, Sample.number, Condition, transcript, Gene, DNA,
         Fetal.genotype, SPRT, MCMC, Z.score)

write.csv(new_bespoke_results_table, 
          file = (paste0("analysis_outputs/new_bespoke_results_table",
                         format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")),
          row.names = FALSE)

######################################
# Analysis results for entire cohort
######################################

total_cases <- nrow(supplementary_table)

sprt_correct <- nrow(supplementary_table[supplementary_table$sprt_outcome == "correct", "sample_id"])

mcmc_correct <- nrow(supplementary_table[supplementary_table$mcmc_outcome == "correct", "sample_id"])

zscore_correct <- nrow(supplementary_table[supplementary_table$zscore_outcome == "correct", "sample_id"])

sprt_inconclusive <- nrow(supplementary_table[supplementary_table$sprt_outcome == "inconclusive", "sample_id"])

mcmc_inconclusive <- nrow(supplementary_table[supplementary_table$mcmc_outcome == "inconclusive", "sample_id"])

zscore_inconclusive <- nrow(supplementary_table[supplementary_table$zscore_outcome == "inconclusive", "sample_id"])

sprt_incorrect <- nrow(supplementary_table[supplementary_table$sprt_outcome == "incorrect", "sample_id"])

mcmc_incorrect <- nrow(supplementary_table[supplementary_table$mcmc_outcome == "incorrect", "sample_id"])

zscore_incorrect <- nrow(supplementary_table[supplementary_table$zscore_outcome == "incorrect", "sample_id"])

paste0("Overall, across both the SCD and bespoke cohorts, SPRT, MCMC and z-score analysis correctly classified ",
       round((sprt_correct/total_cases)*100, 0), "%, ",
       round((mcmc_correct/total_cases)*100, 0), "% and ",
       round((zscore_correct/total_cases)*100, 0), "% of fetal samples in total, with ",
       round((sprt_inconclusive/total_cases)*100, 0), "%, ",
       round((mcmc_inconclusive/total_cases)*100, 0), "% and ",
       round((zscore_inconclusive/total_cases)*100, 0), "% inconclusive results and ",
       round((sprt_incorrect/total_cases)*100, 0), "%, ",
       round((mcmc_incorrect/total_cases)*100, 0), "% and ",
       round((zscore_incorrect/total_cases)*100, 0), "% discordant results, respectively"
       )

# Inconclusive results removed.

sprt_calls <- nrow(supplementary_table %>% filter(sprt_outcome != "inconclusive"))

mcmc_calls <- nrow(supplementary_table %>% filter(mcmc_outcome != "inconclusive"))

zscore_calls <- nrow(supplementary_table %>% filter(zscore_outcome != "inconclusive"))

paste0("By removing the inconclusive results and thus looking at only the total number of 
       predictions made by each analysis method, the correct classification rates increase to ",
       round((sprt_correct/sprt_calls)*100, 0), "% (",
       sprt_correct, "/", sprt_calls, ") for SPRT, ",
       round((mcmc_correct/mcmc_calls)*100, 0), "% (",
       mcmc_correct, "/", mcmc_calls, ") for MCMC, and ",
       round((zscore_correct/zscore_calls)*100, 0), "% (",
       zscore_correct, "/", zscore_calls, ") for Z score.")

################################################################################