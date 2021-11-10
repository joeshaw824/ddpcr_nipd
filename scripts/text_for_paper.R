##############################
# Text for paper
##############################

################
# Abstract
################

paste0(
  nrow(supplementary_table),
  " cfDNA samples were tested in total: ",
  nrow(supplementary_table %>%
         filter(vf_assay == "HBB c.20A>T")),
  " cases for sickle cell disease, ",
  nrow(supplementary_table %>%
         filter(inheritance_chromosomal == "autosomal" &
                  inheritance_pattern == "dominant")),
  " for autosomal dominant conditions, ",
  nrow(supplementary_table %>%
         filter(inheritance_chromosomal == "autosomal" &
                  inheritance_pattern == "recessive" &
                  vf_assay != "HBB c.20A>T")),
  " for rare recessive conditions, and ",
  nrow(supplementary_table %>%
         filter(inheritance_chromosomal == "x_linked")),
  " for x-linked conditions.")

paste0(
  "Overall, SPRT analysis classified ",
  nrow(supplementary_table %>%
         filter(outcome_sprt == "correct")),
  " fetal genotypes correctly out of ",
  nrow(supplementary_table %>%
         filter(outcome_sprt %in% c("correct", "incorrect"))), 
  " predictions, whilst MCMC classified ", 
  nrow(supplementary_table %>%
         filter(outcome_mcmc == "correct")),
  " fetal genotypes correctly out of ",
  nrow(supplementary_table %>%
         filter(outcome_mcmc %in% c("correct", "incorrect"))),
  "predictions.")

# However, both approaches gave erroneous results when applied to gDNA data, 
paste0("and misclassified ",
       nrow(supplementary_table %>%
              filter(outcome_mcmc == "incorrect" |
                       outcome_sprt == "incorrect")),
       " cfDNA samples in total,")

# indicating that incorrect results came from interpretation of technical 
# variation as a fetal signal. In contrast, z score analysis accounted for 
# technical variation and 
paste0("correctly classified ",
       nrow(supplementary_table %>%
              filter(outcome_zscore == "correct")),
       " fetal genotypes out of ",
       nrow(supplementary_table %>%
              filter(outcome_zscore %in% c("correct", "incorrect"))),
       " predictions, but also generated ", 
       nrow(supplementary_table %>%
              filter(outcome_zscore == "incorrect")),
       " incorrect fetal genotype predictions.")

################

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
# "For the bespoke cohort, the two methods again performed similarly well,"
paste0("with the SPRT analysis correctly predicting ",
       nrow(all_samples_unblinded %>%
              filter(vf_assay != "HBB c.20A>T" &
                       outcome_sprt == "correct")),
       " fetal genotypes and the MCMC analysis correctly predicting ",
       nrow(all_samples_unblinded %>%
              filter(vf_assay != "HBB c.20A>T" &
                       outcome_mcmc == "correct")),
       ".") 

# Correct fetal genotypes were predicted by both the SPRT and MCMC methods
paste0("for 5 pregnancies with a risk of Aicardi-Goutières syndrome (",
       supplementary_table[supplementary_table$vf_assay == "RNASEH2C c.205C>T",
                           "sample_id"],
       " and ",
       supplementary_table[supplementary_table$vf_assay == "ADAR c.2997G>T",
                           "sample_id"],
       "), vitamin B12-responsive methylmalonic aciduria (",
       supplementary_table[supplementary_table$vf_assay == "MMAA c.733+1G>A",
                           "sample_id"],
       "), congenital disorder of glycosylation type 1a (CDG1a) (",
       supplementary_table[supplementary_table$vf_assay == "PMM2 c.691G>A" &
                             supplementary_table$outcome_sprt == "correct",
                           "sample_id"],
       ") and severe combined immunodeficiency (",
       supplementary_table[supplementary_table$vf_assay == "ADA c.556G>A",
                           "sample_id"],
       ").")

paste0("A sample for a later pregnancy of the CDG1a couple (",
       supplementary_table[supplementary_table$vf_assay == "PMM2 c.691G>A" &
                             supplementary_table$outcome_sprt == "inconclusive",
                           "sample_id"],
       ") was deemed inconclusive by both methods,")

# "which was likely a combination of"

paste0("a low fetal fraction (",
       round(supplementary_table[supplementary_table$vf_assay == "PMM2 c.691G>A" &
                                   supplementary_table$outcome_sprt == "inconclusive",
                                 "fetal_percent"], 1),
       " %) and low cfDNA concentration (",
       round(supplementary_table[supplementary_table$vf_assay == "PMM2 c.691G>A" &
                                   supplementary_table$outcome_sprt == "inconclusive",
                                 "totalGE_ml_plasma"], 0),
       " GE/ml).")


#The variant fraction of the heterozygous gDNA controls 
paste0("showed substantial variation (",
       round(min(het_gdna$variant_percent), 1),
       "-",
       round(max(het_gdna$variant_percent), 1),
       "%) when fewer than 2000 haploid genome equivalents (GE) were measured.")  

# The additional classification rules generated a higher number of 
# inconclusive results for the z score approach, 
paste0("although ",
       nrow(all_samples_unblinded %>%
              filter(vf_assay == "HBB c.20A>T" &
                       outcome_zscore == "correct")),
       " correct fetal genotype predictions were made for the sickle cell disease assay, and ",
       nrow(all_samples_unblinded %>%
              filter(vf_assay != "HBB c.20A>T" &
                       outcome_zscore == "correct")),
       " for the bespoke design cohort.")

# Figure 2 legend
paste0("The ddPCR variant fraction results for ",
       length(unique(het_gdna$worksheet_well_sample)),
       " replicates of ",
       length(unique(het_gdna$r_number)),
       " heterozygous parental gDNA samples at varying DNA inputs.")


#"For the HBB c.20A>T assay, a sample with a 

paste0("variant fraction of ",
       round(supplementary_table[supplementary_table$vf_assay == "HBB c.20A>T" &
                                   supplementary_table$outcome_zscore == "incorrect",
                                 "variant_percent"], 1),
       "% and a fetal fraction of ",
       round(supplementary_table[supplementary_table$vf_assay == "HBB c.20A>T" &
                                   supplementary_table$outcome_zscore == "incorrect",
                                 "fetal_percent"], 1),
       "% (",
       supplementary_table[supplementary_table$vf_assay == "HBB c.20A>T" &
                             supplementary_table$outcome_zscore == "incorrect",
                           "sample_id"],
       ") was incorrectly predicted to have a heterozygous fetus"
)

# In the bespoke design cohort, a false positive result was generated for a 
# cfDNA sample from a woman who was heterozygous for 
# the IDS c.182_189del variant, 
paste0("which causes mucopolysaccharidosis type 2 (",
       supplementary_table[supplementary_table$vf_assay != "HBB c.20A>T" &
                             supplementary_table$outcome_zscore == "incorrect",
                           "sample_id"],
       ").")

# The ZFX/ZFY ddPCR assay detected a fetal fraction 
paste0("of ",
       round(supplementary_table[supplementary_table$vf_assay != "HBB c.20A>T" &
                                   supplementary_table$outcome_zscore == "incorrect",
                                 "fetal_percent"], 1),
       "%, and there was a significant increase in the fraction of the pathogenic variant in the cfDNA (",
       round(supplementary_table[supplementary_table$vf_assay != "HBB c.20A>T" &
                                   supplementary_table$outcome_zscore == "incorrect",
                                 "variant_percent"], 1),
       "%) when ",
       round(supplementary_table[supplementary_table$vf_assay != "HBB c.20A>T" &
                                   supplementary_table$outcome_zscore == "incorrect",
                                 "vf_assay_molecules"], 1),
       " molecules were measured.")


#In addition, our cohort had a lower average 
paste0("number of molecules tested per case (median: ",
       round(median(supplementary_table$total_molecules), 0),
       ").")

###################