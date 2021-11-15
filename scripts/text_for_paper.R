##############################
# Text for paper
##############################

# Load tidyverse and pipeline exports if not performing this 
# in the same R session as running the pipeline.

#library(tidyverse)

#supplementary_table <- read.csv("W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/Digital PCR Paper Drafts/ddPCR cohort paper/ddPCR cohort paper v2/Supplementary_data_20211110_170750.csv")

#vf_assay_count_table <- read.csv("W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/Digital PCR Paper Drafts/ddPCR cohort paper/ddPCR cohort paper v2/Table 1 20211110_171532.csv")

######################
# Abstract
######################

paste0(
  nrow(all_samples_blinded),
  " cfDNA samples were tested in total: ",
  nrow(all_samples_blinded %>%
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
         filter(sprt_outcome == "correct")),
  " fetal genotypes correctly out of ",
  nrow(supplementary_table %>%
         filter(sprt_outcome %in% c("correct", "incorrect"))), 
  " predictions, whilst MCMC classified ", 
  nrow(supplementary_table %>%
         filter(mcmc_outcome == "correct")),
  " fetal genotypes correctly out of ",
  nrow(supplementary_table %>%
         filter(mcmc_outcome %in% c("correct", "incorrect"))),
  " predictions.")

# In addition, 
paste0(nrow(supplementary_table %>%
              filter(mcmc_outcome == "incorrect" |
                       sprt_outcome == "incorrect")),
       " cfDNA samples were misclassified in total,")

# indicating that incorrect results came from interpretation of technical 
# variation as a fetal signal. In contrast, z score analysis accounted for 
# technical variation and 
paste0("correctly classified ",
       nrow(supplementary_table %>%
              filter(zscore_outcome == "correct")),
       " fetal genotypes out of ",
       nrow(supplementary_table %>%
              filter(zscore_outcome %in% c("correct", "incorrect"))),
       " predictions, but also generated ", 
       nrow(supplementary_table %>%
              filter(zscore_outcome == "incorrect")),
       " incorrect fetal genotype predictions.")

######################
# Methods
######################

paste0(
  "In total, ",
  nrow(all_samples_blinded),
  " samples were analysed during this study: ",
  nrow(all_samples_blinded %>%
         filter(vf_assay == "HBB c.20A>T")),
  " were for SCD (HBB c.20A>T) and ",
  nrow(all_samples_blinded %>%
         filter(vf_assay != "HBB c.20A>T")),
  " for a range of variants in different disease genes")


paste0("Variant discrimination assays were designed for ",
       length(unique(supplementary_table$vf_assay)),
       " pathogenic variants: ",
       nrow(vf_assay_count_table %>%
              filter(Inheritance == "XLR")),
       " for X-linked recessive conditions, ",
       nrow(vf_assay_count_table %>%
              filter(Inheritance == "XLD")),
       " for an X-linked dominant condition, ",
       nrow(vf_assay_count_table %>%
              filter(Inheritance == "AD")),
       " for autosomal dominant conditions, and ",
       nrow(vf_assay_count_table %>%
              filter(Inheritance == "AR")),
       " for autosomal recessive conditions.")

######################
# Results
######################

paste0("Overall, we analysed ", nrow(all_samples_blinded), 
       " samples from ", 
       nrow(families),
       " families, for ", 
       length(unique(all_samples_blinded$vf_assay)), 
       " pathogenic variants in ",
       length(unique(vf_assay_count_table$Gene)),
       " genes.")

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
       nrow(supplementary_table %>%
              filter(vf_assay == "HBB c.20A>T" &
                       sprt_outcome == "correct")),
       " and ",
       nrow(supplementary_table %>%
              filter(vf_assay == "HBB c.20A>T" &
                       mcmc_outcome == "correct")),
       " correct fetal genotype predictions, respectively, including ",
       nrow(supplementary_table %>%
              filter(vf_assay == "HBB c.20A>T" &
                       fetal_genotype == "homozygous variant",
                     sprt_outcome == "correct")),
       " and ",
       nrow(supplementary_table %>%
              filter(vf_assay == "HBB c.20A>T" &
                       fetal_genotype == "homozygous variant",
                     mcmc_outcome == "correct")),
       " fetuses affected with SCD. However, each method also generated ",
       nrow(supplementary_table %>%
              filter(vf_assay == "HBB c.20A>T" &
                       sprt_outcome == "incorrect")),
       " and ",
       nrow(supplementary_table %>%
              filter(vf_assay == "HBB c.20A>T" &
                       mcmc_outcome == "incorrect")),
       " incorrect fetal gentoype predictions.")

# Example of illogical predictions: samples 30065 and 20763
paste0("For example, sample ",
       supplementary_table[supplementary_table$r_number == "30065",
                           "sample_id"],
       " which had a variant fraction of ",
       round(supplementary_table[supplementary_table$r_number == "30065",
                           "variant_percent"],1),
       "%, was incorrectly predicted")
paste0("By contrast, sample ",
       supplementary_table[supplementary_table$r_number == "20763",
                           "sample_id"],
       " which had a lower variant fraction of ",
       round(supplementary_table[supplementary_table$r_number == "20763",
                           "variant_percent"], 1),
       "%, was correctly called")

##############
# Bespoke cohort
##############

# "For the bespoke cohort, the two methods again performed similarly well,"
paste0("with the SPRT analysis correctly predicting ",
       nrow(supplementary_table %>%
              filter(vf_assay != "HBB c.20A>T" &
                       sprt_outcome == "correct")),
       " fetal genotypes and the MCMC analysis correctly predicting ",
       nrow(supplementary_table %>%
              filter(vf_assay != "HBB c.20A>T" &
                       mcmc_outcome == "correct")),
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
                             supplementary_table$sprt_outcome == "correct",
                           "sample_id"],
       ") and severe combined immunodeficiency (",
       supplementary_table[supplementary_table$vf_assay == "ADA c.556G>A",
                           "sample_id"],
       ").")

paste0(
  "The CDG1a family (family ",
  supplementary_table[supplementary_table$vf_assay == "PMM2 c.691G>A"&
                        supplementary_table$sprt_outcome == "correct",
                      "family_number"],
  ") who were known to be consanguineous also donated a sample for a later pregnancy (",
  supplementary_table[supplementary_table$vf_assay == "PMM2 c.691G>A" &
                        supplementary_table$sprt_outcome == "inconclusive",
                      "sample_id"],
  ") which was deemed inconclusive")

# "which was likely a combination of"
paste0("a low fetal fraction (",
       round(supplementary_table[supplementary_table$vf_assay == "PMM2 c.691G>A" &
                                   supplementary_table$sprt_outcome == "inconclusive",
                                 "fetal_percent"], 1),
       "%) and low cfDNA concentration (",
       round(supplementary_table[supplementary_table$vf_assay == "PMM2 c.691G>A" &
                                   supplementary_table$sprt_outcome == "inconclusive",
                                 "totalGE_ml_plasma"], 0),
       " GE/ml).")

# However, both analysis methods generated false positive
# and false negative results for two X-linked recessive
# variants in the IDS and ABCD1 genes
paste0("(",
       supplementary_table[supplementary_table$vf_assay == "IDS c.182_189del" &
                             supplementary_table$sprt_outcome == "incorrect" &
                             supplementary_table$mcmc_outcome == "incorrect",
                           "sample_id"],
       " and ",
       supplementary_table[supplementary_table$vf_assay == "ABCD1 c.3G>A" &
                             supplementary_table$sprt_outcome == "incorrect" &
                             supplementary_table$mcmc_outcome == "incorrect",
                           "sample_id"],
       ")")

##############
# Z score analysis
##############

#The variant fraction of the heterozygous gDNA controls 
paste0("showed substantial variation (",
       round(min(het_gdna$variant_percent), 1),
       "-",
       round(max(het_gdna$variant_percent), 1),
       "%) when fewer than ",
       vf_assay_molecules_limit,
       " haploid genome equivalents (GE) were measured.")  

# The additional classification rules generated a higher number of 
# inconclusive results for the z score approach, 
paste0("although ",
       nrow(supplementary_table %>%
              filter(vf_assay == "HBB c.20A>T" &
                       zscore_outcome == "correct")),
       " correct fetal genotype predictions were made for the sickle cell disease assay, and ",
       nrow(supplementary_table %>%
              filter(vf_assay != "HBB c.20A>T" &
                       zscore_outcome == "correct")),
       " for the bespoke design cohort.")

# For the 6 samples with incorrect predictions from SPRT and MCMC analysis,
paste0("z score analysis gave ",
       nrow(supplementary_table %>%
              filter(zscore_outcome == "incorrect")),
       " incorrect predictions")

##############
# Sensitivity and incorrect predictions
##############

#"For the HBB c.20A>T assay, a sample with a 
paste0("variant fraction of ",
       round(supplementary_table[supplementary_table$vf_assay == "HBB c.20A>T" &
                                   supplementary_table$zscore_outcome == "incorrect",
                                 "variant_percent"], 1),
       "% and a fetal fraction of ",
       round(supplementary_table[supplementary_table$vf_assay == "HBB c.20A>T" &
                                   supplementary_table$zscore_outcome == "incorrect",
                                 "fetal_percent"], 1),
       "% (",
       supplementary_table[supplementary_table$vf_assay == "HBB c.20A>T" &
                             supplementary_table$zscore_outcome == "incorrect",
                           "sample_id"],
       ") was incorrectly predicted to have a heterozygous fetus")

# In the bespoke design cohort, a false positive result was generated for a 
# cfDNA sample from a woman who was heterozygous for 
# the IDS c.182_189del variant, 
paste0("which causes mucopolysaccharidosis type 2 (",
       supplementary_table[supplementary_table$vf_assay != "HBB c.20A>T" &
                             supplementary_table$zscore_outcome == "incorrect",
                           "sample_id"],
       ").")

# The ZFX/ZFY ddPCR assay detected a fetal fraction 
paste0("of ",
       round(supplementary_table[supplementary_table$vf_assay != "HBB c.20A>T" &
                                   supplementary_table$zscore_outcome == "incorrect",
                                 "fetal_percent"], 1),
       "%, and there was a significant increase in the fraction of the pathogenic variant in the cfDNA (",
       round(supplementary_table[supplementary_table$vf_assay != "HBB c.20A>T" &
                                   supplementary_table$zscore_outcome == "incorrect",
                                 "variant_percent"], 1),
       "%) when ",
       round(supplementary_table[supplementary_table$vf_assay != "HBB c.20A>T" &
                                   supplementary_table$zscore_outcome == "incorrect",
                                 "vf_assay_molecules"], 1),
       " molecules were measured.")

######################
# Conclusions
######################

#In addition, our cohort had a lower average 
paste0("number of molecules tested per case (median: ",
       round(median(supplementary_table$total_molecules), 0),
       ").")

######################
# Figures
######################

# Figure 1 sample numbers
#Sickle cell disease
nrow(all_samples_blinded %>%
       filter(vf_assay == "HBB c.20A>T"))

# Bespoke cohort
nrow(all_samples_blinded %>%
       filter(vf_assay != "HBB c.20A>T"))

# X-linked samples
nrow(all_samples_blinded %>%
       filter(inheritance_chromosomal == "x_linked"))

# Autosomal dominant samples
nrow(all_samples_blinded %>%
       filter(inheritance_chromosomal == "autosomal" &
                inheritance_pattern == "dominant"))

# Autosomal recessive samples (not sickle cell disease)
nrow(all_samples_blinded %>%
       filter(inheritance_chromosomal == "autosomal" &
                inheritance_pattern == "recessive" &
                vf_assay != "HBB c.20A>T"))

unique(supplementary_table$ff_determination)

# Supplemental Figure 1
nrow(supplementary_table %>%
       filter(ff_determination == "NGS SNP panel - workflow 1"))

nrow(supplementary_table %>%
       filter(ff_determination == "ddPCR SNP panel - workflow 2"))

# Figure 2 legend
paste0("The ddPCR variant fraction results for ",
       length(unique(het_gdna$worksheet_well_sample)),
       " replicates of ",
       length(unique(het_gdna$r_number)),
       " heterozygous parental gDNA samples at varying DNA inputs.")

######################