################################################################################
## Enzyme cut site analysis
## July 2021
## Joseph.Shaw@gosh.nhs.uk
## Analysis of occurence of DNASE1L3 end motifs from Serpas et al 
## (PMID: 30593563) in ddPCR assay amplicons
################################################################################

#########################
# Load resources
#########################

# Load packages
library(tidyverse)

# Read in amplicons
amplicons <- read_excel("W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/ddPCR Assay Design/ddPCR_Assay_Ordering.xlsx",
                        sheet = "amplicons") %>%
  filter(reference != "Not provided") %>%
  pivot_longer(
    cols = -assay_name,
    names_to = "amplicon_class",
    values_to = "sequence") %>%
  mutate(length = nchar(sequence))

#########################
# Count motifs
#########################

# Specify DNASE1L3 end motifs from Serpas et al PMID: 30593563 in 
# both forward and reverse complement.
dnase1l3_end_motifs <- c("CCCA", "CCAG", "CCTG", "CCAA", "CCCT", "CCAT", 
                         "TGGG", "CTGG", "CAGG", "TTGG", "AGGG", "ATGG") 

dnase1l3_motif_count <- c()

for (i in amplicons$sequence) {
  
  motif_count <- sum(str_count(i,dnase1l3_end_motifs))
  
  dnase1l3_motif_count <- c(dnase1l3_motif_count, motif_count)
  
  rm(motif_count)
}

amplicons_count <- cbind(amplicons, dnase1l3_motif_count) %>%
  mutate(bp_per_site = round(length/dnase1l3_motif_count,1))

#########################
# Plots
#########################

amplicons_count %>%
  filter(amplicon_class == "reference") %>%
  mutate(assay_name = fct_reorder(assay_name, bp_per_site)) %>%
  ggplot(aes(x = bp_per_site, y = assay_name)) +
  geom_col() +
  theme_bw() +
  labs( x = "Number of DNAse1L3 CC motifs in target amplicon",
        y = "ddPCR assay",
        title = "Frequency of CC motifs in target reference amplicons")

amplicons_count %>%
  pivot_wider(
    id_cols = assay_name,
    names_from = amplicon_class,
    values_from = c(length, dnase1l3_motif_count)) %>%
  mutate(cut_site_diff = abs(dnase1l3_motif_count_reference - 
                               dnase1l3_motif_count_variant),
         assay_name = fct_reorder(assay_name, cut_site_diff)) %>%
  ggplot(aes(x = cut_site_diff, y = assay_name)) +
  geom_point() +
  theme_bw()
