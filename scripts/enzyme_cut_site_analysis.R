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
                        sheet = "amplicons")

#########################
# Count motifs
#########################

# Specify DNASE1L3 end motifs from Serpas et al PMID: 30593563 in 
# both forward and reverse complement.
dnase1l3_end_motifs <- c("CCCA", "CCAG", "CCTG", "CCAA", "CCCT", "CCAT", 
                         "TGGG", "CTGG", "CAGG", "TTGG", "AGGG", "ATGG") 


dnase1l3_motif_count <- c()

for (i in amplicons$amplicon_reference) {
  
  motif_count <- sum(str_count(i,dnase1l3_end_motifs))
  
  dnase1l3_motif_count <- c(dnase1l3_motif_count, motif_count)
  
  rm(motif_count)
}

amplicons_count <- cbind(amplicons, dnase1l3_motif_count)

#########################
# Plots
#########################

ggplot(amplicons_count %>%
         mutate(assay_name = fct_reorder(assay_name, dnase1l3_motif_count)), 
       aes(x = dnase1l3_motif_count, y = assay_name)) +
  geom_col() +
  theme_bw() +
  labs( x = "Number of DNAse1L3 CC motifs in target amplicon",
        y = "ddPCR assay",
        title = "Frequency of CC motifs in target reference amplicons")
