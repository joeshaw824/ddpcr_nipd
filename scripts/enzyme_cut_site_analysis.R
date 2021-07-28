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
  filter(reference_amplicon != "Not provided") %>%
  mutate(reference_amplicon_length = nchar(reference_amplicon),
         variant_amplicon_length = nchar(variant_amplicon))

primer_probes <- read_excel("W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/ddPCR Assay Design/ddPCR_Assay_Ordering.xlsx",
                        sheet = "sequences") %>%
  filter(forward_primer != "Not provided") %>%
  # Remove "+" from locked nucleic acid probes
  mutate(variant_probe_no_lnas = gsub("+","",variant_probe, fixed = TRUE),
         reference_probe_no_lnas = gsub("+","",reference_probe, fixed = TRUE),
         
         # Create reverse complement primers
         reverse_primer_rc = reverse_complement(reverse_primer),
         forward_primer_rc = reverse_complement(forward_primer))
         
amplicons_primers_probes <- full_join(
  amplicons,
  primer_probes,
  by = "assay_name") %>%
  # Checks
  mutate(
    var_in_var_amp = str_detect(variant_amplicon, variant_probe_no_lnas),
    ref_in_ref_amp = str_detect(reference_amplicon, reference_probe_no_lnas),
    var_in_ref_amp = str_detect(variant_amplicon, reference_probe_no_lnas),
    ref_in_var_amp = str_detect(reference_amplicon, variant_probe_no_lnas),
    # Check forward primer
    forward_in_ref_amp = str_detect(reference_amplicon, forward_primer),
    forward_rc_in_ref_amp = str_detect(reference_amplicon, forward_primer_rc),
    forward_in_var_amp = str_detect(variant_amplicon, forward_primer),
    forward_rc_in_var_amp = str_detect(variant_amplicon, forward_primer_rc),
    
    reverse_in_ref_amp = str_detect(reference_amplicon, reverse_primer),
    reverse_rc_in_ref_amp = str_detect(reference_amplicon, reverse_primer_rc),
    reverse_in_var_amp = str_detect(variant_amplicon, reverse_primer),
    reverse_rc_in_var_amp = str_detect(variant_amplicon, reverse_primer_rc),
    
    # Sense check
    sense_check = (case_when(
      
      # Option 1
      var_in_var_amp == "TRUE" &
      ref_in_ref_amp ==  "TRUE" &
      var_in_ref_amp == "FALSE" &
      ref_in_var_amp == "FALSE" &
      forward_in_ref_amp == "TRUE" &
      forward_rc_in_ref_amp == "FALSE" &
      forward_in_var_amp == "TRUE" &
      forward_rc_in_var_amp == "FALSE" &
      reverse_in_ref_amp == "FALSE" &
      reverse_rc_in_ref_amp == "TRUE" &
      reverse_in_var_amp == "FALSE" &
      reverse_rc_in_var_amp == "TRUE"
      ~ "pass",
      
      # Option 2
      var_in_var_amp == "TRUE" &
      ref_in_ref_amp ==  "TRUE" &
      var_in_ref_amp == "FALSE" &
      ref_in_var_amp == "FALSE" &
      forward_in_ref_amp == "FALSE" &
      forward_rc_in_ref_amp == "TRUE" &
      forward_in_var_amp == "FALSE" &
      forward_rc_in_var_amp == "TRUE" &
      reverse_in_ref_amp == "TRUE" &
      reverse_rc_in_ref_amp == "FALSE" &
      reverse_in_var_amp == "TRUE" &
      reverse_rc_in_var_amp == "FALSE"
      ~ "pass",
      
      TRUE ~"fail")))
    

good_assays <- c("ZFXY",
                 "VHL c.583C>T", "F8 c.6046C>T",
                 "FGFR3 c.1138G>A", "OTC c.905A>G",
                 "PMM2 c.691G>A", "MAGED2 c.1426C>T",
                 "HBB c.20A>T")

bad_assays <- c("IDS c.182_189del")

for_rossa <- amplicons_primers_probes %>%
  select(assay_name, annealing_temperature, 
         forward_primer, reverse_primer, 
         reference_probe, variant_probe, fluorophore, 
         modifications, manufacturer, 
         reference_amplicon, variant_amplicon) %>%
  mutate(category = case_when(
    
    assay_name %in% good_assays ~"good",
    assay_name %in% bad_assays ~"bad",
    TRUE ~"unsure")) %>%
  dplyr::rename(ddpcr_annealing_temperature = annealing_temperature)


write.csv(for_rossa, "analysis_outputs/ddpcr_assays.csv",
          row.names = FALSE)

#########################
# Count motifs
#########################

# Specify DNASE1L3 end motifs from Serpas et al PMID: 30593563 in 
# both forward and reverse complement.
dnase1l3_end_motifs <- c("CCCA", "CCAG", "CCTG", "CCAA", "CCCT", "CCAT", 
                         "TGGG", "CTGG", "CAGG", "TTGG", "AGGG", "ATGG") 

amplicons_longer <- for_rossa %>%
  select(assay_name, category, reference_amplicon, variant_amplicon) %>%
  pivot_longer(
    cols = c(reference_amplicon, variant_amplicon),
    names_to = "amplicon_class",
    values_to = "sequence") %>%
    mutate(length = nchar(sequence))
  
  
dnase1l3_motif_count <- c()

for (i in amplicons_longer$sequence) {
  
  motif_count <- sum(str_count(i,dnase1l3_end_motifs))
  
  dnase1l3_motif_count <- c(dnase1l3_motif_count, motif_count)
  
  rm(motif_count)
}

amplicons_count <- cbind(amplicons_longer, dnase1l3_motif_count) %>%
  mutate(bp_per_site = round(length/dnase1l3_motif_count,1))

#########################
# Plots
#########################

ggplot(amplicons_count %>%
         filter(category != "unsure"), aes(x = assay_name, 
                                           y = dnase1l3_motif_count,
                                           colour = category))+
  geom_jitter() +
  ylim(0, 15)


amplicons_count %>%
  filter(amplicon_class == "reference_amplicon") %>%
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

