################################################################################
## HbAS ddPCR limit of detection experiment
## September 2021
## Joseph.Shaw@gosh.nhs.uk
################################################################################

library(tidyverse)

#########################
# Preparing the worksheet
#########################

# The structure of this experiment is to have three tube levels:
#  1) The stock tube of sonicated genomic DNA
#  2) An "intermediate" tube of diluted genomic DNA at a concentration
#     to make pipetting easiest.
#  3) The ddPCR mastermix.

# Specify the genomic equivalents (ge) levels and spike percents.

lod_input <- data.frame(
  ge_level = rep(c(2000, 4000, 8000, 12000), each = 7),
  hom_dna_genotype = rep(rep(c("0", "AA", "SS"), c(1, 3, 3)), 4),
  hom_dna_percent = rep(c(0, 2, 4, 8, 2, 4, 8), 4),
  het_dna_percent = rep(c(100, 98, 96, 92, 98, 96, 92), 4))

# Determine intermediate tubes (level 2 above)
# The gDNA samples were sonicated to approximately 150bp fragments, 
# and quantified with ddPCR (21-3019)
# 19RG-220G0190 HbAS –105ul, 300 copies/ul
# 19RG-220G0193 HbAA – 112ul, 294 copies/ul
# 19RG-220G0191 HbSS – 110ul, 290 copies/ul

int_dilutions <- data.frame(
  episode = rep(c("219RG-220G0190", "19RG-220G0193", 
                  "19RG-220G0191"), each = 4),
  genotype = rep(c("AS", "AA", "SS"), each = 4),
  c1 = rep(c((300*22), (294*22), (290*22)), each = 4),
  v2 = rep(c(75, 300, 300), each = 4),
  c2= c(640, 1280, 2560, 3840, 56, 112, 224, 336, 56, 112, 224, 336)) %>%
  mutate(dna = round((v2*c2) / c1, 1),
         water_to_add = v2 - dna,
         tube_name = paste(genotype, c2, sep = "_"))

write.csv(
  int_dilutions,
  "W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/ddPCR Worksheets/2021/21-3233/dilutions_21_3233.csv",
  row.names = FALSE)

# Number of ddPCR replicates for each mastermix
replicates <- 3

# Calculate the number of molecules needed for each mastermix

lod_calculations <- lod_input %>%
  mutate(
    # The number of molecules to go into a single ddPCR well
    het_dna_single_well = ge_level,
    hom_dna_single_well = (ge_level/het_dna_percent)*hom_dna_percent,
    
    # The number of molecules to go into each mastermix (for 3 ddPCR 
    # replicates)
    het_dna_mastermix = het_dna_single_well * (replicates+0.2),
    hom_dna_mastermix = hom_dna_single_well * (replicates+0.2),
    
    # The concentration of DNA from a stock tube to be pipetted.
    # I chose these concentrations to keep the pipetting with manageable
    # volumes and minimise pipetting error (2.3ul to 9.9ul).
    het_dil_tube_cpm = case_when(
      ge_level == 2000 ~640,
      ge_level == 4000 ~1280,
      ge_level == 8000 ~2560,
      ge_level == 12000 ~3840),
    hom_dil_tube_cpm = case_when(
      ge_level == 2000 ~56,
      ge_level == 4000 ~112,
      ge_level == 8000 ~224,
      ge_level == 12000 ~336),
    
    # The volume of each DNA to pipette into the mastermix
    het_dna_ul = het_dna_mastermix/het_dil_tube_cpm,
    hom_dna_ul = round(hom_dna_mastermix/hom_dil_tube_cpm, 1),
    
    # The name of the mastermix
    mastermix = paste(ge_level, hom_dna_percent, hom_dna_genotype, sep = "_"),
    
    # The volumes of other reagents
    supermix = 11*(replicates+0.2),
    taqman_assay = (replicates+0.2),
    water = (22*(replicates+0.2))-(supermix + taqman_assay +
                                     het_dna_ul +
                                     hom_dna_ul),
    total = supermix + taqman_assay + water + het_dna_ul +
      hom_dna_ul,
    
    # Specify which intermediate tube to use
    hom_dna = paste(hom_dna_genotype, hom_dil_tube_cpm, sep = "_"),
    het_dna = paste("AS", het_dil_tube_cpm, sep = "_"),
    
    # Volume of homozygous DNA to add.
    hom_dna_mastermix = round(hom_dna_mastermix)) %>%
  
  # Arrange in a logical way for use in the lab
  select(ge_level, het_dna_percent, hom_dna_percent, 
         hom_dna_genotype, het_dna_single_well, hom_dna_single_well,
         het_dna_mastermix, hom_dna_mastermix, het_dil_tube_cpm,
         hom_dil_tube_cpm, 
         mastermix, het_dna, het_dna_ul, hom_dna, hom_dna_ul,
         supermix, taqman_assay, water, total)

write.csv(
  lod_calculations,
  "W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/ddPCR Worksheets/2021/21-3233/mastermixes_21_3233.csv",
  row.names = FALSE)

#########################
# Analysing the results
#########################

lod_data <- read_csv("data/21-3233.csv", col_names = TRUE) %>%
  janitor::clean_names() %>%
  dplyr::rename(droplets = accepted_droplets) %>%
  filter(substr(well, 1, 1) != "M" & !(sample == "NTC")) %>%
  mutate(sample = factor(sample, levels = 
                           c("SS 8%", "SS 4%", "SS 2%",
                             "0%",
                             "AA 2%", "AA 4%", "AA 8%")))

lod_data_concentration <- lod_data %>%
  pivot_wider(id_cols = c(well, sample),
              names_from = target,
              values_from = concentration,
              names_glue = "{target}_{.value}") %>%
  # Try using the concentration values and multiplying by 22
  mutate(HbS_total = HbS_concentration *22,
         HbA_total = HbA_concentration *22,
         vf_assay_molecules = HbS_total + HbA_total,
         variant_percent = (HbS_total/vf_assay_molecules) *100)

#########################
# Plotting the results
#########################

# This bit isn't finished yet

ggplot(lod_data_concentration, 
       aes(x = vf_assay_molecules, 
           y = variant_percent)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "none", 
        plot.title = element_text(size = 11),
        legend.title = element_blank(),
        legend.text = element_text(size = 9)) +
  ylim(42, 58)+
  xlim(0, 30000) +
  scale_shape_manual(values = c(24, 24, 24, 21, 
                                25,25, 25)) +
  geom_point(size = 2, aes(shape = sample)) +
  
  labs(x = "",
       y = "Variant fraction (%)",
       title = "ddPCR for gDNA mixtures") +
  upper_line +
  lower_line


upper_line <-   geom_function(fun = "calc_SS_boundary",
                              aes(x = vf_assay_molecules, y =
                                    calc_SS_boundary(vf_assay_molecules, ff, 
                                                     lr)),
                              colour = "black",
                              args = c(ff, lr))

lower_line <- geom_function(fun = "calc_AA_boundary",
                            aes(x = vf_assay_molecules, y =
                                  calc_AA_boundary(vf_assay_molecules, ff, 
                                                   lr)),
                            colour = "black",
                            args = c(ff, lr))


ff <- 0.04
lr <- 100


