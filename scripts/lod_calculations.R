################################################################################
## HbAS ddPCR limit of detection experiment
## September 2021
## Joseph.Shaw@gosh.nhs.uk
################################################################################

library(tidyverse)

# Worksheet 21-3019. HbAS = 150 copies per ul for each allele.
# 300 copies per microlitre for all DNA.
# 300 * 22 = 6600

HbAS_conc <- (150*2) * 22
HbSS_conc <- 2000
HbAA_conc <- 2000
replicates <- 3

# Normalise AA and SS to 2000 GE per ul


# 4 spike-in total DNA levels
# 6 spike-in levels
# 3 replicates each
# Plus 3 NTC wells
# 75 wells total

# There are 125ul of sheared DNA to play with.

lod_input <- data.frame(
  ge_level = rep(c(2000, 4000, 8000, 16000), each = 7),
  hom_spike = rep(rep(c("null", "AA", "SS"), c(1, 3, 3)), 4),
  hom_dna_percent = rep(c(0, 2, 4, 8, 2, 4, 8), 4),
  het_dna_percent = rep(c(100, 98, 96, 92, 98, 96, 92), 4))

lod_calculations <- lod_input %>%
  mutate(
    het_dna_single = (het_dna_percent/100) * ge_level,
    hom_dna_single = (hom_dna_percent/100) * ge_level,
    het_dna_total = het_dna_single * (replicates+0.5),
    hom_dna_total = hom_dna_single * (replicates+0.5),
    het_dna_ul_input = round(het_dna_total / HbAS_conc, 1),
    hom_dna_ul_input = round(hom_dna_total / 2000, 1))

sum(lod_calculations$het_dna_ul_input)

lod_mastermixes <- lod_calculations %>%
  mutate(
    mastermix = paste(ge_level, hom_dna_percent, hom_spike, sep = "_"),
    supermix = 11*(replicates+0.5),
    taqman_assay = (replicates+0.5),
    water = (22*(replicates+0.5))-(supermix + taqman_assay +
                                          het_dna_ul_input +
                                          hom_dna_ul_input)) %>%
  select(mastermix, het_dna_ul_input, hom_dna_ul_input,
         supermix, taqman_assay, water) %>%
  pivot_longer(cols = -c(mastermix), 
               names_to = "reagent",
               values_to = "volume")


# 12ul input into each mastermix