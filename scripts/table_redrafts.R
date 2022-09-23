################################################################################
## Tables for Paper Redraft
################################################################################

# The original analysis was run on CGEN-D0W15R2 with:
# R v4.1.0
# Cmdstanr v0.4.0
# Cmdstan v2.26.1

# For table redrafts I am going to use the original results generated on 15/11/2021.

library(tidyverse)
library(readxl)

setwd("W:/MolecularGenetics/NIPD translational data/NIPD Droplet Digital PCR/ddPCR_R_Analysis/ddpcr_nipd")

supplementary_table <- read_csv("analysis_outputs/supplementary_table20211124_142542.csv")
  
################################################################################