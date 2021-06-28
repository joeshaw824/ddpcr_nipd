#############################################################
## Functions for ddPCR NIPD
## March 2021
## Joseph.Shaw@gosh.nhs.uk
#############################################################

#############################################################
# Load libraries
#############################################################

# Load packages
library(tidyverse)

#############################################################
# ddPCR functions
#############################################################

# Function for Poisson correction. See Barrett et al 2012 supplementary information (PMID: 22451622)

Poisson_correct <- function(N, P) {
  num_molecules <- as.integer(-log((N-P)/N)*N)
  return(num_molecules)}

# Calculate the fetal fraction from paternal allele copies.

calc_ff <- function(Maternal_copies, Paternal_copies) {
  fetal_fraction <- (Paternal_copies*2) / (Paternal_copies + Maternal_copies)
  return(fetal_fraction)
}

# Calculates the likelihood ratio for a ddPCR test with X-linked inheritance.

calc_LR_X_linked <- function(fetal_fraction, overrep_fraction, total_copies) {
  q0 <- (1 - fetal_fraction) / (2 - fetal_fraction)
  q1 <- 1 / (2 - fetal_fraction)
  Delta <- (1- q1)/(1-q0)
  Gamma <- ((q1 * (1-q0))/ (q0*(1-q1)))
  LR <- exp((((overrep_fraction*log(Gamma)) + log(Delta))*total_copies))
  return(LR)
}

# Calculates the likelihood ratio for a ddPCR test when the variant is on an autosome (recessive or dominant).

calc_LR_autosomal <- function(fetal_fraction, overrep_fraction, total_copies) {
  q0 = 0.5
  # I modified the q1 expression to make it easier to use. Fetal fraction should be in the right format
  # I.e. 0.05 not 5
  q1 <- 0.5+(fetal_fraction/2)
  Delta <- (1- q1)/(1-q0)
  Gamma <- ((q1 * (1-q0))/ (q0*(1-q1)))
  LR <- exp((((overrep_fraction*log(Gamma)) + log(Delta))*total_copies))
  return(LR)
}

# Calculates the 95% Poisson confidence intervals.
# Poisson intervals are lamda +/- 1.96 (sqrt(lamda/N)), when lambda is copies per droplet and N is total droplets.

Poisson_max <- function(Copies_per_droplet, AcceptedDroplets) {
  Cpd_poisson_max = (Copies_per_droplet + 1.96*(sqrt(Copies_per_droplet / AcceptedDroplets))) * AcceptedDroplets
  return(Cpd_poisson_max)
}

Poisson_min <- function(Copies_per_droplet, AcceptedDroplets) {
  Cpd_poisson_min = (Copies_per_droplet - 1.96*(sqrt(Copies_per_droplet / AcceptedDroplets))) * AcceptedDroplets
  return(Cpd_poisson_min)
}

# Calculates the 95% Poisson fractional abundances based on the 95% max and min values for 
# one allele, but not the other allele.

Poisson_fraction_max <- function(copies_allele_max, copies_allele_other) {
  Allele_fraction_max = copies_allele_max / (copies_allele_max + copies_allele_other)
  return(Allele_fraction_max)
}

Poisson_fraction_min <- function(copies_allele_min, copies_allele_other) {
  Allele_fraction_min = copies_allele_min / (copies_allele_min + copies_allele_other)
  return(Allele_fraction_min)
}

#############################################################
# ggplot2 plot theme
#############################################################

# Theme for plots
ddPCR_plot_theme <- theme_bw()+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18), 
        plot.title = element_text(size=20),
        legend.position = c(0.85, 0.85), legend.background = element_rect(fill="white"), 
        legend.title = element_text(size= 18), legend.text = element_text(size= 15))

#Add a horizontal line to show 50%
fifty_percent_line <- geom_hline(yintercept=50, linetype="dashed", size = 1)







