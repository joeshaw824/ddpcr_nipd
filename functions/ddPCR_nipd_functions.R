################################################################################
## Functions for ddPCR NIPD
## October 2021
## Joseph.Shaw@gosh.nhs.uk
################################################################################

#########################
# Load packages
#########################

library(tidyverse)

#########################
# ddPCR functions
#########################

# Function for Poisson correction. See Barrett et al 2012 
# supplementary information (PMID: 22451622)

poisson_correct <- function(N, P) {
  num_molecules <- as.integer(-log((N-P)/N)*N)
  return(num_molecules)}

# Function for reversing the Poisson correction to calculate
# number of positive partitions from number of molecules.
# See sprt_calculations.rmd for full explanation.

reverse_poisson <- function(num_molecules, N) {
  
  P <- as.integer(N - (N / exp(num_molecules/N)))
  
  return(P)
}

# Calculate the fetal fraction from paternal allele copies.

calc_ff <- function(maternal_copies, paternal_copies) {
  fetal_fraction <- (paternal_copies*2) / (paternal_copies + maternal_copies)
  return(fetal_fraction)
}

# Calculates the 95% Poisson confidence intervals.
# Poisson intervals are lamda +/- 1.96 (sqrt(lamda/N)), 
# when lambda is copies per droplet (cpd) and N is total droplets.

poisson_max <- function(copies_per_droplet, droplets) {
  copies_poisson_max = (copies_per_droplet + 
                       1.96*(sqrt(copies_per_droplet / droplets))) * droplets
  return(copies_poisson_max)
}

poisson_min <- function(copies_per_droplet, droplets) {
  copies_poisson_min = (copies_per_droplet - 
                       1.96*(sqrt(copies_per_droplet / droplets))) * droplets
  return(copies_poisson_min)
}

# Calculates the 95% Poisson fractional abundances based on the 95% max 
# and min values for one allele, but not the other allele.

poisson_fraction_max <- function(copies_allele_max, copies_allele_other) {
  allele_fraction_max = copies_allele_max / 
    (copies_allele_max + copies_allele_other)
  return(allele_fraction_max)
}

poisson_fraction_min <- function(copies_allele_min, copies_allele_other) {
  allele_fraction_min = copies_allele_min / 
    (copies_allele_min + copies_allele_other)
  return(allele_fraction_min)
}

#########################
# SPRT functions
#########################

# These functions all should use inputs in decimal format.
# I.e. fetal fraction formatted as "0.04" rather than "4%"

# Smaller functions to prevent repetition of code:

calc_q0_x_linked <- function(ff) {
  q0 <-   (1 - ff) / (2 - ff)
  return(q0)
}

calc_q1_autosomal <- function(ff) {
  # I modified the q1 expression to make it easier to use. 
  # Fetal fraction should be in the right format
  # I.e. 0.05 not 5
  q1 <-  0.5+(ff/2)
  return(q1)
}

calc_q1_x_linked <- function(ff) {
  q1 <-   1 / (2 - ff)
  return(q1)
}

calc_delta <- function(q0, q1) {
  delta <- (1- q1)/(1-q0)
  return(delta)
}

calc_gamma <- function(q0, q1) {
  gamma <- ((q1 * (1-q0))/ (q0*(1-q1)))
  return(gamma)
}

# Calculates the likelihood ratio (lr) for a ddPCR test with 
# X-linked inheritance.

calc_lr_x_linked <- function(ff, overrep_fraction, total_copies) {
  q0 <- calc_q0_x_linked(ff)
  q1 <- calc_q1_x_linked(ff)
  delta <- calc_delta(q0, q1)
  gamma <- calc_gamma(q0, q1)
  lr <- exp((((overrep_fraction*log(gamma)) + log(delta))*total_copies))
  return(lr)
}

# Calculates the likelihood ratio for a ddPCR test when the variant is 
# on an autosome (recessive or dominant).

calc_lr_autosomal <- function(ff, overrep_fraction, total_copies) {
  q0 = 0.5
  q1 <- calc_q1_autosomal(ff)
  delta <- calc_delta(q0, q1)
  gamma <- calc_gamma(q0, q1)
  lr <- exp((((overrep_fraction*log(gamma)) + log(delta))*total_copies))
  return(lr)
}

# These functions calculate the SPRT thresholds with likelihood ratio 
# supplied and fetal fraction supplied as a decimal.
calc_hom_var_boundary <- function(total_copies, ff, lr) {
  q0 <- 0.5
  q1 <- calc_q1_autosomal(ff)
  delta <- calc_delta(q0, q1)
  gamma <- calc_gamma(q0, q1)
  hom_var_boundary <- ((log(lr)/total_copies) - log(delta))/log(gamma)*100
  # Convert to a percentage for output
  return(hom_var_boundary)
}

calc_het_upper_boundary <- function(total_copies, ff, lr) {
  q0 <- 0.5
  q1 <- calc_q1_autosomal(ff)
  delta <- calc_delta(q0, q1)
  gamma <- calc_gamma(q0, q1)
  het_upper_boundary <- ((log(1/lr)/total_copies) - log(delta))/log(gamma)*100
  # Convert to a percentage for output
  return(het_upper_boundary)
}

calc_het_lower_boundary <- function(total_copies, ff, lr) {
  het_upper_boundary <- calc_het_upper_boundary(total_copies, ff, lr)
  het_lower_boundary <- 50-(het_upper_boundary-50)
  return(het_lower_boundary)
}

calc_hom_ref_boundary <- function(total_copies, ff, lr) {
  hom_var_boundary <- calc_hom_var_boundary(total_copies, ff, lr)
  hom_ref_boundary <- 50-(hom_var_boundary-50)
  # Convert to a percentage for output
  return(hom_ref_boundary)
}

calc_hemi_var_boundary <- function(total_copies, ff, lr) {
  q0 <- calc_q0_x_linked(ff)
  q1 <- calc_q1_x_linked(ff)
  delta <- calc_delta(q0, q1)
  gamma <- calc_gamma(q0, q1)
  hemi_var_boundary <- (((log(lr)/total_copies) - log(delta))/log(gamma))*100
  return(hemi_var_boundary)
}

calc_hemi_ref_boundary <- function(total_copies, ff, lr) {
  hemi_var_boundary <- calc_hemi_var_boundary(total_copies, ff, lr)
  hemi_ref_boundary <- 50 -(hemi_var_boundary -50)
  return(hemi_ref_boundary)
}

#########################
# MCMC functions
#########################

# This function runs Tristan's MCMC pipeline for 3 inheritance patterns.

# n_K	= number of droplets tested for variant assay
# K_M	= number of droplets positive for variant (mutant) allele
# K_N	= number of droplets positive for normal (reference) allele
# n_Z	= number of droplets tested for fetal fraction assay
# Z_X	= number of droplets positive for maternal homozygous allele
# Z_Y	= number of droplets positive for paternal allele

# Autosomal dominant:
# p_G1: probability fetus is heterozygous
# p_G2: probability fetus is homozygous reference

# Autosomal recessive:
# p_G1: probability fetus is homozygous reference
# p_G2: probability fetus is heterozygous
# p_G3: probability fetus is homozygous variant

# X-linked:
# p_G0: probability fetus is hemizygous reference
# p_G1: probability fetus is hemizygous variant

# Threshold should be 0.95

run_mcmc <- function(data_input, threshold) {
  
  # Input data must have the required columns
  stopifnot(c("inheritance_chromosomal", "inheritance_pattern",
            "n_K", "K_M", "K_N", "n_Z", "Z_X", "Z_Y") %in% colnames(data_input))
  
  # Compile the models
  dominant_model <- cmdstan_model("models/nipt_dominant.stan")

  x_linked_model <- cmdstan_model("models/nipt_x_linked.stan")
  
  recessive_model <- cmdstan_model("models/nipt_recessive.stan")
  
  # Initialise the chains

  initialise_chains_dominant <- function() list(rho = runif(1, 0.1, 0.5), 
                                              M_K = runif(1, 0.1, 0.5), 
                                              M_Z = runif(1, 0.1, 0.5))

  initialise_chains_xlinked <- function() list(rho = rbeta(1, 4, 32),
                                             M_K = abs(rnorm(1, sd = 0.05)),
                                             M_Z = abs(rnorm(1, sd = 0.05)))

  initialise_chains_recessive <- function() list(rho = runif(1, 0.1, 0.5),
                                               M_K = runif(1, 0.1, 0.5),
                                               M_Z = runif(1, 0.1, 0.5))
  
  # Set probability threshold for accepting fetal genotype predictions
  mcmc_threshold <- threshold
  
  # Generate the probabilities for the fetal genotype based on 
  # the inheritance pattern of the variant.
  
  data_with_fits <- data_input %>%
    nest(data = n_K:Z_Y) %>%
    mutate(
      
      data = map(data, as.list),
      
      fit = case_when(
        inheritance_chromosomal == "autosomal" &
        inheritance_pattern == "dominant" ~map(
          data, ~ dominant_model$sample(data = .,
                                        init = initialise_chains_dominant,
                                        step_size = 0.2,
                                        parallel_chains = parallel::detectCores())),
      
      inheritance_chromosomal == "autosomal" & 
        inheritance_pattern == "recessive" ~ map(
          data, ~ recessive_model$sample(data = .,
                                         init = initialise_chains_recessive,
                                         step_size = 0.2,
                                         parallel_chains = parallel::detectCores())),
      
      inheritance_chromosomal == "x_linked" ~map(
        data, ~ x_linked_model$sample(data = .,
                                      init = initialise_chains_xlinked,
                                      step_size = 0.2,
                                      parallel_chains = parallel::detectCores()))),
    
    results = case_when(
      inheritance_chromosomal == "autosomal" &
        inheritance_pattern == "dominant" ~map(
          fit,  ~ setNames(pivot_wider(.$summary(c("pG", "rho"), "mean"),
                                       names_from = "variable",
                                       values_from = "mean"),
                           c("p_G1", "p_G2", "rho_est"))),
      
      inheritance_chromosomal == "autosomal" & 
        inheritance_pattern == "recessive" ~map(
          fit,  ~ setNames(pivot_wider(.$summary(c("pG", "rho"), "mean"),
                                       names_from = "variable",
                                       values_from = "mean"),
                           c("p_G1", "p_G2", "p_G3", "rho_est"))),
      
      inheritance_chromosomal == "x_linked" ~map(
        fit,  ~ setNames(pivot_wider(.$summary(c("pG", "rho"), "mean"),
                                     names_from = "variable",
                                     values_from = "mean"),
                         c("p_G0", "p_G1", "rho_est"))))) %>%
  
  unnest_wider(results)
  
  # Add on fetal genotype predictions based on the previously set threshold.
  
  data_with_predictions <- data_with_fits %>%
    select(-c(data, fit)) %>%
    rename(fetal_fraction = rho_est) %>%
    mutate(mcmc_prediction = case_when(
      
      # Dominant predictions
      inheritance_chromosomal == "autosomal" & 
        inheritance_pattern == "dominant" & 
        p_G1 > mcmc_threshold ~"heterozygous",
      
      inheritance_chromosomal == "autosomal" &
        inheritance_pattern == "dominant" & 
        p_G2 > mcmc_threshold ~"homozygous reference",
      
      inheritance_chromosomal == "autosomal" &
        inheritance_pattern == "dominant" & 
        p_G1 < mcmc_threshold & 
        p_G2 < mcmc_threshold ~"inconclusive",
      
      # Recessive predictions
      inheritance_chromosomal == "autosomal" &
        inheritance_pattern == "recessive" & 
        p_G1 > mcmc_threshold ~"homozygous reference",
      
      inheritance_chromosomal == "autosomal" &
        inheritance_pattern == "recessive" &
        p_G2 > mcmc_threshold ~"heterozygous",
      
      inheritance_chromosomal == "autosomal" &
        inheritance_pattern == "recessive" &
        p_G3 > mcmc_threshold ~"homozygous variant",
      
      inheritance_chromosomal == "autosomal" &
        inheritance_pattern == "recessive" & 
        p_G1 < mcmc_threshold & 
        p_G2 < mcmc_threshold &
        p_G3 < mcmc_threshold~"inconclusive",
      
      # X-linked predictions
      inheritance_chromosomal == "x_linked" &
        p_G0 > mcmc_threshold ~"hemizygous reference",
      
      inheritance_chromosomal == "x_linked" &
        p_G1 > mcmc_threshold ~"hemizygous variant",
      
      inheritance_chromosomal == "x_linked" &
        p_G0 < mcmc_threshold &
        p_G1 < mcmc_threshold ~"inconclusive"))

  return(data_with_predictions)
  
}

#########################
# Grouped functions
#########################

# These functions calculate the number of molecules and fractional
# abundance for each target, including 95% Poisson confidence intervals.
# This is to prevent massive duplication of code when performing these
# calculations for gDNA, cfDNA and limit of detection datasets.

var_ref_calculations <- function(data_input) {
  
  # The data_input needs these variables for this function to work:
  # vf_assay_droplets, variant_positives, reference_positives
  
  stopifnot(c("vf_assay_droplets", "variant_positives", 
              "reference_positives")
            %in% colnames(data_input))
  
  data_output <- data_input %>%
    
    mutate(
      # Perform Poisson correction to determine the total number of molecules 
      # detected for each target.
      variant_molecules = poisson_correct(vf_assay_droplets,
                                          variant_positives), 
      reference_molecules = poisson_correct(vf_assay_droplets,
                                            reference_positives),   
      
      # Calculate the total molecules for the variant fraction (vf) assay
      vf_assay_molecules = variant_molecules + reference_molecules,
      
      # Calculate the fractional abundance of each allele
      reference_fraction = reference_molecules / vf_assay_molecules,
      variant_fraction = variant_molecules / vf_assay_molecules,
      variant_percent  = variant_fraction * 100,
      reference_percent  = reference_fraction * 100,
      
      # Select the major and minor alleles
      major_allele = case_when(
        # If both fractions are exactly 50, then the variant is arbitrarily 
        # chosen as the major allele
        variant_percent >= reference_percent ~ "variant allele",
        reference_percent > variant_percent ~ "reference allele"),
      
      minor_allele = case_when(
        major_allele == "variant allele" ~ "reference allele",
        major_allele == "reference allele" ~ "variant allele"),
      
      # Choose the percentage of the major and minor alleles
      major_allele_percent = case_when(
        major_allele == "variant allele" ~variant_percent,
        major_allele == "reference allele" ~reference_percent),
      
      # Choose the molecules for each allele
      major_allele_molecules = case_when(
        major_allele == "variant allele" ~variant_molecules,
        major_allele == "reference allele" ~reference_molecules),
      
      minor_allele_molecules = case_when(
        minor_allele == "variant allele" ~variant_molecules,
        minor_allele == "reference allele" ~reference_molecules),
      
      # Find the difference between the numbers of molecules for each allele
      difference_molecules = abs(major_allele_molecules - 
                                   minor_allele_molecules),
      
      # Calculate the 95% confidence intervals of the numbers of molecules
      # for each allele
      variant_molecules_max = poisson_max((variant_molecules/ 
                                             vf_assay_droplets), 
                                          vf_assay_droplets),
      
      variant_molecules_min = poisson_min((variant_molecules/
                                             vf_assay_droplets),
                                          vf_assay_droplets),
      
      reference_molecules_max = poisson_max((reference_molecules/
                                               vf_assay_droplets),
                                            vf_assay_droplets),
      
      reference_molecules_min = poisson_min((reference_molecules/
                                               vf_assay_droplets), 
                                            vf_assay_droplets),
      
      major_allele_molecules_max = case_when(
        major_allele == "variant allele" ~variant_molecules_max,
        major_allele == "reference allele" ~reference_molecules_max),
      
      major_allele_molecules_min = case_when(
        major_allele == "variant allele" ~variant_molecules_min,
        major_allele == "reference allele" ~reference_molecules_min),
      
      minor_allele_molecules_max = case_when(
        minor_allele == "variant allele" ~variant_molecules_max,
        minor_allele == "reference allele" ~reference_molecules_max),
      
      minor_allele_molecules_min = case_when(
        minor_allele == "variant allele" ~variant_molecules_min,
        minor_allele == "reference allele" ~reference_molecules_min),
      
      difference_molecules_max = major_allele_molecules_max - 
        minor_allele_molecules_min,
      
      difference_molecules_min = major_allele_molecules_min - 
        minor_allele_molecules_max,
      
      # Calculate the 95% confidence intervals of the fractional abundances
      variant_percent_max = (poisson_fraction_max(
        variant_molecules_max, reference_molecules))*100,
      
      variant_percent_min = (poisson_fraction_min(
        variant_molecules_min, reference_molecules))*100,
      
      reference_percent_max = (poisson_fraction_max(
        reference_molecules_max, variant_molecules))*100,
      
      reference_percent_min = (poisson_fraction_min(
        reference_molecules_min, variant_molecules))*100,
      
      major_allele_percent_max = case_when(
        major_allele == "variant allele" ~variant_percent_max,
        major_allele == "reference allele" ~reference_percent_max),
      
      major_allele_percent_min = case_when(
        major_allele == "variant allele" ~variant_percent_min,
        major_allele == "reference allele" ~reference_percent_min),
      
      # Calculate the 95% confidence intervals of the numbers of molecules 
      # detected by each assay
      vf_assay_molecules_max = poisson_max((
        vf_assay_molecules/vf_assay_droplets), 
        vf_assay_droplets),
      
      vf_assay_molecules_min = poisson_min((
        vf_assay_molecules/vf_assay_droplets),
        vf_assay_droplets))
  
  return(data_output)
}

ff_calculations <- function(data_input) {
  
  # The data_input needs these variables for this function to work:
  # maternal_positives, paternal_positives, ff_assay_droplets
  
  stopifnot(c("maternal_positives", "paternal_positives",
              "ff_assay_droplets") %in% colnames(data_input))
  
  data_output <- data_input %>%
    mutate(
      
      # Calculate the same variables for the fetal fraction assay
      
      maternal_molecules = poisson_correct(
        ff_assay_droplets,maternal_positives),
      
      paternal_molecules = poisson_correct(
        ff_assay_droplets,paternal_positives),
      
      # Calculate the fetal fraction
      fetal_fraction = calc_ff(maternal_molecules, paternal_molecules),
      
      fetal_percent = fetal_fraction*100,
      
      ff_assay_molecules = maternal_molecules + paternal_molecules,
      
      paternal_molecules_max = poisson_max((
        paternal_molecules / ff_assay_droplets), 
        ff_assay_droplets),
      
      paternal_molecules_min = poisson_min((
        paternal_molecules / ff_assay_droplets), 
        ff_assay_droplets),
      
      maternal_molecules_max = poisson_max((
        maternal_molecules / ff_assay_droplets), 
        ff_assay_droplets),
      
      maternal_molecules_min = poisson_min((
        maternal_molecules / ff_assay_droplets), 
        ff_assay_droplets),
      
      ff_assay_molecules_max = poisson_max((
        ff_assay_molecules / ff_assay_droplets), 
        ff_assay_droplets),
      
      ff_assay_molecules_min = poisson_min((
        ff_assay_molecules / ff_assay_droplets), 
        ff_assay_droplets),
      
      # When calculating for the fetal fraction, the paternal fraction 
      # must be multiplied by 2.
      fetal_percent_max = 200* (poisson_fraction_max(
        paternal_molecules_max, maternal_molecules)),
      
      fetal_percent_min = 200* (poisson_fraction_min(
        paternal_molecules_min, maternal_molecules)))
  
  return(data_output)
}

#########################
# RMD plot function
#########################

# This is the function for plotting relative mutation dosage (RMD) 
# results for ddPCR, including parental gDNA controls. The amount of 
# information on this plot can be modified to suit user preference.

# Maternal sample only example
# draw_rmd_plot("30065", "21-1863.csv_M07", "21-1862.csv_M10")

# Maternal and paternal samples example
# draw_rmd_plot("20238", 
              # c("21-1413.csv_M02", "21-1413.csv_M08"), 
              # c("21-1412.csv_M08",	"21-1412.csv_M10"))

draw_rmd_plot <- function(cfdna_sample, parent_vf_wells, parent_ff_wells) {
  
  # Get cfDNA data
  cfDNA_rmd <- ddpcr_sprt_unblinded %>%
    filter(r_number %in% cfdna_sample) %>%
    dplyr::rename(sample = r_number) %>%
    mutate(identity = "cfDNA") %>%
    select(sample, identity, sprt_prediction, fetal_percent,
           variant_percent, 
           vf_assay, ff_assay, mutation_genetic_info_fetus,
           variant_molecules, reference_molecules,
           variant_molecules_max, variant_molecules_min, 
           reference_molecules_max, reference_molecules_min, 
           maternal_molecules, paternal_molecules,
           paternal_molecules_max, paternal_molecules_min,
           maternal_molecules_max, maternal_molecules_min)
  
  # Get parent data and cfDNA data for a single NIPT case and pivot into
  # the format for plotting
  parents_rmd <- left_join(
    parent_gDNA_var_ref %>%
      # First table
      filter(worksheet_well_sample %in% parent_vf_wells) %>%
      select(sample, identity, vf_assay, variant_molecules, 
             reference_molecules, variant_molecules_max, variant_molecules_min, 
             reference_molecules_max, reference_molecules_min),
    # Second table
    parent_gDNA_ff %>%
      filter(worksheet_well_sample %in% parent_ff_wells) %>%
      select(sample, ff_assay, maternal_molecules, paternal_molecules,
             paternal_molecules_max, paternal_molecules_min,
             maternal_molecules_max, maternal_molecules_min),
    by = "sample") %>%
    select(sample, identity, vf_assay, ff_assay,
           variant_molecules, reference_molecules,
           variant_molecules_max, variant_molecules_min, 
           reference_molecules_max, reference_molecules_min, 
           maternal_molecules, paternal_molecules,
           paternal_molecules_max, paternal_molecules_min,
           maternal_molecules_max, maternal_molecules_min)
  
  # Bind with cfDNA data
  case_rmd <- rbind(parents_rmd,
                    cfDNA_rmd %>%
                      select(sample, identity, vf_assay, 
                             ff_assay, variant_molecules, 
                             reference_molecules, variant_molecules_max, 
                             variant_molecules_min, reference_molecules_max, 
                             reference_molecules_min, maternal_molecules, 
                             paternal_molecules, paternal_molecules_max, 
                             paternal_molecules_min, maternal_molecules_max, 
                             maternal_molecules_min))
  
  # The fetal fraction and variant assay for all samples in a case 
  # should be the same
  # na.omit used because some cases do not have both assays tested for 
  # parental samples
  stopifnot(length(unique(na.omit(case_rmd$ff_assay)))==1)
  stopifnot(length(unique(na.omit(case_rmd$vf_assay)))==1)
  
  # This section can probably be simplified with a clever pivot
  # to get the max and min values in separate columns
  
  case_data_variant <- case_rmd %>%
    select(sample, identity, variant_molecules,
           variant_molecules_max, variant_molecules_min) %>%
    mutate(target_type = "Variant") %>%
    rename(molecules = variant_molecules,
           molecules_max = variant_molecules_max,
           molecules_min = variant_molecules_min)
  
  case_data_reference <- case_rmd %>%
    select(sample, identity, reference_molecules,
           reference_molecules_max, reference_molecules_min) %>%
    mutate(target_type = "Reference") %>%
    rename(molecules = reference_molecules,
           molecules_max = reference_molecules_max,
           molecules_min = reference_molecules_min)
  
  case_data_maternal <- case_rmd %>%
    select(sample, identity, maternal_molecules,
           maternal_molecules_max, maternal_molecules_min) %>%
    mutate(target_type = "Shared allele") %>%
    rename(molecules = maternal_molecules,
           molecules_max = maternal_molecules_max,
           molecules_min = maternal_molecules_min)
  
  case_data_paternal <- case_rmd %>%
    select(sample, identity, paternal_molecules,
           paternal_molecules_max, paternal_molecules_min) %>%
    mutate(target_type = "Fetal-specific allele") %>%
    rename(molecules = paternal_molecules,
           molecules_max = paternal_molecules_max,
           molecules_min = paternal_molecules_min)
  
  case_data_long <- rbind(case_data_variant,
                          case_data_reference,
                          case_data_maternal,
                          case_data_paternal)%>%
    
    # Control factor order for plot
    mutate(identity = factor(identity, levels = c("maternal gDNA",
                                                  "paternal gDNA",
                                                  "cfDNA")),
           target_type = factor(target_type, levels = c("Shared allele",
                                                        "Fetal-specific allele",
                                                        "Reference",
                                                        "Variant")))
  
  rmd_plot <- ggplot(case_data_long, 
                     aes(x = identity, y = molecules, fill = target_type)) +
    geom_col(position = position_dodge(width = 0.9), 
             colour="black", alpha = 0.6)+
    geom_errorbar(aes(ymin = molecules_min, ymax = molecules_max, 
                      width = 0.3), position = position_dodge(width = 0.9))+
    scale_fill_manual(values = c("#99FFFF", "#FFCC99", "#3366FF", "#FF0000"))+
    theme_bw()+
    theme(axis.text=element_text(size=18), axis.title = element_text(size=18),
          legend.position = "bottom", legend.title = element_blank(),
          legend.text = element_text(size= 14),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    labs(x = "", y = "Molecules", 
         title = paste("cfDNA:", cfdna_sample),
         subtitle = paste("Variant fraction assay:", cfDNA_rmd[1,6], "  ",
                          "Fetal fraction assay:", cfDNA_rmd[1,7],
                          "  ", "SPRT prediction: ",
                          cfDNA_rmd[1,3],
                          "  ", "Fetal fraction: ",
                          round(cfDNA_rmd[1,4], 1),"%  ",
                          "Variant fraction: ",
                          round(cfDNA_rmd[1,5], 1),"%  ",
                          "Invasive result: ",
                          cfDNA_rmd[1,8]))+
    geom_text(aes(x = identity, y = molecules_max, label = molecules), 
              position = position_dodge(width = 0.9), vjust = -1)
  
  return(rmd_plot)
}

#########################
# Sequence functions
#########################

reverse_complement <- function(input_sequence){
  
  rev_comp <- stri_reverse(chartr("ATGC","TACG",input_sequence))
  
  return(rev_comp)
  
}

#########################
# Prediction functions
#########################

# Function to make fetal genotype predictions following SPRT analysis

predict_sprt_genotypes <- function(df, lr_threshold) {
  
  stopifnot(c("inheritance_chromosomal", "inheritance_pattern",
              "likelihood_ratio")
            %in% colnames(df))
  
  predictions <- df %>%
    mutate(
      # Classify based on likelihood ratio threshold supplied
      # Fetal genotype predictions are named consistently as 
      # "inconclusive", "heterozygous", "homozygous/hemizygous reference/variant"
      sprt_prediction = case_when(
        inheritance_chromosomal == "autosomal" &
          likelihood_ratio > lr_threshold &
          major_allele == "reference allele" 
        ~ "homozygous reference",
        
        inheritance_chromosomal == "autosomal" &
          likelihood_ratio > lr_threshold &
          major_allele == "variant allele" 
        ~ "homozygous variant",
        
        inheritance_chromosomal == "autosomal" &
          likelihood_ratio < (1/lr_threshold)
        ~ "heterozygous",
        
        inheritance_chromosomal == "autosomal" &
          likelihood_ratio < lr_threshold &
          likelihood_ratio > (1/lr_threshold) 
        ~ "inconclusive",
        
        inheritance_chromosomal == "x_linked" &
          likelihood_ratio > lr_threshold &
          major_allele == "reference allele" 
        ~ "hemizygous reference",
        
        inheritance_chromosomal == "x_linked" &
          likelihood_ratio > lr_threshold &
          major_allele == "variant allele" 
        ~ "hemizygous variant",
        
        inheritance_chromosomal == "x_linked" &
          likelihood_ratio < lr_threshold 
        ~ "inconclusive"),
      
      # Factorise for plotting
      sprt_prediction = factor(sprt_prediction, levels = 
                                 c("hemizygous variant",
                                   "homozygous variant",
                                   "heterozygous",
                                   "homozygous reference",
                                   "hemizygous reference",
                                   "inconclusive")))
  
  return(predictions)
  
}

# This function converts fetal genotype predictions for multiple inheritance 
# patterns into binary "positive" and "negative" classifiers, to help
# with sensitivity calculations.

binary_predictions <- function(df, prediction) {
  
  stopifnot(c("inheritance_chromosomal", "inheritance_pattern")
            %in% colnames(df))
  
  new_df <- df %>%
    mutate(binary_call = case_when(
      # Autosomal dominant inheritance
      inheritance_chromosomal == "autosomal" &
        inheritance_pattern == "dominant" &
        !!prediction == "heterozygous" ~"positive",
      inheritance_chromosomal == "autosomal" &
        inheritance_pattern == "dominant" &
        !!prediction == "homozygous reference" ~"negative",
      
      # Autosomal recessive inheritance
      inheritance_chromosomal == "autosomal" &
        inheritance_pattern == "recessive" &
        !!prediction %in% c("homozygous variant", 
                            "homozygous reference")  ~"positive",
      inheritance_chromosomal == "autosomal" &
        inheritance_pattern == "recessive" &
        !!prediction == "heterozygous" ~"negative",
      
      # X linked inheritance
      inheritance_chromosomal == "x_linked" &
        !!prediction == "hemizygous variant"  ~"positive",
      inheritance_chromosomal == "x_linked" &
        !!prediction == "hemizygous reference"  ~"negative",
      TRUE ~"inconclusive"))
  
  return(new_df)
}

sensitivity_metrics <- function(df, prediction_binary, outcome, cohort_input,
                          cohort_name) {
  
  # Example
  # sensitivity_metrics(df = supplementary_table, 
      # prediction_binary = quo(sprt_binary), 
      # outcome = quo(outcome_sprt), 
      # cohort_input = c("sickle cell disease", "bespoke design"),
      # cohort_name = "all") 
  
  # Filter by cohort
  df_cohort <- df %>%
    filter(cohort %in% cohort_input)
  
  true_positives <- nrow(df_cohort %>%
                           filter(!!prediction_binary == "positive" &
                                    !!outcome == "correct"))
  
  true_negatives <- nrow(df_cohort %>%
                           filter(!!prediction_binary == "negative" &
                                    !!outcome == "correct"))
  
  false_positives <- nrow(df_cohort %>%
                            filter(!!prediction_binary %in% "positive" &
                                     !!outcome == "incorrect"))
  
  false_negatives <- nrow(df_cohort %>%
                            filter(!!prediction_binary == "negative" &
                                     !!outcome == "incorrect"))
  
  inconclusives <- nrow(df_cohort %>%
                          filter(!!prediction_binary == "inconclusive"))
  
  # Generate an input table for the epiR::epi.tests function
  # Order must go: TP, FP, FN, TN
  data_table <- as.table(matrix(c(true_positives, false_positives, 
                                  false_negatives, true_negatives), 
                                nrow = 2, byrow = TRUE))
  
  # Use epiR to calculate sensitivity and specificity
  data_results <- epiR::epi.tests(data_table, conf.level = 0.95)
  
  # Sensitivity
  
  sensitivity_paste <- paste0(round(data_results$detail$se$est, 3)*100, 
         " (",
         round(data_results$detail$se$lower, 3)*100,
         "-",
         round(data_results$detail$se$upper, 3)*100,
         ")")
  
  # Specificity
  
  specificity_paste <- paste0(round(data_results$detail$sp$est, 3)*100, 
         " (",
         round(data_results$detail$sp$lower, 3)*100,
         "-",
         round(data_results$detail$sp$upper, 3)*100,
         ")")
  
  output <- data.frame(
    "group" = cohort_name,
    true_positive = c(as.character(true_positives)),
    true_negative = c(as.character(true_negatives)),
    false_positive = c(as.character(false_positives)),
    false_negative = c(as.character(false_negatives)),
    inconclusive = c(as.character(inconclusives)),
    sensitivity = c(sensitivity_paste),
    specificity = c(specificity_paste)) %>%
    
    pivot_longer(cols = c(-group),
                 names_to = "category",
                 values_to = "analysis_method")
  
  return(output)
  
}

#########################