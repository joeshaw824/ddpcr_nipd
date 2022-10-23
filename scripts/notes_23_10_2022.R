## Done
# Reliance on RAPID biobank removed
# Old sickle cell results tabled removed
# Need to streamline supp info columns
# Removed gene information table




# Add in quality filtering step
dplyr::rename(sprt_prediction_pre_qc = sprt_prediction) %>%
  mutate(,
    # Factorise for plotting
    )

# Convert SPRT predictions into a binary format, for assistance with
# sensitivity calculations later on
all_samples_sprt <- binary_predictions(df = all_samples_sprt, 
                                       prediction = quo(sprt_prediction)) %>%
  dplyr::rename(sprt_binary = binary_call)

all_samples_mcmc <- binary_predictions(
  df = all_samples_mcmc,
  prediction = quo(mcmc_prediction)) %>%
  dplyr::rename(mcmc_binary = binary_call)


all_samples_zscore <- binary_predictions(
  df = all_samples_zscore,
  prediction = quo(zscore_prediction)) %>%
  dplyr::rename(zscore_binary = binary_call)

###################
# Sensitivity and specificity table - Table 2
###################

sprt_scd <- sensitivity_metrics(df = supplementary_table, 
                                prediction_binary = quo(sprt_binary), 
                                outcome = quo(sprt_outcome), 
                                cohort_input = "sickle cell disease",
                                cohort_name = "sickle cell disease") %>%
  dplyr::rename(sprt = analysis_method)

sprt_bespoke <- sensitivity_metrics(df = supplementary_table, 
                                    prediction_binary = quo(sprt_binary), 
                                    outcome = quo(sprt_outcome), 
                                    cohort_input = "bespoke design",
                                    cohort_name = "bespoke design") %>%
  dplyr::rename(sprt = analysis_method)

sprt_all <- sensitivity_metrics(df = supplementary_table, 
                                prediction_binary = quo(sprt_binary), 
                                outcome = quo(sprt_outcome), 
                                cohort_input = c("bespoke design", 
                                                 "sickle cell disease"),
                                cohort_name = "all") %>%
  dplyr::rename(sprt = analysis_method)

mcmc_scd <- sensitivity_metrics(df = supplementary_table, 
                                prediction_binary = quo(mcmc_binary), 
                                outcome = quo(mcmc_outcome), 
                                cohort_input = "sickle cell disease",
                                cohort_name = "sickle cell disease") %>%
  dplyr::rename(mcmc = analysis_method)

mcmc_bespoke <- sensitivity_metrics(df = supplementary_table, 
                                    prediction_binary = quo(mcmc_binary), 
                                    outcome = quo(mcmc_outcome), 
                                    cohort_input = "bespoke design",
                                    cohort_name = "bespoke design") %>%
  dplyr::rename(mcmc = analysis_method)

mcmc_all <- sensitivity_metrics(df = supplementary_table, 
                                prediction_binary = quo(mcmc_binary), 
                                outcome = quo(mcmc_outcome), 
                                cohort_input = c("bespoke design", 
                                                 "sickle cell disease"),
                                cohort_name = "all") %>%
  dplyr::rename(mcmc = analysis_method)

zscore_scd <- sensitivity_metrics(df = supplementary_table, 
                                  prediction_binary = quo(zscore_binary), 
                                  outcome = quo(zscore_outcome), 
                                  cohort_input = "sickle cell disease",
                                  cohort_name = "sickle cell disease") %>%
  dplyr::rename(zscore = analysis_method)

zscore_bespoke <- sensitivity_metrics(df = supplementary_table, 
                                      prediction_binary = quo(zscore_binary), 
                                      outcome = quo(zscore_outcome), 
                                      cohort_input = "bespoke design",
                                      cohort_name = "bespoke design") %>%
  dplyr::rename(zscore = analysis_method)

zscore_all <- sensitivity_metrics(df = supplementary_table, 
                                  prediction_binary = quo(zscore_binary), 
                                  outcome = quo(zscore_outcome), 
                                  cohort_input = c("bespoke design", 
                                                   "sickle cell disease"),
                                  cohort_name = "all") %>%
  dplyr::rename(zscore = analysis_method)


# Row  bind analysis methods together

sprt_metrics <- rbind(sprt_scd, sprt_bespoke, sprt_all)
mcmc_metrics <- rbind(mcmc_scd, mcmc_bespoke, mcmc_all)
zscore_metrics <- rbind(zscore_scd, zscore_bespoke, zscore_all)

analysis_metrics <- sprt_metrics %>%
  left_join(mcmc_metrics,
            by = c("group", "category")) %>%
  left_join(zscore_metrics,
            by = c("group", "category")) %>%
  # Tidy up names for paper
  dplyr::rename(
    `SPRT` = sprt,
    `MCMC` = mcmc,
    `Z score` = zscore) %>%
  mutate(`Result` = case_when(
    category == "true_positive" ~ "True positive",
    category == "true_negative" ~ "True negative",
    category == "false_positive" ~ "False positive",
    category == "false_negative" ~ "False negative",
    category == "sensitivity" ~ "Sensitivity (%)",
    category == "specificity" ~ "Specificity (%)",
    category == "inconclusive" ~ "Inconclusive"),
    `Cohort` = case_when(
      group == "sickle cell disease" ~ "Sickle cell disease",
      group == "bespoke design" ~ "Bespoke design",
      group == "all" ~"All samples")) %>%
  select(`Cohort`, `Result`, `SPRT`, `MCMC`, `Z score`)


spe_sen_table <- analysis_metrics %>%
  filter(Cohort == "All samples" &
           Result %in% c("Sensitivity (%)", "Specificity (%)")) %>%
  select(`Result`, `SPRT`, `MCMC`, `Z score`)

write.csv(spe_sen_table,
          file = paste0("analysis_outputs/spe_sen_table",
                        format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
          row.names = FALSE)

###################

###################
# Plot: SPRT analysis
###################

cfdna_sprt_plot_title <- paste0("ddPCR for ",
                                nrow(all_samples_unblinded),
                                " cfDNA samples with SPRT classification")

plot_2a <- ggplot(all_samples_unblinded, aes(x = vf_assay_molecules, 
                                             y = variant_percent)) +
  theme_bw() +
  multiplot_theme + 
  multiplot_y +
  multiplot_x +
  cfdna_fill +
  cfdna_alpha +
  cfdna_shape +
  geom_point(size = 2, aes(fill = sprt_outcome,
                           alpha = sprt_outcome,
                           shape = fetal_genotype),
             colour = "black") +
  
  labs(y = "Variant fraction (%)", x = "", 
       title = cfdna_sprt_plot_title)

###################
# Plot: MCMC analysis
###################

cfdna_mcmc_plot_title <- paste0("ddPCR for ",
                                nrow(all_samples_unblinded),
                                " cfDNA samples with MCMC classification")

plot_2b <- ggplot(all_samples_unblinded, aes(x = vf_assay_molecules, 
                                             y = variant_percent)) +
  theme_bw() +
  multiplot_theme + 
  multiplot_y +
  multiplot_x +
  cfdna_fill +
  cfdna_alpha +
  cfdna_shape +
  geom_point(size = 2, aes(fill = mcmc_outcome,
                           alpha = mcmc_outcome,
                           shape = fetal_genotype),
             colour = "black") +
  
  labs(y = "", x = "", 
       title = cfdna_mcmc_plot_title)

###################
# Plot: heterozygous gDNA controls
###################

het_gdna_plot_title <- paste0("ddPCR for ",
                              length(unique(het_gdna$r_number)),
                              " heterozygous gDNA controls")

plot_2c <- ggplot(het_gdna %>%
                    mutate(sample_type = "het gDNA"), 
                  aes(x = vf_assay_molecules, 
                      y = variant_percent)) +
  z3_line +
  zminus3_line +
  z2_line +
  zminus2_line +
  vertical_line +
  multiplot_y +
  scale_x_continuous(limits = c(0,32000),
                     breaks = c(0, 2000, 10000, 20000 ,30000)) +
  scale_shape_manual(values = c(21)) +
  geom_point(size = 2, aes(shape = sample_type), fill= "white",
             colour = "black") +
  theme_bw() +
  multiplot_theme +
  labs(y = "Variant fraction (%)", x = "Genome equivalents (GE)",
       title = het_gdna_plot_title)

###################
# Plot: z score analysis
###################

cfdna_zscore_plot_title <- paste0("ddPCR for ",
                                  nrow(all_samples_unblinded),
                                  " cfDNA samples with Z score classification")

plot_2d <- ggplot(all_samples_unblinded, aes(x = vf_assay_molecules, 
                                             y = variant_percent)) +
  theme_bw() +
  multiplot_theme + 
  z3_line +
  zminus3_line +
  z2_line +
  zminus2_line +
  vertical_line +
  multiplot_y +
  scale_x_continuous(limits = c(0,32000),
                     breaks = c(0, 2000, 10000, 20000 ,30000)) +
  cfdna_fill +
  cfdna_alpha +
  cfdna_shape +
  geom_point(size = 2, aes(fill = zscore_outcome,
                           alpha = zscore_outcome,
                           shape = fetal_genotype),
             colour = "black") +
  
  labs(y = "", x = "Genome equivalents (GE)", 
       title = cfdna_zscore_plot_title)

###################
# All plots together
###################

ddpcr_cohort <- ggpubr::ggarrange(plot_2a, plot_2b,
                                  plot_2c, plot_2d,
                                  ncol = 2, nrow = 2, align = "v",
                                  labels = c("A", "B", "C", "D"))

ggsave(plot = ddpcr_cohort, 
       filename = paste0("figure_2_",
                         format(Sys.time(), "%Y%m%d_%H%M%S"),
                         ".tiff"),
       path = "plots/", device='tiff', dpi=600,
       units = "in",
       width = 12.5,
       height = 7)

#######################


duplicated(ddpcr_target_panel$assay)




###################
# Sickle cell disease results table
################### 





scd_results <- supplementary_table %>%
  filter(vf_assay == "HBB c.20A>T")

genotype_count <- count(scd_results, fetal_genotype)

get_scd_counts <- function(analysis_method, analysis_outcome, df) {
  
  stopifnot("fetal_genotype" %in% colnames(df))
  
  output <- count(df, fetal_genotype,
                  !!analysis_outcome) %>%
    mutate(analysis = analysis_method) %>%
    pivot_wider(
      id_cols = fetal_genotype,
      names_from = c(analysis, !!analysis_outcome),
      values_from = n)
  return(output)
}

sprt_scd_count <- get_scd_counts(analysis_method = "sprt", 
                                 analysis_outcome = quo(sprt_outcome),
                                 df = scd_results)

mcmc_scd_count <- get_scd_counts(analysis_method = "mcmc", 
                                 analysis_outcome = quo(mcmc_outcome),
                                 df = scd_results)

zscore_scd_count <- get_scd_counts(analysis_method = "zscore", 
                                   analysis_outcome = quo(zscore_outcome),
                                   df = scd_results)

scd_analysis_results <- cbind(
  genotype_count,
  sprt_scd_count %>%
    select(-fetal_genotype),
  mcmc_scd_count%>%
    select(-fetal_genotype),
  zscore_scd_count%>%
    select(-fetal_genotype)) %>%
  dplyr::rename(Samples = n) %>%
  mutate(fetal_genotype = as.character(fetal_genotype))

# Replace NAs with 0s
scd_analysis_results[is.na(scd_analysis_results)] = 0
scd_analysis_results[scd_analysis_results == "homozygous variant"] <- "HbSS"
scd_analysis_results[scd_analysis_results == "heterozygous"] <- "HbAS"
scd_analysis_results[scd_analysis_results == "homozygous reference"] <- "HbAA"

write.csv(scd_analysis_results,
          file = paste0("analysis_outputs/scd_analysis_results",
                        format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
          row.names = FALSE)

###################
# Gene information table
###################

xl_count_table <- all_samples_blinded %>%
  filter(inheritance_chromosomal == "x_linked") %>%
  count(vf_assay)

ad_count_table <- all_samples_blinded %>%
  filter(inheritance_chromosomal == "autosomal" &
           inheritance_pattern == "dominant") %>%
  count(vf_assay)

ar_count_table <- all_samples_blinded %>%
  filter(inheritance_chromosomal == "autosomal" &
           inheritance_pattern == "recessive") %>%
  count(vf_assay)

vf_assay_count_table <- rbind(xl_count_table, 
                              ad_count_table,
                              ar_count_table) %>%
  dplyr::rename(samples = n) %>%
  
  # Bind to gene information
  left_join(gene_info,
            by = "vf_assay") %>%
  
  # Add inheritance information
  left_join(
    ddpcr_target_panel %>%
      select(assay, inheritance_chromosomal, inheritance_pattern) %>%
      #Remove duplicate rows
      distinct(.keep_all = TRUE) %>%
      dplyr::rename(vf_assay = assay),
    by = "vf_assay") %>%
  
  # Add on abbreviation of inheritance patterns
  mutate(inheritance_abbreviation = case_when(
    inheritance_chromosomal == "autosomal" &
      inheritance_pattern == "recessive" ~"AR",
    inheritance_chromosomal == "autosomal" &
      inheritance_pattern == "dominant" ~"AD",
    inheritance_chromosomal == "x_linked" &
      inheritance_pattern == "dominant" ~"XLD",
    inheritance_chromosomal == "x_linked" &
      inheritance_pattern == "recessive" ~"XLR"),
    
    # Factorise for ordering
    inheritance_abbreviation = factor(inheritance_abbreviation,
                                      levels = c("AR", "AD", "XLR", "XLD"))) %>%
  
  arrange(inheritance_abbreviation, vf_assay) %>%
  
  select(inheritance_abbreviation, gene, transcript, variant_dna, 
         variant_protein, condition,
         samples) %>%
  
  # Rename without underscores
  dplyr::rename(Inheritance = inheritance_abbreviation,
                Gene = gene,
                Transcript = transcript,
                DNA = variant_dna,
                Protein = variant_protein,
                Condition = condition,
                Samples = samples)

write.csv(vf_assay_count_table, 
          file = (paste0("analysis_outputs/Table 1 ",
                         format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")),
          row.names = FALSE)






###################
# SPRT
###################
sprt_hbss_correct <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous variant" 
                                              & sickle_cell_results$sprt_outcome == "correct", "sample_id"])

sprt_hbss_inconclusive <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous variant" 
                                                   & sickle_cell_results$sprt_outcome == "inconclusive", "sample_id"])                          

sprt_hbss_incorrect <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous variant" 
                                                & sickle_cell_results$sprt_outcome == "incorrect", "sample_id"])                          

sprt_hbas_correct <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "heterozygous" 
                                              & sickle_cell_results$sprt_outcome == "correct", "sample_id"])

sprt_hbas_inconclusive <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "heterozygous" 
                                                   & sickle_cell_results$sprt_outcome == "inconclusive", "sample_id"])                          

sprt_hbas_incorrect <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "heterozygous" 
                                                & sickle_cell_results$sprt_outcome == "incorrect", "sample_id"])  

sprt_hbaa_correct <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous reference" 
                                              & sickle_cell_results$sprt_outcome == "correct", "sample_id"])

sprt_hbaa_inconclusive <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous reference" 
                                                   & sickle_cell_results$sprt_outcome == "inconclusive", "sample_id"])                          

sprt_hbaa_incorrect <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous reference" 
                                                & sickle_cell_results$sprt_outcome == "incorrect", "sample_id"])  

sprt_correct_total <- nrow(sickle_cell_results[sickle_cell_results$sprt_outcome == "correct", "sample_id"])  

sprt_inconclusive_total <- nrow(sickle_cell_results[sickle_cell_results$sprt_outcome == "inconclusive", "sample_id"])  

sprt_incorrect_total <- nrow(sickle_cell_results[sickle_cell_results$sprt_outcome == "incorrect", "sample_id"])  

###################
# MCMC
###################
mcmc_hbss_correct <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous variant" 
                                              & sickle_cell_results$mcmc_outcome == "correct", "sample_id"])

mcmc_hbss_inconclusive <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous variant" 
                                                   & sickle_cell_results$mcmc_outcome == "inconclusive", "sample_id"])

mcmc_hbss_incorrect <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous variant" 
                                                & sickle_cell_results$mcmc_outcome == "incorrect", "sample_id"])

mcmc_hbas_correct <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "heterozygous" 
                                              & sickle_cell_results$mcmc_outcome == "correct", "sample_id"])

mcmc_hbas_inconclusive <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "heterozygous" 
                                                   & sickle_cell_results$mcmc_outcome == "inconclusive", "sample_id"])

mcmc_hbas_incorrect <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "heterozygous" 
                                                & sickle_cell_results$mcmc_outcome == "incorrect", "sample_id"])

mcmc_hbaa_correct <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous reference" 
                                              & sickle_cell_results$mcmc_outcome == "correct", "sample_id"])

mcmc_hbaa_inconclusive <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous reference" 
                                                   & sickle_cell_results$mcmc_outcome == "inconclusive", "sample_id"])

mcmc_hbaa_incorrect <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous reference" 
                                                & sickle_cell_results$mcmc_outcome == "incorrect", "sample_id"])

mcmc_correct_total <- nrow(sickle_cell_results[sickle_cell_results$mcmc_outcome == "correct", "sample_id"])  

mcmc_inconclusive_total <- nrow(sickle_cell_results[sickle_cell_results$mcmc_outcome == "inconclusive", "sample_id"])  

mcmc_incorrect_total <- nrow(sickle_cell_results[sickle_cell_results$mcmc_outcome == "incorrect", "sample_id"])  

###################
# Z score
###################

zscore_hbss_correct <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous variant" 
                                                & sickle_cell_results$zscore_outcome == "correct", "sample_id"])

zscore_hbss_inconclusive <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous variant" 
                                                     & sickle_cell_results$zscore_outcome == "inconclusive", "sample_id"])

zscore_hbss_incorrect <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous variant" 
                                                  & sickle_cell_results$zscore_outcome == "incorrect", "sample_id"])

zscore_hbas_correct <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "heterozygous" 
                                                & sickle_cell_results$zscore_outcome == "correct", "sample_id"])

zscore_hbas_inconclusive <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "heterozygous" 
                                                     & sickle_cell_results$zscore_outcome == "inconclusive", "sample_id"])

zscore_hbas_incorrect <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "heterozygous" 
                                                  & sickle_cell_results$zscore_outcome == "incorrect", "sample_id"])

zscore_hbaa_correct <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous reference" 
                                                & sickle_cell_results$zscore_outcome == "correct", "sample_id"])

zscore_hbaa_inconclusive <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous reference" 
                                                     & sickle_cell_results$zscore_outcome == "inconclusive", "sample_id"])

zscore_hbaa_incorrect <- nrow(sickle_cell_results[sickle_cell_results$fetal_genotype == "homozygous reference" 
                                                  & sickle_cell_results$zscore_outcome == "incorrect", "sample_id"])

zscore_correct_total <- nrow(sickle_cell_results[sickle_cell_results$zscore_outcome == "correct", "sample_id"])  

zscore_inconclusive_total <- nrow(sickle_cell_results[sickle_cell_results$zscore_outcome == "inconclusive", "sample_id"])  

zscore_incorrect_total <- nrow(sickle_cell_results[sickle_cell_results$zscore_outcome == "incorrect", "sample_id"])  

###################
# New sickle cell disease table
###################

sickle_table_new <- data.frame(
  Outcome = c("Correct", "Inconclusive", "Incorrect"),
  HbSS = c(sprt_hbss_correct, sprt_hbss_inconclusive, sprt_hbss_incorrect),
  HbAS = c(sprt_hbas_correct, sprt_hbas_inconclusive, sprt_hbas_incorrect),
  HbAA = c(sprt_hbaa_correct, sprt_hbaa_inconclusive, sprt_hbaa_incorrect),
  total = c(sprt_correct_total, sprt_inconclusive_total, sprt_incorrect_total),
  HbSS = c(mcmc_hbss_correct, mcmc_hbss_inconclusive, mcmc_hbss_incorrect),
  HbAS = c(mcmc_hbas_correct, mcmc_hbas_inconclusive, mcmc_hbas_incorrect),
  HbAA = c(mcmc_hbaa_correct, mcmc_hbaa_inconclusive, mcmc_hbaa_incorrect),
  total = c(mcmc_correct_total, mcmc_inconclusive_total, mcmc_incorrect_total),
  HbSS = c(zscore_hbss_correct, zscore_hbss_inconclusive, zscore_hbss_incorrect),
  HbAS = c(zscore_hbas_correct, zscore_hbas_inconclusive, zscore_hbas_incorrect),
  HbAA = c(zscore_hbaa_correct, zscore_hbaa_inconclusive, zscore_hbaa_incorrect),
  total = c(zscore_correct_total, zscore_inconclusive_total, zscore_incorrect_total))


results_by_genotype <- function(cohort_input, analysis_type, outcome_input) {
  
  results <- supplementary_table %>%
    filter(cohort == cohort_input) %>%
    filter(!!analysis_type == outcome_input) %>%
    count(!!quo(fetal_genotype)) 
  
  return(results)
}

results_by_genotype("sickle cell disease", quo(sprt_outcome),
                    "incorrect")







table_2 <- table_1 %>%
  filter(!!analysis_type != "inconclusive")%>%
  select(-percent) %>%
  mutate(percent_2 = round((total_count/count_incon_excluded)*100,0)) %>%
  select(!!analysis_type, percent_2)

table_joined <- table_1 %>%
  left_join(table_2, by = c(table_1[,1] = table_2[,1])) %>%
  dplyr::rename(HbSS = `homozygous variant`,
                HbAS = heterozygous,
                HbAA = `homozygous reference`) %>%
  replace(is.na(.),0)






























scd_mcmc_results <- supplementary_table %>%
  filter(cohort == "sickle cell disease") %>%
  group_by(mcmc_outcome, fetal_genotype) %>%
  summarise(total = n()) %>%
  pivot_wider(id_cols = mcmc_outcome,
              names_from = fetal_genotype,
              values_from = total) %>%
  replace(is.na(.), 0) %>%
  mutate(total_count = `homozygous variant` + heterozygous +
           `homozygous reference`,
         percent = round((total_count/genotype_count)*100,0),
         percent_2 = round((total_count/count_qc_fail_excluded)*100,0))






genotype_count <- nrow(scd_results)

results_by_genotype <- function(cohort_input, analysis_type, outcome_input) {
  
  results <- supplementary_table %>%
    filter(cohort == cohort_input) %>%
    filter(!!analysis_type == outcome_input) %>%
    count(!!quo(fetal_genotype)) 
  
  return(results)
}





total <- sum(results$n)

output_table <- results %>%
  mutate(outcome = outcome_input) %>%
  pivot_wider(
    id_cols = outcome,
    names_from = fetal_genotype,
    values_from = n) %>%
  mutate(total = total_correct,
         percent = round((total/genotype_count)*100,0))



scd_sprt <- rbind(results_by_genotype("sickle cell disease", quo(sprt_outcome),
                                      "correct"),
                  results_by_genotype("sickle cell disease", quo(sprt_outcome),
                                      "incorrect"),
                  results_by_genotype("sickle cell disease", quo(sprt_outcome),
                                      "inconclusive"))








results_by_analysis <- function(cohort_input, analysis_type) {
  
  correct_row <- results_by_genotype(cohort_input, !!analysis_type,
                                     "correct")
  
  incorrect_row <- results_by_genotype(cohort_input, !!analysis_type,
                                       "incorrect")
  
  inconclusive_row <- results_by_genotype(cohort_input, !!analysis_type,
                                          "inconclusive")
  
  output <- rbind(correct_row, incorrect_row, inconclusive_row)
  
}

results_by_genotype("sickle cell disease", quo(sprt_outcome),
                    "correct")

results_by_analysis("sickle cell disease", quo(sprt_outcome))



sickle_cell_results <- supplementary_table %>%
  filter(vf_assay == "HBB c.20A>T") %>%
  filter(sprt_outcome == "correct")






sickle_cell_results %>%
  filter(sprt_outcome == outcome) %>%
  count(!!quo(fetal_genotype))



correct_results_wide <- correct_results %>%
  mutate(outcome = "correct") %>%
  pivot_wider(
    id_cols = outcome,
    names_from = fetal_genotype,
    values_from = n) %>%
  
  
  
  
  
  
  
  mutate(sprt_outcome = factor(sprt_outcome, levels = 
                                 c("correct", "inconclusive",
                                   "incorrect"))) %>%
  arrange(sprt_outcome) %>%
  mutate(total = )



get_scd_counts <- function(analysis_method, analysis_outcome, df) {
  
  stopifnot("fetal_genotype" %in% colnames(df))
  
  output <- count(df, fetal_genotype,
                  !!analysis_outcome) %>%
    mutate(analysis = analysis_method) %>%
    pivot_wider(
      id_cols = fetal_genotype,
      names_from = c(analysis, !!analysis_outcome),
      values_from = n)
  return(output)
}



get_bespoke_table(quo(sprt_outcome))


%>%
  pivot_wider(id_cols = !!analysis_type,
              names_from = fetal_genotype,
              values_from = total) %>%
  replace(is.na(.), 0) %>%
  mutate(total_count = `homozygous variant` + heterozygous +
           `homozygous reference`,
         percent = round((total_count/genotype_count)*100,0),
         percent_2 = round((total_count/count_incon_excluded)*100,0)) %>%
  dplyr::rename(HbSS = `homozygous variant`,
                HbAS = heterozygous,
                HbAA = `homozygous reference`)



bespoke_genotype_count <- nrow(supplementary_table %>%
                                 filter(cohort != "sickle cell disease"))

bespoke_count_incon_excluded <- nrow(supplementary_table %>%
                                       filter(cohort != "sickle cell disease"
                                              & !!analysis_type != "inconclusive"))

bespoke_sprt <- supplementary_table %>%
  filter(vf_assay != "HBB c.20A>T") %>%
  group_by(sprt_outcome) %>%
  summarise(sprt = n()) %>%
  dplyr::rename(outcome = sprt_outcome)

bespoke_mcmc <- supplementary_table %>%
  filter(vf_assay != "HBB c.20A>T") %>%
  group_by(mcmc_outcome) %>%
  summarise(mcmc = n()) %>%
  select(-mcmc_outcome)

bespoke_zscore <- supplementary_table %>%
  filter(vf_assay != "HBB c.20A>T") %>%
  group_by(zscore_outcome) %>%
  summarise("Z score" = n()) %>%
  select(-zscore_outcome)

bespoke_table_new <- cbind(bespoke_sprt, bespoke_mcmc, bespoke_zscore)



new_bespoke_results_table <- bespoke_results %>%
  left_join(vf_assay_info %>%
              dplyr::rename(DNA = variant_dna) %>%
              select(transcript, DNA), by = "DNA") %>%
  select(Inheritance, Sample.number, Condition, transcript, Gene, DNA,
         Fetal.genotype, SPRT, MCMC, Z.score)



