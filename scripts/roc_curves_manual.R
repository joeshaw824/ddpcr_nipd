###################
# Making ROC curves using tidyverse
###################

roc_binary_calls <- all_samples_unblinded %>%
  # Convert invasive results to binary outcomes
  mutate(
    fetal_genotype_binary = case_when(
      # Autosomal dominant inheritance
      inheritance_chromosomal == "autosomal" &
        inheritance_pattern == "dominant" &
        fetal_genotype == "heterozygous" ~"TRUE",
      inheritance_chromosomal == "autosomal" &
        inheritance_pattern == "dominant" &
        fetal_genotype == "homozygous reference" ~"FALSE",
      
      # Autosomal recessive inheritance
      inheritance_chromosomal == "autosomal" &
        inheritance_pattern == "recessive" &
        fetal_genotype %in% c("homozygous variant", 
                              "homozygous reference")  ~"TRUE",
      inheritance_chromosomal == "autosomal" &
        inheritance_pattern == "recessive" &
        fetal_genotype == "heterozygous" ~"FALSE",
      
      # X linked inheritance
      inheritance_chromosomal == "x_linked" &
        fetal_genotype == "hemizygous variant"  ~"TRUE",
      inheritance_chromosomal == "x_linked" &
        fetal_genotype == "hemizygous reference"  ~"FALSE"),
    
    # Convert to a Boolean vector
    fetal_genotype_binary = as.logical(fetal_genotype_binary),
    
    # Convert the MCMC calls to a binary outcome
    # Annoyingly, pG2 has different meanings for recessive and dominant
    # cases.
    mcmc_unbalanced_call = case_when(
      inheritance_chromosomal == "autosomal" &
        inheritance_pattern == "recessive" ~pmax(p_G1, p_G3),
      inheritance_chromosomal == "autosomal" &
        inheritance_pattern == "dominant" ~p_G2,
      inheritance_chromosomal == "x_linked" ~pmax(p_G0, p_G1)),
    # Remove minus signs from z score
    zscore_unbalanced_call = abs(zscore))

sprt_roc2 <- roc_binary_calls %>%
  arrange(desc(likelihood_ratio)) %>%
  # tpr is "true positive rate" and fpr is "false positive rate"
  mutate(sprt_tpr = cumsum(fetal_genotype_binary)/sum(fetal_genotype_binary)) %>%
  mutate(sprt_fpr = cumsum(!fetal_genotype_binary)/sum(!fetal_genotype_binary)) %>%
  select(r_number, likelihood_ratio, sprt_tpr, sprt_fpr, fetal_genotype_binary)

mcmc_roc <- roc_binary_calls %>%
  arrange(desc(mcmc_unbalanced_call)) %>%
  mutate(mcmc_tpr = cumsum(unbalanced)/sum(unbalanced)) %>%
  mutate(mcmc_fpr = cumsum(!unbalanced)/sum(!unbalanced)) %>%
  select(r_number, mcmc_unbalanced_call, mcmc_tpr, mcmc_fpr)

zscore_roc <- roc_binary_calls %>%
  arrange(desc(zscore_unbalanced_call)) %>%
  mutate(zscore_tpr = cumsum(unbalanced)/sum(unbalanced)) %>%
  mutate(zscore_fpr = cumsum(!unbalanced)/sum(!unbalanced)) %>%
  select(r_number, zscore_tpr, zscore_fpr)

sprt_roc_plot <- ggplot(sprt_roc2, aes(x = sprt_fpr, y = sprt_tpr))+
  geom_line(size = 2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(linetype = "dashed")+
  ylim(0, 1)+
  xlim(0,1)+
  labs(x = "",
       y = "True positive rate", 
       title = "SPRT analysis") +
  annotate(geom = "text", x = 0.25, y = 0.2,
           label = paste0("Sensitivity = ", as.character(sprt_sen))) +
  annotate(geom = "text", x = 0.25, y = 0.1,
           label = paste0("Specificity = ", as.character(sprt_spe))) +
  geom_vline(xintercept = 0.939, linetype = "dashed") +
  geom_hline(yintercept = 0.957, linetype = "dashed")
  

mcmc_roc_plot <- ggplot(mcmc_roc, aes(x = mcmc_fpr, y = mcmc_tpr))+
  geom_line(size = 2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(linetype = "dashed")+
  ylim(0, 1)+
  xlim(0,1)+
  labs(x = "False positive rate",
       y = "", 
       title = "MCMC analysis")

zscore_roc_plot <- ggplot(zscore_roc, aes(x = zscore_fpr, y = zscore_tpr))+
  geom_line(size = 2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(linetype = "dashed")+
  ylim(0, 1)+
  xlim(0,1)+
  labs(x = "",
       y = "", 
       title = "Z score analysis")

# WE FORGOT TO FILTER THE INCONCLUSIVE RESULTS
# AND WE FORGOT THAT autosomal dominant cases will have "negative" predictions
# for "positive" affected status

# All plots together
roc_plot <- ggpubr::ggarrange(sprt_roc_plot, mcmc_roc_plot, 
                              zscore_roc_plot, 
                              ncol = 3, nrow = 1)

ggsave(plot = roc_plot, 
       filename = "roc_plot.tiff",
       path = "plots/", device='tiff')



write.csv(sprt_roc2, "analysis_outputs/sprt_roc2.csv",
          row.names = FALSE)

