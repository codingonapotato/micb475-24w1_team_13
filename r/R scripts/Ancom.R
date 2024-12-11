library(tidyverse)
library(phyloseq) 
library(ggpubr)
library(randomcoloR)
library(ANCOMBC)
library(microbiome)

set.seed(711)

load("depression_phyloseq.Rdata")

# Aggretate data to family level for analysis
family = tax_glom(depression,'Family')
ntaxa(depression); ntaxa(family)

# Optional filtering
calculate_relative_abundance <- function(x) x / sum(x)
total_counts <- taxa_sums(family)
relative_abundance <- calculate_relative_abundance(total_counts)
abundant <- relative_abundance > 0.001
family <- prune_taxa(abundant, family)
family

# Ensure the grouping variable is set correctly
sample_data(family)$BDI_category_antidepressant_use <- factor(
  sample_data(family)$BDI_category_antidepressant_use, 
  levels = c("low+no", "low+yes", "high+no", "high+yes")
)

# ANCOM-BC
ancom.family = ancombc(phyloseq = family,
                       formula = 'BDI_category_antidepressant_use',
                       p_adj_method = "fdr",
                       prv_cut=0.10,
                       lib_cut = 1000,
                       group = 'BDI_category_antidepressant_use',
                       struc_zero = T)

# Update column names for each part of the results
colnames(ancom.family$res$lfc) = paste(colnames(ancom.family$res$lfc),'_beta',sep='')
colnames(ancom.family$res$se) = paste(colnames(ancom.family$res$se),'_se',sep='')
colnames(ancom.family$res$W) = paste(colnames(ancom.family$res$W),'_W',sep='')
colnames(ancom.family$res$p_val) = paste(colnames(ancom.family$res$p_val),'_p_val',sep='')
colnames(ancom.family$res$q_val) = paste(colnames(ancom.family$res$q_val),'_q_val',sep='')
colnames(ancom.family$res$diff_abn) = paste(colnames(ancom.family$res$diff_abn),'_diff_abn',sep='')

results = lapply(ancom.family$res,function(x) rownames_to_column(x,'Family')) %>% 
  lapply(as_tibble) %>% reduce(full_join)

# Remove (Intercept) terms
results = results %>% dplyr::select(-contains('Intercept'))

# Plotting the results with significant taxa
q_val_columns <- grep("_q_val$", colnames(results), value = TRUE)
results.sig = results %>%
  filter(if_any(all_of(q_val_columns), ~ . < 0.05))

hits = results.sig$Family
hits = as.numeric(hits)

family_tss = family %>% transform('compositional')
selected_taxa <- taxa_names(family)[hits]
family_tss <- prune_taxa(selected_taxa, family_tss)
family_tss_melt = family_tss %>% psmelt()

bug = unique(family_tss_melt$Family)
family_tss_melt = family_tss_melt %>% 
  dplyr::select(-c(OTU,Domain:Order)) %>%
  pivot_wider(names_from = Family, values_from = Abundance)

for(b in bug){
    p = family_tss_melt %>% 
    mutate(BDI_Category_Antidepressant_Use = str_to_title(BDI_category_antidepressant_use))
    
    colors <- randomColor(4)
  
  ggplot(p,aes(x = BDI_category_antidepressant_use, y = .data[[b]],fill=BDI_category_antidepressant_use)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(height = 0, width = 0.2) +
    theme_classic(base_size = 16) +
    theme(legend.position='none') +
    scale_fill_manual(values = colors) +
    xlab('BDI_category_antidepressant_use') + ylab(paste(b,'(% Ab.)',sep=' '))
  print(p)
  ggsave(paste('ancom3_',b,'.jpeg',sep=''),height = 5,width = 5)
}

# Save ANCOM results
saveRDS(ancom.family,'ancom_results_family.rds')

# Save reseults as .csv
write.csv(results,'ancom_results_family.csv',row.names = F)

#Save results as excel, with significant results on separate sheet
library(writexl)
write_xlsx(list('all_results' = results,'sig_results' = results.sig),
           'ancom_results_family.xlsx')
