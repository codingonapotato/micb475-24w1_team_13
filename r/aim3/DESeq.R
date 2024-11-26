library(tidyverse)
library(phyloseq)
library(DESeq2)

load("depression_phyloseq.Rdata")

#### Comparing between BDI + antidepressant use
depression_plus1 <- transform_sample_counts(depression, function(x) x+1)
deseq_object <- phyloseq_to_deseq2(depression_plus1, ~BDI_category_antidepressant_use)
DESeq_output <- DESeq(deseq_object)

#### high/no vs high/yes
res_hnhy <- results(DESeq_output, tidy=TRUE,
                    contrast = c("BDI_category_antidepressant_use","high+yes", "high+no"))
# Plot results as volcano plot
vol_plot_hnhy <- res_hnhy %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
vol_plot_hnhy
# To get table of results (Significant ASVs)
sigASVs_hnhy <- res_hnhy %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
# Get only asv names
sigASVs_hnhy_vec <- sigASVs_hnhy %>%
  pull(ASV)
# Prune phyloseq file
depression_DESeq_hnhy <- prune_taxa(sigASVs_hnhy_vec,depression_plus1)
sigASVs_hnhy <- tax_table(depression_DESeq_hnhy) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_hnhy) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))
# Plot results as bar plot
hnhy <- ggplot(sigASVs_hnhy) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("High/No vs High/Yes")
hnhy

#### high/no vs low/yes
res_hnly <- results(DESeq_output, tidy=TRUE,
                    contrast = c("BDI_category_antidepressant_use","low+yes", "high+no"))
# Plot results as volcano plot
vol_plot_hnly <- res_hnly %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
# To get table of results (Significant ASVs)
sigASVs_hnly <- res_hnly %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
# Get only asv names
sigASVs_hnly_vec <- sigASVs_hnly %>%
  pull(ASV)
# Prune phyloseq file
depression_DESeq_hnly <- prune_taxa(sigASVs_hnly_vec,depression_plus1)
sigASVs_hnly <- tax_table(depression_DESeq_hnly) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_hnly) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))
# Plot results as bar plot
hnly <- ggplot(sigASVs_hnly) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("High/No vs Low/Yes")
hnly

#### high/no vs low/no
res_hnln <- results(DESeq_output, tidy=TRUE,
                    contrast = c("BDI_category_antidepressant_use","low+no", "high+no"))
# Plot results as volcano plot
vol_plot_hnln <- res_hnln %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
# To get table of results (Significant ASVs)
sigASVs_hnln <- res_hnln %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
# Get only asv names
sigASVs_hnln_vec <- sigASVs_hnln %>%
  pull(ASV)
# Prune phyloseq file
depression_DESeq_hnln <- prune_taxa(sigASVs_hnln_vec,depression_plus1)
sigASVs_hnln <- tax_table(depression_DESeq_hnln) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_hnln) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))
# Plot results as bar plot
hnln <- ggplot(sigASVs_hnln) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  ggtitle("High/No vs Low/No")
hnln

#### high/yes vs low/no
res_hyln <- results(DESeq_output, tidy=TRUE,
                    contrast = c("BDI_category_antidepressant_use","low+no", "high+yes"))
# Plot results as volcano plot
vol_plot_hyln <- res_hyln %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
# To get table of results (Significant ASVs)
sigASVs_hyln <- res_hyln %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
# Get only asv names
sigASVs_hyln_vec <- sigASVs_hyln %>%
  pull(ASV)
# Prune phyloseq file
depression_DESeq_hyln <- prune_taxa(sigASVs_hyln_vec,depression_plus1)
sigASVs_hyln <- tax_table(depression_DESeq_hyln) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_hyln) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))
# Plot results as bar plot
hyln <- ggplot(sigASVs_hyln) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  ggtitle("High/Yes vs Low/No")
hyln

#### high/yes vs low/yes
res_hyly <- results(DESeq_output, tidy=TRUE,
                    contrast = c("BDI_category_antidepressant_use","low+yes", "high+yes"))
# Plot results as volcano plot
vol_plot_hyly <- res_hyly %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
# To get table of results (Significant ASVs)
sigASVs_hyly <- res_hyly %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
# Get only asv names
sigASVs_hyly_vec <- sigASVs_hyly %>%
  pull(ASV)
# Prune phyloseq file
depression_DESeq_hyly <- prune_taxa(sigASVs_hyly_vec,depression_plus1)
sigASVs_hyly <- tax_table(depression_DESeq_hyly) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_hyln) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))
# Plot results as bar plot
hyly <- ggplot(sigASVs_hyly) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  ggtitle("High/Yes vs Low/Yes")
hyly

#### low/yes vs low/no
res_lyln <- results(DESeq_output, tidy=TRUE,
                    contrast = c("BDI_category_antidepressant_use","low+no", "low+yes"))
# Plot results as volcano plot
vol_plot_lyln <- res_lyln %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
# To get table of results (Significant ASVs)
sigASVs_lyln <- res_lyln %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
# Get only asv names
sigASVs_lyln_vec <- sigASVs_lyln %>%
  pull(ASV)
# Prune phyloseq file
depression_DESeq_lyln <- prune_taxa(sigASVs_lyln_vec,depression_plus1)
sigASVs_lyln <- tax_table(depression_DESeq_lyln) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_lyln) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))
# Plot results as bar plot
lyln <- ggplot(sigASVs_lyln) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  ggtitle("Low/Yes vs Low/No")
lyln


#### Comparing between high vs low BDI
bdi_deseq_object <- phyloseq_to_deseq2(depression_plus1, ~BDI_category)
bdi_DESeq_output <- DESeq(bdi_deseq_object)

res_bdi <- results(bdi_DESeq_output, tidy=TRUE,
                    contrast = c("BDI_category","high", "low"))
# Plot results as volcano plot
vol_plot_bdi <- res_bdi %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
vol_plot_bdi
ggsave(filename="vol_plot_bdi.png",vol_plot_bdi)
# To get table of results (Significant ASVs)
sigASVs_bdi <- res_bdi %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
# Get only asv names
sigASVs_bdi_vec <- sigASVs_bdi %>%
  pull(ASV)
# Prune phyloseq file
depression_DESeq_bdi <- prune_taxa(sigASVs_bdi_vec,depression_plus1)
sigASVs_bdi <- tax_table(depression_DESeq_bdi) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_bdi) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))
# Plot results as bar plot
bdi <- ggplot(sigASVs_bdi) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("High vs Low BDI")
bdi
ggsave(filename="DESeq_bdi.png",bdi)


#### Comparing between no vs yes antidepressant use
atdp_deseq_object <- phyloseq_to_deseq2(depression_plus1, ~Antidepressant_use)
atdp_DESeq_output <- DESeq(atdp_deseq_object)

res_atdp <- results(atdp_DESeq_output, tidy=TRUE,
                   contrast = c("Antidepressant_use","No", "Yes"))
# Plot results as volcano plot
vol_plot_atdp <- res_atdp %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
vol_plot_atdp
ggsave(filename="vol_plot_atdp.png",vol_plot_atdp)
# To get table of results (Significant ASVs)
sigASVs_atdp <- res_atdp %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
# Get only asv names
sigASVs_atdp_vec <- sigASVs_atdp %>%
  pull(ASV)
# Prune phyloseq file
depression_DESeq_atdp <- prune_taxa(sigASVs_atdp_vec,depression_plus1)
sigASVs_atdp <- tax_table(depression_DESeq_atdp) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_atdp) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))
# Plot results as bar plot
atdp <- ggplot(sigASVs_atdp) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("No vs Yes Antidepressant Use")
atdp
ggsave(filename="DESeq_atdp.png",atdp)
                                            
