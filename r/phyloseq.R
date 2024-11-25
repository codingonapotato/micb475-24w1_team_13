library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)

#### Load data ####
# Change file paths as necessary
metafp <- "data/pd_filtered_final_metadata.txt"
meta <- read_delim(metafp, delim = "\t")

otufp <- "data/feature-table.txt"
otu <- read_delim(file = otufp, delim = "\t", skip = 1)

taxfp <- "data/taxonomy.tsv"
tax <- read_delim(taxfp, delim = "\t")

phylotreefp <- "data/tree.nwk"
phylotree <- read.tree(phylotreefp)

#### Format OTU table ####
# OTU tables should be a matrix
# with rownames and colnames as OTUs and sampleIDs, respectively
# Note: tibbles do not allow rownames so if you imported with read_delim, change back

# save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[, -1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
class(OTU)

#### Format sample metadata ####
# Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[, -1])
# Make sampleids the rownames
rownames(samp_df) <- meta$"#SampleID"
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP)

#### Formatting taxonomy ####
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>%
  select(-Confidence) %>%
  separate(
    col = Taxon, sep = "; ",
    into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  ) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs
tax_mat <- tax_mat[, -1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX)

#### Create phyloseq object ####
# Merge all into a phyloseq object
depression <- phyloseq(OTU, SAMP, TAX, phylotree)

#### Looking at phyloseq object #####
# View components of phyloseq object with the following commands
otu_table(depression)
sample_data(depression)
tax_table(depression)
phy_tree(depression)

#### Quality Control ####
# TODO(?)

#### Save phyloseq object ####
save(depression, file = "data/depression_phyloseq.RData")

#### Conduct PERMANOVA and create distance matrices ####
# Weighted Unifrac
dm_wunifrac <- UniFrac(depression, weighted=TRUE) 
adonis2(dm_wunifrac ~ Antidepressant_use*BDI_category, data = meta)

# Unweighted Unifrac
dm_unifrac <- UniFrac(depression, weighted=FALSE)
adonis2(dm_unifrac ~ Antidepressant_use*BDI_category, data = meta)

# Bray Curtis
dm_braycurtis <- vegdist(t(otu_table(depression)), method="bray")
adonis2(dm_braycurtis ~ Antidepressant_use*BDI_category, data = meta)

#### Plot as an ordination to a PCoA plot ####
ord.unifrac <- ordinate(depression, method="PCoA", distance="unifrac")
plot_ordination(depression, ord.unifrac, color="BDI_category", shape = "Antidepressant_use")

ord.wunifrac <- ordinate(depression, method="PCoA", distance="wunifrac")
plot_ordination(depression, ord.wunifrac, color="BDI_category", shape = "Antidepressant_use")

ord.bray <- ordinate(depression, method="PCoA", distance="bray")
plot_ordination(depression, ord.bray, color="BDI_category", shape = "Antidepressant_use")

#### DESeq ####

depression_plus1 <- transform_sample_counts(depression, function(x) x+1)
deseq_object <- phyloseq_to_deseq2(depression_plus1, ~BDI_category_antidepressant_use)
DESeq_output <- DESeq(deseq_object)

## high/no vs high/yes
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
View(sigASVs_hnhy)
# Get only asv names
sigASVs_hnhy_vec <- sigASVs_hnhy %>%
  pull(ASV)
# Prune phyloseq file
depression_DESeq_hnhy <- prune_taxa(sigASVs_hnhy_vec,depression_plus1)
sigASVs_hnhy <- tax_table(depression_DESeq_hnhy) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))
# Plot results as bar plot
hnhy <- ggplot(sigASVs_hnhy) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
hnhy

## high/no vs low/yes
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
View(sigASVs_hnly)
# Get only asv names
sigASVs_hnly_vec <- sigASVs_hnly %>%
  pull(ASV)
# Prune phyloseq file
depression_DESeq_hnly <- prune_taxa(sigASVs_hnly_vec,depression_plus1)
sigASVs_hnly <- tax_table(depression_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_hnly) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))
# Plot results as bar plot
hnly <- ggplot(sigASVs_hnly) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
hnly

## high/no vs low/no
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
View(sigASVs_hnln)
# Get only asv names
sigASVs_hnln_vec <- sigASVs_hnln %>%
  pull(ASV)
# Prune phyloseq file
depression_DESeq_hnln <- prune_taxa(sigASVs_hnln_vec,depression_plus1)
sigASVs_hnln <- tax_table(depression_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_hnln) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))
# Plot results as bar plot
hnln <- ggplot(sigASVs_hnln) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
hnln

