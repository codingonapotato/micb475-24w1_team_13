library(phyloseq)
library(vegan)
library(ape)
library(tidyverse)

#### Load Data Files ####

# Load OTU table
otu <- read_delim(file="data/feature-table.txt", delim="\t", skip=1)

# Load Metadata
metadata <- read_delim(file="data/pd_metadata_bdi_categorized.txt", "\t")

# Load taxonomy table
taxonomy <- read_delim(file="data/taxonomy.tsv", delim="\t")

# Load rooted phylogenetic tree
tree <- read.tree(file="data/tree.nwk")

#### Format data ####

# Replace the index column with OTU_ID as the row names
otu_matrix <- as.matrix(otu[, -1])
rownames(otu_matrix) <- otu$`#OTU ID`

# Convert metadata into dataframe + make sample_name rownames
metadata_df <- as.data.frame(metadata[, -1])
rownames(metadata_df) <- metadata$`#SampleID`

# Format taxonomy table
taxonomy_matrix <- taxonomy %>% 
  select(-Confidence) %>%
  separate(col=Taxon, sep="; ", 
           into=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  as.matrix()

# Make feature id into rownames for taxonomy_matrix
taxonomy_matrix <- taxonomy_matrix[, -1]
rownames(taxonomy_matrix) <- taxonomy$"Feature ID"

#### Generate phyloseq objects ####
OTU <- otu_table(otu_matrix, taxa_are_rows=TRUE)
SAMP <- sample_data(metadata_df)
TAX <- tax_table(taxonomy_matrix)

# Aggregate phyloseq objects (+phylogenic tree) into one object
pd_phyloseq <- phyloseq(OTU, SAMP, TAX, tree)

# TODO: Additional QC? (e.g.) removing ASV's w/ < 5 counts or samples w/ < 100 reads?

#### Rarefy data ####
rarecurve(t(as.data.frame(otu_table(pd_phyloseq))), cex=0.1)
pd_rare <- rarefy_even_depth(pd_phyloseq, rngseed=5, sample.size=5696) # sample.size determined from Qiime2 viewer

#### Save phyloseq object (before and after rarefaction) ####
save(pd_phyloseq, file="data/pd_phyloseq.RData")
save(pd_rare, file="data/pd_rare.RData")