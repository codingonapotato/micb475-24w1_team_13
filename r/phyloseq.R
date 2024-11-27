library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)

#### Load data ####
# Change file paths as necessary
metafp <- "data/pd_filtered_metadata.txt"
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

#### rarefy for diversity metrics ####

rarecurve(t(as.data.frame(otu_table(depression))), cex=0.1)
depression_rare <- rarefy_even_depth(depression, rngseed = 1, sample.size = 4907)

#### Quality Control ####
# TODO(?)

#### Save phyloseq object ####
save(depression, file = "data/depression_phyloseq.RData")
save(depression_rare, file = "data/phyloseq_rare.RData")

