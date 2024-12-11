library(tidyverse)
library(phyloseq)
library(vegan)

load("depression_phyloseq.Rdata")

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
