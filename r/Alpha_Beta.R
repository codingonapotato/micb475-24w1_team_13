library(phyloseq)
library(ape)
library(tidyverse)
library(picante)

load("depression_phyloseq.Rdata")

#### Beta Diversity ####
# Plot as an ordination to a PCoA plot ####
ord.unifrac <- ordinate(depression, method="PCoA", distance="unifrac")
plot_ordination(depression, ord.unifrac, color="BDI_category", shape = "Antidepressant_use")

ord.wunifrac <- ordinate(depression, method="PCoA", distance="wunifrac")
plot_ordination(depression, ord.wunifrac, color="BDI_category", shape = "Antidepressant_use")

ord.bray <- ordinate(depression, method="PCoA", distance="bray")
plot_ordination(depression, ord.bray, color="BDI_category", shape = "Antidepressant_use")
