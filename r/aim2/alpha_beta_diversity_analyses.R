library(phyloseq)
library(tidyverse)
library(picante)

#### Import RData ####
load("data/pd_rare.RData")

#### Alpha Diversity ####
# Calculate Faith's Phylogenetic Diversity
faiths_pd <- pd(t(otu_table(pd_rare)), phy_tree(pd_rare), include.root=F)
# Add the Faith's PD results as to metadata table of phyloseq object
sample_data(pd_rare)$PD <- faiths_pd$PD
# Box plot for Shannon Diveristy
gg_richness <- plot_richness(pd_rare, x = "BDI_category", measures = c("Shannon")) +
  geom_boxplot()
# Box plot for Faith's PD
gg_faiths_pd <- ggplot(sample_data(pd_rare), aes(BDI_category, PD)) +
  geom_boxplot()
# TODO: Pielou's Evenness
# Save diversity figures
ggsave(filename="figures/shannon.png", gg_richness,
       height=5, width=3)
ggsave(filename="figures/faiths_pd.png", gg_faiths_pd, 
       height=5, width=3)

#### Beta Diversity ####