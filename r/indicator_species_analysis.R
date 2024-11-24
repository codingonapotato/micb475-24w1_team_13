library(tidyverse)
library(phyloseq)
library(indicspecies)
library(ggpubr)

#### Load Data ####
load("data/depression_phyloseq.RData")

#### Indicator Species/Taxa Analysis ####
# Glom to Species
depression_species <- tax_glom(depression, "Species", NArm = TRUE)
depression_species_RA <- transform_sample_counts(depression_species, fun = function(x) x / sum(x))

# ISA
isa_depression <- multipatt(t(otu_table(depression_species_RA)), cluster = sample_data(depression_species_RA)$BDI_category)
summary(isa_depression)
taxtable <- tax_table(depression_species) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV")

isa_final <- isa_depression$sign %>%
  rownames_to_column(var = "ASV") %>%
  left_join(taxtable) %>%
  filter(p.value < 0.05)

isa_final %>% View()

#### Save data ####
save(isa_final, file = "data/isa_final.RData")
