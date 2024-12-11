library(tidyverse)
library(phyloseq)
library(indicspecies)

#### Load Data ####
load("data/depression_phyloseq.RData")
bdi_anti_high_no_unique_ASVs <- as.vector(read.csv("data/bdi_anti_high_no_unique_members.csv"))

#### Indicator Species/Taxa Analysis ####
# Glom to Species
depression_species <- tax_glom(depression, "Species", NArm = TRUE)
depression_species_RA <- transform_sample_counts(depression_species, fun = function(x) x / sum(x))

# Modify taxnomy table
taxtable <- tax_table(depression_species) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV")

# ISA BDI Category
isa_bdi_category <- multipatt(t(otu_table(depression_species_RA)), 
                              cluster = sample_data(depression_species_RA)$BDI_category)
summary(isa_bdi_category)

isa_final_bdi_category <- isa_bdi_category$sign %>%
  rownames_to_column(var = "ASV") %>%
  left_join(taxtable) %>%
  filter(p.value < 0.05)

# ISA Antidepressant_use
isa_antidepressant <- multipatt(t(otu_table(depression_species_RA)), 
                              cluster = sample_data(depression_species_RA)$Antidepressant_use)
summary(isa_antidepressant)

isa_final_antidepressant <- isa_antidepressant$sign %>%
  rownames_to_column(var = "ASV") %>%
  left_join(taxtable) %>%
  filter(p.value < 0.05)

# ISA Antidepressant_use and BDI_category
isa_mixed <- multipatt(t(otu_table(depression_species_RA)), 
                       cluster = sample_data(depression_species_RA)$BDI_category_antidepressant_use)
summary(isa_mixed)

isa_final_mixed <- isa_mixed$sign %>%
  rownames_to_column(var = "ASV") %>%
  left_join(taxtable) %>%
  filter(p.value < 0.05)

# Filter and extract ASVs from isa_final_mixed with at least indicator species value of 0.5 and
# only the ASVs that are unique to the high BDI, no antidepressant use condition
isa_final_mixed_filtered <- filter(isa_final_mixed, stat >= 0.5 & 
                                          `s.high+no` == 1 &
                                          `s.high+yes` + `s.low+yes` + `s.low+no` == 0)
isa_final_mixed_ASVs_list <- as.list(isa_final_mixed_filtered$ASV)

# Find list of common ASVs shared between core microbiome group + ISA groups
common_ASVs <- intersect(isa_final_mixed_ASVs_list, bdi_anti_high_no_unique_ASVs)
# Filter table with this intersection group
common_isa <- filter(isa_final_mixed, isa_final_mixed$ASV %in% common_ASVs)

#### Save data ####
save(isa_final_bdi_category, file = "data/isa_final_bdi_category.RData")
save(isa_final_antidepressant, file = "data/isa_final_antidepressant.RData")
save(isa_final_mixed, file = "data/isa_final_mixed.RData")
save(isa_final_mixed_filtered, file = "data/isa_bdi_anti_high_no.RData")
save(common_isa, file = "data/common_isa.RData")