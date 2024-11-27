library(phyloseq)
library(microbiome)
library(tidyverse)
library(ggVennDiagram)

#### Load Data ####
load("data/depression_phyloseq.RData")

#### Core Microbiome Analysis ####
#### Transform relative abundances (RA) ####
depression_RA <- transform_sample_counts(depression, fun=function(x) x/sum(x))

#### Filtering Dataset ####
# Filter by Antidepressant_use
depression_anti_yes <- subset_samples(depression_RA, Antidepressant_use == "Yes")
depression_anti_no <- subset_samples(depression_RA, Antidepressant_use =="No")

# Filter by BDI_category
depression_bdi_high <- subset_samples(depression_RA, BDI_category == "high")
depression_bdi_low <- subset_samples(depression_RA, BDI_category =="low")

# Filter by BDI_category_antidepressant_use
depression_bdi_anti_low_no <- subset_samples(depression_RA, 
                                             BDI_category_antidepressant_use == "low+no")
depression_bdi_anti_low_yes <- subset_samples(depression_RA, 
                                             BDI_category_antidepressant_use == "low+yes")
depression_bdi_anti_high_no <- subset_samples(depression_RA, 
                                             BDI_category_antidepressant_use == "high+no")
depression_bdi_anti_high_yes <- subset_samples(depression_RA, 
                                             BDI_category_antidepressant_use == "high+yes")

#### Evaluate ASV's found in more than 'X'% of samples for each metadata category ####
# Constant for prevalence parameter passed to core_members function
PREVALENCE <- 0.7
# Constant for detection parameter passed to core_members function
DETECTION <- 0

# Antidepressant_use ASV's
anti_yes_ASVs <- core_members(depression_anti_yes, detection=DETECTION, prevalence = PREVALENCE) 
anti_no_ASVs <- core_members(depression_anti_no, detection=DETECTION, prevalence = PREVALENCE)
anti_list_full <- list(antidepressant_yes = anti_yes_ASVs, antidepressant_no = anti_no_ASVs)

# BDI_category ASV's
bdi_high_ASVs <- core_members(depression_bdi_high, detection=DETECTION, prevalence = PREVALENCE)
bdi_low_ASVs <- core_members(depression_bdi_low, detection=DETECTION, prevalence = PREVALENCE)
bdi_list_full <- list(bdi_high = bdi_high_ASVs, bdi_low = bdi_low_ASVs)

# BDI_category_antidepressant_use ASV's
bdi_anti_high_no_ASVs <- core_members(depression_bdi_anti_high_no, detection=DETECTION, prevalence = PREVALENCE)
bdi_anti_high_yes_ASVs <- core_members(depression_bdi_anti_high_yes, detection=DETECTION, prevalence = PREVALENCE)
bdi_anti_low_no_ASVs <- core_members(depression_bdi_anti_low_no, detection=DETECTION, prevalence = PREVALENCE)
bdi_anti_low_yes_ASVs <- core_members(depression_bdi_anti_low_yes, detection=DETECTION, prevalence = PREVALENCE)
bdi_anti_list_full <- list(bdi_anti_high_no = bdi_anti_high_no_ASVs,
                           bdi_anti_high_yes = bdi_anti_high_yes_ASVs,
                           bdi_anti_low_no = bdi_anti_low_no_ASVs,
                           bdi_anti_low_yes = bdi_anti_low_yes_ASVs)

#### Generate Venn Diagrams by metadata group ####
# Antidepressant_use
anti_venn <- ggVennDiagram(x = anti_list_full) +
  scale_x_continuous(expand = expansion(mult = .7))

# BDI_category
bdi_venn <- ggVennDiagram(x = bdi_list_full) +
  scale_x_continuous(expand = expansion(mult = .3))

# BDI_category_antidepressant_use
bdi_anti_venn <- ggVennDiagram(x = bdi_anti_list_full) +
  scale_x_continuous(expand = expansion(mult = .2))
bdi_anti_venn

#### Save plots to file ####
# Constant for dimensions
WIDTH <- 7
HEIGHT <- 5
ggsave(anti_venn, filename = "figures/antidepressant_use_venn_diagram.png", 
       width=WIDTH, height=HEIGHT)
ggsave(bdi_venn, filename = "figures/bdi_category_venn_diagram.png", 
       width=WIDTH, height=HEIGHT)
ggsave(bdi_anti_venn, 
       filename = "figures/combined_venn_diagram.png", 
       width=WIDTH, height=HEIGHT)

