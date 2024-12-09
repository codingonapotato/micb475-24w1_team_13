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
anti_list_full <- list(yes = anti_yes_ASVs, no = anti_no_ASVs)

# BDI_category ASV's
bdi_high_ASVs <- core_members(depression_bdi_high, detection=DETECTION, prevalence = PREVALENCE)
bdi_low_ASVs <- core_members(depression_bdi_low, detection=DETECTION, prevalence = PREVALENCE)
bdi_list_full <- list(high = bdi_high_ASVs, low = bdi_low_ASVs)

# BDI_category_antidepressant_use ASV's
bdi_anti_high_no_ASVs <- core_members(depression_bdi_anti_high_no, detection=DETECTION, prevalence = PREVALENCE)
bdi_anti_high_yes_ASVs <- core_members(depression_bdi_anti_high_yes, detection=DETECTION, prevalence = PREVALENCE)
bdi_anti_low_no_ASVs <- core_members(depression_bdi_anti_low_no, detection=DETECTION, prevalence = PREVALENCE)
bdi_anti_low_yes_ASVs <- core_members(depression_bdi_anti_low_yes, detection=DETECTION, prevalence = PREVALENCE)
bdi_anti_list_full <- list(HN = bdi_anti_high_no_ASVs,
                           HY = bdi_anti_high_yes_ASVs,
                           LN = bdi_anti_low_no_ASVs,
                           LY = bdi_anti_low_yes_ASVs)

# Unpack the ASVs into separate elements for downstream analysis
unpacked_bdi_anti_high_no_ASVs <- lapply(bdi_anti_high_no_ASVs, function(x) strsplit(x, "\t")[[1]])
unpacked_bdi_anti_high_yes_ASVs <- lapply(bdi_anti_high_yes_ASVs, function(x) strsplit(x, "\t")[[1]])
unpacked_bdi_anti_low_no_ASVs <- lapply(bdi_anti_low_no_ASVs, function(x) strsplit(x, "\t")[[1]])
unpacked_bdi_anti_low_yes_ASVs <- lapply(bdi_anti_low_yes_ASVs, function(x) strsplit(x, "\t")[[1]])

#### Generate Venn Diagrams by metadata group ####
# Antidepressant_use
anti_venn <- ggVennDiagram(x = anti_list_full) +
  scale_x_continuous(expand = expansion(mult = .5)) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(hjust = 0.5, margin=margin(10,0,30,0))) +
  ggtitle("Antidepressant Use")
anti_venn

# BDI_category
bdi_venn <- ggVennDiagram(x = bdi_list_full) +
  scale_x_continuous(expand = expansion(mult = .3)) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(hjust = 0.5, margin=margin(10,0,30,0))) +
  ggtitle("BDI Category")

# BDI_category_antidepressant_use
bdi_anti_venn <- ggVennDiagram(x = bdi_anti_list_full) +
  scale_x_continuous(expand = expansion(mult = .2)) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(hjust = 0.5, margin=margin(10,0,30,0))) +
  ggtitle("BDI Category & Antidepressant Use")

#### Extract the unique members from the BDI_category = High and Antidepressant_use = No group
# Calculate set difference
unique_ASVs <- setdiff(unpacked_bdi_anti_high_no_ASVs, unpacked_bdi_anti_high_yes_ASVs) %>%
    setdiff(unpacked_bdi_anti_low_yes_ASVs) %>%
    setdiff(unpacked_bdi_anti_low_no_ASVs)

# Coerce into dataframe
df_unique_ASVs <- as.data.frame(unique_ASVs)

# Write dataframe to file
write.csv(df_unique_ASVs, "data/bdi_anti_high_no_unique_members.csv")

#### Save plots to file ####
# Constant for dimensions
WIDTH <- 9
HEIGHT <- 5
ggsave(anti_venn, filename = "figures/antidepressant_use_venn_diagram.png", 
       width=WIDTH, height=HEIGHT) 
ggsave(bdi_venn, filename = "figures/bdi_category_venn_diagram.png", 
       width=WIDTH, height=HEIGHT)
ggsave(bdi_anti_venn, 
       filename = "figures/combined_venn_diagram.png", 
       width=WIDTH, height=HEIGHT)

