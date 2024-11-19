library("tidyverse")

#### Load PD metadata ####
file_path <- "raw/parkinsons_metadata.txt"
metadata <- read_tsv(file_path)

#### Filter for Disease status == "PD" ####
pd_metadata <- metadata %>% filter(Disease == "PD")

#### Check if BDI_depression_score column is tidy ####
na_count <- pd_metadata %>%
  filter(is.na(BDI_depression_score)) %>%
  nrow()

# BDI Category Legend:
# Low: 0 - 3
# Medium: 4 -5
# High: >= 16
LOW_BDI_UPPER_THRESHOLD <- 3
HIGH_BDI_LOWER_THRESHOLD <- 16

#### Generate new column of metadata to categorize by BDI score ####
updated_metadata <- pd_metadata %>%
  mutate(BDI_category = BDI_depression_score %>%
    map_chr(function(x) {
      if (x <= LOW_BDI_UPPER_THRESHOLD) {
        "low"
      } else if (x >= HIGH_BDI_LOWER_THRESHOLD) {
        "high"
      } else {
        "medium"
      }
    }))

#### Filter for only rows with a BDI_category of "low" or "high" ####
high_low_metadata <- updated_metadata %>% filter(BDI_category != "medium")

#### Generate a new column of metadata containing combinations of BDI_category and Antidepressant_use values (delimited with "+") ####
# Combination legend:
# 1. "low+yes"
# 2. "low+no"
# 3. "high+yes"
# 4. "high+no"
final_metadata <- high_low_metadata %>% mutate(
  BDI_category_antidepressant_use =
    tolower(paste(as.character(BDI_category), as.character(Antidepressant_use), sep = "+"))
)

#### Export BDI categorized metadata as a new file ####
export_file_path <- "pd_filtered_final_metadata.tsv"
write_tsv(updated_metadata, export_file_path)
