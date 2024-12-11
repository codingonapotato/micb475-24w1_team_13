# Purpose: Generates an updated manifest file that is a subset of the original manifest file based on the sample ID's in the supplied metadata file
library(tidyverse)

#### Load manifest and final filtered metadata ####
pd_manifest <- read_tsv("raw/parkinsons_manifest.txt")
filtered_final_metadata <- read_tsv("pd_filtered_final_metadata.txt")

#### Subset pd_manifest such that only rows with sample ID in filtered_final_metadata are included ####
# Get vector of all sample ID's present in the filtered_final_metadata
samples <- filtered_final_metadata$`#SampleID`
# Subset pd_manifest
subset_pd_manifest <- pd_manifest %>% filter(`sample-id` %in% samples)

#### Save to file ####
export_file_path <- "filtered_parkinsons_manifest.txt"
write_tsv(subset_pd_manifest, export_file_path)
