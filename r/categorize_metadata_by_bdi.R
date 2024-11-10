library('tidyverse')

# Load PD metadata
file_path <- './parkinsons_metadata.txt'
metadata <- read_tsv(file_path)

# Check if BDI_depression_score column is tidy
na_count <- metadata %>%
  filter(is.na(BDI_depression_score)) %>%
  nrow()

# BDI Category Legend: 
# Low: 0 - 3
# High: >= 16
LOW_BDI_UPPER_THRESHOLD <- 3
HIGH_BDI_LOWER_THRESHOLD <- 16

# Filter metadata for samples that are within our range of interest (low: 0 - 3, high: >= 16 )
filtered_metadata <- metadata %>%
    filter(BDI_depression_score <= 3 | BDI_depression_score >= 16)

# Generate new column of metadata to categorize by BDI score:
updated_metadata <- filtered_metadata %>%
    mutate(BDI_category = BDI_depression_score %>%
             map_chr(function(x) if (x <= LOW_BDI_UPPER_THRESHOLD) 'low' else 'high'))

# Export BDI categorized metadata as a new file
export_file_path <- 'pd_metadata_bdi_categorized.txt'
write_tsv(updated_metadata, export_file_path)