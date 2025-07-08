# ---- Libraries ----
library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(data.table)


# ---- Settings ----
input_dir <- "data/analysis_results/probatch_diagnostics/batch_correction"
output_dir <- "data/processed/proteomics_unified"
data_types <- c("gene", "iso_log", "iso_frac")

# Create output directory if not exists
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- Merge Function ----
merged_matrices <- list()

for (data_type in data_types) {
  files <- list.files(input_dir, pattern = paste0("_", data_type, "_batch_corrected.csv$"), full.names = TRUE)
  
  all_samples <- lapply(files, function(file) {
    mat <- read.csv(file, row.names = 1, check.names = FALSE)
    cancer_type <- str_extract(basename(file), "PDC[0-9]+")
    colnames(mat) <- paste0(colnames(mat))
    return(mat)
  })
  
  all_samples_cleaned <- lapply(all_samples, function(df) {
    df <- df[, !colnames(df) %in% "feature"]  # Remove existing 'feature' column if present
    rownames_to_column(df, var = "feature")
  })
  
  # Full outer join across all
  all_samples_full <- Reduce(function(x, y) full_join(x, y, by = "feature"), all_samples_cleaned)
  
  # Convert back to matrix with feature as rownames
  merged_matrix <- all_samples_full %>%
    column_to_rownames("feature") %>%
    as.matrix()
  
  merged_matrices[[data_type]] <- merged_matrix
  
  # Save for use
  fwrite(
    data.table::data.table(feature = rownames(merged_matrix), merged_matrix),
    file = file.path(output_dir, paste0("proteomics_", data_type, "_all_cancers.csv")),
    quote = FALSE
  )
}


