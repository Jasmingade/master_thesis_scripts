library(data.table)
library(readr)
library(vroom)

meta <- read_csv("data/processed/metaData_multiomics.csv")
rna_dir <- "data/processed/RNA"
out_dir <- "data/processed/transcriptomics"
dir.create(out_dir, showWarnings = FALSE)

rna_files <- list(
  gene = "transcriptomics_gene_tpm.csv",
  iso_log = "transcriptomics_iso_tpm.csv",
  iso_frac = "transcriptomics_iso_frac.csv"
)

cancer_types <- unique(meta$tumor_code)

# Loop through RNA data types
for (block_name in names(rna_files)) {
  full_path <- file.path(rna_dir, rna_files[[block_name]])
  
  # Read only header to get all sample IDs
  header <- vroom(full_path, n_max = 0)
  all_cols <- names(header)
  clean_ids <- sub("(^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}).*", "\\1", all_cols)
  clean_ids[1] <- "feature_id"  # Fix first column name
  names(header) <- clean_ids
  
  for (cancer in cancer_types) {
    case_ids <- meta$case_id[meta$tumor_code == cancer]
    matched_cols <- which(clean_ids %in% case_ids)
    
    if (length(matched_cols) > 0) {
      sel_cols <- c(1, matched_cols)
      df <- fread(full_path, select = sel_cols)
      colnames(df)[1] <- "feature_id"
      out_file <- file.path(out_dir, paste0("RNA_", cancer, "_", block_name, ".csv"))
      fwrite(df, out_file)
    }
  }
}
