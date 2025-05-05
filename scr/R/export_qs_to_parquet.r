#!/usr/bin/env Rscript
library(qs)
library(arrow)
library(fs)

export_cleaned <- function(qs_path, prefix) {
  message("▶ Loading ", qs_path)
  cleaned <- qread(qs_path)
  
  for (study in names(cleaned)) {
    message("  • Study: ", study)
    # cleaned[[study]] is a tibble with cols: study_folder, data (list of run‐tibbles)
    study_tbl <- cleaned[[study]]
    
    # iterate over rows
    for (i in seq_len(nrow(study_tbl))) {
      run_info   <- study_tbl[i, ]             # one‐row tibble
      folder     <- run_info$study_folder      # the folder/run name
      run_tbl    <- run_info$data[[1]]         # unpack the list column → tibble with cols file_name, filepath, aggregated
      file_name  <- run_tbl$file_name[[1]]     # the original .psm file name (no extension)
      # aggregated itself is a list of length 1 containing the cleaned tibble
      df         <- run_tbl$aggregated[[1]]
      
      outdir     <- path("data/processed/python_input", prefix, study)
      dir_create(outdir, recurse = TRUE)
      out_file   <- path(outdir, paste0(file_name, ".parquet"))
      
      write_parquet(df, out_file, compression = "snappy")
      message("    – wrote ", out_file)
    }
  }
}

# point to your actual .qs files
export_cleaned("data/processed/cleaned_TMT10.qs", "TMT10")
export_cleaned("data/processed/cleaned_TMT11.qs", "TMT11")
message("✅ All runs exported to data/processed/python_input/")
