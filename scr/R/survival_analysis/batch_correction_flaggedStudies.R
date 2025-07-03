library(readr)
library(dplyr)
library(sva)
library(stringr)

# === USER SETTINGS ===
summary_file <- "data/analysis_results/figures/batch_effect_testing_pre_batch_correction/summary/batch_coxph_summary_all.csv"
data_dir <- "data/processed/proteomics"
output_dir <- "data/processed/proteomics_batch_corrected"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Thresholds ----
PC1_VAR_THRESHOLD <- 0.10
PC1_P_THRESHOLD   <- 0.01
MIN_SAMPLES       <- 30
MIN_FEATURES      <- 500

summary_df <- read_csv(summary_file, show_col_types = FALSE)
flagged <- summary_df %>% filter(batch_effect_flag == TRUE)

n_flagged <- nrow(flagged)
pb <- txtProgressBar(min = 0, max = n_flagged, style = 3)

corrected_files <- c()
skipped_files   <- c()

for (i in seq_len(n_flagged)) {
  setTxtProgressBar(pb, i)  # Progress bar FIRST
  
  label <- flagged$label[i]
  pdc_id <- flagged$pdc_id[i]
  type <- flagged$type[i]
  n_samples <- as.numeric(flagged$n_samples[i])
  n_features <- as.numeric(flagged$n_features[i])
  batch_count_vector <- as.numeric(unlist(strsplit(gsub(" ", "", flagged$batches[i]), ",")))
  pc1_var <- if ("PC1Var" %in% names(flagged)) as.numeric(flagged$PC1Var[i]) else NA
  pc1_p   <- if ("PC1Pval" %in% names(flagged)) as.numeric(flagged$PC1Pval[i]) else NA
  
  input_pattern <- paste0("TMT\\d+_", type, "_", pdc_id, "_combined\\.csv$")
  files <- list.files(file.path(data_dir, type), pattern = input_pattern, full.names = TRUE)
  if (length(files) == 0) {
    cat("No file found for:", label, "\n")
    skipped_files <- c(skipped_files, label)
    next
  }
  input_file <- files[1]
  file_name <- basename(input_file)
  output_subdir <- file.path(output_dir, type)
  dir.create(output_subdir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(output_subdir, gsub("_combined\\.csv$", "_batch_corrected.csv", file_name))
  
  cat("\nProcessing:", label, "|", file_name, "\n")
  expr_df <- read_csv(input_file, show_col_types = FALSE)
  sample_names <- colnames(expr_df)[-1]
  batch <- str_extract(sample_names, "^\\d{2}CPTAC")
  batch <- ifelse(is.na(batch), "Unknown", batch)
  batch <- factor(batch)
  expr_matrix <- expr_df[, -1] %>% as.matrix()
  rownames(expr_matrix) <- expr_df$feature
  valid_samples <- which(batch != "Unknown")
  expr_matrix <- expr_matrix[, valid_samples, drop = FALSE]
  batch <- droplevels(batch[valid_samples])
  sample_names <- sample_names[valid_samples]
  n_batches <- length(unique(batch))
  batch_sizes <- table(batch)
  is_constant <- apply(expr_matrix, 1, function(x) var(x, na.rm = TRUE) == 0)
  if (any(is_constant)) expr_matrix <- expr_matrix[!is_constant, , drop = FALSE]
  has_na_inf <- apply(expr_matrix, 1, function(x) any(is.na(x) | is.infinite(x)))
  if (any(has_na_inf)) expr_matrix <- expr_matrix[!has_na_inf, , drop = FALSE]
  batch_sizes <- table(batch)
  
  # --- Checks and individualized logic ---
  skip_batch_correction <- FALSE
  combat_par.prior <- TRUE
  
  if (any(batch_sizes < 2) || nrow(expr_matrix) < 10) {
    warning("Too few samples/features for reliable batch correction. Skipping batch correction.")
    write_csv(as.data.frame(expr_df), output_file)
    skipped_files <- c(skipped_files, output_file)
    next
  }
  
  if (n_batches < 2) {
    message("Only one batch detected. Skipping batch correction.")
    skip_batch_correction <- TRUE
  } else if (any(batch_count_vector < 5)) {
    cat("Using non-parametric ComBat: at least one batch has <5 samples\n")
    combat_par.prior <- FALSE
  }
  if (n_samples < MIN_SAMPLES || n_features < MIN_FEATURES) {
    warning("Low n_samples or n_features: batch correction may be unreliable.")
  }
  if (!is.na(pc1_var) && !is.na(pc1_p)) {
    if (pc1_var < PC1_VAR_THRESHOLD || pc1_p > PC1_P_THRESHOLD) {
      cat("PC1 variance or p-value below batch effect threshold. Skipping batch correction.\n")
      skip_batch_correction <- TRUE
    } else {
      cat("Strong batch effect: PC1 Var =", pc1_var, ", PC1 p =", pc1_p, "\n")
    }
  }
  
  # --- Apply or Skip ---
  if (skip_batch_correction) {
    write_csv(as.data.frame(expr_df), output_file)
    skipped_files <- c(skipped_files, output_file)
    cat("Skipped batch correction for:", output_file, "\n")
    next
  } else {
    result <- try({
      corrected <- ComBat(dat = expr_matrix, batch = batch, par.prior = combat_par.prior, prior.plots = FALSE)
      corrected_df <- cbind(feature = rownames(corrected), as.data.frame(corrected))
      write_csv(corrected_df, output_file)
      corrected_files <- c(corrected_files, output_file)
      cat("Batch-corrected file:", output_file, "\n")
    }, silent = TRUE)
    if (inherits(result, "try-error")) {
      cat("ComBat failed for", label, ":", conditionMessage(result), "\n")
      write_csv(as.data.frame(expr_df), output_file)
      skipped_files <- c(skipped_files, output_file)
    }
  }
  rm(expr_matrix); gc()
}

close(pb)
cat("\nBatch correction complete.\n")
cat("\n--- Batch-corrected files ---\n")
print(corrected_files)
cat("\n--- Skipped files (saved as-is) ---\n")
print(skipped_files)

# Optionally write to files
#writeLines(corrected_files, file.path(output_dir, "batch_corrected_files.txt"))
#writeLines(skipped_files,   file.path(output_dir, "skipped_files.txt"))
