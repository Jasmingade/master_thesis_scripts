library(Coxmos)
library(readr)
library(dplyr)
library(stringr)
library(tibble)
library(data.table)
library(readxl)
library(purrr)
library(impute)
library(tidyr)



# Setup
data_types <- c("gene", "iso_log", "iso_frac")
clinical_file <- "data/processed/metaData_multiomics.csv"
input_dir <- "data/analysis_results/probatch_diagnostics/batch_correction"
output_dir <- "data/coxmos_results"

# Load multiomics clinical metadata
clin_df_all <- read_csv(clinical_file, show_col_types = FALSE)



# 1) find all of your batch-corrected files
files <- list.files(input_dir,
                    pattern = "^PDC[0-9]+_(gene|iso_log|iso_frac)_batch_corrected\\.csv$",
                    full.names = TRUE)

# 2) read each into a matrix
omics_blocks <- lapply(files, function(f) {
  df <- fread(f, data.table = FALSE)
  rownames(df) <- df[[1]]
  df <- df[ , -1, drop = FALSE]
  as.matrix(df)
})

# 3) build the names by extracting PDC and data_type
pdc_ids    <- str_extract(basename(files), "^PDC[0-9]+")
data_types <- str_extract(basename(files), "(gene|iso_log|iso_frac)")
names(omics_blocks) <- paste0(pdc_ids, "_", data_types)



for (type in data_types) {
  cat("\n=== Running Coxmos for:", type, "===\n")
  
  # Filter omics_blocks for the current data_type
  selected_blocks <- omics_blocks[str_detect(names(omics_blocks), paste0("_", type, "$"))]
  
  # Preprocess each block individually
  selected_blocks <- lapply(selected_blocks, function(X_mat) {
    sample_ids <- colnames(X_mat)
    case_ids <- str_extract(sample_ids, "(?<=_matrix_)([0-9]{2}[A-Z]{2}[0-9]{3}|C3[NL]-[0-9]{5}|NX[0-9]+)")
    colnames(X_mat) <- case_ids
    X_mat <- X_mat[complete.cases(X_mat), ]
    X_mat
  })
  
  # Extract case IDs from sample names (e.g., "TCGA-BH-A18S-01A|PDC001" → "TCGA-BH-A18S-01A")
  sample_ids <- names(X_mat)
  case_ids <- str_extract(sample_ids, "(?<=_matrix_)([0-9]{2}[A-Z]{2}[0-9]{3}|C3[NL]-[0-9]{5}|NX[0-9]+)")
  colnames(X_mat) <- case_ids
  X_mat <- t(X_mat)
  X_mat <- X_mat[complete.cases(X_mat), ]
  
  # Align with clinical data
  sample_annotation <- tibble(
    case_id = rownames(X_mat),
    FullRunName = sample_ids,
    batch = str_extract(sample_ids, "^\\d{2}CPTAC")) %>%
    left_join(clin_df_all, by = "case_id") %>%
    filter(!is.na(OS_time) & !is.na(OS_event))
  
  # Check for duplicates
  dups <- duplicated(sample_annotation$case_id)
  
  # Inspect duplicates
  sample_annotation[dups, ]
  
  # Remove duplicates
  sample_annotation <- sample_annotation[!dups, ]
  
  X_mat <- X_mat[rownames(X_mat) %in% sample_annotation$case_id, ]
  X_mat <- X_mat[!dups,]
  sample_annotation <- sample_annotation %>% filter(case_id %in% rownames(X_mat))
  sample_annotation <- sample_annotation[match(rownames(X_mat), sample_annotation$case_id), ]
  
  Y <- as.data.frame(sample_annotation) %>% select(time = OS_time, event = OS_event)
  rownames(Y) <- sample_annotation$case_id
  
  # Optionally include clinical covariates
  X_clinical <- sample_annotation %>%
    select(case_id, age_group, sex, tumor_stage_clean, race, histologic_grade) %>%
    column_to_rownames("case_id") %>%
    factorToBinary(all = TRUE, sep = "_")
  
  X_blocks <- list(omics = X_mat, clinical = X_clinical)
  X_blocks$clinical <- factorToBinary(X = X_blocks$clinical, all = TRUE, sep = "_")

  # Clean rownames of X_mat (omics block)
  rownames(X_blocks$omics) <- rownames(X_blocks$omics) %>%
    str_remove("^X") %>%
    str_replace_all("\\.", "-")
  
  # EPV Check:
  EPV <- getEPV.mb(X_blocks, Y)
  print(EPV)
  
  # Train/test split
  split_data <- getTrainTest(X_blocks, Y, p = 0.7, seed = 1234)
  X_train <- split_data$X_train
  Y_train <- split_data$Y_train
  X_test <- split_data$X_test
  Y_test <- split_data$Y_test
  
  # Fit Coxmos models
  # ---- SB.sPLS-ICOX ----
  cv.sb.splsicox_res <- cv.mb.coxmos(method = "sb.splsicox",
                                     X = X_train, Y = Y_train,
                                     max.ncomp = 2, penalty.list = c(0.5, 0.9),
                                     n_run = 1, k_folds = 5)
  
  sb.splsicox_model <- mb.coxmos(method = "sb.splsicox",
                                 X = X_train, Y = Y_train,
                                 n.comp = cv.sb.splsicox_res$opt.comp,
                                 penalty = cv.sb.splsicox_res$opt.penalty,
                                 remove_non_significant = TRUE)
  
  # ---- SB.sPLS-DRCOX ----
  cv.sb.splsdrcox_res <- cv.mb.coxmos(method = "sb.splsdrcox",
                                      X = X_train, Y = Y_train,
                                      max.ncomp = 2, vector = NULL,
                                      n_run = 1, k_folds = 5)
  
  sb.splsdrcox_model <- mb.coxmos(method = "sb.splsdrcox",
                                  X = X_train, Y = Y_train,
                                  n.comp = cv.sb.splsdrcox_res$opt.comp,
                                  vector = cv.sb.splsdrcox_res$opt.nvar,
                                  remove_non_significant = TRUE)
  
  # ---- SB.sPLS-DACOX ----
  cv.sb.splsdacox_res <- cv.mb.coxmos(method = "sb.splsdacox",
                                      X = X_train, Y = Y_train,
                                      max.ncomp = 2, vector = NULL,
                                      n_run = 1, k_folds = 5)
  
  sb.splsdacox_model <- mb.coxmos(method = "sb.splsdacox",
                                  X = X_train, Y = Y_train,
                                  n.comp = cv.sb.splsdacox_res$opt.comp,
                                  vector = cv.sb.splsdacox_res$opt.nvar,
                                  remove_non_significant = TRUE)
  
  # ---- MB.sPLS-DRCOX ----
  cv.mb.splsdrcox_res <- cv.mb.coxmos(method = "mb.splsdrcox",
                                      X = X_train, Y = Y_train,
                                      max.ncomp = 2, vector = NULL,
                                      MIN_NVAR = 10, MAX_NVAR = NULL, n.cut_points = 10,
                                      EVAL_METHOD = "AUC", n_run = 1, k_folds = 5)
  
  mb.splsdrcox_model <- mb.coxmos(method = "mb.splsdrcox",
                                  X = X_train, Y = Y_train,
                                  n.comp = cv.mb.splsdrcox_res$opt.comp,
                                  vector = cv.mb.splsdrcox_res$opt.nvar)
  
  # ---- MB.sPLS-DACOX ----
  cv.mb.splsdacox_res <- cv.mb.coxmos(method = "mb.splsdacox",
                                      X = X_train, Y = Y_train,
                                      max.ncomp = 2, vector = NULL,
                                      n_run = 1, k_folds = 5)
  
  mb.splsdacox_model <- mb.coxmos(method = "mb.splsdacox",
                                  X = X_train, Y = Y_train,
                                  n.comp = cv.mb.splsdacox_res$opt.comp,
                                  vector = cv.mb.splsdacox_res$opt.nvar)
  
  # ---- Evaluation ----
  lst_models <- list(
    "SB.sPLS-ICOX" = sb.splsicox_model,
    "SB.sPLS-DRCOX" = sb.splsdrcox_model,
    "SB.sPLS-DACOX" = sb.splsdacox_model,
    "MB.sPLS-DRCOX" = mb.splsdrcox_model,
    "MB.sPLS-DACOX" = mb.splsdacox_model
  )
  
  eval_results <- eval_Coxmos_models(lst_models = lst_models,
                                     X_test = X_test, Y_test = Y_test,
                                     times = NULL)
  
  plot_list <- plot_evaluation(eval_results, evaluation = "AUC", pred.attr = "mean")
  print(plot_list$lst_plots$lineplot.mean)
  
}