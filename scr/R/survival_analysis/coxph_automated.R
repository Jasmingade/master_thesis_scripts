# =====================================================
# Efficient CoxPH Survival Analysis for Multiple PDC Studies
# =====================================================

# Load required packages
library(survival)
library(dplyr)
library(readr)
library(stringr)
library(pbapply)
library(ggplot2)
library(qvalue)
library(survminer)
library(sva)
library(grid)
library(gridGraphics)

# === USER SETTINGS ===
pdc_ids <- c("PDC000110", "PDC000116", "PDC000120", "PDC000125", "PDC000127", "PDC000153", 
             "PDC000204", "PDC000221", "PDC000234", "PDC000270")
data_types <- c("iso_log", "iso_frac", "gene")
base_dir <- "data"
input_dir <- file.path(base_dir, "processed", "proteomics_batch_corrected")
clinical_file <- file.path(base_dir, "processed", "all_cancers_clinical_core.csv")
log_base_dir <- file.path(base_dir, "analysis_results", "proteomics", "log_files")
dir.create(log_base_dir, recursive = TRUE, showWarnings = FALSE)

# === Load and clean clinical data ===
clin_df_all <- read_csv(clinical_file, show_col_types = FALSE) %>%
  rename(
    case_id = Patient_ID,
    age = `consent/age`,
    sex = `consent/sex`,
    race = `consent/race`,
    ethnicity = `consent/ethnicity`,
    tumor_size_cm = `baseline/tumor_size_cm`,
    bmi = `medical_history/bmi`,
    alcohol_consumption = `medical_history/alcohol_consumption`,
    tobacco_smoking_history = `medical_history/tobacco_smoking_history`,
    vital_status = `follow-up/vital_status_at_date_of_last_contact`,
    days_to_last_contact = `follow-up/number_of_days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact`,
    days_to_death = `follow-up/number_of_days_from_date_of_initial_pathologic_diagnosis_to_date_of_death`,
    tumor_stage = `baseline/tumor_stage_pathological`,
    histologic_subtype = `baseline/histologic_type`,
    histologic_grade = `cptac_path/histologic_grade`,
    OS_time = `Overall survival, days`,
    OS_event = `Survival status (1, dead; 0, alive)`
  ) %>%
  mutate(
    age = suppressWarnings(as.numeric(age)),
    age_log2 = ifelse(!is.na(age) & age > 0, log2(age), NA),
    age_group = cut(age,
                    breaks = c(0, 40, 50, 60, 70, 80, Inf),
                    right = FALSE,
                    labels = c("<40", "40–49", "50–59", "60–69", "70–79", "80+"))
  )
clin_df_all$case_id <- toupper(clin_df_all$case_id)

write_csv(clin_df_all, file.path(input_dir, paste0("clin_df_all", ".csv")))


# === Function: CoxPH for a single matrix ===
run_coxph <- function(expr_file, clin_df_all, out_file, log_file, pdc_id, tmt, type) {
  cat("Reading:", expr_file, "\n")
  expr_df <- read_csv(expr_file, show_col_types = FALSE)
  if (!"feature" %in% colnames(expr_df)) stop("Missing 'feature' column")
  
  sample_cols <- colnames(expr_df)[-1]
  batch <- str_extract(sample_cols, "^\\d{2}CPTAC")
  batch <- ifelse(is.na(batch), "Unknown", batch)
  batch <- factor(batch)
  
  expr_df <- expr_df %>% filter(!is.na(feature) & feature != "")
  expr_df <- as.data.frame(expr_df)
  rownames(expr_df) <- expr_df$feature
  expr_df <- expr_df[, !(colnames(expr_df) == "feature")]
  
  withdrawn_mask <- grepl("Withdrawn", sample_cols)
  expr_df <- expr_df[, !withdrawn_mask]
  batch <- batch[!withdrawn_mask]
  sample_cols <- sample_cols[!withdrawn_mask]
  
  sample_case_ids <- str_extract(sample_cols, "(?<=_matrix_)([0-9]{2}[A-Z]{2}[0-9]{3}|C3[NL]-[0-9]{5}|NX[0-9]+)")
  keep_cols <- !is.na(sample_case_ids)
  expr_df <- expr_df[, keep_cols]
  batch <- batch[keep_cols]
  sample_case_ids <- sample_case_ids[keep_cols]
  colnames(expr_df) <- sample_case_ids
  
  matched_case_ids <- intersect(sample_case_ids, clin_df_all$case_id)
  if (length(matched_case_ids) < 25) return(NULL)
  
  expr_df <- expr_df[, matched_case_ids]
  batch <- batch[match(matched_case_ids, sample_case_ids)]
  
  clin_df_filtered <- clin_df_all %>%
    filter(case_id %in% matched_case_ids) %>%
    mutate(tumor_stage_clean = factor(str_extract(tumor_stage, "Stage [IVX]+"),
                                      levels = c("Stage I", "Stage II", "Stage III", "Stage IV"),
                                      ordered = TRUE)) %>%
    arrange(match(case_id, matched_case_ids))
  expr_df <- expr_df[, clin_df_filtered$case_id]
  
  cat("Matched cases:", length(matched_case_ids), "\n")
  
  results <- pblapply(seq_len(nrow(expr_df)), function(i) {
    feature <- rownames(expr_df)[i]
    x <- as.numeric(expr_df[i, ])
    df <- data.frame(expr = x,
                     time = clin_df_filtered$OS_time,
                     event = clin_df_filtered$OS_event,
                     stage = clin_df_filtered$tumor_stage_clean,
                     age = clin_df_filtered$age_log2,
                     batch = batch)
    df <- df[complete.cases(df), ]
    if (nrow(df) < 25) return(NULL)
    fit <- try(coxph(Surv(time, event) ~ expr + stage, data = df), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    coef <- summary(fit)$coefficients["expr", ]
    data.frame(feature = feature,
               coef = coef["coef"],
               HR = coef["exp(coef)"],
               p = coef["Pr(>|z|)"])
  })
  
  results_df <- bind_rows(results)
  results_df <- results_df[!is.na(results_df$p), ]
  if (nrow(results_df) > 0) {
    results_df$adj_p <- qvalue(results_df$p)$qvalues
  } else {
    results_df$adj_p <- numeric(0)
  }
  
  if ("feature" %in% colnames(results_df)) {
    results_df <- results_df %>% dplyr::relocate(feature)
  }
  write_csv(results_df, out_file)
  
  # Plot p-value histogram
  if (nrow(results_df) > 0 && "p" %in% colnames(results_df)) {
    pvalues <- results_df$p
    qobj <- qvalue(p = pvalues)
    results_df$adj_p <- qobj$qvalues
    
    summary_path <- file.path(dirname(out_file), paste0("qval_summary_", type, "_", pdc_id, ".txt"))
    capture.output(summary(qobj), file = summary_path)
    
    png(filename = file.path(dirname(out_file), paste0("pval_hist_", type, "_", pdc_id, ".png")),
        width = 6, height = 5, units = "in", res = 300)
    hist(pvalues,
         breaks = 25,
         col = "grey80",
         border = "black",
         main = paste("P-value distribution for", pdc_id, "-", type),
         xlab = "p-value",
         ylab = "Frequency")
    abline(v = 0.05, col = "red", lty = 2, lwd = 2)
    dev.off()
    
    png(filename = file.path(dirname(out_file), paste0("qval_diagnostics_", type, "_", pdc_id, ".png")),
        width = 6, height = 5, units = "in", res = 300)
    plot(qobj)
    dev.off()
  }
  
  return(results_df)
}

# === Run CoxPH for all combinations ===
for (pdc_id in pdc_ids) {
  cat("\n===== Processing:", pdc_id, "=====\n")
  output_dir <- file.path(base_dir, "analysis_results", "proteomics", pdc_id)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (type in data_types) {
    # Find any TMT version for this type and pdc_id
    pattern <- paste0("TMT\\d+_", type, "_", pdc_id, "_batch_corrected\\.csv$")
    files <- list.files(file.path(input_dir, type), pattern = pattern, full.names = TRUE)
    if (length(files) == 0) {
      cat("Skipping missing file for:", type, pdc_id, "\n")
      next
    }
    for (expr_file in files) {
      tmt <- stringr::str_extract(basename(expr_file), "TMT\\d+")
      out_file <- file.path(output_dir, paste0("coxph_", type, "_", pdc_id, "_", tmt, ".csv"))
      log_file <- file.path(log_base_dir, paste0("coxph_", type, "_", pdc_id, "_", tmt, ".log"))
      
      sink(log_file, split = TRUE)
      cat("=== CoxPH LOG ===\n")
      cat("PDC Study:", pdc_id, "\n")
      cat("Data Type:", type, "\n")
      cat("TMT Set:", tmt, "\n")
      cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
      
      results_df <- run_coxph(expr_file, clin_df_all, out_file, log_file, pdc_id, tmt, type)
      
      if (!is.null(results_df)) {
        sig_count <- sum(results_df$p < 0.05, na.rm = TRUE)
        cat("Significant features (p < 0.05):", sig_count, "\n")
      } else {
        cat("No results returned.\n")
      }
      
      sink()
    }
  }
}
