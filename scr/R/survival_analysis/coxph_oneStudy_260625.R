# ================================================
# Efficient CoxPH Survival Analysis one PDC Study with Batch Correction
# ================================================

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
library(grid)         # Required for grid graphics
library(gridGraphics) # Needed to capture base plots as grid objects

# === USER SETTINGS ===
pdc_id <- "PDC000125"
base_dir <- "data"
input_dir <- file.path(base_dir, "processed", "proteomics")
clinical_file <- file.path(base_dir, "processed", "all_cancers_clinical_core.csv")
output_dir <- file.path(base_dir, "analysis_results", "proteomics", pdc_id)
log_dir <- file.path(base_dir, "analysis_results", "proteomics", "log_files")
data_types <- c("iso_log", "iso_frac", "gene")

# Create output and log directories if they don't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

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
    age = suppressWarnings(as.numeric(age)),  # suppress warning if you want
    age_log2 = ifelse(!is.na(age) & age > 0, log2(age), NA),
    age_group = cut(age,
                    breaks = c(0, 40, 50, 60, 70, 80, Inf),
                    right = FALSE,
                    labels = c("<40", "40–49", "50–59", "60–69", "70–79", "80+"))
  )



clin_df_all$case_id <- toupper(clin_df_all$case_id)

required_cols <- c("case_id", "OS_time", "OS_event")
if (!all(required_cols %in% colnames(clin_df_all))) {
  stop("Your clinical file must contain: case_id, OS_time, OS_event")
}

# === Function: CoxPH for a single matrix ===
run_coxph <- function(expr_file, clin_df_all, out_file, log_file, pdc_id, tmt, type) {
  cat("Reading:", expr_file, "\n")
  expr_df <- read_csv(expr_file, show_col_types = FALSE)
  if (!"feature" %in% colnames(expr_df)) {
    stop("Expression file must contain a column named 'feature'")
  }
  
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
  
  cat("Number of matched_case_ids:", length(matched_case_ids), "\n")
  cat("Number of features in expr_df:", nrow(expr_df), "\n")
  cat("Batch levels:", paste(levels(batch), collapse = ", "), "\n")
  
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
    fit <- try(coxph(Surv(time, event) ~ expr + stage + batch, data = df), silent = TRUE)
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
    warning("No valid p-values for qvalue calculation")
  }
  
  if ("feature" %in% colnames(results_df)) {
    results_df <- results_df %>% dplyr::relocate(feature)
  }
  write_csv(results_df, out_file)
  
  # Plotting of p-value histograms
  if (nrow(results_df) > 0 && "p" %in% colnames(results_df)) {
    pvalues <- results_df$p 
    qobj <- qvalue(p = pvalues)
    results_df$adj_p <- qobj$qvalues
    
    # Save qvalue summary to text file
    summary_path <- file.path(dirname(out_file), paste0("qval_summary_", type, "_", pdc_id, ".txt"))
    capture.output(summary(qobj), file = summary_path)
    
    # 1. Save histogram of raw p-values
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
    
    # 2. Save q-value diagnostic plot
    png(filename = file.path(dirname(out_file), paste0("qval_diagnostics_", type, "_", pdc_id, ".png")),
        width = 6, height = 5, units = "in", res = 300)
    plot(qobj)
    dev.off()
  }
  
  return(results_df)
}

# === Run Cox analysis per data type ===
for (type in data_types) {
  expr_file <- file.path(input_dir, type, paste0("TMT10_", type, "_", pdc_id, "_batch_corrected.csv"))
  out_file <- file.path(output_dir, paste0("coxph_", type, "_", pdc_id, ".csv"))
  log_file <- file.path(log_dir, paste0("coxph_", type, "_", pdc_id, ".log"))
  if (!file.exists(expr_file)) next
  
  sink(log_file, split = TRUE)
  cat("=== CoxPH LOG ===\n")
  cat("PDC Study:", pdc_id, "\n")
  cat("Data Type:", type, "\n")
  cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  results_df <- run_coxph(expr_file, clin_df_all, out_file, log_file, pdc_id, "TMT10", type)
  
  if (!is.null(results_df)) {
    sig_count <- sum(results_df$p < 0.05, na.rm = TRUE)
    cat("Significant features (p < 0.05):", sig_count, "\n")
  } else {
    cat("No results returned.\n")
  }
  
  sink()
}





