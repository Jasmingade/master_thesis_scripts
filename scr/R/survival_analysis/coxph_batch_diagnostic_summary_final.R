# ========================================================
# CoxPH + Batch Effect Diagnostics Across All PDC Studies
# With TMT11/10 support and batch count summaries
# ========================================================

library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(patchwork)
library(rlang)

# === Settings ===
pdc_ids <- c("PDC000110", "PDC000116", "PDC000120", "PDC000125", "PDC000127", "PDC000153", 
             "PDC000204", "PDC000221", "PDC000234", "PDC000270")
data_types <- c("gene", "iso_log", "iso_frac")
base_dir <- "data"
# For pre-correction
#input_dir <- file.path(base_dir, "processed", "proteomics")
# For post-correction
input_dir <- file.path(base_dir, "processed", "proteomics_batch_corrected")
# For pre-correction
#output_base_dir <- "data/analysis_results/figures/batch_effect_testing_pre_batch_correction"
# For post-correction
output_base_dir <- "data/analysis_results/figures/batch_effect_testing_post_batch_correction"
clinical_file <- file.path(base_dir, "processed", "all_cancers_clinical_core.csv")

# Create summary folder
dir.create(file.path(output_base_dir, "summary"), recursive = TRUE, showWarnings = FALSE)

# === Function: Flag problematic batch effect ===
is_problematic <- function(anova_csv_path) {
  if (!file.exists(anova_csv_path)) return(FALSE)
  df <- read_csv(anova_csv_path, show_col_types = FALSE)
  
  # Ensure correct order, just in case
  df <- df[order(as.numeric(gsub("PC", "", df$PC))), ]
  
  # Flag 1: sum of explained variance by PCs 1–5 with p < 0.01
  n_pcs_to_sum <- 5
  flag1_var_sum <- sum(df$PctVar[1:n_pcs_to_sum][df$Pval[1:n_pcs_to_sum] < 0.01])
  
  # Flag 2: number of first 10 PCs with p < 0.01
  n_pcs_to_check <- 10
  n_signif_pcs <- sum(df$Pval[1:n_pcs_to_check] < 0.01)
  
  # Final decision (adjust thresholds as needed)
  return(flag1_var_sum > 0.15 | n_signif_pcs >= 5)
}

# === Load clinical data ===
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
    age_group = cut(age, breaks = c(0, 40, 50, 60, 70, 80, Inf), right = FALSE,
                    labels = c("<40", "40–49", "50–59", "60–69", "70–79", "80+"))
  )

# === Function: Analyze Batch Effects + Diagnostics ===
analyze_batch_effects <- function(expr_file, label, title_prefix, clin_df_all, pdc_id, type) {
  expr_df <- read_csv(expr_file, show_col_types = FALSE)
  sample_names <- colnames(expr_df)[-1]
  batch_labels <- str_extract(sample_names, "^\\d{2}CPTAC")
  batch_labels <- ifelse(is.na(batch_labels), "Unknown", batch_labels)
  
  # Create output directory per study and data type
  study_dir <- file.path(output_base_dir, pdc_id, type)
  dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save batch counts
  batch_count_df <- data.frame(batch = batch_labels) %>%
    group_by(batch) %>%
    summarise(n_samples = n()) %>%
    arrange(desc(n_samples))
  write_csv(batch_count_df, file.path(study_dir, paste0("batch_counts_", label, ".csv")))
  
  # Matrix setup
  expr_matrix <- t(as.matrix(expr_df[, -1]))
  colnames(expr_matrix) <- expr_df$feature
  rownames(expr_matrix) <- sample_names
  
  # Clean
  expr_matrix <- expr_matrix[, colSums(is.na(expr_matrix)) < 0.5 * nrow(expr_matrix)]
  expr_matrix <- expr_matrix[rowSums(is.na(expr_matrix)) < 0.3 * ncol(expr_matrix), ]
  expr_matrix <- apply(expr_matrix, 2, function(x) { x[is.na(x)] <- median(x, na.rm = TRUE); x })
  expr_matrix <- expr_matrix[, apply(expr_matrix, 2, var) > 0]
  
  n_features <- ncol(expr_matrix)
  sample_case_ids <- str_extract(sample_names, "(?<=_matrix_)([0-9]{2}[A-Z]{2}[0-9]{3}|C3[NL]-[0-9]{5}|NX[0-9]+)")
  matched_case_ids <- intersect(sample_case_ids, clin_df_all$case_id)
  if (length(matched_case_ids) < 25) return(NULL)
  n_samples <- length(matched_case_ids)
  n_lowvar <- sum(apply(expr_matrix, 2, var) < 1e-5)
  pct_lowvar <- round(100 * n_lowvar / n_features, 2)
  
  # PCA
  pca_res <- prcomp(expr_matrix, center = TRUE, scale. = TRUE)
  pca_df <- as.data.frame(pca_res$x)[, 1:15]
  pca_df$batch <- batch_labels
  pc_var <- summary(pca_res)$importance[2, ]
  
  # Scree plot
  png(file.path(study_dir, paste0("scree_plot_", label, ".png")), width = 800, height = 600)
  plot(pc_var[1:20], type = "b", ylab = "Variance Explained", xlab = "Principal Component",
       main = paste(title_prefix, "Scree Plot"))
  dev.off()
  
  # PCA grid
  grid <- expand.grid(x = 1:5, y = 1:5) %>% filter(x < y)
  get_pc_label <- function(pc_num) paste0("PC", pc_num, " (", round(100 * pc_var[pc_num], 1), "%)")
  
  plots <- lapply(1:nrow(grid), function(i) {
    x <- grid$x[i]; y <- grid$y[i]
    ggplot(pca_df, aes(x = !!sym(paste0("PC", x)), y = !!sym(paste0("PC", y)), color = batch)) +
      geom_point(size = 1.2, alpha = 0.7) +
      labs(title = paste0("PC", x, " vs PC", y),
           x = get_pc_label(x), y = get_pc_label(y)) +
      theme_minimal() +
      theme(plot.title = element_text(size = 10),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 6),
            legend.position = ifelse(i == 1, "right", "none"))
  })
  
  ggsave(file.path(study_dir, paste0("pca_grid_", label, ".png")),
         wrap_plots(plots, ncol = 4, guides = "collect"), width = 12, height = 9, dpi = 300)
  
  # ANOVA
  # ANOVA (already present)
  anova_summary <- data.frame(
    PC = paste0("PC", 1:15),
    PctVar = pc_var[1:15],
    Pval = sapply(1:15, function(i) {
      summary(aov(pca_res$x[, i] ~ batch_labels))[[1]]$`Pr(>F)`[1]
    })
  )
  anova_path <- file.path(study_dir, paste0("anova_summary_", label, ".csv"))
  write_csv(anova_summary, anova_path)
  
  # Extract PC1 stats
  PC1Var <- pc_var[1]     # Variance explained by PC1 (as decimal, e.g. 0.12)
  PC1Pval <- anova_summary$Pval[1]  # ANOVA p-value for PC1
  
  # Boxplots
  pca_long <- pivot_longer(pca_df, starts_with("PC"), names_to = "PC", values_to = "value")
  pca_long$PC <- factor(pca_long$PC, levels = paste0("PC", 1:15))
  
  ggsave(file.path(study_dir, paste0("pc_boxplots_", label, ".png")),
         ggplot(pca_long, aes(x = batch, y = value, fill = batch)) +
           geom_boxplot() +
           facet_wrap(~PC, scales = "free_y", ncol = 5) +
           theme_minimal() +
           theme(axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.x = element_blank(),
                 legend.position = "right") +
           ggtitle(paste(title_prefix, "PC Scores by Batch")),
         width = 12, height = 7, dpi = 300)
  
  # Return diagnostic summary
  censoring <- table(clin_df_all$OS_event[clin_df_all$case_id %in% colnames(expr_matrix)])
  list(
    label = label,
    pdc_id = pdc_id,
    type = type,
    n_features = n_features,
    n_samples = n_samples,
    n_lowvar_features = n_lowvar,
    pct_lowvar = pct_lowvar,
    PC1Var = PC1Var,
    PC1Pval = PC1Pval,
    batches = toString(table(batch_labels)),
    censoring = toString(censoring),
    batch_effect_flag = is_problematic(anova_path)
    
  )
}

# === Run all combinations ===
summary_list <- list()

for (pdc_id in pdc_ids) {
  for (type in data_types) {
    # For pre-correction
    #pattern <- paste0("TMT\\d+_", type, "_", pdc_id, "_combined\\.csv$")
    # For post-correction
    pattern <- paste0("TMT\\d+_", type, "_", pdc_id, "_batch_corrected\\.csv$")
    files <- list.files(file.path(input_dir, type), pattern = pattern, full.names = TRUE)
    if (length(files) == 0) next
    expr_file <- files[1]
    label <- paste0(type, "_", pdc_id)
    title <- paste0(toupper(type), " – ", pdc_id)
    cat("Processing:", label, "\n")
    res <- analyze_batch_effects(expr_file, label, title, clin_df_all, pdc_id, type)
    summary_list[[label]] <- res
  }
}

summary_df <- bind_rows(summary_list)
write_csv(summary_df, file.path(output_base_dir, "summary", "batch_coxph_summary_all.csv"))
print(summary_df)
