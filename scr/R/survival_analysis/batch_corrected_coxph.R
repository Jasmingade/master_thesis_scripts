# ---- Libraries ----
library(readr)
library(dplyr)
library(stringr)
library(survival)
library(tibble)
library(ggplot2)
library(patchwork)
library(qvalue)
library(pbapply)

# ---- Settings ----
pdc_ids <- c("PDC000153")
data_types <- c("gene", "iso_log", "iso_frac")
clinical_file <- "data/processed/clin_df_all.csv"
input_base_dir <- "data/analysis_results/probatch_diagnostics/batch_correction"
output_base_dir <- "data/analysis_results//coxph"

clin_df_all <- read_csv(clinical_file, show_col_types = FALSE)

for (pdc_id in pdc_ids) {
  cat("CoxPH for", pdc_id, "\n")
  p_hist_list <- list()  # To hold ggplot objects for patchwork
  for (type in data_types) {
    # ---- Load batch corrected matrix ----
    infile <- file.path(input_base_dir, paste0(pdc_id, "_", type, "_batch_corrected.csv"))
    if (!file.exists(infile)) {
      warning("No batch corrected file: ", infile)
      next
    }
    batch_matrix <- as.matrix(read.csv(infile, row.names = 1, check.names = FALSE))
    sample_names <- colnames(batch_matrix)
    cat("Dimensions of batch corrected input", dim(batch_matrix), "\n")
    cat("Length of sample names", length(sample_names), "\n")
    
    # ---- Prepare sample annotation ----
    case_ids <- stringr::str_extract(sample_names, "(?<=_matrix_)([0-9]{2}[A-Z]{2}[0-9]{3}|C3[NL]-[0-9]{5}|NX[0-9]+)")
    sample_annotation <- tibble(
      FullRunName = sample_names,
      case_id = case_ids
    ) %>%
      left_join(
        clin_df_all %>%
          select(case_id, age, sex, tumor_stage, race, tumor_code, histologic_subtype, histologic_grade, age_group, OS_time, OS_event),
        by = "case_id"
      )
    
    # tumor_stage as ordered factor
    sample_annotation <- sample_annotation %>%
      mutate(tumor_stage_clean = factor(str_extract(tumor_stage, "Stage [IVX]+"),
                                        levels = c("Stage I", "Stage II", "Stage III", "Stage IV"),
                                        ordered = TRUE))
    
    # Filter: must have OS_time, OS_event, and match
    valid <- !is.na(sample_annotation$OS_time) & !is.na(sample_annotation$OS_event)
    batch_matrix <- batch_matrix[, valid, drop = FALSE]
    sample_annotation <- sample_annotation[valid, ]
    cat("Valid for survival analysis", type, ":", ncol(batch_matrix), "samples; ", nrow(batch_matrix), "features\n")
    cat("Number of case id's in sample annotion", length(sample_annotation$FullRunName), "\n")

    
    # ---- Run CoxPH for each feature ----
    results <- pblapply(seq_len(nrow(batch_matrix)), function(i) {
      feature <- rownames(batch_matrix)[i]
      expr <- as.numeric(batch_matrix[i, ])
      
      df <- data.frame(
        expr = expr,
        time = sample_annotation$OS_time,
        event = sample_annotation$OS_event,
        stage = sample_annotation$tumor_stage_clean,
        age = sample_annotation$age_group
      )
      
      df <- df[complete.cases(df), ]
      # Skip if not enough samples or no variance
      if (nrow(df) < 25) return(NULL)
      if (length(unique(df$expr)) < 3) return(NULL)
      
      fit <- try(coxph(Surv(time, event) ~ expr + stage, data = df), silent = TRUE)
      if (inherits(fit, "try-error")) return(NULL)
      s <- summary(fit)
      data.frame(
        feature = feature,
        coef = s$coefficients["expr", "coef"],
        HR = s$coefficients["expr", "exp(coef)"],
        se = s$coefficients["expr", "se(coef)"],
        z = s$coefficients["expr", "z"],
        p = s$coefficients["expr", "Pr(>|z|)"],
        n = s$n
      )
    })
    results <- bind_rows(results)
    if (nrow(results) > 0 && "p" %in% colnames(results)) {
      results$adj_p <- qvalue(results$p)$qvalues
    }
    
    # ---- Save results ----
    cox_outfile <- file.path(output_base_dir, paste0(pdc_id, "_", type, "_coxph.csv"))
    write.csv(results, cox_outfile, row.names = FALSE)
    cat("  Saved: ", cox_outfile, "\n")
    
    # ---- P-value histogram plot ----
    if (nrow(results) > 0 && "p" %in% colnames(results)) {
      pvalues <- results$p
      adj_p <- results$adj_p
      qobj <- qvalue(p = pvalues)
      results$adj_p <- qobj$qvalues
      
      summary_path <- file.path(dirname(cox_outfile), paste0("qval_summary_", type, "_", pdc_id, ".txt"))
      capture.output(summary(qobj), file = summary_path)
      
      
      p_hist <- hist(pvalues,
                     breaks = 25,
                     col = "grey80",
                     border = "black",
                     main = paste("P-value distribution for", pdc_id, "-", type),
                     xlab = "p-value",
                     ylab = "Frequency")
      abline(v = 0.05, col = "red", lty = 2, lwd = 2)
      
      
      q_hist <- plot(qobj)
    }
  }
  
  
  
  # ---- Combine and save grid for this pdc_id ----
  if (length(p_hist_list) > 0) {
    # patchwork: gene / iso_log / iso_frac vertically
    hist_grid <- (p_hist_list[["gene"]] / p_hist_list[["iso_log"]] / p_hist_list[["iso_frac"]]) +
      plot_annotation(title = paste0("CoxPH p-value histograms: ", pdc_id))
    png(file.path(output_base_dir, paste0(pdc_id, "_coxph_pvalue_histograms.png")),
        width = 1000, height = 1800, res = 150)
    print(hist_grid)
    dev.off()
    cat("  Saved p-value histogram grid for", pdc_id, "\n")
  }
}





