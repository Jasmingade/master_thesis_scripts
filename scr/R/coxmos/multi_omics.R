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
library(vroom)
library(mixOmics)
library(utils)

# –– 0) helper to save ggplots or lists of ggplots ––––––––––––––––––––––
save_plot <- function(plot_obj, filename, width=8, height=6, dpi=300) {
  # if a list, find the first ggplot
  if (is.list(plot_obj) && !inherits(plot_obj, "ggplot")) {
    plot_obj <- plot_obj[[ which(vapply(plot_obj, inherits, logical(1), "ggplot"))[1] ]]
  }
  ggsave(filename, plot = plot_obj, width=width, height=height, dpi=dpi)
}

# Setup
data_types <- c("gene", "iso_log", "iso_frac")
clinical_file <- "data/processed/clin_df_all_cleaned.csv"
input_dir <- "data/analysis_results/probatch_diagnostics/batch_correction"
transcriptomics_dir <- "data/processed/transcriptomics"
output_dir <- "data/coxmos_results"

# Load clinical metadata
clin_df_all <- read_csv(clinical_file, show_col_types = FALSE)


cat("\n=== Creating proteomics block ===\n")
# Load proteomics blocks
prot_files <- list.files(input_dir, pattern = "^PDC[0-9]+_(gene|iso_log|iso_frac)_batch_corrected\\.csv$", full.names = TRUE)

# Mapping PDC IDs to TCGA cancer type codes
pdc_to_tcga <- c(
  PDC000110 = "OV",
  PDC000116 = "COAD",
  PDC000120 = "BRCA",
  PDC000125 = "UCEC",
  PDC000127 = "KIRC",
  PDC000153 = "LUAD",
  PDC000204 = "GBM",
  PDC000221 = "HNSC",
  PDC000234 = "LUSC",
  PDC000270 = "PAAD"
)

# Load and name blocks with proteomics prefix
prot_blocks <- setNames(
  lapply(prot_files, function(f) {
    df <- fread(f, data.table = FALSE)
    rownames(df) <- df[[1]]
    df <- df[, -1, drop = FALSE]
    as.matrix(df)
  }),
  sapply(basename(prot_files), function(f) {
    pdc <- str_extract(f, "^PDC[0-9]+")
    type <- str_extract(f, "(gene|iso_log|iso_frac)")
    tcga <- pdc_to_tcga[[pdc]]
    if (is.null(tcga)) stop(paste("Unmapped PDC ID:", pdc))
    paste0("PROT_", tcga, "_", type)
  })
)

cat("\n=== Creating transcriptomic block ===\n")

# Select subset of cancer types to test
selected_ids <- c("LUAD", "BRCA", "OV")

# Load only relevant transcriptomics blocks
trans_files <- list.files(transcriptomics_dir, pattern = "^RNA_.*_(gene|iso_log|iso_frac)\\.csv$", full.names = TRUE)

# Filter files by selected cancer types
trans_files <- trans_files[
  str_detect(basename(trans_files), paste0("^RNA_(", paste(selected_ids, collapse = "|"), ")_"))
]

# Load filtered transcriptomics data
trans_blocks <- setNames(
  lapply(trans_files, function(f) {
    df <- fread(f, data.table = FALSE)
    rownames(df) <- df[[1]]
    df <- df[, -1, drop = FALSE]
    as.matrix(df)
  }),
  str_replace_all(basename(trans_files), c("^RNA_" = "TRANS_", "\\.csv$" = ""))
)

# Merge proteomics and transcriptomics
omics_blocks <- c(prot_blocks, trans_blocks)


all_block_evals <- list()

# Iterate through selected cancer types
for (cancer in selected_ids) {
  cat("\n=== Coxmos for", cancer, "===\n")
  
  # pick just the PROT_* and TRANS_* blocks for this cancer
  this_blocks <- omics_blocks[
    str_detect(names(omics_blocks), paste0("^(PROT|TRANS)_", cancer, "_"))
  ]
  
  # figure out how many total “jobs” we have:
  all_jobs <- unlist(
    lapply(selected_ids, function(cancer) {
      names(this_blocks)[str_detect(names(this_blocks),
                                     paste0("^(PROT|TRANS)_", cancer, "_"))]
    })
  )
  n_jobs <- length(all_jobs)
  
  # create a progress bar
  pb <- txtProgressBar(min = 0, max = n_jobs, style = 2)
  for(i in job_idx){Sys.sleep(0.5); setTxtProgressBar(pb, i)}
  Sys.sleep(1)
  job_idx <- 0
  
  # for each block, run Coxmos independently
  for (block_name in names(this_blocks)) {
    # advance progress bar
    job_idx <- job_idx + 1
    setTxtProgressBar(pb, job_idx)
    for (i in job_idx) { print(i); flush.console() }
    
    
    cat("\n→ [", job_idx, "/", n_jobs, "] Block:", block_name, "\n")
    
    cat("→ Block:", block_name, "\n")
    X_raw <- this_blocks[[block_name]]
    
    # transpose, filter columns with >80% missing, impute
    sample_ids <- colnames(X_raw)
    case_ids <- str_extract(sample_ids, "(?<=_matrix_)([0-9]{2}[A-Z]{2}[0-9]{3}|C3[NL]-[0-9]{5}|NX[0-9]+)")
    colnames(X_raw) <- case_ids
    X <- t(X_raw)
    X <- X[, colMeans(is.na(X)) < 0.8, drop = FALSE]
    X <- impute.knn(X)$data
    
    # assemble Y & clinical just for these samples
    samp_anno <- tibble(
      case_id = case_ids) %>%
      left_join(clin_df_all, by = "case_id") %>%
      dplyr::filter(!is.na(OS_time) & !is.na(OS_event)) %>%
      distinct(case_id, .keep_all = TRUE)
    
    # re‐subset X to those with clinical
    # Check for duplicates
    dups <- duplicated(samp_anno$case_id)

    
    keep <- samp_anno$case_id[!dups]
    X <- X[keep, , drop=FALSE]
    samp_anno <- samp_anno[!dups, ]
    
    # 1) Factor your columns
    samp_anno <- samp_anno %>%
      mutate(
        age_group  = factor(age_group,  levels = c("<40","40–49","50–59","60–69","70–79","80+")),
        sex        = factor(sex,        levels = c("female","male")),
        tumor_stage = factor(tumor_stage_clean, levels = c("stage i","stage ii","stage iii","stage iv", "unknown")),
        hist_grade = factor(histologic_grade,
                            levels = c("G1",
                                       "G2",
                                       "G3",
                                       "G4",
                                       "GX",
                                       "[unknown]"))
      )
    
    # 2) Select and one-hot encode
    X_clin_df <- samp_anno %>%
      dplyr::select(case_id, age_group, sex, tumor_stage, hist_grade) %>%
      column_to_rownames("case_id")
    
    X_clin_onehot <- factorToBinary(X = X_clin_df, all = TRUE, sep = "_")
    
    # 3) Force numeric
    X_clin_num <- X_clin_onehot %>%
      mutate(across(everything(), ~ as.integer(.))) %>%    # turn "0"/"1" → 0/1 integers
      as.matrix()
    
    # 4) Build your multiblock list
    X_blocks <- list(
      molecular = X,         # your imputed omics matrix
      clinical  = X_clin_num # now a true numeric matrix
    )
    
    X_blocks$clinical <- factorToBinary(X = X_blocks$clinical, all = TRUE, sep = "_")
    
    
    Y <- as.data.frame(samp_anno) %>% dplyr::select(time = OS_time, event = OS_event)
    rownames(Y) <- samp_anno$case_id
    
    
    # quick EPV
    cat(" EPV:\n")
    print(getEPV.mb(X_blocks, Y))
    
    cat("    • Splitting into train/test… \n")
    split_data <- getTrainTest(X_blocks, Y, p = 0.7, seed = 1234)
    X_train <- split_data$X_train
    Y_train <- split_data$Y_train
    X_test <- split_data$X_test
    Y_test <- split_data$Y_test
  
  
    # Fit Coxmos models
    # ---- SB.sPLS-ICOX ----
    cat("    • CV & fitting SB.sPLS-ICOX… \n")
    cv.sb.splsicox_res <- cv.mb.coxmos(method = "sb.splsicox",
                                       X = X_train,
                                       Y = Y_train,
                                       MIN_NVAR = 5, MAX_NVAR = 20,
                                       max.ncomp = 1,
                                       n_run = 2,
                                       k_folds = 5,
                                       MIN_EPV = 0.01,
                                       remove_zero_variance = TRUE,
                                       remove_near_zero_variance = TRUE,
                                       remove_variance_at_fold_level = TRUE,
                                       remove_non_significant = FALSE)
  
    sb.splsicox_model <- mb.coxmos(method = "sb.splsicox",
                                   X = X_train, Y = Y_train,
                                   n.comp = cv.sb.splsicox_res$opt.comp,
                                   penalty = cv.sb.splsicox_res$opt.penalty,
                                   remove_non_significant = TRUE)
    cat("      ✔︎ SB.sPLS-ICOX done \n.")
  
    # ---- SB.sPLS-DRCOX ----
    cat("    • CV & fitting SB.sPLS-DRCOX… \n")
    cv.sb.splsdrcox_res <- cv.mb.coxmos(method = "sb.splsdrcox",
                                        X = X_train,
                                        Y = Y_train,
                                        MIN_NVAR = 5, MAX_NVAR = 20,
                                        max.ncomp = 1,
                                        n_run = 2,
                                        k_folds = 5,
                                        MIN_EPV = 0.01,
                                        remove_zero_variance = TRUE,
                                        remove_near_zero_variance = TRUE,
                                        remove_variance_at_fold_level = TRUE,
                                        remove_non_significant = FALSE)
  
    sb.splsdrcox_model <- mb.coxmos(method = "sb.splsdrcox",
                                    X = X_train, Y = Y_train,
                                    n.comp = cv.sb.splsdrcox_res$opt.comp,
                                    vector = cv.sb.splsdrcox_res$opt.nvar,
                                    remove_non_significant = TRUE)
    cat("      ✔︎ SB.sPLS-DRCOX done\n.")
  
    # ---- SB.sPLS-DACOX ----
    cat("    • CV & fitting SB.sPLS-DACOX… \n")
    cv.sb.splsdacox_res <- cv.mb.coxmos(method = "sb.splsdacox",
                                        X = X_train,
                                        Y = Y_train,
                                        MIN_NVAR = 5, MAX_NVAR = 20,
                                        max.ncomp = 1,
                                        n_run = 2,
                                        k_folds = 5,
                                        MIN_EPV = 0.01,
                                        remove_zero_variance = TRUE,
                                        remove_near_zero_variance = TRUE,
                                        remove_variance_at_fold_level = TRUE,
                                        remove_non_significant = FALSE)
  
    sb.splsdacox_model <- mb.coxmos(method = "sb.splsdacox",
                                    X = X_train, Y = Y_train,
                                    n.comp = cv.sb.splsdacox_res$opt.comp,
                                    vector = cv.sb.splsdacox_res$opt.nvar,
                                    remove_non_significant = TRUE)
    cat("      ✔︎ SB.sPLS-DACOX done\n.")
  
    # ---- MB.sPLS-DRCOX ----
    cat("    • CV & fitting MB.sPLS-DRCOX…\n")
    cv.mb.splsdrcox_res <- cv.mb.coxmos(method = "mb.splsdrcox",
                                        X = X_train,
                                        Y = Y_train,
                                        MIN_NVAR = 5, MAX_NVAR = 20,
                                        max.ncomp = 1,
                                        n_run = 2,
                                        k_folds = 5,
                                        MIN_EPV = 0.01,
                                        remove_zero_variance = TRUE,
                                        remove_near_zero_variance = TRUE,
                                        remove_variance_at_fold_level = TRUE,
                                        remove_non_significant = FALSE)
  
    mb.splsdrcox_model <- mb.coxmos(method = "mb.splsdrcox",
                                    X = X_train, Y = Y_train,
                                    n.comp = cv.mb.splsdrcox_res$opt.comp,
                                    vector = cv.mb.splsdrcox_res$opt.nvar)
    cat("      ✔︎ MB.sPLS-DRCOX done\n.")
  
    # ---- MB.sPLS-DACOX ----
    cat("    • CV & fitting MB.sPLS-DACOX… \n")
    cv.mb.splsdacox_res <- cv.mb.coxmos(method = "mb.splsdacox",
                                        X = X_train,
                                        Y = Y_train,
                                        MIN_NVAR = 5, MAX_NVAR = 20,
                                        max.ncomp = 1,
                                        n_run = 2,
                                        k_folds = 5,
                                        MIN_EPV = 0.01,
                                        remove_zero_variance = TRUE,
                                        remove_near_zero_variance = TRUE,
                                        remove_variance_at_fold_level = TRUE,
                                        remove_non_significant = FALSE)
  
    mb.splsdacox_model <- mb.coxmos(method = "mb.splsdacox",
                                    X = X_train, Y = Y_train,
                                    n.comp = cv.mb.splsdacox_res$opt.comp,
                                    vector = cv.mb.splsdacox_res$opt.nvar)
    cat("      ✔︎ MB.sPLS-DACOX done\n.")
  
    # ---- Evaluation ----
    cat("    • Evaluating all models…\n")
    lst_models <- list("SB.sPLS-ICOX" = sb.splsicox_model,
                     #"iSB.sPLS-ICOX" = isb.splsicox_model,
                     "SB.sPLS-DRCOX-Dynamic" = sb.splsdrcox_model,
                     #"iSB.sPLS-DRCOX-Dynamic" = isb.splsdrcox_model,
                     "SB.sPLS-DACOX-Dynamic" = sb.splsdacox_model,
                     #"iSB.sPLS-DACOX-Dynamic" = isb.splsdacox_model,
                     "MB.sPLS-DRCOX" = mb.splsdrcox_model,
                     "MB.sPLS-DACOX" = mb.splsdacox_model)
  
    eval_results <- eval_Coxmos_models(lst_models = lst_models,
                                     X_test = X_test, Y_test = Y_test, 
                                     times = NULL,
                                     verbose = TRUE,
                                     PARALLEL = TRUE)
    
    # store it keyed by block_name
    eval_results <- eval_Coxmos_models(
      lst_models   = lst_models,
      X_test       = X_test,
      Y_test       = Y_test,
      verbose      = FALSE,
      progress_bar = FALSE
    )
    all_block_evals[[block_name]] <- eval_results
  
    # –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
    # now plot & save **all** of your validation diagnostics right here:
    dir.create(file.path(output_dir, cancer), showWarnings=FALSE, recursive=TRUE)
    odir <- file.path(output_dir, cancer)
  
    # 1) time-dependent AUC & IBS
    res_auc    <- plot_evaluation(eval_results, evaluation="AUC",  pred.attr="mean")
    save_plot(res_auc$lst_plots$lineplot.mean,   file.path(odir,"eval_AUC_times.png"))
    save_plot(res_auc$lst_plot_comparisons$anova, file.path(odir,"eval_AUC_ANOVA.png"))
  
    res_ibs <- plot_evaluation(eval_results, evaluation="IBS", pred.attr="mean")
    save_plot(res_ibs$lst_plots$lineplot.mean, file.path(odir,"eval_IBS_times.png"))
  
    # 2) computation times
    lst_time_models <- c(
      sb.splsicox_model,
      sb.splsdrcox_model,
      sb.splsdacox_model,
      mb.splsdrcox_model,
      mb.splsdacox_model,
      list(eval=eval_results)
      )
    ggp_time <- plot_time.list(lst_time_models, txt.x.angle=90)
    save_plot(ggp_time, file.path(odir,"model_timing.png"))
  
    # 3) PH diagnostics (list of plots)
    ph_list <- plot_proportionalHazard(sb.splsicox_model)
    for (nm in names(ph_list)) {
      save_plot(ph_list[[nm]], file.path(odir, paste0("PH_", nm, ".png")))
      }
  
    # 4) forest
    forest_list <- plot_forest(sb.splsicox_model)
    for (nm in names(forest_list)) {
      save_plot(forest_list[[nm]], file.path(odir, paste0("forest_", nm, ".png")))
      }
  
    # 5) linear-predictor density + histogram
    dlist <- plot_cox.event(sb.splsicox_model, type="lp")
    save_plot(dlist$plot.density,   file.path(odir,"lp_density.png"))
    save_plot(dlist$plot.histogram, file.path(odir,"lp_histogram.png"))
  
    # 6) per-variable AUC
    var_auc <- eval_Coxmos_model_per_variable(
      model  = sb.splsicox_model,
      X_test = sb.splsicox_model$X_input,
      Y_test = sb.splsicox_model$Y_input
      )
    pvar  <- plot_evaluation(var_auc, evaluation="AUC")$lst_plots$lineplot.mean
    save_plot(pvar, file.path(odir,"AUC_by_variable.png"))
  
    # 7) sPLS plotting
    sc <- plot_sPLS_Coxmos(sb.splsicox_model, comp=1:2, mode="scores")
    ld <- plot_sPLS_Coxmos(sb.splsicox_model, comp=1:2, mode="loadings", top=10)
    bp <- plot_sPLS_Coxmos(sb.splsicox_model, comp=1:2, mode="biplot", top=15,
                         only_top=TRUE, overlaps=20)
    save_plot(sc, file.path(odir,"spls_scores.png"))
    save_plot(ld, file.path(odir,"spls_loadings.png"))
    save_plot(bp, file.path(odir,"spls_biplot.png"))
  
    # 8) pseudobeta
    pb <- plot_pseudobeta(sb.splsicox_model,
                        error.bar=TRUE, onlySig=FALSE,
                        alpha=0.05, zero.rm=TRUE,
                        auto.limits=TRUE, top=20,
                        show_percentage=TRUE, size_percentage=2)
    save_plot(pb$plot,          file.path(odir,"pseudobeta_main.png"))
    save_plot(pb$mb_plot$plot,  file.path(odir,"pseudobeta_molecular.png"))
  
    # 9) Kaplan–Meier splits (LP + comps)
    km_lp   <- getAutoKM(type="LP",   model=sb.splsicox_model)
    km_cmp  <- getAutoKM(type="COMP", model=sb.splsicox_model, comp=1:2)
    all_km  <- c(km_lp$LST_PLOTS, km_cmp$LST_PLOTS$molecular, km_cmp$LST_PLOTS$clinical)
    for (nm in names(all_km)) {
      save_plot(all_km[[nm]], file.path(odir,paste0("KM_", nm, ".png")), width=7, height=5)
    }
  }
}

combined_time_df <- bind_rows(
  lapply(names(all_block_evals), function(block) {
    # this is the tibble with one row per (model, time_point)
    tdf <- all_block_evals[[block]]$df
    tdf$block <- block
    tdf
  })
)

close(pb)
cat("\nAll blocks finished!\n")
