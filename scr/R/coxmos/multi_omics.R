library(Coxmos)
library(impute)
library(dplyr)
library(stringr)
library(purrr)



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
  PDC000110 = "OV",PDC000116 = "COAD",PDC000120 = "BRCA",PDC000125 = "UCEC",PDC000127 = "KIRC",
  PDC000153 = "LUAD",PDC000204 = "GBM",PDC000221 = "HNSC",PDC000234 = "LUSC",PDC000270 = "PAAD"
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
selected_ids <- c("LUAD")

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

for(cancer in selected_ids){
  cat("=== Coxmos for", cancer, "===\n")
  
  # 1) grab all your PROT_/TRANS_ blocks for this cancer
  omics_for_cancer <- prot_blocks[
    str_detect(names(prot_blocks), paste0("^(PROT|TRANS)_", cancer, "_"))
  ]
  
  # 2) Preprocess each matrix:
  molecular_blocks <- lapply(omics_for_cancer, function(raw_mat) {
    # extract case IDs from the original column names
    samp_ids <- colnames(raw_mat)
    case_ids <- str_extract(
      samp_ids,
      "(TCGA-[0-9]{2}-[0-9]{4}-[0-9]{2}|[0-9]{2}[A-Z]{2}[0-9]{3}|C3[NL]-[0-9]{5}|NX[0-9]+)"
    )
    if (any(is.na(case_ids))) {
      stop("Failed to parse case_id for: ",
           paste(samp_ids[is.na(case_ids)], collapse = ", "))
    }
    colnames(raw_mat) <- case_ids
    
    # transpose & drop >80% missing
    m    <- t(raw_mat)
    keep <- colMeans(is.na(m)) < 0.8
    m2   <- m[, keep, drop = FALSE]
    
    # remember the rownames so we can restore them after imputation
    original_rn <- rownames(m2)
    
    # impute
    imp <- impute.knn(m2)$data
    
    # restore row and col names
    rownames(imp) <- original_rn
    colnames(imp) <- colnames(m2)
    
    # drop any accidental "NA." columns
    imp[, colnames(imp) != "NA.", drop = FALSE]
  })
  cat("  • Processed", length(molecular_blocks), "molecular blocks\n")
  
  # 3) build your clinical annotation once, based on the first block’s case_ids
  sample_ids <- rownames(molecular_blocks[[1]])
  samp_anno  <- tibble(case_id = sample_ids) %>%
    left_join(clin_df_all, by="case_id") %>%
    filter(!is.na(OS_time), !is.na(OS_event)) %>%
    distinct(case_id, .keep_all = TRUE)
  cat("  • after clinical join/filter:", nrow(samp_anno), "samples remain\n")
  
  # 3a) **synchronize** your molecular blocks to exactly these case_ids
  keep_ids <- samp_anno$case_id
  molecular_blocks <- lapply(molecular_blocks, function(m) {
    # check that every case_id is present
    missing <- setdiff(keep_ids, rownames(m))
    if (length(missing) > 0) {
      stop("These case_ids are missing from your block:\n  ", paste(missing, collapse=", "))
    }
    # reorder & subset rows to exactly keep_ids
    m[keep_ids, , drop = FALSE]
  })
  cat("  • synchronized molecular blocks to", length(keep_ids), "samples\n")
  
  # 4) build your clinical block
  samp_anno <- samp_anno %>%
    mutate(
      age_group  = factor(age_group,  levels = c("<40","40–49","50–59","60–69","70–79","80+")),
      sex        = factor(sex,        levels = c("female","male")),
      tumor_stage = factor(tumor_stage_clean, levels = c("stage i","stage ii","stage iii","stage iv", "unknown")),
      hist_grade = factor(histologic_grade,
                          levels = c("G1","G2","G3","G4","GX","[unknown]"))
    )
  
  # Select and one-hot encode
  X_clin_df <- samp_anno %>%
    dplyr::select(case_id, age_group, sex, tumor_stage, hist_grade) %>%
    column_to_rownames("case_id")
  
  X_clin_onehot <- factorToBinary(X = X_clin_df, all = TRUE, sep = "_")
  
  # Force numeric
  X_clin_num <- X_clin_onehot %>%
    mutate(across(everything(), ~ as.integer(.))) %>%    # turn "0"/"1" → 0/1 integers
    as.matrix()
  cat("  • clinical block dims:", dim(X_clin), "\n")
  
  # 5) combine into one multiblock list
  X_blocks <- c(molecular_blocks, list(clinical = X_clin_num))
  
  # 6) response
  Y <- samp_anno %>%
    transmute(time = OS_time, event = OS_event) %>%
    as.data.frame()
  rownames(Y) <- samp_anno$case_id
  
  # 6) split once
  sp        <- getTrainTest(X_blocks, Y, p = 0.7, seed = 1234)
  X_train   <- sp$X_train;  Y_train <- sp$Y_train
  X_test    <- sp$X_test;   Y_test  <- sp$Y_test
  
  # 7) explicitly fit your five models
  cv1 <- cv.mb.coxmos("sb.splsicox", X_train, Y_train,
                      MIN_NVAR=5, MAX_NVAR=20, max.ncomp=1,
                      n_run=2, k_folds=5, MIN_EPV=0.01,
                      remove_zero_variance=TRUE,
                      remove_near_zero_variance=TRUE,
                      remove_variance_at_fold_level=TRUE,
                      remove_non_significant=FALSE,
                      PARALLEL = TRUE)
  m1  <- mb.coxmos("sb.splsicox", X_train, Y_train,
                   n.comp=cv1$opt.comp,
                   penalty=cv1$opt.penalty,
                   remove_non_significant=TRUE)
  
  cv2 <- cv.mb.coxmos("sb.splsdrcox", X_train, Y_train,
                      MIN_NVAR = 5, MAX_NVAR = 20,
                      max.ncomp = 1,
                      n_run = 2,
                      k_folds = 5,
                      MIN_EPV = 0.01,
                      remove_zero_variance = TRUE,
                      remove_near_zero_variance = TRUE,
                      remove_variance_at_fold_level = TRUE,
                      remove_non_significant = FALSE,
                      PARALLEL = TRUE) 
  m2  <- mb.coxmos(method = "sb.splsdrcox",
                   X = X_train, Y = Y_train,
                   n.comp = cv2$opt.comp,
                   vector = cv2$opt.nvar,
                   remove_non_significant = TRUE)
  
  cv3 <- cv.mb.coxmos(method = "sb.splsdacox",
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
                      remove_non_significant = FALSE,
                      PARALLEL = TRUE)
  m3  <- mb.coxmos(method = "sb.splsdacox",
                   X = X_train, Y = Y_train,
                   n.comp = cv3$opt.comp,
                   vector = cv3$opt.nvar,
                   remove_non_significant = TRUE)
  
  cv4 <- cv.mb.coxmos(method = "mb.splsdrcox",
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
                      remove_non_significant = FALSE,
                      PARALLEL = TRUE)
  m4  <- mb.coxmos(method = "mb.splsdrcox",
                   X = X_train, Y = Y_train,
                   n.comp = cv4$opt.comp,
                   vector = cv4$opt.nvar)
  
  cv5 <- cv.mb.coxmos(method = "mb.splsdacox",
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
                      remove_non_significant = FALSE,
                      PARALLEL = TRUE)
  m5  <- mb.coxmos(method = "mb.splsdacox",
                   X = X_train, Y = Y_train,
                   n.comp = cv5$opt.comp,
                   vector = cv5$opt.nvar)
  
  lst_models <- list(
    "SB.sPLS-ICOX"          = m1,
    "SB.sPLS-DRCOX-Dynamic" = m2,
    "SB.sPLS-DACOX-Dynamic" = m3,
    "MB.sPLS-DRCOX"         = m4,
    "MB.sPLS-DACOX"         = m5
  )
  
  # 8) evaluate *once*
  eval_res <- eval_Coxmos_models(
    lst_models   = lst_models,
    X_test       = X_test,
    Y_test       = Y_test,
    times        = NULL,
    progress_bar = TRUE
  )
  
  all_block_evals[[cancer]] <- eval_res
  
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


combined_time_df <- bind_rows(
  lapply(names(all_block_evals), function(block) {
    # this is the tibble with one row per (model, time_point)
    tdf <- all_block_evals[[block]]$df
    tdf$block <- block
    tdf
  })

cat("All cancers done!\n")
