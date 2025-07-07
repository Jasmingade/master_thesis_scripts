# ---- Load libraries ----
library(readr)
library(dplyr)
library(stringr)
library(proBatch)
library(tibble)
library(ggplot2)
library(patchwork)
library(grid)
library(gridExtra)
library(png)
library(pals)
library(limma)

# ---- Settings ----
pdc_ids <- c("PDC000110")
data_types <- c("gene", "iso_log", "iso_frac")
base_dir <- "data"
input_dir <- file.path(base_dir, "processed", "proteomics")
clinical_file <- file.path(base_dir, "processed", "clin_df_all.csv")
output_base_dir <- "data/analysis_results/probatch_diagnostics/batch_correction"
dir.create(output_base_dir, recursive = TRUE, showWarnings = FALSE)
clin_df_all <- read_csv(clinical_file, show_col_types = FALSE)

# Custom color palettes
batch_colors      <- pals::alphabet(25)
tumor_stage_lvls  <- c("Stage I", "Stage II", "Stage III", "Stage IV", "Staging is not applicable or unknown", "NA")
tumor_stage_colors <- c(RColorBrewer::brewer.pal(4, "Set1"), "#999999", "#CCCCCC")
age_group_lvls    <- c("<40", "40–49", "50–59", "60–69", "70–79", "80+", "NA")
age_group_colors  <- c(RColorBrewer::brewer.pal(6, "Dark2"), "#CCCCCC")
sex_lvls          <- c("Female", "Male", "NA")
sex_colors        <- c("lightcoral", "lightgreen", "grey")


for (pdc_id in pdc_ids) {
  cat("## -------------------- ##
## Processing PDC000110 ##
## -------------------- ##", pdc_id, "\n")
  pca_panel_before <- list(); pca_panel_after <- list()
  
  for (type in data_types) {
    # ---- Load data ----
    pattern <- paste0("TMT\\d+_", type, "_", pdc_id, "_combined\\.csv$")
    files <- list.files(file.path(input_dir, type), pattern = pattern, full.names = TRUE)
    if (length(files) == 0) next
    expr_file <- files[1]
    expr_df <- read_csv(expr_file, show_col_types = FALSE)
    expr_matrix <- as.matrix(expr_df[, -1])
    rownames(expr_matrix) <- expr_df$feature
    sample_names <- colnames(expr_df)[-1]
    cat("Processing", pdc_id, type, "from file:", expr_file, "\n")
    
    # ---- Annotation ----
    sample_annotation <- tibble(
      FullRunName = sample_names,
      case_id = str_extract(sample_names, "(?<=_matrix_)([0-9]{2}[A-Z]{2}[0-9]{3}|C3[NL]-[0-9]{5}|NX[0-9]+)"),
      batch = str_extract(sample_names, "^\\d{2}CPTAC"),
      injection_date = str_extract(sample_names, "\\d{8}"),
      order = seq_along(sample_names),
      data_type = type
    ) %>%
      left_join(
        clin_df_all %>%
          select(case_id, tumor_code, age, sex, race, ethnicity, age_group, bmi,
                 tumor_stage, histologic_subtype, histologic_grade, tumor_size_cm,
                 alcohol_consumption, tobacco_smoking_history, vital_status, OS_time, OS_event),
        by = "case_id"
      ) %>%
      rename(cancer_type = tumor_code)
    
    # Remove all-NA features/samples and sync annotation
    expr_matrix <- expr_matrix[rowSums(is.na(expr_matrix)) < ncol(expr_matrix), , drop = FALSE]
    expr_matrix <- expr_matrix[, colSums(is.na(expr_matrix)) < nrow(expr_matrix), drop = FALSE]
    good_samples <- intersect(colnames(expr_matrix), sample_annotation$FullRunName)
    expr_matrix <- expr_matrix[, good_samples, drop = FALSE]
    sample_annotation <- sample_annotation[sample_annotation$FullRunName %in% good_samples, ]
    
    # Color list for factors
    color_list <- sample_annotation_to_colors(
      sample_annotation,
      factor_columns = c('batch', 'cancer_type', 'sex', 'tumor_stage', 'age_group'),
      numeric_columns = c('order')
    )
    
    # ---- 1. Cyclic Loess normalization ----
    norm_matrix <- normalizeCyclicLoess(expr_matrix, method = "fast")
    
    # Save normalized
    norm_matrix_outfile <- file.path(output_base_dir, paste0(pdc_id, "_", type, "_cyclic_loess.csv"))
    write.csv(norm_matrix, norm_matrix_outfile, quote = FALSE, row.names = TRUE)
    cat("Dimensions of norm_matrix before:", dim(norm_matrix), "\n")
    
    # Remove features with NA rownames
    na_rownames <- is.na(rownames(norm_matrix))
    if (any(na_rownames)) {
      cat("Removing", sum(na_rownames), "features with NA rownames\n")
      norm_matrix <- norm_matrix[!na_rownames, , drop = FALSE]
    }
    stopifnot(!any(is.na(rownames(norm_matrix))))  # This should now pass
    
    # i. Sync annotation and matrix columns
    good_samples <- intersect(colnames(norm_matrix), sample_annotation$FullRunName)
    norm_matrix <- norm_matrix[, good_samples, drop = FALSE]
    sample_annotation <- sample_annotation[sample_annotation$FullRunName %in% good_samples, ]
    sample_annotation <- sample_annotation[match(colnames(norm_matrix), sample_annotation$FullRunName), ]
    stopifnot(all(sample_annotation$FullRunName == colnames(norm_matrix)))
    cat("All sample names match and are ordered!\n")
    
    # --- Debug checks ---
    cat("dim(norm_matrix):", dim(norm_matrix), "\n")
    stopifnot(ncol(norm_matrix) == nrow(sample_annotation))
    stopifnot(identical(colnames(norm_matrix), sample_annotation$FullRunName))
    stopifnot(length(rownames(norm_matrix)) == length(unique(rownames(norm_matrix))))
    stopifnot(!any(is.na(rownames(norm_matrix))))
    
    
    # ii. Build batch index list
    by_batch <- split(seq_along(colnames(norm_matrix)), sample_annotation$batch)
    cat("Batch sample sizes:", sapply(by_batch, length), "\n")
    cat("max idx in by_batch:", max(unlist(by_batch)), "\n")
    
    # iii. Filter features
    bad_features <- sapply(rownames(norm_matrix), function(f) {
      any(sapply(by_batch, function(idx) {
        # Defensive: Check for out-of-bounds
        if (any(idx > ncol(norm_matrix))) {
          cat(sprintf("Batch index out of bounds: max idx %d, ncol %d\n", max(idx), ncol(norm_matrix)))
          return(TRUE)
        }
        # The actual two conditions:
        all(is.na(norm_matrix[f, idx])) | sum(!is.na(norm_matrix[f, idx])) < 2
      }))
    })
    
    cat("Number of bad features:", sum(bad_features), "\n")
    cat("Length of norm_matrix rows:", nrow(norm_matrix), "\n")
    stopifnot(length(bad_features) == nrow(norm_matrix))  # Must match!
    
    # iiii. Subset
    norm_matrix <- norm_matrix[!bad_features, , drop = FALSE]
    cat("Dimensions of norm_matrix after:", dim(norm_matrix), "\n")
    
    
    # ---- 2. Batch correction ----
    df_long <- matrix_to_long(norm_matrix,
                              sample_annotation = sample_annotation,
                              feature_id_col = "feature",
                              measure_col = "Intensity",
                              sample_id_col = "FullRunName")
    
    batch_matrix_df <- correct_batch_effects_df(
      df_long = df_long,
      sample_annotation = sample_annotation,
      discrete_func = "ComBat",
      sample_id_col = "FullRunName",
      batch_col = "batch",
      feature_id_col = "feature",
      measure_col = "Intensity",
      
    )
    # Convert batch-corrected long df back to wide matrix
    batch_matrix <- long_to_matrix(
      batch_matrix_df,
      feature_id_col = "feature",
      sample_id_col = "FullRunName",
      measure_col = "Intensity"
    )
    batch_matrix <- as.matrix(batch_matrix)
    
    # Save batch-corrected
    batch_matrix_outfile <- file.path(output_base_dir, paste0(pdc_id, "_", type, "_batch_corrected.csv"))
    write.csv(batch_matrix, batch_matrix_outfile, quote = FALSE, row.names = TRUE)
    
    # ---- 3. PCA Plots (Before/After Correction, 4 variables each) ----
    pca_vars <- c(batch = "Batch", tumor_stage = "Tumor Stage", age_group = "Age Group", sex = "Sex")
    color_scheme_list <- list(
      batch       = setNames(batch_colors, sort(unique(sample_annotation$batch))),
      tumor_stage = setNames(tumor_stage_colors, tumor_stage_lvls),
      age_group   = setNames(age_group_colors, age_group_lvls),
      sex         = setNames(sex_colors, sex_lvls)
    )
    # Before
    panel_before <- lapply(names(pca_vars), function(v) {
      plot_PCA(
        norm_matrix, sample_annotation,
        color_by = v,
        plot_title = paste0(pca_vars[[v]], " (before correction)"),
        base_size = 10,
        color_scheme = color_scheme_list[[v]]
      ) + coord_fixed()
    })
    pca_panel_before[[type]] <- (panel_before[[1]] + panel_before[[2]]) /
      (panel_before[[3]] + panel_before[[4]])
    # After
    panel_after <- lapply(names(pca_vars), function(v) {
      plot_PCA(
        batch_matrix, sample_annotation,
        color_by = v,
        plot_title = paste0(pca_vars[[v]], " (after correction)"),
        base_size = 10,
        color_scheme = color_scheme_list[[v]]
      ) + coord_fixed()
    })
    pca_panel_after[[type]] <- (panel_after[[1]] + panel_after[[2]]) /
      (panel_after[[3]] + panel_after[[4]])
    
    # ---- 4. Save HCLUST PNGs before and after ----
    hclust_before_file <- file.path(output_base_dir, paste0(pdc_id, "_", type, "_hclust_before.png"))
    hclust_after_file  <- file.path(output_base_dir, paste0(pdc_id, "_", type, "_hclust_after.png"))
    # Before correction
    png(hclust_before_file, width = 1200, height = 900, res = 120)
    plot_hierarchical_clustering(
      norm_matrix,
      sample_annotation = sample_annotation,
      color_list = color_list,
      factors_to_plot = c("batch", "tumor_stage"),
      distance = 'euclidean',
      agglomeration = 'complete',
      label_samples = FALSE,
      plot_title = paste(type, "(before correction)")
    )
    dev.off()
    # After correction
    png(hclust_after_file, width = 1200, height = 900, res = 120)
    plot_hierarchical_clustering(
      batch_matrix,
      sample_annotation = sample_annotation,
      color_list = color_list,
      factors_to_plot = c("batch", "tumor_stage"),
      distance = 'euclidean',
      agglomeration = 'complete',
      label_samples = FALSE,
      plot_title = paste(type, "(after correction)")
    )
    dev.off()
  }
  
  ## ---- 5. Save PCA side-by-side grid (all types) ----
  titles <- wrap_elements(grid::textGrob("Gene", gp=grid::gpar(fontsize=15, fontface="bold"))) +
    wrap_elements(grid::textGrob("Iso_log", gp=grid::gpar(fontsize=15, fontface="bold"))) +
    wrap_elements(grid::textGrob("Iso_frac", gp=grid::gpar(fontsize=15, fontface="bold")))
  
  pca_grid_before <- (pca_panel_before[["gene"]] | pca_panel_before[["iso_log"]] | pca_panel_before[["iso_frac"]]) +
    plot_layout(widths = c(1, 1, 1), guides = "collect") +
    plot_annotation(title = paste0("PCA Before Batch Correction: ", pdc_id))
  
  pca_grid_after <- (pca_panel_after[["gene"]] | pca_panel_after[["iso_log"]] | pca_panel_after[["iso_frac"]]) +
    plot_layout(widths = c(1, 1, 1), guides = "collect") +
    plot_annotation(title = paste0("PCA After Batch Correction: ", pdc_id))
  
  final_pca_side_by_side <- pca_grid_before / pca_grid_after + plot_layout(heights = c(1, 1))
  
  ggsave(file.path(output_base_dir, paste0(pdc_id, "_PCA_side_by_side.png")),
         final_pca_side_by_side, width = 24, height = 14)
  
  ## ---- 6. HCLUST grid panel ----
  hclust_pngs_before <- file.path(output_base_dir, paste0(pdc_id, "_", data_types, "_hclust_before.png"))
  hclust_pngs_after  <- file.path(output_base_dir, paste0(pdc_id, "_", data_types, "_hclust_after.png"))
  grobs_before <- lapply(hclust_pngs_before, function(f) {
    if (!file.exists(f)) grid::textGrob("Missing file") else grid::rasterGrob(png::readPNG(f), interpolate = TRUE)
  })
  grobs_after <- lapply(hclust_pngs_after, function(f) {
    if (!file.exists(f)) grid::textGrob("Missing file") else grid::rasterGrob(png::readPNG(f), interpolate = TRUE)
  })
  outfile <- file.path(output_base_dir, paste0(pdc_id, "_hclust_side_by_side.png"))
  png(outfile, width = 3600, height = 1800)
  grid.arrange(
    gridExtra::arrangeGrob(grobs_before[[1]], grobs_before[[2]], grobs_before[[3]],
                           nrow = 1,
                           top = grid::textGrob("Hierarchical Clustering Before Correction", gp = grid::gpar(fontsize = 24, fontface = "bold"))
    ),
    gridExtra::arrangeGrob(grobs_after[[1]], grobs_after[[2]], grobs_after[[3]],
                           nrow = 1,
                           top = grid::textGrob("Hierarchical Clustering After Correction", gp = grid::gpar(fontsize = 24, fontface = "bold"))
    ),
    nrow = 2
  )
  dev.off()
  cat("Saved PCA and HCLUST grids for", pdc_id, "\n")
}







