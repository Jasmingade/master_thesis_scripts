library(readr)
library(dplyr)
library(stringr)
library(proBatch)
library(tibble)
library(limma)
library(ggpubr)
library(viridis)
library(Polychrome)
library(gridExtra)
library(grid)
library(cowplot)
library(ggpubr)
library(rlang)
library(ggplot2)
library(patchwork)

# ---- Settings ----
pdc_id <- "PDC000270"
data_types <- c("gene", "iso_log", "iso_frac")
base_dir <- "data"
input_dir <- file.path(base_dir, "processed", "proteomics")
clinical_file <- file.path(base_dir, "processed", "clin_df_all.csv")
output_base_dir <- "data/analysis_results/probatch_diagnostics/panel_layouts"
clin_df_all <- read_csv(clinical_file, show_col_types = FALSE)

# For holding plots
plot_list <- list(
  sample_mean_raw = list(),
  sample_boxplot_raw = list(),
  sample_mean_norm = list(),
  sample_boxplot_norm = list(),
  pca = list(),
  heatmap = list(),
  hclust = list()
)

# ---- Plot per data type ----
for (type in data_types) {
  # ---- File finding and input ----
  pattern <- paste0("TMT\\d+_", type, "_", pdc_id, "_combined\\.csv$")
  files <- list.files(file.path(input_dir, type), pattern = pattern, full.names = TRUE)
  if (length(files) == 0) next
  expr_file <- files[1]
  expr_df <- read_csv(expr_file, show_col_types = FALSE)
  expr_matrix <- as.matrix(expr_df[, -1])
  rownames(expr_matrix) <- expr_df$feature
  sample_names <- colnames(expr_df)[-1]
  
  # ---- Sample annotation ----
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
  
  # Remove all-NA features/samples
  expr_matrix <- expr_matrix[rowSums(is.na(expr_matrix)) < ncol(expr_matrix), , drop = FALSE]
  expr_matrix <- expr_matrix[, colSums(is.na(expr_matrix)) < nrow(expr_matrix), drop = FALSE]
  good_samples <- intersect(colnames(expr_matrix), sample_annotation$FullRunName)
  expr_matrix <- expr_matrix[, good_samples, drop = FALSE]
  sample_annotation <- sample_annotation[sample_annotation$FullRunName %in% good_samples, ]
  
  color_list <- sample_annotation_to_colors(
    sample_annotation,
    factor_columns = c('batch', 'cancer_type', 'sex', 'tumor_stage', 'age_group'),
    numeric_columns = c('order', 'injection_date')
  )
  batch_col <- "batch"
  
  # --- RAW (Mean/Boxplot) ---
  df_long <- matrix_to_long(expr_matrix, sample_annotation = NULL, feature_id_col = "feature", measure_col = "Intensity", sample_id_col = "FullRunName")
  
  main_title <- paste0("PDC: ", pdc_id, " | ", type)
  
  # --- NORMALIZATION ---
  cyclic_loess_matrix <- normalizeCyclicLoess(expr_matrix, method = "fast")
  cyclic_loess_long <- matrix_to_long(
    cyclic_loess_matrix,
    feature_id_col = "feature",
    sample_id_col = "FullRunName",
    measure_col = "Intensity"
  )
  
  # ----- HERE: Calculate shared ylim for mean plots -----
  raw_means  <- apply(expr_matrix, 2, mean, na.rm = TRUE)
  norm_means <- apply(cyclic_loess_matrix, 2, mean, na.rm = TRUE)
  mean_ylim <- range(c(raw_means, norm_means), na.rm = TRUE)
  
  # ----- Plots -----
  # Raw sample mean
  plot_list$sample_mean_raw[[type]] <- plot_sample_mean(
    expr_matrix, sample_annotation,
    sample_id_col = "FullRunName",
    order_col = "order",
    batch_col = batch_col,
    color_by_batch = TRUE,
    base_size = 8,
    plot_title = paste0(type),
    color_scheme = color_list[[batch_col]]
  ) + 
    coord_cartesian(ylim = mean_ylim)
  
  # Raw boxplot
  plot_list$sample_boxplot_raw[[type]] <- plot_boxplot(
    df_long, sample_annotation,
    sample_id_col = "FullRunName",
    batch_col = "batch",
    base_size = 8,
    plot_title = paste0(type),
    order_col = "order",
    width = 0.8  # <-- this controls box width!
  )
  
  # Norm sample mean
  plot_list$sample_mean_norm[[type]] <- plot_sample_mean(
    cyclic_loess_matrix, sample_annotation,
    order_col = "order",
    batch_col = batch_col,
    color_by_batch = TRUE,
    base_size = 8,
    plot_title = paste0(type),
    color_scheme = color_list[[batch_col]]
  ) +
    coord_cartesian(ylim = mean_ylim)
  # Norm boxplot
  plot_list$sample_boxplot_norm[[type]] <- plot_boxplot(
    cyclic_loess_long, sample_annotation,
    sample_id_col = "FullRunName",
    batch_col = "batch",
    base_size = 8,
    plot_title = paste0(type),
    order_col = "order",
    width = 0.8  # <-- this controls box width!
  )
}

# ---- Arrange final layouts ----

mean_plots <- list(
  plot_list$sample_mean_raw[["gene"]],
  plot_list$sample_mean_norm[["gene"]],
  plot_list$sample_mean_raw[["iso_log"]],
  plot_list$sample_mean_norm[["iso_log"]],
  plot_list$sample_mean_raw[["iso_frac"]],
  plot_list$sample_mean_norm[["iso_frac"]]
)

box_plots <- list(
  plot_list$sample_boxplot_raw[["gene"]],
  plot_list$sample_boxplot_norm[["gene"]],
  plot_list$sample_boxplot_raw[["iso_log"]],
  plot_list$sample_boxplot_norm[["iso_log"]],
  plot_list$sample_boxplot_raw[["iso_frac"]],
  plot_list$sample_boxplot_norm[["iso_frac"]]
)

# Build your 3x2 plot grid first (NO legend, NO titles)
plot_grid_means <- (mean_plots[[1]] + mean_plots[[2]]) /
  (mean_plots[[3]] + mean_plots[[4]]) /
  (mean_plots[[5]] + mean_plots[[6]]) +
  plot_layout(guides = "collect") & 
  theme(legend.direction= 'vertical')

plot_grid_box <- (box_plots[[1]] + box_plots[[2]]) /
  (box_plots[[3]] + box_plots[[4]]) /
  (box_plots[[5]] + box_plots[[6]]) +
  plot_layout(guides = "collect") & 
  theme(legend.direction= 'vertical')

# Add column headers with patchwork
headers <- wrap_elements(grid::textGrob("RAW", gp=grid::gpar(fontsize=15, fontface="bold"))) +
  wrap_elements(grid::textGrob("CYCLIC LOESS", gp=grid::gpar(fontsize=15, fontface="bold")))


# Add super-titles (example for mean and box grids)
final_mean <- (
  headers / plot_grid_means +
    plot_layout(heights = c(0.08, 1))
) +
  plot_annotation(
    title = paste0("Sample Means: Raw vs. Cyclic Loess for ", pdc_id),
    theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))
  )

final_box <- (
  headers / plot_grid_box +
    plot_layout(heights = c(0.08, 1))
) +
  plot_annotation(
    title = paste0("Sample Boxplots: Raw vs. Cyclic Loess for ", pdc_id),
    theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))
  )



# Save it
ggsave(
  file.path(output_base_dir, paste0(pdc_id, "_mean_grid.png")),
  final_mean, width = 15, height = 10
)
ggsave(
  file.path(output_base_dir, paste0(pdc_id, "_box_grid.png")),
  final_box, width = 15, height = 10
)
