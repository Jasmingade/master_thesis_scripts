# --- 1. Hierarchical clustering, saving each as PNG ---
library(readr)
library(dplyr)
library(stringr)
library(proBatch)
library(tibble)
library(limma)
library(viridis)
library(RColorBrewer)
library(Polychrome)
library(gridExtra)
library(grid)
library(pals)

pdc_id <- "PDC000270"
data_types <- c("gene", "iso_log", "iso_frac")
base_dir <- "data"
input_dir <- file.path(base_dir, "processed", "proteomics")
clinical_file <- file.path(base_dir, "processed", "clin_df_all.csv")
output_base_dir <- "data/analysis_results/probatch_diagnostics/panel_layouts"
hclust_dir <- file.path(output_base_dir, "hclust_panels")
dir.create(output_base_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(hclust_dir, recursive = TRUE, showWarnings = FALSE)
clin_df_all <- read_csv(clinical_file, show_col_types = FALSE)

batch_colors      <- Polychrome::palette36.colors(25)[1:25]
tumor_stage_lvls  <- c("Stage I", "Stage II", "Stage III", "Stage IV", "Staging is not applicable or unknown", "NA")
tumor_stage_colors <- c(brewer.pal(4, "Pastel1"), "#999999", "#CCCCCC")

for (type in data_types) {
  pattern <- paste0("TMT\\d+_", type, "_", pdc_id, "_combined\\.csv$")
  files <- list.files(file.path(input_dir, type), pattern = pattern, full.names = TRUE)
  if (length(files) == 0) next
  expr_file <- files[1]
  expr_df <- read_csv(expr_file, show_col_types = FALSE)
  expr_matrix <- as.matrix(expr_df[, -1])
  rownames(expr_matrix) <- expr_df$feature
  sample_names <- colnames(expr_df)[-1]
  
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
  
  expr_matrix <- expr_matrix[rowSums(is.na(expr_matrix)) < ncol(expr_matrix), , drop = FALSE]
  expr_matrix <- expr_matrix[, colSums(is.na(expr_matrix)) < nrow(expr_matrix), drop = FALSE]
  good_samples <- intersect(colnames(expr_matrix), sample_annotation$FullRunName)
  expr_matrix <- expr_matrix[, good_samples, drop = FALSE]
  sample_annotation <- sample_annotation[sample_annotation$FullRunName %in% good_samples, ]
  
  run_date_lvls <- sort(unique(sample_annotation$injection_date))
  color_list <- sample_annotation_to_colors(sample_annotation, 
                                            factor_columns = c("batch","tumor_stage"), 
                                            numeric_columns = c("injection_date", "order"))
  
  cyclic_loess_matrix <- normalizeCyclicLoess(expr_matrix, method = "fast")
  selected_annotations <- c("batch", "injection_date", "tumor_stage")
  file_out <- file.path(hclust_dir, paste0(pdc_id, "_", type, "_hclust.png"))
  plot_hierarchical_clustering(
    cyclic_loess_matrix,
    sample_annotation = sample_annotation,
    color_list = color_list,
    factors_to_plot = selected_annotations,
    distance = 'euclidean',
    agglomeration = 'complete',
    label_samples = FALSE,
    base_size = 8,
    plot_title = NULL,
    filename = file_out
  )
}

# --- 2. Combine individual PNGs in a grid ---
library(png)
library(grid)
library(gridExtra)

types <- c("gene", "iso_log", "iso_frac")
png_files <- file.path(hclust_dir, paste0(pdc_id, "_", types, "_hclust.png"))
subtitles <- c("Gene", "Iso_log", "Iso_frac")  # edit as needed

# Read PNGs as grobs and add subtitles ABOVE each
grobs_with_titles <- mapply(function(pngfile, subtitle) {
  if (!file.exists(pngfile)) {
    warning(sprintf("File does not exist: %s", pngfile))
    raster <- grid::textGrob("Missing file")
  } else {
    img <- png::readPNG(pngfile)
    raster <- grid::rasterGrob(img, interpolate = TRUE)
  }
  # Subtitle above, then image
  arrangeGrob(
    grid::textGrob(subtitle, gp=grid::gpar(fontsize=16, fontface="bold")),
    raster,
    ncol = 1,
    heights = c(1, 10)
  )
}, png_files, subtitles, SIMPLIFY = FALSE)

# Save the grid
outfile <- file.path(output_base_dir, paste0(pdc_id, "_hclust_grid.png"))
png(outfile, width = 2300, height = 500)
gridExtra::grid.arrange(
  grobs = grobs_with_titles,
  nrow = 1,
  top = grid::textGrob(
    label = paste0("Hierarchical Clustering Panels for ", pdc_id),
    gp = grid::gpar(fontsize = 30, fontface = "bold")
  )
)
dev.off()
cat("Saved PNG grid to", outfile, "\n")



