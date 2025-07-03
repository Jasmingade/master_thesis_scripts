library(readr)
library(dplyr)
library(stringr)
library(proBatch)
library(tibble)
library(limma)
library(ggpubr)
library(viridis)
library(RColorBrewer)
library(Polychrome)
library(gridExtra)
library(grid)
library(patchwork)

# ---- Settings ----
pdc_id <- "PDC000270"
data_types <- c("gene", "iso_log", "iso_frac")
base_dir <- "data"
input_dir <- file.path(base_dir, "processed", "proteomics")
clinical_file <- file.path(base_dir, "processed", "clin_df_all.csv")
output_base_dir <- "data/analysis_results/probatch_diagnostics/panel_layouts"
dir.create(output_base_dir, recursive = TRUE, showWarnings = FALSE)
clin_df_all <- read_csv(clinical_file, show_col_types = FALSE)

# --- Custom color palettes for categorical PCA variables ---
batch_colors      <- Polychrome::palette36.colors(25)[1:25] # Up to 13 batches
tumor_stage_lvls  <- c("Stage I", "Stage II", "Stage III", "Stage IV", "Staging is not applicable or unknown", "NA")
tumor_stage_colors <- c(brewer.pal(4, "Set1"), "#999999", "#CCCCCC") # For up to 6 levels
age_group_lvls    <- c("<40", "40–49", "50–59", "60–69", "70–79", "80+", "NA")
age_group_colors  <- c(brewer.pal(6, "Dark2"), "#CCCCCC") # 6 + 1 levels
sex_lvls          <- c("Female", "Male", "NA")
sex_colors        <- c("lightcoral", "lightgreen", "grey") # blue, gray

# Store panels here
pca_panel_list <- list()

for (type in data_types) {
  # ---- Input file ----
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
  
  # Remove all-NA features/samples and sync annotation
  expr_matrix <- expr_matrix[rowSums(is.na(expr_matrix)) < ncol(expr_matrix), , drop = FALSE]
  expr_matrix <- expr_matrix[, colSums(is.na(expr_matrix)) < nrow(expr_matrix), drop = FALSE]
  good_samples <- intersect(colnames(expr_matrix), sample_annotation$FullRunName)
  expr_matrix <- expr_matrix[, good_samples, drop = FALSE]
  sample_annotation <- sample_annotation[sample_annotation$FullRunName %in% good_samples, ]
  
  # ---- Color scheme ----
  # Build a color_list as usual, but we will override for PCA plots below
  color_list <- sample_annotation_to_colors(
    sample_annotation,
    factor_columns = c('batch', 'cancer_type', 'sex', 'tumor_stage', 'age_group'),
    numeric_columns = c('order', 'injection_date')
  )
  
  # ---- Normalize (cyclic loess) ----
  cyclic_loess_matrix <- normalizeCyclicLoess(expr_matrix, method = "fast")
  
  # ---- Only plot if enough samples for PCA ----
  if (ncol(cyclic_loess_matrix) > 2) {
    main_title <- paste0("PDC: ", pdc_id, " | ", type)
    pca_vars <- c(batch = "Batch", tumor_stage = "Tumor Stage", age_group = "Age Group", sex = "Sex")
    
    # --- Custom color schemes for each PCA variable
    color_scheme_list <- list(
      batch       = setNames(batch_colors, sort(unique(sample_annotation$batch))),
      tumor_stage = setNames(tumor_stage_colors, tumor_stage_lvls),
      age_group   = setNames(age_group_colors, age_group_lvls),
      sex         = setNames(sex_colors, sex_lvls)
    )
    
    # --- Build PCA plots with consistent coloring
    pca_plots <- lapply(names(pca_vars), function(v) {
      plot_PCA(
        cyclic_loess_matrix, sample_annotation,
        color_by = v, 
        plot_title = paste0(pca_vars[[v]]),
        base_size = 10,
        color_scheme = color_scheme_list[[v]]
      ) + 
        coord_fixed()  # <-- Forces a 1:1 aspect ratio, so all panels are square
    })
    
    # 2x2 patchwork grid panel for this data type
    panel <- (pca_plots[[1]] + pca_plots[[2]]) /
      (pca_plots[[3]] + pca_plots[[4]]) +
      plot_annotation(title = paste0(type),
                      theme = theme(plot.title = element_text(size=16, hjust=0.5))
      )
    pca_panel_list[[type]] <- panel
  }
}

# ---- Combine all data type PCA panels into one 1x3 grid ----
# Add a super-title
titles <- wrap_elements(grid::textGrob("Gene", gp=grid::gpar(fontsize=15, fontface="bold"))) +
  wrap_elements(grid::textGrob("Iso_log", gp=grid::gpar(fontsize=15, fontface="bold"))) +
  wrap_elements(grid::textGrob("Iso_frac", gp=grid::gpar(fontsize=15, fontface="bold")))
final_pca_grid <- titles / 
  (pca_panel_list[["gene"]] | pca_panel_list[["iso_log"]] | pca_panel_list[["iso_frac"]]) +
  plot_layout(heights = c(0.09, 1), guides = "collect") +
  plot_annotation(
    title = paste0("PCA Panel Comparison for ", pdc_id),
    theme = theme(plot.title = element_text(size=20, hjust=0.5, face="bold"))
  )

final_pca_grid <- final_pca_grid + plot_layout(widths = c(1, 1, 1))

# Save it
ggsave(
  file.path(output_base_dir, paste0(pdc_id, "_PCA_grid.png")),
  final_pca_grid, width = 30, height = 10
)
