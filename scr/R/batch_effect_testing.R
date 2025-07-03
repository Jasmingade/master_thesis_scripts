library(readr)
library(dplyr)
library(stringr)
library(patchwork)
library(ggplot2)
library(tidyr)
library(rlang)


analyze_batch_effects <- function(expr_file, title_prefix = "", label = "dataset") {
  expr_df <- read_csv(expr_file, show_col_types = FALSE)
  
  sample_names <- colnames(expr_df)[-1]
  batch_labels <- str_extract(sample_names, "JHU_\\d{8}")
  batch_labels <- ifelse(is.na(batch_labels), "Unknown", batch_labels)
  
  # Transpose expression matrix
  expr_matrix <- t(as.matrix(expr_df[, -1]))
  colnames(expr_matrix) <- expr_df$feature
  rownames(expr_matrix) <- sample_names
  
  # Remove features with too many NAs
  na_thresh <- 0.5 * nrow(expr_matrix)
  expr_matrix <- expr_matrix[, colSums(is.na(expr_matrix)) < na_thresh]
  
  # Optionally remove samples with too many NAs
  expr_matrix <- expr_matrix[rowSums(is.na(expr_matrix)) < 0.3 * ncol(expr_matrix), ]
  
  # Impute remaining NAs with median of each feature
  expr_matrix <- apply(expr_matrix, 2, function(x) {
    x[is.na(x)] <- median(x, na.rm = TRUE)
    x
  })
  
  # Remove constant (zero-variance) features
  zero_var_cols <- apply(expr_matrix, 2, function(x) var(x) == 0)
  expr_matrix <- expr_matrix[, !zero_var_cols]
  
  # Check for missing or infinite values
  sum_na <- sum(is.na(expr_matrix))
  sum_inf <- sum(is.infinite(expr_matrix))
  cat("Missing values (NA):", sum_na, "\n")
  cat("Infinite values (Inf):", sum_inf, "\n")
  
  cat("Features before:", nrow(expr_df), "\n")
  cat("Remaining features:", ncol(expr_matrix), "\n")
  
  cat("Samples before:", ncol(expr_df), "\n")
  cat("Remaining samples:", nrow(expr_matrix), "\n")
  
  
  pca_res <- prcomp(expr_matrix, center = TRUE, scale. = TRUE)
  pca_df <- as.data.frame(pca_res$x)[, 1:15]
  pca_df$batch <- batch_labels
  
  # Scree plot
  scree <- summary(pca_res)$importance[2, 1:20]
  scree_file <- paste0("figures/batch_effect_testing/scree_plot_", label, ".png")
  png(scree_file, width = 800, height = 600)
  plot(scree, type = "b", ylab = "Variance Explained", xlab = "Principal Component",
       main = paste(title_prefix, "Scree Plot"))
  dev.off()
  
  # PCA grid
  grid <- expand.grid(x = 1:5, y = 1:5)
  grid <- grid[grid$x < grid$y, ]
  pc_var <- summary(pca_res)$importance[2, ]
  get_pc_label <- function(pc_num) {
    var <- round(100 * pc_var[pc_num], 1)
    paste0("PC", pc_num, " (", var, "%)")
  }
  
  plots <- lapply(1:nrow(grid), function(i) {
    x_num <- grid$x[i]
    y_num <- grid$y[i]
    x_pc <- paste0("PC", x_num)
    y_pc <- paste0("PC", y_num)
    
    x_label <- get_pc_label(x_num)
    y_label <- get_pc_label(y_num)
    
    p <- ggplot(pca_df, aes(x = !!sym(x_pc), y = !!sym(y_pc), color = batch)) +
      geom_point(size = 1.2, alpha = 0.7) +
      labs(title = paste0(x_pc, " vs ", y_pc), x = x_label, y = y_label) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6)
      )
    if (i != 1) {
      p <- p + guides(color = "none")
    }
    return(p)
  })
  
  g <- wrap_plots(plots, ncol = 4, guides = "collect") & theme(legend.position = "right")
  ggsave(paste0("figures/batch_effect_testing/pca_grid_", label, ".png"),
         g, width = 12, height = 9, dpi = 300)
  
  # ANOVA summary
  summary_df <- data.frame(
    PC = paste0("PC", 1:15),
    PctVar = pc_var[1:15],
    Pval = sapply(1:15, function(i) {
      pc_scores <- pca_res$x[, i]
      anova_res <- aov(pc_scores ~ batch_labels)
      summary(anova_res)[[1]]$`Pr(>F)`[1]
    })
  )
  write.csv(summary_df, paste0("figures/batch_effect_testing/anova_summary_", label, ".csv"),
            row.names = FALSE)
  
  # Boxplots
  pca_long <- pca_df %>%
    pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "value") %>%
    mutate(PC = factor(PC, levels = paste0("PC", 1:15)))
  
  boxplot_file <- paste0("figures/batch_effect_testing/pc_boxplots_", label, ".png")
  p <- ggplot(pca_long, aes(x = batch, y = value, fill = batch)) +
    geom_boxplot() +
    facet_wrap(~PC, scales = "free_y", ncol = 5) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "right"
    ) +
    ggtitle(paste(title_prefix, "PC Scores by Batch"))
  ggsave(boxplot_file, p, width = 12, height = 7, dpi = 300)
}


analyze_batch_effects(
  expr_file = "data/processed/proteomics/gene/TMT10_gene_PDC000110_combined.csv",
  title_prefix = "Gene – ",
  label = "gene"
)

analyze_batch_effects(
  expr_file = "data/processed/proteomics/iso_log/TMT10_iso_log_PDC000110_combined.csv",
  title_prefix = "Isoform (log) – ",
  label = "iso_log"
)

analyze_batch_effects(
  expr_file = "data/processed/proteomics/iso_frac/TMT10_iso_frac_PDC000110_combined.csv",
  title_prefix = "Isoform (fraction) – ",
  label = "iso_frac"
)
