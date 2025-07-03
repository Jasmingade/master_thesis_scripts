# === CoxPH Visualization Script ===

library(ggplot2)
library(dplyr)
library(readr)

# === SETTINGS ===
pdc_id <- "PDC000110"
data_types <- c("iso_log", "iso_frac", "gene")
results_dir <- file.path("data", "analysis_results", "proteomics", pdc_id)
plot_dir <- file.path(results_dir, "plots", data_types)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# === Function to visualize results ===
visualize_results <- function(type) {
  infile <- file.path(results_dir, paste0("coxph_", type, "_", pdc_id, ".csv"))
  if (!file.exists(infile)) return(NULL)
  
  df <- read_csv(infile, show_col_types = FALSE)
  if (!all(c("HR", "p", "adj_p") %in% colnames(df))) return(NULL)
  
  df <- df %>%
    mutate(logHR = log2(HR),
           neglogP = -log10(p),
           significant = adj_p < 0.05)
  
  ## 1. Boxplot with colored points
  p1 <- ggplot(df, aes(x = "", y = logHR)) +
    geom_boxplot(outlier.shape = NA, fill = "grey95") +
    geom_jitter(aes(color = significant), width = 0.2, alpha = 0.5, size = 1) +
    scale_color_manual(values = c("black", "red")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = paste("log2(HR) distribution –", type),
         y = "log2(Hazard Ratio)", x = "") +
    theme_minimal()
  
  ggsave(file.path(plot_dir, paste0("boxplot_logHR_", type, "_", pdc_id, ".png")), p1, width = 6, height = 5)
  
  ## 2. Volcano plot
  p2 <- ggplot(df, aes(x = logHR, y = neglogP, color = significant)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("grey60", "firebrick")) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
    labs(title = paste("Volcano plot –", type),
         x = "log2(HR)", y = "-log10(p-value)") +
    theme_minimal()
  
  ggsave(file.path(plot_dir, paste0("volcano_", type, "_", pdc_id, ".png")), p2, width = 6, height = 5)
  
  ## 3. (Optional) Barplot of top 20 features
  top_df <- df %>% filter(significant) %>% arrange(p) %>% head(20)
  if (nrow(top_df) > 0) {
    p3 <- ggplot(top_df, aes(x = reorder(feature, HR), y = HR)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      geom_hline(yintercept = 1, linetype = "dashed") +
      coord_flip() +
      labs(title = paste("Top 20 features –", type),
           x = "Feature", y = "Hazard Ratio") +
      theme_minimal()
    
    ggsave(file.path(plot_dir, paste0("top_features_", type, "_", pdc_id, ".png")), p3, width = 6, height = 5)
  }
}

# === Run for each type ===
invisible(lapply(data_types, visualize_results))
