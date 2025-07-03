library(tidyverse)
library(biomaRt)
library(pheatmap)
library(ggridges)

# ----------------------------------------
# Set paths
# ----------------------------------------
gene_dir <- "gene"
pa_dir <- "protein_assembly_summary"
out_dir <- "qc_outputs"
dir.create(out_dir, showWarnings = FALSE)

# ----------------------------------------
# Prepare biomart connection (only once)
# ----------------------------------------
mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

# ----------------------------------------
# List gene-level TMT files
# ----------------------------------------
gene_files <- list.files(gene_dir, pattern = "TMT(10|11)_gene_.*\\.csv$", full.names = TRUE)

# ----------------------------------------
# Collect correlation results
# ----------------------------------------
all_corrs <- list()

for (file in gene_files) {
  fname <- basename(file)
  study_id <- str_extract(fname, "PDC\\d+")
  message("Processing: ", study_id)
  
  # ---- Load gene-level TMT matrix ----
  tmt_gene <- read.csv(file, check.names = FALSE) %>%
    rename(ensembl_gene_id = 1) %>%
    filter(!is.na(ensembl_gene_id) & ensembl_gene_id != "") %>%
    mutate(ensembl_gene_id = sub("\\..*", "", ensembl_gene_id))
  
  gene_map <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = tmt_gene$ensembl_gene_id,
    mart = mart
  )
  
  tmt_gene <- tmt_gene %>%
    as_tibble() %>%
    left_join(gene_map, by = "ensembl_gene_id") %>%
    filter(hgnc_symbol != "") %>%
    dplyr::select(-ensembl_gene_id) %>%
    group_by(hgnc_symbol) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") %>%
    rename(Gene = hgnc_symbol)
  
  # ---- Aggregate TMT log2 values per study ----
  tmt_long <- tmt_gene %>%
    pivot_longer(-Gene, names_to = "Sample", values_to = "Log2") %>%
    mutate(Study = str_extract(Sample, "^\\d{2}CPTAC_\\w+?_\\w+?_\\w+"))
  
  tmt_study_avg <- tmt_long %>%
    group_by(Gene, Study) %>%
    summarise(Log2 = mean(Log2, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Study, values_from = Log2)
  
  # ---- Load matching Protein Assembly summary ----
  pa_file <- list.files(pa_dir, pattern = paste0("^", study_id, ".*summary.*\\.tsv$"), full.names = TRUE)
  if (length(pa_file) == 0) {
    warning("No PA summary for ", study_id)
    next
  }
  
  pa_sum <- read.delim(pa_file[1], check.names = FALSE)
  colnames(pa_sum)[1] <- "Gene"
  
  # Use Distinct Peptides instead of Spectral Counts:
  distinct_cols <- grep("Distinct Peptides$", colnames(pa_sum), value = TRUE)
  pa_counts <- pa_sum %>%
    dplyr::select(Gene, all_of(distinct_cols)) %>%
    filter(!is.na(Gene) & Gene != "") %>%
    distinct(Gene, .keep_all = TRUE) %>%
    mutate(across(-Gene, as.numeric)) %>%
    rowwise() %>%
    mutate(TotalDistinctPeptides = sum(c_across(-Gene), na.rm = TRUE)) %>%
    ungroup() %>%
    dplyr::select(Gene, TotalDistinctPeptides)
  
  # ---- Align on gene symbols ----
  common_genes <- intersect(tmt_study_avg$Gene, pa_counts$Gene)
  if (length(common_genes) < 50) {
    warning("Too few common genes in ", study_id)
    next
  }
  
  tmt_common <- tmt_study_avg %>% filter(Gene %in% common_genes)
  pa_common  <- pa_counts %>% filter(Gene %in% common_genes)
  
  combined <- inner_join(tmt_common, pa_common, by = "Gene")
  
  # ---- Compute Spearman correlation per TMT Study ----
  study_cols <- setdiff(names(combined), c("Gene", "TotalDistinctPeptides"))
  cor_vals <- map_dbl(study_cols, function(col) {
    cor(combined[[col]], combined$TotalDistinctPeptides, method = "spearman", use = "pairwise.complete.obs")
  })
  
  all_corrs[[study_id]] <- tibble(
    Study = study_id,
    TMT_Study = study_cols,
    Spearman = cor_vals
  )
}

# ----------------------------------------
# Export results
# ----------------------------------------
qc_df <- bind_rows(all_corrs)
write.csv(qc_df, file.path(out_dir, "study_level_spearman_qc_DistinctPeptides.csv"), row.names = FALSE)

# ----------------------------------------
# Visualization
# ----------------------------------------

# Boxplot across studies
ggplot(qc_df, aes(x = Study, y = Spearman)) +
  geom_boxplot(outlier.shape = NA, fill = "lightblue") +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Boxplot of Spearman correlations (TMT vs Distinct Peptides)", 
       x = "Study", y = "Spearman correlation") +
  theme_minimal(base_size = 13)
ggsave(file.path(out_dir, "boxplot_spearman_TMT_vs_DistinctPeptides.png"), dpi = 300, width = 8, height = 6)

# Ridge density plot
ggplot(qc_df, aes(x = Spearman, y = Study, fill = Study)) +
  geom_density_ridges(alpha = 0.7) +
  labs(title = "Distribution of Spearman correlations per study", 
       x = "Spearman correlation", y = "") +
  theme_ridges() + 
  theme(legend.position = "none")
ggsave(file.path(out_dir, "ridgeplot_spearman_TMT_vs_DistinctPeptides.png"), dpi = 300, width = 8, height = 6)

# Faceted barplots per study (looped plots)
for (s in unique(qc_df$Study)) {
  p <- ggplot(filter(qc_df, Study == s), aes(x = Spearman, y = reorder(TMT_Study, Spearman))) +
    geom_col(fill = "steelblue") +
    labs(title = paste("Spearman correlation (TMT vs Distinct Peptides):", s),
         x = "Spearman correlation", y = "TMT Study") +
    theme_minimal(base_size = 10)
  
  ggsave(file.path(out_dir, paste0("per_study_spearman_", s, ".png")), plot = p, dpi = 300, width = 7, height = 6)
}