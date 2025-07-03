# ---------------------------------------------------
# Proteomics QC pipeline for thesis (clean version with log2 transform)
# ---------------------------------------------------

library(tidyverse)
library(ggplot2)
library(pheatmap)
library(ggfortify)
library(reshape2)
library(impute)

# ---------------------------------------------------
# Load data
# ---------------------------------------------------

# Adjust path as needed
tmt_file <- "gene/TMT10_gene_PDC000116_combined.csv"

# Load TMT data
tmt <- read_csv(tmt_file, show_col_types = FALSE) %>%
  rename(Gene = 1) %>%
  filter(!is.na(Gene), Gene != "")

# ---------------------------------------------------
# Step 1 - Sample-level QC
# ---------------------------------------------------

# Compute sample-wise total abundance (raw, not log2 transformed here)
sample_sums <- tmt %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Abundance") %>%
  group_by(Sample) %>%
  summarize(Total_Abundance = sum(Abundance, na.rm = TRUE))

# Boxplot of total signal per sample with red outliers
ggplot(sample_sums, aes(y = Total_Abundance)) +
  geom_boxplot(fill = "skyblue", outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  labs(title = "Total TMT abundance per sample", y = "Total abundance") +
  theme_minimal(base_size = 14)

# ---------------------------------------------------
# Prepare data for downstream analyses with log2 transform (+1 pseudocount)
# ---------------------------------------------------

tmt_log <- tmt %>%
  mutate(across(-Gene, ~ log2(.x + 1)))

# ---------------------------------------------------
# Step 2 - PCA to check sample similarity
# ---------------------------------------------------

# Prepare matrix for PCA (log2 transformed, imputed)
expr_mat <- tmt_log %>%
  column_to_rownames("Gene") %>%
  as.matrix()

expr_mat_imputed <- impute.knn(expr_mat)$data

# Transpose so samples are rows for PCA
pca_res <- prcomp(t(expr_mat_imputed), scale. = TRUE)

# Import metadata
metadata <- read_delim("~/Desktop/master_thesis_local/master_thesis_scripts/data/raw/clinical_Pan-cancer.May2022.tsv",
                       delim = "\t", escape_double = FALSE,
                       trim_ws = TRUE)

# Extract sample short names (adjust regex to your naming convention)
sample_annot <- tibble(
  Sample = colnames(expr_mat),
  Sample_short = str_extract(colnames(expr_mat), "[0-9]{2}CO[0-9]{3}")
)

# Join with clinical metadata on Sample_short vs case_id
sample_annot <- left_join(sample_annot, metadata, by = c("Sample_short" = "case_id"))

# Check exact column name in metadata; adjust if needed
# For example, if column name is "Tumor Stage (Pathological)", you can rename it for easy access:
sample_annot <- sample_annot %>%
  rename(tumor_stage_pathological = `baseline/tumor_stage_pathological`)

# PCA plot colored by tumor stage
autoplot(pca_res, data = sample_annot, colour = "tumor_stage_pathological") +
  labs(title = "PCA of samples by tumor stage (log2 transformed data)") +
  theme_minimal(base_size = 14)



# Grouped by study
# Prepare matrix for PCA (log2 transformed)
expr_mat <- tmt_log %>%
  column_to_rownames("Gene") %>%
  as.matrix()

# Transpose so samples are rows for PCA
pca_res <- prcomp(t(expr_mat), scale. = TRUE)

# Extract clean study ID and short sample name for PCA
sample_annot <- tibble(
  Sample = colnames(expr_mat),
  Study = str_extract(colnames(expr_mat), "^\\d{2}CPTAC"),
  PatientID = str_extract(colnames(expr_mat), "[0-9]{2}CO[0-9]{3}")
)

# PCA plot (colored by study)
autoplot(pca_res, data = sample_annot, colour = "Study") +
  labs(title = "PCA of samples by study (log2 transformed data)") +
  theme_minimal(base_size = 14)


# ---------------------------------------------------
# Step 3 - Hierarchical clustering heatmap (cleaned version)
# ---------------------------------------------------

# Compute Spearman correlation matrix (log2 transformed data)
cor_mat <- cor(expr_mat, method = "spearman", use = "pairwise.complete.obs")

# Clean sample names for heatmap labels
clean_names <- colnames(expr_mat) %>%
  str_extract("[0-9]{2}CO[0-9]{3}")

# Add annotation for heatmap
annot <- tibble(Sample = colnames(expr_mat),
                Study = str_extract(colnames(expr_mat), "^\\d{2}CPTAC")) %>%
  column_to_rownames("Sample")

# Heatmap with cleaned labels and annotation
pheatmap(cor_mat,
         annotation_col = annot,
         labels_col = clean_names,
         labels_row = clean_names,
         fontsize = 4,
         main = "Sample clustering heatmap (log2 transformed data)",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         color = colorRampPalette(c("green", "white", "red"))(100))

# ---------------------------------------------------
# Step 4 - Protein-level QC (low abundance proteins)
# ---------------------------------------------------

# Calculate missingness (% missing per protein) — no transformation needed
protein_missing <- tmt %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Abundance") %>%
  group_by(Gene) %>%
  summarize(Missingness = mean(is.na(Abundance))) %>%
  arrange(desc(Missingness))

# Median abundance per protein (log2 transformed)
protein_abundance <- tmt_log %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Abundance") %>%
  group_by(Gene) %>%
  summarize(Median_Abundance = median(Abundance, na.rm = TRUE))

# Plot abundance distribution
ggplot(protein_abundance, aes(x = Median_Abundance)) +
  geom_histogram(bins = 50, fill = "skyblue") +
  labs(title = "Protein median abundance distribution (log2 transformed)", x = "Log2 abundance", y = "Count") +
  theme_minimal(base_size = 14)

# ---------------------------------------------------
# Filter proteins for further analysis (optional)
# ---------------------------------------------------

filtered_tmt <- tmt %>%
  filter(Gene %in% protein_missing$Gene[protein_missing$Missingness < 0.3])

# ---------------------------------------------------
# Step 5 - Variance vs abundance (log2 transformed)
# ---------------------------------------------------

protein_stats <- tmt_log %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Abundance") %>%
  group_by(Gene) %>%
  summarize(Mean_Abundance = mean(Abundance, na.rm = TRUE),
            Variance = var(Abundance, na.rm = TRUE))

ggplot(protein_stats, aes(x = Mean_Abundance, y = Variance)) +
  geom_point(alpha = 0.5) +
  labs(title = "Variance vs mean abundance (log2 transformed)", x = "Mean abundance", y = "Variance") +
  theme_minimal(base_size = 14)

# ---------------------------------------------------
# Step 6 - Dynamic range sanity check
# ---------------------------------------------------

# Prepare long data for dynamic range check (log2 transformed)
tmt_long <- tmt_log %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Abundance") %>%
  mutate(Sample_short = str_extract(Sample, "[0-9]{2}OV[0-9]{3}"))

# Calculate boxplot stats manually
bp_stats <- tmt_long %>%
  group_by(Sample_short) %>%
  summarize(Q1 = quantile(Abundance, 0.25, na.rm = TRUE),
            Q3 = quantile(Abundance, 0.75, na.rm = TRUE),
            IQR = Q3 - Q1,
            Lower = Q1 - 1.5 * IQR,
            Upper = Q3 + 1.5 * IQR)

# Join back to tag outliers
tmt_long <- tmt_long %>%
  left_join(bp_stats, by = "Sample_short") %>%
  mutate(IsOutlier = Abundance < Lower | Abundance > Upper)

# Plot with outliers
ggplot(tmt_long, aes(x = Sample_short, y = Abundance)) +
  geom_boxplot(fill = "lightblue", outlier.shape = NA) +  # boxplot without outliers
  geom_point(data = filter(tmt_long, IsOutlier), color = "red", size = 1, alpha = 0.6) +  # outliers in red
  coord_flip() +
  labs(title = "Abundance distributions across samples (log2 transformed)", y = "Abundance (log2 scale)", x = "Sample") +
  theme_minimal(base_size = 10)

# ---------------------------------------------------
# Step 7 - Export important summary QC tables
# ---------------------------------------------------

write_csv(sample_sums, file = "qc_outputs/sample_total_abundance.csv")
write_csv(protein_missing, file = "qc_outputs/protein_missingness.csv")
write_csv(protein_abundance, file = "qc_outputs/protein_abundance.csv")
write_csv(protein_stats, file = "qc_outputs/protein_variance.csv")

# ---------------------------------------------------
# Optional: Save all plots automatically to PDF (commented out)
# ---------------------------------------------------

# pdf("qc_outputs/proteomics_QC_report.pdf", width = 10, height = 8)
# print(...)
# dev.off()

# ---------------------------------------------------
cat("QC report completed! ✅\n")
