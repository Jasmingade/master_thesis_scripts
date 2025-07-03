library(tidyverse)
library(readxl)
library(biomaRt)

# ---------------------------------------------------
# File paths (modify if needed)
# ---------------------------------------------------

tmt_file <- "gene/TMT10_gene_PDC000110_combined.csv"
ibaq_file <- "PDC000110_CPTAC2_OVarian_iBAQ.xlsx"
out_dir <- "qc_outputs"
dir.create(out_dir, showWarnings = FALSE)

# ---------------------------------------------------
# Step 1 — Load TMT data
# ---------------------------------------------------

tmt <- read_csv(tmt_file, show_col_types = FALSE) %>%
  rename(ensembl_gene_id = 1) %>%
  filter(!is.na(ensembl_gene_id), ensembl_gene_id != "") %>%
  mutate(ensembl_gene_id = sub("\\..*", "", ensembl_gene_id))

# Map Ensembl -> HGNC symbol
mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = tmt$ensembl_gene_id,
  mart = mart
)

# Merge gene symbols into TMT data
tmt_mapped <- tmt %>%
  left_join(gene_map, by = "ensembl_gene_id") %>%
  filter(hgnc_symbol != "") %>%
  dplyr::select(-ensembl_gene_id) %>%
  rename(Gene = hgnc_symbol) %>%
  relocate(Gene, .before = everything())

# Long format for TMT
tmt_long <- tmt_mapped %>%
  pivot_longer(-Gene, names_to = "TMT_Sample", values_to = "TMT_Value") %>%
  mutate(Sample_ID = str_extract(TMT_Sample, "[0-9]{2}OV[0-9]{3}"))

# ---------------------------------------------------
# Step 2 — Load iBAQ data
# ---------------------------------------------------

ibaq <- read_excel(ibaq_file, sheet = 2)
colnames(ibaq)[2] <- "Gene"

# Long format for iBAQ
ibaq_long <- ibaq %>%
  pivot_longer(-c(refseq_peptide, Gene), names_to = "Sample_ID", values_to = "iBAQ_Value") %>%
  mutate(Sample_ID = sub("^C", "", Sample_ID))

# ---------------------------------------------------
# Step 3 — Merge datasets
# ---------------------------------------------------

merged <- inner_join(tmt_long, ibaq_long, by = c("Gene", "Sample_ID")) %>%
  filter(!is.na(TMT_Value), !is.na(iBAQ_Value)) %>%
  dplyr::select(Gene, Sample_ID, TMT_Value, iBAQ_Value, refseq_peptide, TMT_Sample)

# How many samples overlap?
cat("Number of shared samples:", length(unique(merged$Sample_ID)), "\n")

# ---------------------------------------------------
# Step 4 — Calculate correlations
# ---------------------------------------------------

# Global Spearman
overall_cor <- cor(merged$TMT_Value, merged$iBAQ_Value, method = "spearman")
cat("Overall Spearman correlation:", round(overall_cor, 3), "\n")

# Per-sample Spearman
sample_corr <- merged %>%
  group_by(Sample_ID) %>%
  summarize(Spearman = cor(TMT_Value, iBAQ_Value, method = "spearman"))

# ---------------------------------------------------
# Step 5 — Plotting & saving
# ---------------------------------------------------

# Boxplot
p <- ggplot(sample_corr, aes(x = reorder(Sample_ID, Spearman), y = Spearman)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Per-sample Spearman correlation (TMT vs iBAQ)",
       x = "Sample", y = "Spearman correlation") +
  theme_minimal(base_size = 12)
# Save to PNG
ggsave(filename = file.path(out_dir, "per_sample_spearman_TMT_vs_iBAQ.png"), plot = p, width = 10, height = 8, dpi = 300)

# Boxplot
p_box <- ggplot(sample_corr, aes(x = "", y = Spearman)) +
  geom_boxplot(fill = "lightblue") +
  geom_jitter(width = 0.1, alpha = 0.6, size = 2) +
  labs(title = "Distribution of Spearman correlations (TMT vs iBAQ)",
       x = "", y = "Spearman correlation") +
  theme_minimal(base_size = 14)

# Save boxplot as PNG
ggsave(filename = file.path(out_dir, "spearman_boxplot_TMT_vs_iBAQ.png"), plot = p_box, width = 6, height = 6, dpi = 300)


# Also save numeric results
write.csv(sample_corr, file.path(out_dir, "spearman_per_sample_TMT_vs_iBAQ.csv"), row.names = FALSE)
