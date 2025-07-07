# Prepare overview storage
dim_overview <- tibble::tibble(
  pdc_id = character(),
  data_type = character(),
  n_features_before = integer(),
  n_samples_before = integer(),
  n_features_after = integer(),
  n_samples_after = integer()
)

for (pdc_id in pdc_ids) {
  cat("Checking", pdc_id, "\n")
  for (type in data_types) {
    pattern <- paste0("TMT\\d+_", type, "_", pdc_id, "_combined\\.csv$")
    files <- list.files(file.path(input_dir, type), pattern = pattern, full.names = TRUE)
    if (length(files) == 0) next
    expr_file <- files[1]
    expr_df <- readr::read_csv(expr_file, show_col_types = FALSE)
    expr_matrix <- as.matrix(expr_df[, -1])
    rownames(expr_matrix) <- expr_df$feature
    
    # Normalize
    norm_matrix <- normalizeCyclicLoess(expr_matrix, method = "fast")
    dim_before <- dim(norm_matrix)
    
    # Remove NA rownames
    na_rownames <- is.na(rownames(norm_matrix))
    if (any(na_rownames)) norm_matrix <- norm_matrix[!na_rownames, , drop = FALSE]
    
    # Sync annotation
    good_samples <- intersect(colnames(norm_matrix), sample_annotation$FullRunName)
    norm_matrix <- norm_matrix[, good_samples, drop = FALSE]
    sample_annotation_sub <- sample_annotation[sample_annotation$FullRunName %in% good_samples, ]
    sample_annotation_sub <- sample_annotation_sub[match(colnames(norm_matrix), sample_annotation_sub$FullRunName), ]
    stopifnot(all(sample_annotation_sub$FullRunName == colnames(norm_matrix)))
    
    # Batch feature filtering
    by_batch <- split(seq_along(colnames(norm_matrix)), sample_annotation_sub$batch)
    bad_features <- sapply(rownames(norm_matrix), function(f) {
      any(sapply(by_batch, function(idx) {
        all(is.na(norm_matrix[f, idx])) | sum(!is.na(norm_matrix[f, idx])) < 2
      }))
    })
    norm_matrix <- norm_matrix[!bad_features, , drop = FALSE]
    dim_after <- dim(norm_matrix)
    
    # Save row to overview
    dim_overview <- dplyr::add_row(
      dim_overview,
      pdc_id = pdc_id,
      data_type = type,
      n_features_before = dim_before[1],
      n_samples_before = dim_before[2],
      n_features_after = dim_after[1],
      n_samples_after = dim_after[2]
    )
  }
}

# Print as table
print(dim_overview)
# Optionally save to CSV
write.csv(dim_overview, "matrix_dimension_overview.csv", row.names = FALSE)
