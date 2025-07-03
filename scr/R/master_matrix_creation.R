# Load data
load("iso_gene_mtx_fraction_TMT11_log2.RData")

# List of data types
data_types <- list(
  iso_log = iso_mtx_TMT11_log,
  iso_frac = iso_fraction_TMT11,
  gene = gene_mtx_TMT11
)

for (type_name in names(data_types)) {
  data_type <- data_types[[type_name]]
  
  for (study_name in names(data_type)) {
    study_matrices <- data_type[[study_name]]
    
    # Check all are matrices and not NULL
    if (!all(sapply(study_matrices, is.matrix))) {
      message(paste("Skipping", type_name, study_name, "- non-matrix found"))
      next
    }
    
    # Remove REF columns and rename columns to avoid duplicate sample names
    cleaned <- lapply(names(study_matrices), function(name) {
      mat <- study_matrices[[name]]
      ref_cols <- grep("^REF_", colnames(mat), value = TRUE)
      mat <- mat[, !colnames(mat) %in% ref_cols, drop = FALSE]
      # Add prefix to sample names to make them unique
      colnames(mat) <- paste(name, colnames(mat), sep = "_")
      return(mat)
    })
    names(cleaned) <- names(study_matrices)
    
    # Get common isoform IDs
    common_rows <- Reduce(intersect, lapply(cleaned, rownames))
    if (length(common_rows) == 0) {
      message(paste("No common isoforms in", type_name, study_name))
      next
    }
    
    # Align and combine
    aligned <- lapply(cleaned, function(mat) mat[common_rows, , drop = FALSE])
    combined <- tryCatch({
      do.call(cbind, aligned)
    }, error = function(e) {
      message(paste("Still failed combining", type_name, study_name))
      return(NULL)
    })
    
    # Save
    if (!is.null(combined)) {
      file_name <- paste0("TMT11_", type_name, "_", study_name, "_combined.csv")
      write.csv(
        cbind(feature = rownames(combined), combined),
        file = file_name,
        row.names = FALSE
      )
      message(paste("Saved:", file_name))
    }
  }
}
