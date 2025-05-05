#!/usr/bin/env Rscript
#
# create_peptide_intensity_matrices_listFormat.R
#  - Stagewise: flatten at fraction level → nest by patient → build per-patient
#    fraction-level matrices
#  - Per-stage checkpoints so you can resume anywhere
#  - Progress bars via progress::progress_bar
#

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(furrr)
  library(progress)    # <-- for progress_bar
  library(qs)
})

# -----------------------------------------------------------------------------
# 0) Parallel backend
# -----------------------------------------------------------------------------
ncpus <- as.integer(
  Sys.getenv("SLURM_CPUS_PER_TASK", parallel::detectCores(logical = FALSE))
)
plan(multisession, workers = ncpus)

# -----------------------------------------------------------------------------
# 00) checkpoint helper
# -----------------------------------------------------------------------------
dir.create("checkpoints", showWarnings = FALSE)
chkpt <- function(path, code) {
  if (file.exists(path)) {
    message("⚡️  skipping (found) ", path)
    qread(path)
  } else {
    message("⚡️  running and saving ", path)
    out <- code()
    qsave(out, path)
    out
  }
}

# -----------------------------------------------------------------------------
# helpers
# -----------------------------------------------------------------------------
flatten_fractional <- function(study_tbl, study_id) {
  study_tbl %>%
    unnest(data) %>%
    unnest(aggregated) %>%
    transmute(
      study_id,
      case_id        = `Participant ID`,
      run_folder     = `Folder Name`,
      fraction       = file_name,
      PeptideSequence,
      total_intensity
    )
}

make_frac_matrix <- function(df) {
  df %>%
    unite(sample, study_id, run_folder, fraction, sep = "_", remove = FALSE) %>%
    select(PeptideSequence, sample, total_intensity) %>%
    pivot_wider(
      names_from   = sample,
      values_from  = total_intensity,
      values_fill  = list(total_intensity = 0L)
    )
}

# -----------------------------------------------------------------------------
# 1) Load cleaned nested TMT objects (with checkpoint)
# -----------------------------------------------------------------------------
ck1_10 <- "checkpoints/cleaned_TMT10.qs"
cleaned_TMT10 <- chkpt(ck1_10, function() {
    message("\n--- STAGE 1: Loading cleaned_TMT10 from workspace ---")
    ct10 <- qread("data_proteomics/workspaces/cleaned_TMT10.qs")
    message(sprintf("  • cleaned_TMT10: %d studies loaded", length(ct10)))
    ct10
  })

ck1_11 <- "checkpoints/cleaned_TMT11.qs"
cleaned_TMT11 <- chkpt(ck1_11, function() {
    message("\n--- STAGE 1: Loading cleaned_TMT11 from workspace ---")
    ct11 <- qread("data_proteomics/workspaces/cleaned_TMT11.qs")
    message(sprintf("  • cleaned_TMT11: %d studies loaded", length(ct11)))
    ct11
  })

# -----------------------------------------------------------------------------
# 2) Flatten TMT10 at fraction level
# -----------------------------------------------------------------------------
ck2 <- "checkpoints/flat_frac_TMT10.qs"
flat_TMT10 <- chkpt(ck2, function() {
  message("\n--- STAGE 2: Flattening TMT10 at fraction level ---")
  pb <- progress_bar$new(
    total = length(cleaned_TMT10),
    format = "  TMT10 [:bar] :percent  eta: :eta"
  )
  imap_dfr(cleaned_TMT10, function(st, id) {
    pb$tick(tokens = list())
    flatten_fractional(st, id)
  })
})
message("→ flat_frac_TMT10 rows = ", nrow(flat_TMT10))
rm(flat_TMT10); gc()

# -----------------------------------------------------------------------------
# 3) Flatten TMT11 at fraction level
# -----------------------------------------------------------------------------
ck3 <- "checkpoints/flat_frac_TMT11.qs"
flat_TMT11 <- chkpt(ck3, function() {
  message("\n--- STAGE 3: Flattening TMT11 at fraction level ---")
  pb <- progress_bar$new(
    total = length(cleaned_TMT11),
    format = "  TMT11 [:bar] :percent  eta: :eta"
  )
  imap_dfr(cleaned_TMT11, function(st, id) {
    pb$tick(tokens = list())
    flatten_fractional(st, id)
  })
})
message("→ flat_frac_TMT11 rows = ", nrow(flat_TMT11))
rm(flat_TMT11); gc()

# -----------------------------------------------------------------------------
# 4) Nest by patient for TMT10
# -----------------------------------------------------------------------------
ck4 <- "checkpoints/patient_list_frac_TMT10.qs"
patient_list_TMT10 <- chkpt(ck4, function() {
  message("\n--- STAGE 4: Nesting TMT10 by patient ---")
  flat <- qread("checkpoints/flat_frac_TMT10.qs")
  ids  <- unique(flat$case_id)
  pb   <- progress_bar$new(
    total = length(ids),
    format = "  TMT10 patients [:bar] :percent  eta: :eta"
  )
  lst <- map(ids, function(cid) {
    pb$tick()
    filter(flat, case_id == cid)
  })
  names(lst) <- ids
  lst
})
message("→ patients in TMT10 = ", length(patient_list_TMT10))
rm(patient_list_TMT10); gc()

# -----------------------------------------------------------------------------
# 5) Nest by patient for TMT11
# -----------------------------------------------------------------------------
ck5 <- "checkpoints/patient_list_frac_TMT11.qs"
patient_list_TMT11 <- chkpt(ck5, function() {
  message("\n--- STAGE 5: Nesting TMT11 by patient ---")
  flat <- qread("checkpoints/flat_frac_TMT11.qs")
  ids  <- unique(flat$case_id)
  pb   <- progress_bar$new(
    total = length(ids),
    format = "  TMT11 patients [:bar] :percent  eta: :eta"
  )
  lst <- map(ids, function(cid) {
    pb$tick()
    filter(flat, case_id == cid)
  })
  names(lst) <- ids
  lst
})
message("→ patients in TMT11 = ", length(patient_list_TMT11))
rm(patient_list_TMT11); gc()

# -----------------------------------------------------------------------------
# 6) Build per-patient fraction matrices for TMT10
# -----------------------------------------------------------------------------
ck6 <- "checkpoints/patient_fracmat_TMT10.qs"
patient_matrices_TMT10 <- chkpt(ck6, function() {
  message("\n--- STAGE 6: Building TMT10 fraction matrices ---")
  lst <- qread("checkpoints/patient_list_frac_TMT10.qs")
  pb  <- progress_bar$new(
    total = length(lst),
    format = "  TMT10 matrices [:bar] :percent  eta: :eta"
  )
  imap(lst, function(df, cid) {
    pb$tick()
    make_frac_matrix(df)
  })
})
message("→ built ", length(patient_matrices_TMT10), " TMT10 matrices")
rm(patient_matrices_TMT10); gc()

# -----------------------------------------------------------------------------
# 7) Build per-patient fraction matrices for TMT11
# -----------------------------------------------------------------------------
ck7 <- "checkpoints/patient_fracmat_TMT11.qs"
patient_matrices_TMT11 <- chkpt(ck7, function() {
  message("\n--- STAGE 7: Building TMT11 fraction matrices ---")
  lst <- qread("checkpoints/patient_list_frac_TMT11.qs")
  pb  <- progress_bar$new(
    total = length(lst),
    format = "  TMT11 matrices [:bar] :percent  eta: :eta"
  )
  imap(lst, function(df, cid) {
    pb$tick()
    make_frac_matrix(df)
  })
})
message("→ built ", length(patient_matrices_TMT11), " TMT11 matrices")
rm(patient_matrices_TMT11); gc()

message("\n✅ All stages complete!")
