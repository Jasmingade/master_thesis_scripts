#!/usr/bin/env Rscript
#
# create_peptide_intensity_matrices_listFormat.R
#
#  - Stagewise: flatten at the fraction level → nest by patient → build per-patient
#    fraction-level matrices
#  - qsave() checkpoints so you can resume
#  - txt progress bars via progressr
#

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(furrr)
  library(progressr)
  library(qs)
  library(tibble)
  library(tidyverse)
})

# -----------------------------------------------------------------------------
# 0) parallel back-end & progress bar style
# -----------------------------------------------------------------------------
ncpus <- as.integer(
  Sys.getenv("SLURM_CPUS_PER_TASK", parallel::detectCores(logical = FALSE))
)
plan(multisession, workers = ncpus)
handlers("txtprogressbar")

# -----------------------------------------------------------------------------
# 00) checkpoint helper
# -----------------------------------------------------------------------------
dir.create("checkpoints", showWarnings = FALSE)
chkpt <- function(path, code) {
  if (file.exists(path)) {
    message("🔥 loading checkpoint: ", path)
    qread(path)
  } else {
    message("⚡️ generating checkpoint: ", path)
    out <- code()
    qsave(out, path)
    out
  }
}

# -----------------------------------------------------------------------------
# helper: flatten one study’s nested runs into a long table, keeping fraction
# -----------------------------------------------------------------------------
flatten_fractional <- function(study_tbl, study_id) {
  study_tbl %>%
    unnest(data) %>%              # get file_name & aggregated list
    unnest(aggregated) %>%        # get one row per peptide per fraction
    transmute(
      study_id,
      case_id        = `Participant ID`,
      run_folder     = `Folder Name`,
      fraction       = file_name,         # <-- preserve each .psm
      PeptideSequence,
      total_intensity
    )
}

# -----------------------------------------------------------------------------
# 1) Load your cleaned nested TMT10 & TMT11 objects
# -----------------------------------------------------------------------------
#message("\n--- STAGE 1: Loading cleaned nested TMT objects ---")
#cleaned_TMT10 <- qread("data_proteomics/workspaces/cleaned_TMT10.qs")
#cleaned_TMT11 <- qread("data_proteomics/workspaces/cleaned_TMT11.qs")
#message(sprintf("  • cleaned_TMT10: %d studies", length(cleaned_TMT10)))
#message(sprintf("  • cleaned_TMT11: %d studies", length(cleaned_TMT11)))

# -----------------------------------------------------------------------------
# 2) Flatten ALL TMT10 studies at the **fraction** level
# -----------------------------------------------------------------------------
message("\n--- STAGE 2: Flattening TMT10 at fraction level ---")
flat_TMT10 <- chkpt("checkpoints/flat_frac_TMT10.qs", function() {
  with_progress({
    p <- progressor(along = cleaned_TMT10)
    imap_dfr(cleaned_TMT10, function(st, id) {
      p(sprintf("TMT10 %s", id))
      flatten_fractional(st, id)
    })
  })
})
message("  → flat_TMT10 rows =", nrow(flat_TMT10))
rm(cleaned_TMT10); gc()

# -----------------------------------------------------------------------------
# 3) Flatten ALL TMT11 studies at the **fraction** level
# -----------------------------------------------------------------------------
message("\n--- STAGE 3: Flattening TMT11 at fraction level ---")
flat_TMT11 <- chkpt("checkpoints/flat_frac_TMT11.qs", function() {
  with_progress({
    p <- progressor(along = cleaned_TMT11)
    imap_dfr(cleaned_TMT11, function(st, id) {
      p(sprintf("TMT11 %s", id))
      flatten_fractional(st, id)
    })
  })
})
message("  → flat_TMT11 rows =", nrow(flat_TMT11))
rm(cleaned_TMT11); gc()

# -----------------------------------------------------------------------------
# 4) Nest by patient (case_id) — now each tibble has a 'fraction' column
# -----------------------------------------------------------------------------
message("\n--- STAGE 4: Nesting by patient (fraction-aware) ---")
patient_list_TMT10 <- chkpt("checkpoints/patient_list_frac_TMT10.qs", function() {
  ids <- unique(flat_TMT10$case_id)
  with_progress({
    p <- progressor(steps = length(ids))
    lst <- map(ids, function(cid) {
      p(sprintf("TMT10 patient %s", cid))
      filter(flat_TMT10, case_id == cid)
    })
    names(lst) <- ids
    lst
  })
})
message("  → TMT10 patients =", length(patient_list_TMT10))
rm(flat_TMT10); gc()

patient_list_TMT11 <- chkpt("checkpoints/patient_list_frac_TMT11.qs", function() {
  ids <- unique(flat_TMT11$case_id)
  with_progress({
    p <- progressor(steps = length(ids))
    lst <- map(ids, function(cid) {
      p(sprintf("TMT11 patient %s", cid))
      filter(flat_TMT11, case_id == cid)
    })
    names(lst) <- ids
    lst
  })
})
message("  → TMT11 patients =", length(patient_list_TMT11))
rm(flat_TMT11); gc()

# -----------------------------------------------------------------------------
# helper: build a per-patient fraction-level wide matrix
# -----------------------------------------------------------------------------
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
# 5) Build per-patient **fraction-level** matrices for TMT10
# -----------------------------------------------------------------------------
message("\n--- STAGE 5: Building per-patient fraction matrices for TMT10 ---")
patient_mat_TMT10 <- chkpt("checkpoints/patient_fracmat_TMT10.qs", function() {
  with_progress({
    p <- progressor(along = patient_list_TMT10)
    imap(patient_list_TMT10, function(df, cid) {
      p(sprintf("Patient %s", cid))
      make_frac_matrix(df)
    })
  })
})
message("  → built matrices for", length(patient_mat_TMT10), "TMT10 patients")
rm(patient_list_TMT10); gc()

# -----------------------------------------------------------------------------
# 6) Build per-patient **fraction-level** matrices for TMT11
# -----------------------------------------------------------------------------
message("\n--- STAGE 6: Building per-patient fraction matrices for TMT11 ---")
patient_mat_TMT11 <- chkpt("checkpoints/patient_fracmat_TMT11.qs", function() {
  with_progress({
    p <- progressor(along = patient_list_TMT11)
    imap(patient_list_TMT11, function(df, cid) {
      p(sprintf("Patient %s", cid))
      make_frac_matrix(df)
    })
  })
})
message("  → built matrices for", length(patient_mat_TMT11), "TMT11 patients")
rm(patient_list_TMT11); gc()

# -----------------------------------------------------------------------------
# 7) Save final objects
# -----------------------------------------------------------------------------
message("\n--- STAGE 7: Saving final fraction-level matrices ---")
qsave(patient_mat_TMT10, "data_proteomics/workspaces/patient_fracmat_TMT10.qs")
qsave(patient_mat_TMT11, "data_proteomics/workspaces/patient_fracmat_TMT11.qs")
message("✅  All done!")
