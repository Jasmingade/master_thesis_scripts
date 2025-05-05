#!/usr/bin/env Rscript
#
# create_peptide_intensity_matrices_listFormat.R
#
#  - Stagewise flatten → nest → matrix creation, with qsave checkpoints
#  - Live txt progress bars via progressr
#

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(furrr)
  library(progressr)
  library(qs)
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
chkpt <- function(path, expr) {
  if (file.exists(path)) {
    message("🔥 loading checkpoint: ", path)
    qread(path)
  } else {
    message("⚡️ generating checkpoint: ", path)
    out <- expr()
    qsave(out, path)
    out
  }
}

# -----------------------------------------------------------------------------
# helper: flatten one study’s nested runs into a long table
# -----------------------------------------------------------------------------
flatten_study <- function(study_tbl, study_id) {
  study_tbl %>%
    unnest(data) %>%              # expand runs
    unnest(aggregated) %>%        # expand peptides
    transmute(
      study_id,
      case_id        = `Participant ID`,
      run_folder     = `Folder Name`,
      PeptideSequence,
      total_intensity
    )
}

# -----------------------------------------------------------------------------
# 1) Load cleaned nested TMT objects
# -----------------------------------------------------------------------------
message("\n-----------------------------------------------------------------------------
⏳  STAGE 1: Load cleaned nested TMT objects
-----------------------------------------------------------------------------")
cleaned_TMT10 <- qread("data_proteomics/workspaces/cleaned_TMT10.qs")
message(sprintf("✅  cleaned_TMT10: %d studies loaded", length(cleaned_TMT10)))
cleaned_TMT11 <- qread("data_proteomics/workspaces/cleaned_TMT11.qs")
message(sprintf("✅  cleaned_TMT11: %d studies loaded", length(cleaned_TMT11)))

# -----------------------------------------------------------------------------
# 2) Flatten ALL TMT10 studies with progress
# -----------------------------------------------------------------------------
message("\n-----------------------------------------------------------------------------
⏳  STAGE 2: Flattening TMT10 studies
-----------------------------------------------------------------------------")
flat_TMT10 <- chkpt("checkpoints/flat_TMT10.qs", function() {
  with_progress({
    p <- progressor(along = cleaned_TMT10)
    imap_dfr(cleaned_TMT10, function(study_tbl, study_id) {
      p(sprintf("TMT10 %s", study_id))
      flatten_study(study_tbl, study_id)
    })
  })
})
message("✅  STAGE 2: flat_TMT10 rows = ", nrow(flat_TMT10))
rm(cleaned_TMT10); gc()

# -----------------------------------------------------------------------------
# 3) Flatten ALL TMT11 studies with progress
# -----------------------------------------------------------------------------
message("\n-----------------------------------------------------------------------------
⏳  STAGE 3: Flattening TMT11 studies
-----------------------------------------------------------------------------")
flat_TMT11 <- chkpt("checkpoints/flat_TMT11.qs", function() {
  with_progress({
    p <- progressor(along = cleaned_TMT11)
    imap_dfr(cleaned_TMT11, function(study_tbl, study_id) {
      p(sprintf("TMT11 %s", study_id))
      flatten_study(study_tbl, study_id)
    })
  })
})
message("✅  STAGE 3: flat_TMT11 rows = ", nrow(flat_TMT11))
rm(cleaned_TMT11); gc()

# -----------------------------------------------------------------------------
# 4) Nest by patient for TMT10
# -----------------------------------------------------------------------------
message("\n-----------------------------------------------------------------------------
⏳  STAGE 4: Nesting TMT10 by patient
-----------------------------------------------------------------------------")
patient_list_TMT10 <- chkpt("checkpoints/patient_list_TMT10.qs", function() {
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
message("✅  STAGE 4: patients in TMT10 = ", length(patient_list_TMT10))
rm(flat_TMT10); gc()

# -----------------------------------------------------------------------------
# 5) Nest by patient for TMT11
# -----------------------------------------------------------------------------
message("\n-----------------------------------------------------------------------------
⏳  STAGE 5: Nesting TMT11 by patient
-----------------------------------------------------------------------------")
patient_list_TMT11 <- chkpt("checkpoints/patient_list_TMT11.qs", function() {
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
message("✅  STAGE 5: patients in TMT11 = ", length(patient_list_TMT11))
rm(flat_TMT11); gc()

# -----------------------------------------------------------------------------
# helper: build a wide peptide×run matrix for one patient
# -----------------------------------------------------------------------------
make_matrix <- function(df) {
  df %>%
    unite(sample, study_id, run_folder, remove = FALSE) %>%
    pivot_wider(
      names_from   = sample,
      values_from  = total_intensity,
      values_fill  = list(total_intensity = 0L)   # ← must be a named list
    )
}

# -----------------------------------------------------------------------------
# 6) Build per-patient matrices for TMT10
# -----------------------------------------------------------------------------
message("\n-----------------------------------------------------------------------------
⏳  STAGE 6: Building per-patient matrices for TMT10
-----------------------------------------------------------------------------")
patient_matrices_TMT10 <- chkpt("checkpoints/patient_matrices_TMT10.qs", function() {
  with_progress({
    p <- progressor(along = patient_list_TMT10)
    imap(patient_list_TMT10, function(df, cid) {
      p(sprintf("Matrix TMT10 %s", cid))
      make_matrix(df)
    })
  })
})
message("✅  STAGE 6: built matrices for TMT10 patients")
rm(patient_list_TMT10); gc()

# -----------------------------------------------------------------------------
# 7) Build per-patient matrices for TMT11
# -----------------------------------------------------------------------------
message("\n-----------------------------------------------------------------------------
⏳  STAGE 7: Building per-patient matrices for TMT11
-----------------------------------------------------------------------------")
patient_matrices_TMT11 <- chkpt("checkpoints/patient_matrices_TMT11.qs", function() {
  with_progress({
    p <- progressor(along = patient_list_TMT11)
    imap(patient_list_TMT11, function(df, cid) {
      p(sprintf("Matrix TMT11 %s", cid))
      make_matrix(df)
    })
  })
})
message("✅  STAGE 7: built matrices for TMT11 patients")
rm(patient_list_TMT11); gc()

# -----------------------------------------------------------------------------
# 8) Final save of everything into your workspace
# -----------------------------------------------------------------------------
message("\n-----------------------------------------------------------------------------
⏳  STAGE 8: Saving final objects to workspace
-----------------------------------------------------------------------------")
qsave(patient_list_TMT10,     "data_proteomics/workspaces/patient_list_TMT10.qs")
qsave(patient_list_TMT11,     "data_proteomics/workspaces/patient_list_TMT11.qs")
qsave(patient_matrices_TMT10, "data_proteomics/workspaces/patient_matrices_TMT10.qs")
qsave(patient_matrices_TMT11, "data_proteomics/workspaces/patient_matrices_TMT11.qs")
message("✅  STAGE 8: all done!")
