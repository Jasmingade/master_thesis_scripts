#!/usr/bin/env Rscript
#
# create_peptide_intensity_matrices_listFormat.R
#  - Stagewise flatten → nest → matrix creation, with qsave checkpoints
#  - Uses group_nest() for efficient per‐patient splitting
#  - Live txt progress bars via progressr
#

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(furrr)
  library(progressr)
  library(qs)
  library(tibble)
})

# -----------------------------------------------------------------------------
# 0) parallel back-end & simple progress bars
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
# 1) Load cleaned nested TMT objects
# -----------------------------------------------------------------------------
#message("
#-----------------------------------------------------------------------------
#⏳  STAGE 1: Loading cleaned nested TMT objects
#-----------------------------------------------------------------------------")
#cleaned_TMT10 <- qread("data_proteomics/workspaces/cleaned_TMT10.qs")
#message(sprintf("✅  cleaned_TMT10: %d studies loaded", length(cleaned_TMT10)))
#cleaned_TMT11 <- qread("data_proteomics/workspaces/cleaned_TMT11.qs")
#message(sprintf("✅  cleaned_TMT11: %d studies loaded", length(cleaned_TMT11)))

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
# 2) Flatten ALL TMT10 studies with progress
# -----------------------------------------------------------------------------
message("
-----------------------------------------------------------------------------
⏳  STAGE 2: Flattening TMT10 studies
-----------------------------------------------------------------------------")
flat_TMT10 <- chkpt("checkpoints/flat_TMT10.qs", function(){
  # load nested TMT10 only when rebuilding
  cleaned_TMT10 <- qread("data_proteomics/workspaces/cleaned_TMT10.qs")
  with_progress({
    p <- progressor(along = cleaned_TMT10)
    imap_dfr(cleaned_TMT10, ~ { p(.y); flatten_study(.x, .y) })
  })
})
message(sprintf("✅  STAGE 2: flat_TMT10 rows = %d", nrow(flat_TMT10)))
rm(cleaned_TMT10); gc()


# -----------------------------------------------------------------------------
# 3) Flatten ALL TMT11 studies with progress
# -----------------------------------------------------------------------------
message("
-----------------------------------------------------------------------------
⏳  STAGE 3: Flattening TMT11 studies
-----------------------------------------------------------------------------")
flat_TMT11 <- chkpt("checkpoints/flat_TMT11.qs", function(){
  cleaned_TMT11 <- qread("data_proteomics/workspaces/cleaned_TMT11.qs")
  with_progress({
    p <- progressor(along = cleaned_TMT11)
    imap_dfr(cleaned_TMT11, ~ { p(.y); flatten_study(.x, .y) })
  })
})
message(sprintf("✅  STAGE 3: flat_TMT11 rows = %d", nrow(flat_TMT11)))
rm(cleaned_TMT11); gc()

# -----------------------------------------------------------------------------
# 4) Nest by patient for TMT10 (one-pass grouping)
# -----------------------------------------------------------------------------
message("
-----------------------------------------------------------------------------
⏳  STAGE 4: Nesting TMT10 by patient
-----------------------------------------------------------------------------")
patient_list_TMT10 <- chkpt("checkpoints/patient_list_TMT10.qs", function() {
  flat_TMT10 %>%
    group_nest(case_id) %>%  # split into a tibble(case_id, data)
    deframe()                # convert to named list: names=case_id, values=data
})
message(sprintf("✅  STAGE 4: patients in TMT10 = %d", length(patient_list_TMT10)))
rm(flat_TMT10); gc()

# -----------------------------------------------------------------------------
# 5) Nest by patient for TMT11 (one-pass grouping)
# -----------------------------------------------------------------------------
message("
-----------------------------------------------------------------------------
⏳  STAGE 5: Nesting TMT11 by patient
-----------------------------------------------------------------------------")
patient_list_TMT11 <- chkpt("checkpoints/patient_list_TMT11.qs", function() {
  flat_TMT11 %>%
    group_nest(case_id) %>%
    deframe()
})
message(sprintf("✅  STAGE 5: patients in TMT11 = %d", length(patient_list_TMT11)))
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
      values_fill  = 0L
    )
}

# -----------------------------------------------------------------------------
# 6) Build per-patient matrices for TMT10
# -----------------------------------------------------------------------------
message("
-----------------------------------------------------------------------------
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
message("
-----------------------------------------------------------------------------
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
message("
-----------------------------------------------------------------------------
⏳  STAGE 8: Saving final objects to workspace
-----------------------------------------------------------------------------")
qsave(patient_list_TMT10,     "data_proteomics/workspaces/patient_list_TMT10.qs")
qsave(patient_list_TMT11,     "data_proteomics/workspaces/patient_list_TMT11.qs")
qsave(patient_matrices_TMT10, "data_proteomics/workspaces/patient_matrices_TMT10.qs")
qsave(patient_matrices_TMT11, "data_proteomics/workspaces/patient_matrices_TMT11.qs")
message("✅  STAGE 8: all done!")
