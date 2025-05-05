# ============================================
#  Parallel Clean‐&‐Collapse with future.apply + run-level progress
# ============================================

# 0) Load required libraries -------------------------------------------------
library(tidyverse)      # for dplyr, purrr, stringr, tibble
library(qs)             # for qread/qsave
library(future.apply)   # for future_lapply
library(progressr)      # for with_progress, progressor
library(stringr)        # for str_replace_all

# 1) Configure parallel backend & progress bar --------------------------------
#   Use all physical cores minus one to leave some for the OS
n_workers <- max(1, parallel::detectCores(logical = FALSE) - 34)
plan(multisession, workers = n_workers)

# Tell progressr to use a simple single-line bar
handlers("txtprogressbar")

# 2) Define the “clean & collapse” helper -------------------------------------
#    - strips any non-A–Z from PeptideSequence
#    - groups by Participant ID & Folder & CleanSeq
#    - sums intensities & counts PSMs
clean_and_collapse <- function(agg_df) {
  agg_df %>%
    mutate(
      CleanSeq = str_replace_all(PeptideSequence, "[^A-Z]", "")
    ) %>%
    group_by(`Participant ID`, `Folder Name`, CleanSeq) %>%
    summarize(
      total_intensity = sum(total_intensity, na.rm = TRUE),
      n_psm           = sum(n_psm,           na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(PeptideSequence = CleanSeq)
}

# 3) Read in your nested TMT objects -----------------------------------------
#message("📂 Loading nested TMT10/TMT11 objects…")
#nested_TMT10 <- qread("data_proteomics/workspaces/nested_TMT10.qs")
nested_TMT11 <- qread("data_proteomics/workspaces/nested_TMT11.qs")

# 4) Clean TMT10: one progress tick per run -----------------------------------
#    Count total number of runs (i.e. rows) across all studies
total_runs_T10 <- sum(map_int(nested_TMT10, ~ length(.x$data)))

message("🧹 Cleaning all TMT10 runs…")
with_progress({
  p_run <- progressor(steps = total_runs_T10)
  cleaned_TMT10 <- future_lapply(
    nested_TMT10,
    function(study_tbl) {
      # process each run in this study
      cleaned_runs <- lapply(
        study_tbl$data,
        function(run_tbl) {
          p_run()  # advance progress by one run
          run_tbl %>%
            mutate(
              aggregated = lapply(aggregated, clean_and_collapse)
            )
        }
      )
      # put cleaned runs back into the same study_tbl structure
      study_tbl %>% mutate(data = cleaned_runs)
    },
    future.seed = TRUE
  )
})
# preserve names
names(cleaned_TMT10) <- names(nested_TMT10)

# 5) Clean TMT11: same pattern -----------------------------------------------
total_runs_T11 <- sum(map_int(nested_TMT11, ~ length(.x$data)))

message("🧹 Cleaning all TMT11 runs…")
with_progress({
  p_run <- progressor(steps = total_runs_T11)
  cleaned_TMT11 <- future_lapply(
    nested_TMT11,
    function(study_tbl) {
      cleaned_runs <- lapply(
        study_tbl$data,
        function(run_tbl) {
          p_run()
          run_tbl %>%
            mutate(
              aggregated = lapply(aggregated, clean_and_collapse)
            )
        }
      )
      study_tbl %>% mutate(data = cleaned_runs)
    },
    future.seed = TRUE
  )
})
names(cleaned_TMT11) <- names(nested_TMT11)

# 6) Save out the cleaned objects --------------------------------------------
message("💾 Saving cleaned objects…")
qsave(cleaned_TMT10, "data_proteomics/workspaces/cleaned_TMT10.qs")
qsave(cleaned_TMT11, "data_proteomics/workspaces/cleaned_TMT11.qs")

message("✅ All runs cleaned & saved!")
