library(tidyverse)
library(qs)

# 1) Read in your cleaned, nested TMT objects
#cleaned_TMT10 <- qread("data_proteomics/workspaces/cleaned_TMT10.qs")
#cleaned_TMT11 <- qread("data_proteomics/workspaces/cleaned_TMT11.qs")

# 2) A helper to flatten one study’s nested runs into a long tibble
flatten_study <- function(study_tbl, study_id) {
  study_tbl %>%
    # each row has cols: study_folder, data (a list of run‐tables)
    unnest(cols = c(data)) %>%          # expand runs → file_name, filepath, aggregated (list)
    unnest(cols = c(aggregated)) %>%    # expand aggregated → one row per peptide
    transmute(
      study_id       = study_id,
      case_id        = `Participant ID`,
      run_folder     = `Folder Name`,
      PeptideSequence,
      total_intensity
    )
}

# 3) Flatten ALL studies into one long table
flat_TMT10 <- imap_dfr(cleaned_TMT10, flatten_study)
flat_TMT11 <- imap_dfr(cleaned_TMT11, flatten_study)

# 4) Build a unique sample name and pivot wide
wide_TMT10 <- flat_TMT10 %>%
  unite(sample, study_id, run_folder, remove = FALSE) %>%
  pivot_wider(
    names_from   = sample,
    values_from  = total_intensity,
    values_fill  = list(total_intensity = 0L)
  )



wide_TMT11 <- flat_TMT11 %>%
  unite(sample, study_id, run_folder, remove = FALSE) %>%
  pivot_wider(
    names_from   = sample,
    values_from  = total_intensity,
    values_fill  = list(total_intensity = 0L)
  )


# 5) (Optional) Save out rectangular matrices
qsave(wide_TMT10, "data_proteomics/workspaces/peptide_matrix_TMT10.qs")
qsave(wide_TMT11, "data_proteomics/workspaces/peptide_matrix_TMT11.qs")
