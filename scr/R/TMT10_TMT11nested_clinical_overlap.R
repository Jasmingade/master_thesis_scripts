library(tidyverse)
library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(qs)


#############################################
### 1. Load the Clinical Data  ###
#############################################
clinical_data <- read_delim("data_proteomics/clinical_metadata/clinical_Pan-cancer.May2022.tsv", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)
# filter clinical data columns
clinical_data <- clinical_data %>%
  select(
    case_id,
    tumor_code,
    `consent/age`,
    `consent/sex`,
    `consent/race`,
    `consent/ethnicity`,
    `Inferred ancestry`,
    `baseline/tumor_stage_pathological`,
    `baseline/histologic_type`,
    `cptac_path/histologic_grade`,
    `baseline/tumor_size_cm`,
    `medical_history/bmi`,
    `medical_history/alcohol_consumption`,
    `medical_history/tobacco_smoking_history`,
    `follow-up/vital_status_at_date_of_last_contact`,
    `follow-up/number_of_days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact`,
    `follow-up/number_of_days_from_date_of_initial_pathologic_diagnosis_to_date_of_death`,
    `Survival status (1, dead; 0, alive)`,
    `Overall survival, days`
  )

#############################################
### 2. Open nested tibbles ###
#############################################
nested_TMT10 <- qread("data_proteomics/workspaces/nested_TMT10.qs")
nested_TMT11 <- qread("data_proteomics/workspaces/nested_TMT11.qs")


####################################################################
### 3. Write a Function to Summarize Overlap With clinical_data ###
####################################################################

summarize_overlap <- function(study_data, clinical_data) {
  # 1. Extract unique case IDs from this study’s aggregated data
  agg_ids <- study_data %>%
    distinct(`Participant ID`) %>%
    pull(`Participant ID`)
  
  # 2. Filter clinical_data to keep only rows that match these case IDs
  clin_matched <- clinical_data %>%
    filter(case_id %in% agg_ids)
  
  # 3. Filter further to those with non-missing vital status & time-to-event
  clin_valid <- clin_matched %>%
    filter(
      !is.na(`follow-up/vital_status_at_date_of_last_contact`) & 
        `follow-up/vital_status_at_date_of_last_contact` != "",
      !is.na(`follow-up/number_of_days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact`) &
        `follow-up/number_of_days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact` != ""
    )
  
  # 4. Return the summary as a tibble
  tibble(
    total_agg_cases = length(agg_ids),
    total_matched   = nrow(clin_matched),
    total_valid     = nrow(clin_valid)
  )
}

###################################################
### 4. Apply This Function to Each Study Tibble ###
###################################################
results <- map_dfr(
  names(nested_TMT10),
  function(study_name) {
    # Pull the study's nested tibble
    nested_study <- nested_TMT10[[study_name]]
    
    # Flatten to access all psm files with their aggregated data
    flat_study <- nested_study %>%
      unnest(data) %>%
      unnest(aggregated)
    
    # Run summary stats
    agg_ids <- flat_study %>%
      distinct(`Participant ID`) %>%
      pull(`Participant ID`)
    
    clin_matched <- clinical_data %>%
      filter(case_id %in% agg_ids)
    
    valid_ids <- clin_matched %>%
      filter(
        !is.na(`follow-up/vital_status_at_date_of_last_contact`) &
          `follow-up/vital_status_at_date_of_last_contact` != "",
        !is.na(`follow-up/number_of_days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact`) &
          `follow-up/number_of_days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact` != ""
      ) %>%
      pull(case_id)
    
    # Count unique files associated with those valid case IDs
    valid_files <- flat_study %>%
      filter(`Participant ID` %in% valid_ids) %>%
      distinct(`Folder Name`, file_name) %>%
      nrow()
    
    # Get the most common tumor_code in this matched set (there should usually only be one)
    tumor_code <- clin_matched %>%
      count(tumor_code, sort = TRUE) %>%
      slice(1) %>%
      pull(tumor_code)
    
    # Return full summary
    tibble(
      pdc_study_id   = study_name,
      tumor_code     = tumor_code,
      total_agg_cases = length(agg_ids),
      total_matched   = nrow(clin_matched),
      total_valid     = length(valid_ids),
      valid_files     = valid_files
    )
  }
)
# results
print(results)



################################################
### FILTER CLINICAL DATA AND PROTEOMICS DATA ###
################################################
library(tidyverse)
library(readr)
library(tidyr)
library(dplyr)
library(fs)

# --- Step 1: Define folder to store filtered outputs ---
output_dir <- "nested_TMT10_chunks"
dir_create(output_dir)

# --- Step 2: Get all unique proteomics Participant IDs ---
# We loop one-by-one through nested_TMT10 to extract them without blowing up memory
all_proteomics_ids <- map(nested_TMT10, function(study) {
  study %>%
    unnest(data) %>%
    unnest(aggregated) %>%
    distinct(`Participant ID`)
}) %>%
  bind_rows() %>%
  distinct(`Participant ID`) %>%
  pull(`Participant ID`)

# --- Step 3: Filter clinical_data ---
clinical_data_filtered <- clinical_data %>%
  filter(case_id %in% all_proteomics_ids)

# --- Step 4: filter clinical data to including only those with survival metadata
clinical_data_valid_survival <- clinical_data_filtered %>%
  filter(
    !is.na(`follow-up/vital_status_at_date_of_last_contact`) & 
      `follow-up/vital_status_at_date_of_last_contact` != "",
    !is.na(`follow-up/number_of_days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact`) &
      `follow-up/number_of_days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact` != ""
  )
table(clinical_data_valid_survival$tumor_code)


# --- Step 4: Save filtered clinical data ---
write_csv(clinical_data_valid_survival, file.path(output_dir, "clinical_data_matched_to_proteomics_TMT10.csv"))

