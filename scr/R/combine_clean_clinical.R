library(readr)
library(dplyr)
library(stringr)
library(tibble)
library(data.table)
library(readxl)

# Proteomics clinical
# === Load and clean clinical data ===
proteomics_clinical <- read_csv("data/processed/all_cancers_clinical_core.csv") %>%
  rename(
    case_id = Patient_ID,
    age = `consent/age`,
    sex = `consent/sex`,
    race = `consent/race`,
    ethnicity = `consent/ethnicity`,
    tumor_size_cm = `baseline/tumor_size_cm`,
    bmi = `medical_history/bmi`,
    alcohol_consumption = `medical_history/alcohol_consumption`,
    tobacco_smoking_history = `medical_history/tobacco_smoking_history`,
    vital_status = `follow-up/vital_status_at_date_of_last_contact`,
    days_to_last_contact = `follow-up/number_of_days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact`,
    days_to_death = `follow-up/number_of_days_from_date_of_initial_pathologic_diagnosis_to_date_of_death`,
    tumor_stage = `baseline/tumor_stage_pathological`,
    histologic_subtype = `baseline/histologic_type`,
    histologic_grade = `cptac_path/histologic_grade`,
    OS_time = `Overall survival, days`,
    OS_event = `Survival status (1, dead; 0, alive)`
  ) %>%
  mutate(
    age = suppressWarnings(as.numeric(age)),
    age_log2 = ifelse(!is.na(age) & age > 0, log2(age), NA),
    age_group = cut(age,
                    breaks = c(0, 40, 50, 60, 70, 80, Inf),
                    right = FALSE,
                    labels = c("<40", "40–49", "50–59", "60–69", "70–79", "80+"))
  )

# Transcriptomics clinical
transcriptomics_clinical <- read_excel("data/raw/metadata/transcriptomics_clinical.xlsx", 
                                       sheet = "TCGA-CDR") %>%
  rename(
    case_id = bcr_patient_barcode,
    tumor_code = type,
    age = age_at_initial_pathologic_diagnosis,
    sex = gender,
    race = race,
    vital_status = vital_status,
    days_to_last_contact = last_contact_days_to,
    days_to_death = death_days_to,
    tumor_stage = ajcc_pathologic_tumor_stage,
    histologic_subtype = histological_type,
    histologic_grade = histological_grade,
    OS_time = OS.time,
    OS_event = OS
  ) %>%
  mutate(
    age = suppressWarnings(as.numeric(age)),
    age_log2 = ifelse(!is.na(age) & age > 0, log2(age), NA),
    age_group = cut(age, breaks = c(0, 40, 50, 60, 70, 80, Inf), right = FALSE,
                    labels = c("<40", "40–49", "50–59", "60–69", "70–79", "80+"))
  )


# find common columns
common_cols <- intersect(colnames(proteomics_clinical), colnames(transcriptomics_clinical))

# subset to common columns
prot_clin_common <- proteomics_clinical[, common_cols]
trans_clin_common <- transcriptomics_clinical[, common_cols]

# Row-bind both datasets
metaData_multiomics <- bind_rows(prot_clin_common, trans_clin_common)

# Check for duplicates or NA case IDs
metaData_multiomics <- metaData_multiomics %>% distinct(case_id, .keep_all = TRUE) %>% filter(!is.na(case_id))

# Clean and Group Tumor Stages
metaData_multiomics <- metaData_multiomics %>%
  mutate(
    tumor_stage_clean = case_when(
      tumor_stage %in% c("Stage 0", "IS")                       ~ "Stage 0",
      tumor_stage %in% c("Stage I", "Stage IA", "Stage IB")     ~ "Stage I",
      tumor_stage %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC") ~ "Stage II",
      tumor_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC") ~ "Stage III",
      tumor_stage %in% c("Stage IV", "Stage IVA", "Stage IVB", "Stage IVC") ~ "Stage IV",
      tumor_stage %in% c("[Unknown]", "[Not Applicable]", "[Not Available]", 
                         "[Discrepancy]", "Stage X", 
                         "Staging is not applicable or unknown", 
                         "I/II NOS") ~ "Unknown",
      TRUE ~ "Other"
    )
  )
write_csv(metaData_multiomics, "data/processed/metaData_multiomics.csv")
