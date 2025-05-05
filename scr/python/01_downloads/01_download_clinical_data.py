#!/usr/bin/env python

"""
Download and process clinical datasets from the CPTAC package.

This script performs the following steps:

1. Instantiates CPTAC objects for multiple tumour cohorts.
2. Ensures an output directory exists for saving raw clinical tables.
3. Iterates through each cohort to download (if necessary) and save its clinical table as a CSV.
4. Concatenates all individual clinical tables into a single pan-cancer table.
5. Exports the merged table to disk.
6. Defines a mapping of desired variables to their CPTAC column names.
7. Selects, renames, and orders only those variables present in the merged table.
8. Saves the resulting “core” clinical table and reports any requested variables that were missing.
"""

import cptac
import pathlib
import pandas as pd

# ------------------------------------------------------------------
# 1. Instantiate CPTAC tumour cohort datasets
# ------------------------------------------------------------------
datasets = {
    "ov":    cptac.Ov(),
    "coad":  cptac.Coad(),
    "brca":  cptac.Brca(),
    "ucec":  cptac.Ucec(),
    "ccrcc": cptac.Ccrcc(),
    "luad":  cptac.Luad(),
    "gbm":   cptac.Gbm(),
    "hnscc": cptac.Hnscc(),
    "lscc":  cptac.Lscc(),
    "pdac":  cptac.Pdac(),
}

# ------------------------------------------------------------------
# 2. Prepare the output directory for raw clinical tables
# ------------------------------------------------------------------
out_dir = pathlib.Path("data/raw/clinical_tables")
out_dir.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------
# 3. Download and save each cohort’s clinical table
# ------------------------------------------------------------------
for code, ds in datasets.items():
    clin = ds.get_clinical("mssm")
    csv_path = out_dir / f"{code}_clinical.csv"
    clin.to_csv(csv_path, index=False)
    print(f"✔ {code.upper():5s}: {csv_path} ({clin.shape[0]} patients)")

# ------------------------------------------------------------------
# 3b. Build a pan-cancer clinical table by concatenation
# ------------------------------------------------------------------
frames = []
for code, ds in datasets.items():
    clin = ds.get_clinical("mssm").copy()
    clin["cancer"] = code
    frames.append(clin)

all_clinical = pd.concat(frames, axis=0, join="outer", ignore_index=False)

# ------------------------------------------------------------------
# 4. Export the merged pan-cancer clinical table
# ------------------------------------------------------------------
merged_path = out_dir / "all_cancers_clinical.csv"
all_clinical.to_csv(merged_path, index=True)
print(
    "✔ Pan-cancer clinical table written to:", merged_path,
    f"({all_clinical.shape[0]} patients, {all_clinical.shape[1]} columns)"
)

# ------------------------------------------------------------------
# 5. Define mapping from desired variable names to CPTAC column names
# ------------------------------------------------------------------
rename_map = {
    # Identifiers
    "case_id":      "Patient_ID",
    "tumor_code":   "tumor_code",
    # Consent/demographics
    "consent/age":        "age",
    "consent/sex":        "sex",
    "consent/race":       "race",
    "consent/ethnicity":  "ethnicity",
    "Inferred ancestry":  "Inferred ancestry",
    # Baseline pathology
    "baseline/tumor_stage_pathological": "tumor_stage_pathological",
    "baseline/histologic_type":          "histologic_type",
    "cptac_path/histologic_grade":       "histologic_grade",
    "baseline/tumor_size_cm":            "tumor_size_cm",
    # Medical history
    "medical_history/bmi":                 "bmi",
    "medical_history/alcohol_consumption": "alcohol_consumption",
    "medical_history/tobacco_smoking_history": "tobacco_smoking_history",
    # Follow-up & survival outcomes
    "follow-up/vital_status_at_date_of_last_contact":
        "vital_status_at_date_of_last_contact",
    "follow-up/number_of_days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact":
        "number_of_days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact",
    "follow-up/number_of_days_from_date_of_initial_pathologic_diagnosis_to_date_of_death":
        "number_of_days_from_date_of_initial_pathologic_diagnosis_to_date_of_death",
    "Survival status (1, dead; 0, alive)": "Survival status (1, dead; 0, alive)",
    "Overall survival, days":              "Overall survival, days",
}

# ------------------------------------------------------------------
# 6. Select, rename, and order only the requested variables present
# ------------------------------------------------------------------
present     = {k: v for k, v in rename_map.items() if v in all_clinical.columns}
missing_out = [k for k, v in rename_map.items() if v not in all_clinical.columns]

core = (
    all_clinical
      .rename(columns={v: k for k, v in present.items()})
      .loc[:, list(present.keys())]
)

# ------------------------------------------------------------------
# 7. Save the core clinical table and report missing variables
# ------------------------------------------------------------------
core_path = out_dir / "all_cancers_clinical_core.csv"
core.to_csv(core_path, index=True)

print(f"✔ Core clinical table written to: {core_path}"
      f"  ({core.shape[0]} patients, {core.shape[1]} variables)")

if missing_out:
    print("\n⚠ The following requested variables were not found and were omitted:")
    for var in missing_out:
        print("   •", var)
