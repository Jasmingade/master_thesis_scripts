#!/usr/bin/env python
"""
Summarize CPTAC TMT10/11 proteomics cohorts.

Generates a CSV overview of:
- Number of aggregated cases per study
- How many matched to clinical metadata
- How many valid for survival analysis
- Raw PSM file counts
Includes a footer row summarizing totals.

Steps:
1. Define study dictionaries for TMT10 and TMT11.
2. Parse CLI arguments for clinical CSV, raw PSM root, and output path.
3. Load and filter clinical metadata.
4. Summarize each study:
   a. Count unique Participant IDs from parquet outputs
   b. Match to clinical and count valid survival entries
   c. Count raw .psm files
5. Append a totals footer and write the summary CSV.
6. Print the summary to console.
"""
import argparse
from pathlib import Path

import pandas as pd

# ------------------------------------------------------------------
# 1. Study mappings
# ------------------------------------------------------------------
TMT10_STUDIES = {
    "PDC000110": "Prospective_Ovarian_JHU_Proteome",
    "PDC000116": "Prospective_Colon_PNNL_Proteome_Qeplus",
    "PDC000120": "Prospective_Breast_BI_Proteome",
    "PDC000125": "CPTAC_UCEC_Discovery_Study_-_Proteome",
    "PDC000127": "CPTAC_CCRCC_Discovery_Study_-_Proteome",
    "PDC000153": "CPTAC_LUAD_Discovery_Study_-_Proteome",
}
TMT11_STUDIES = {
    "PDC000204": "CPTAC_GBM_Discovery_Study_-_Proteome",
    "PDC000221": "CPTAC_HNSCC_Discovery_Study_-_Proteome",
    "PDC000234": "CPTAC_LSCC_Discovery_Study_-_Proteome",
    "PDC000270": "CPTAC_PDA_Discovery_Study_-_Proteome",
}

# ------------------------------------------------------------------
# 2. Summarization helper
# ------------------------------------------------------------------
def summarize(studies: dict, clinical: pd.DataFrame, raw_root: Path, plex_label: str) -> pd.DataFrame:
    """
    Build a DataFrame summarizing each study:
    - Total aggregated cases
    - Matched to clinical
    - Valid for survival analysis
    - Count of raw .psm files
    """
    records = []
    for pdc_id, tail in studies.items():
        agg_dir = Path("data/processed/aggregates_withREF") / pdc_id
        ids = set()
        if agg_dir.exists():
            for pq in agg_dir.rglob("*.parquet"):
                ids.update(pd.read_parquet(pq, columns=["Participant ID"])["Participant ID"].unique())
        total_agg = len(ids)

        matched = clinical[clinical["Patient_ID"].isin(ids)]
        total_matched = matched.shape[0]

        valid = matched.loc[
            matched["follow-up/vital_status_at_date_of_last_contact"].notna() &
            (matched["follow-up/vital_status_at_date_of_last_contact"] != "") &
            matched["follow-up/number_of_days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact"].notna() &
            (matched["follow-up/number_of_days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact"] != "")
        ]
        total_valid = valid.shape[0]

        raw_dir = raw_root / f"{plex_label}_files" / f"{pdc_id}_{tail}"
        n_psm = sum(1 for _ in raw_dir.rglob("*.psm")) if raw_dir.exists() else 0

        tumor = matched["tumor_code"].mode()
        tumor = tumor.iloc[0] if not tumor.empty else ""

        records.append({
            "Study ID": pdc_id,
            "Cancer Type": tumor,
            "Total Cases": total_agg,
            "Matched to Clinical": total_matched,
            "Valid for Survival": total_valid,
            "Amount of PSM Files": n_psm,
            "Plex": plex_label,
        })
    return pd.DataFrame(records)

# ------------------------------------------------------------------
# 3. Main entrypoint
# ------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarize CPTAC cohorts for TMT10 and TMT11")
    parser.add_argument(
        "--clinical_csv",
        default="data/raw/clinical_tables/all_cancers_clinical_core.csv",
        help="Path to clinical metadata CSV"
    )
    parser.add_argument(
        "--raw_root",
        default="data/raw/TMT_peptide_spectral_matches",
        help="Root directory for raw PSM files"
    )
    parser.add_argument(
        "--out_csv",
        default="reports/tables/01_cohort_summary.csv",
        help="Output path for the cohort summary CSV"
    )
    args = parser.parse_args()

    # Load and filter clinical data
    clinical = pd.read_csv(args.clinical_csv, dtype=str)
    clinical = clinical[[
        "Patient_ID", "tumor_code",
        "follow-up/vital_status_at_date_of_last_contact",
        "follow-up/number_of_days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact"
    ]]

    raw_root = Path(args.raw_root)

    # Generate summaries for both plexes
    df10 = summarize(TMT10_STUDIES, clinical, raw_root, "TMT10")
    df11 = summarize(TMT11_STUDIES, clinical, raw_root, "TMT11")
    summary = pd.concat([df10, df11], ignore_index=True)

    # Add footer row of totals
    footer = {
        "Study ID": "Sum",
        "Cancer Type": summary.shape[0],
        "Total Cases": summary["Total Cases"].sum(),
        "Matched to Clinical": summary["Matched to Clinical"].sum(),
        "Valid for Survival": summary["Valid for Survival"].sum(),
        "Amount of PSM Files": summary["Amount of PSM Files"].sum(),
        "Plex": ""
    }
    summary = pd.concat([summary, pd.DataFrame([footer])], ignore_index=True)

    # Export and display
    summary.to_csv(args.out_csv, index=False)
    print(summary.to_string(index=False))