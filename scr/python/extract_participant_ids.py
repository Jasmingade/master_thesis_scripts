#!/usr/bin/env python
"""
Extract unique Participant IDs from abundance CSV abundance matrices,
compute how many times each appears, and list which files contain each ID.

Outputs a CSV with columns:
- Participant ID
- Count (number of CSV files containing it)
- Files (semicolon-separated list of file paths)
"""
import pandas as pd
from pathlib import Path


def extract_participant_summary(csv_root: Path) -> pd.DataFrame:
    """
    Traverse all CSV files under csv_root, extract sample column names
    (excluding peptide sequence and REF columns), and record file occurrences.

    Returns:
        DataFrame with columns ['Participant ID', 'Count', 'Files']
    """
    summary = {}
    # collect CSV files
    csv_files = sorted(csv_root.rglob("*.csv"))
    for path in csv_files:
        try:
            # read only header row
            cols = pd.read_csv(path, nrows=0).columns.tolist()
        except Exception as e:
            print(f"⚠️ Could not read {path}: {e}")
            continue

        # identify participant ID columns: skip peptide/REF
        pids = [c for c in cols if c != 'PeptideSequence' and not c.startswith('REF_')]
        for pid in pids:
            summary.setdefault(pid, []).append(str(path))

    # Build summary DataFrame
    records = []
    for pid, files in summary.items():
        records.append({
            "Participant ID": pid,
            "Count": len(files),
            "Files": ";".join(files)
        })

    df = pd.DataFrame(records)
    # If empty, warn and return empty DataFrame with correct columns
    if df.empty:
        print("⚠️  No participant IDs found in CSV files under", csv_root)
        return pd.DataFrame(columns=["Participant ID", "Count", "Files"])

    df = df.sort_values(by=["Count", "Participant ID"], ascending=[False, True])
    return df


if __name__ == "__main__":
    # <-- INSERT YOUR PATHS HERE -->
    csv_root = Path("data/for_analysis/fromPython_data_withREF_csv")  # e.g. "data/csv_matrices"
    out_csv  = Path("reports/csv_inspection/participant_summary.csv")

    # Ensure output directory exists
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    # Run extraction and save
    summary_df = extract_participant_summary(csv_root)
    summary_df.to_csv(out_csv, index=False)
    print(f"Wrote participant summary for {len(summary_df)} IDs to {out_csv}")
