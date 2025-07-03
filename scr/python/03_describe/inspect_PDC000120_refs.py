#!/usr/bin/env python
import pandas as pd
from pathlib import Path

def main():
    # 1) Where the “per‐run” matrices are stored for PDC000120:
    mat_dir = Path("data/processed/aggregates_cleaned/TMT10/PDC000120")
    if not mat_dir.is_dir():
        raise FileNotFoundError(f"{mat_dir} not found")

    # 2) Collect every “REF_…” column into a dict of Series:
    ref_series = {}
    for pq in sorted(mat_dir.glob("*_matrix.parquet")):
        df = pd.read_parquet(pq)
        ref_cols = [c for c in df.columns if c.startswith("REF_")]
        if not ref_cols:
            continue

        run_name = pq.stem.replace("_matrix", "")
        ref_series[run_name] = df[ref_cols[0]].rename(run_name)

    if not ref_series:
        print("⚠️  No REF_ columns found under", mat_dir)
        return

    # 3) Make a single DataFrame (union of all peptide indices):
    ref_df = pd.concat(ref_series, axis=1, join="outer")

    # 4) Inspect its shape and first few rows:
    print("REF‐by‐run DataFrame shape:", ref_df.shape)
    print("\nFirst 5 peptides × each REF run:")
    print(ref_df.head())

    # 5) Count non‐missing peptides per REF run:
    counts = ref_df.notna().sum(axis=0).sort_index()
    print("\nNumber of peptides measured in each REF run:")
    print(counts)

    # 6) Compute summary statistics per REF run:
    print("\nSummary (mean/std/min/max) of REF runs:")
    print(ref_df.describe().T[["count","mean","std","min","max"]])

    # 7) Save these tables for later:
    out_csv_folder = Path("reports/qc/csvs")
    out_csv_folder.mkdir(parents=True, exist_ok=True)

    ref_df.to_csv(out_csv_folder / "PDC000120_raw_REF_peptides_by_run.csv")
    counts.to_csv(out_csv_folder / "PDC000120_REF_nonmissing_counts.csv", 
                  header=["non_missing_count"])

    spearman_corr = ref_df.corr(method="spearman")
    spearman_corr.to_csv(out_csv_folder / "PDC000120_REF_spearman.csv")

    print("\n▶ Saved raw peptide×run, counts, and Spearman to reports/qc/csvs/")

if __name__ == "__main__":
    main()
