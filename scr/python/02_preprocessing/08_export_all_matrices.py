#!/usr/bin/env python
import pandas as pd
from pathlib import Path

# 1) Point this at your directory of per‐run Parquet matrices
INPUT_ROOT = Path("data/for_analysis/abundance_matrices")

# 2) Where to write the CSVs (will mirror the structure under INPUT_ROOT)
OUTPUT_ROOT = Path("data/for_analysis/abundance_matrices_csv")

def main():
    # find every *_matrix.parquet under INPUT_ROOT
    for pq in sorted(INPUT_ROOT.rglob("*_matrix.parquet")):
        # compute where to write the CSV
        rel = pq.relative_to(INPUT_ROOT)              # e.g. TMT10/PDC000110/01XYZ_matrix.parquet
        out_csv = (OUTPUT_ROOT / rel).with_suffix(".csv")
        out_csv.parent.mkdir(parents=True, exist_ok=True)

        # read & write
        df = pd.read_parquet(pq)
        df.to_csv(out_csv, index=True, index_label="PeptideSequence")
        print(f"Wrote {out_csv}")

if __name__ == "__main__":
    main()
