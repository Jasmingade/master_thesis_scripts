#!/usr/bin/env python
import pandas as pd
from pathlib import Path

# Change this to match the exact path on disk where your PDC000120 matrices are
study_folder = Path("data/for_analysis/abundance_matrices/TMT11/PDC000270")

if not study_folder.exists():
    raise RuntimeError(f"Folder not found: {study_folder}")

all_channels = set()

for pq_file in sorted(study_folder.glob("*_matrix.parquet")):
    # Read the Parquet file, which should be a peptide×channel matrix
    df = pd.read_parquet(pq_file)
    # Collect all column names (each column is a Participant ID / channel)
    cols = list(df.columns)
    print(f"\n— {pq_file.name} has channels:")
    for c in cols:
        print("    ", c)
    all_channels.update(cols)

# If you just want the unique set of all channels across all runs:
print("\n\n► Unique channels seen in PDC000120 (all runs combined):")
for c in sorted(all_channels):
    print("   ", c)
