#!/usr/bin/env python
import pandas as pd
from pathlib import Path

# 1) Point to your Parquet file:
pq = Path("data/for_analysis/fromPython_data/TMT10/PDC000110/01CPTAC_OVprospective_Proteome_JHU_20161209_matrix.parquet")

# 2) Read it
df = pd.read_parquet(pq)

# 3) Take just the first 50 rows to inspect
preview = df.head(50)

# 4) Ensure output directory exists
out_dir = Path("reports/csv_inspection")
out_dir.mkdir(parents=True, exist_ok=True)

# 5) Build your CSV path by combining the stem with a new suffix
out_csv = out_dir / f"{pq.stem}.preview.csv"

# 6) Write out
preview.to_csv(out_csv, index=True, index_label="PeptideSequence")
print(f"Wrote preview to {out_csv}")
