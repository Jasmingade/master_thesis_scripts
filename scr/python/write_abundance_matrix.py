#!/usr/bin/env python
import pandas as pd
from pathlib import Path

# 1) point to one of your run-level matrices
pq = Path("data/for_analysis/fromPython_data_withREF/TMT10/PDC000110/01CPTAC_OVprospective_Proteome_JHU_20161209_matrix.parquet")

# 2) read it in
df = pd.read_parquet(pq)

# 3) move the index (PeptideSequence) into a real column
df = df.reset_index()  # now you have a "PeptideSequence" column again

# 4) (optional) take a small preview if you just want the first 50 rows
preview = df.head(50)

# 5) make sure the output folder exists
out = Path("reports/csv_inspection")
out.mkdir(parents=True, exist_ok=True)

# 6) write the full matrix (or the preview) to CSV
preview.to_csv(out / f"{pq.stem}.preview.csv", index=False)

print("Wrote:", out / f"{pq.stem}.preview.csv")
