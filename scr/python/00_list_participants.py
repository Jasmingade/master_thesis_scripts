#!/usr/bin/env python
import pandas as pd
from pathlib import Path

def list_participant_ids(agg_root: Path, out_file: Path):
    ids = set()
    for pq in sorted(agg_root.rglob("*.parquet")):
        # only read the Participant ID column for speed
        df = pd.read_parquet(pq, columns=["Participant ID"])
        ids.update(df["Participant ID"].dropna().unique())
    sorted_ids = sorted(ids)
    # write to a file
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with out_file.open("w") as f:
        for pid in sorted_ids:
            f.write(f"{pid}\n")
    # also print
    print(f"Found {len(sorted_ids)} unique Participant IDs:")
    for pid in sorted_ids:
        print(" ", pid)

if __name__=="__main__":
    agg_root = Path("data/for_analysis/fromPython_data_withREF_csv")  # or wherever per‐run matrices live
    out_txt   = Path("reports/csv_inspection/participant_ids.txt")
    list_participant_ids(agg_root, out_txt)
