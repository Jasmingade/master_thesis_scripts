#!/usr/bin/env python
"""
00_count_and_list_ref_ids.py — count and list distinct REF_… Participant IDs per study
"""

import pandas as pd
from pathlib import Path

def find_ref_ids_in_study(study_dir: Path):
    """
    Scan all .parquet files in study_dir, collect any Participant IDs
    that start with "REF_", and return a sorted list of unique such IDs.
    """
    refs = set()
    for pq in study_dir.glob("*.parquet"):
        try:
            df = pd.read_parquet(pq, columns=["Participant ID"])
        except Exception:
            # if parquet has no such column or is unreadable, skip it
            continue
        for pid in df["Participant ID"].dropna().unique():
            if isinstance(pid, str) and pid.startswith("REF_"):
                refs.add(pid)
    return sorted(refs)

if __name__=="__main__":
    base = Path("data/processed/aggregates_cleaned_withREF")  # adjust to your root
    for plex_dir in sorted(base.iterdir()):
        if not plex_dir.is_dir(): 
            continue
        print(f"\n=== {plex_dir.name} ===")
        for study_dir in sorted(plex_dir.iterdir()):
            if not study_dir.is_dir(): 
                continue
            ref_ids = find_ref_ids_in_study(study_dir)
            print(f"{study_dir.name}: {len(ref_ids)} distinct REF IDs")
            if ref_ids:
                for rid in ref_ids:
                    print("   ", rid)
