#!/usr/bin/env python
"""
Build peptide abundance matrices.

This script:
1. Reads each cleaned Parquet file for a TMT run.
2. Pivots data to a peptide × sample matrix (participants as columns).
3. Writes each matrix as a Parquet file under an analysis directory.
"""
import pandas as pd
from pathlib import Path

def build_run_matrix(run_parquet: Path, out_dir: Path):
    """
    Pivot a run-level Parquet file to a peptide × participant intensity matrix.
    """
    df = pd.read_parquet(run_parquet)
    matrix = df.pivot_table(
        index="PeptideSequence",
        columns="Participant ID",
        values="total_intensity",
        aggfunc="sum",
        fill_value=0
    )
    matrix.columns.name = None
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"{run_parquet.stem}_matrix.parquet"
    matrix.to_parquet(out_file, index=True, compression="snappy")
    print(f"Wrote {out_file}: {matrix.shape[0]} peptides × {matrix.shape[1]} samples")

if __name__ == "__main__":
    root     = Path("data/processed/aggregates_cleaned")
    out_root = Path("data/for_analysis/abundance_matrices")
    
    for plex_dir in sorted(root.iterdir()):
        if not plex_dir.is_dir():
            continue
        for study_dir in sorted(plex_dir.iterdir()):
            if not study_dir.is_dir():
                continue
            print(f"Processing {plex_dir.name}/{study_dir.name}")
            for pq in sorted(study_dir.glob("*.parquet")):
                build_run_matrix(
                    run_parquet=pq,
                    out_dir=out_root / plex_dir.name / study_dir.name
                )


