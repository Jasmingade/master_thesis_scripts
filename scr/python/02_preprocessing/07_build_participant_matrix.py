#!/usr/bin/env python
"""
build_global_peptide_run_matrix.py —  
1) Reads every per‐run Parquet file of cleaned, aggregated peptides.
2) Tags any QC/reference channels by study and run (e.g. REF_PDC000110_01CPTAC_OV…).
3) Concatenates all runs into one table and pivots to a single peptide × run matrix.
4) Writes out a Parquet file with rows=peptides and columns=TMT runs.
"""
import argparse
import pandas as pd
from pathlib import Path

# ──────────────────────────────────────────────────────────────────────────────
# Reference‐sample lookup: keys = PDC study ID, values = list of Participant IDs
# that should be renamed to the study‐level REF_<study> label.
REFERENCE_IDS = {
    "PDC000110": ["pooled sample", "QC1","QC2","QC3","QC4","QC5","QC6","QC7","QC8","QC9","Ref"],
    "PDC000116": ["ColonRef"],
    "PDC000120": ["Internal Reference - Pooled Sample", "RetroIR"],
    "PDC000125": ["Ref"],
    "PDC000127": ["pooled sample","NCI7-1","NCI7-2","NCI7-3","NCI7-4","NCI7-5","QC1","QC2","QC3","QC4","QC5","QC6","QC7","QC8","QC9"],
    "PDC000153": ["Internal Reference - Pooled Sample","Tumor Only IR","Normal Only IR","Taiwanese IR"],
    "PDC000204": ["ref"],
    "PDC000221": ["pooled sample","QC1","QC2","QC3","QC4","QC5","QC6","QC7","QC9","LungTumor1","LungTumor2","LungTumor3"],
    "PDC000234": ["LSCC Global CR","LUAD Global CR (pool 2)","JHU HNSCC CR","LUAD Global CR (pool 1)","LSCC Tumor ONLY CR"],
    "PDC000270": ["pooled sample","QC1","QC2","QC3","QC4","QC5","QC6","KoreanReference1","WU-PDA1","WU-pooled sample"],
}
# ──────────────────────────────────────────────────────────────────────────────

def main(indir: Path, outfile: Path):
    dfs = []
    n_files = 0

    for plex_dir in sorted(indir.iterdir()):
        if not plex_dir.is_dir(): 
            continue
        for study_dir in sorted(plex_dir.iterdir()):
            if not study_dir.is_dir(): 
                continue
            study = study_dir.name
            refs = REFERENCE_IDS.get(study, [])
            for pq in sorted(study_dir.glob("*.parquet")):
                df = pd.read_parquet(pq)
                n_files += 1

                # record the run name (folder_name or file stem)
                run_name = pq.stem
                df["Run"] = run_name

                # tag reference channels with a unique REF_<study>_<run> label
                if refs:
                    mask = df["Participant ID"].isin(refs)
                    if mask.any():
                        tag = f"REF_{study}_{run_name}"
                        df.loc[mask, "Participant ID"] = tag
                        print(f"  › Tagged {mask.sum():,} rows in {study}/{run_name} as {tag}")

                dfs.append(df[["PeptideSequence", "Run", "Participant ID", "total_intensity"]])

    if not dfs:
        raise RuntimeError(f"No Parquet files found under {indir}")

    print(f"Read {n_files} files → concatenating…")
    all_runs = pd.concat(dfs, ignore_index=True)

    print("Pivoting to global peptide × run matrix…")
    matrix = all_runs.pivot_table(
        index="PeptideSequence",
        columns="Run",
        values="total_intensity",
        aggfunc="sum",
        fill_value=0
    )
    matrix.columns.name = None

    print(f"Writing {outfile} ({matrix.shape[0]:,} peptides × {matrix.shape[1]:,} runs)…")
    outfile.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_parquet(outfile, index=True, compression="snappy")
    print("Done.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build a single peptides×runs matrix, tagging QC/reference channels by study & run"
    )
    parser.add_argument(
        "--indir", type=Path,
        default=Path("data/processed/aggregates_cleaned"),
        help="Top-level folder (contains TMT10/, TMT11/ subfolders)"
    )
    parser.add_argument(
        "--outfile", type=Path,
        default=Path("data/for_analysis/peptide_by_run_matrix.parquet"),
        help="Path for the combined peptide×run matrix"
    )
    args = parser.parse_args()
    main(args.indir, args.outfile)
