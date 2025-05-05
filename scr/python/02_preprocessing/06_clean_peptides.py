#!/usr/bin/env python
"""
Clean aggregated peptide tables.

This script performs the following steps:
1. Parse CLI arguments for input/output directories.
2. Configure logging to console and file.
3. Iterate each study folder in the input directory:
   a. Skip files already cleaned.
   b. Load Parquet, strip “+mass” tags from peptide sequences.
   c. Filter to sequences containing only the 20 standard amino acids.
   d. Write cleaned Parquet files to the output directory.
4. Report per‐file and overall statistics via the logger.
"""
import re
import argparse
import logging
from pathlib import Path
from datetime import datetime
import pandas as pd
from tqdm.auto import tqdm

# ────────────────────────────────────────────────────────────────
# 0 ▸ CLI arguments & logging setup
# ────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(
    description=(
        "Clean aggregated Parquet peptide tables by "
        "stripping +mass tags and filtering to standard amino acids"
    )
)
parser.add_argument(
    "--indir",
    default="data/processed/aggregates_withREF",
    help="Directory of raw Parquet outputs"
)
parser.add_argument(
    "--outdir",
    default="aggregates_cleaned_withREF",
    help="Directory for cleaned Parquet outputs"
)
args = parser.parse_args()

LOG_DIR = Path("reports") / "logs"
LOG_DIR.mkdir(parents=True, exist_ok=True)
log_file = LOG_DIR / "02_clean_peptides.log"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-7s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(log_file, mode="w")
    ]
)

_start_time = datetime.now()

# ────────────────────────────────────────────────────────────────
# 1 ▸ Constants
# ────────────────────────────────────────────────────────────────
MOD_TAG_REGEX   = re.compile(r"\+\d+(?:\.\d+)?")  # Matches +123.456-style tags
VALID_SEQ_REGEX = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")  # 20 aa only

# ────────────────────────────────────────────────────────────────
# 2 ▸ Main cleaning loop
# ────────────────────────────────────────────────────────────────
indir = Path(args.indir)
outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

total_files = total_before = total_after = total_skipped = 0

for study_dir in sorted(indir.iterdir()):
    if not study_dir.is_dir():
        continue
    target_dir = outdir / study_dir.name
    target_dir.mkdir(parents=True, exist_ok=True)

    for pq in tqdm(sorted(study_dir.glob("*.parquet")),
                   desc=f"Cleaning {study_dir.name}", unit="file"):
        total_files += 1
        out_path = target_dir / pq.name

        # Skip if already cleaned
        if out_path.exists():
            logging.info("↷ Skipping %s (already cleaned)", pq.name)
            total_skipped += 1
            continue

        # Load aggregated peptides
        df = pd.read_parquet(pq)
        before = len(df)

        # Remove modification tags and filter sequences
        df["PeptideSequence_clean"] = (
            df["PeptideSequence"]
              .astype(str)
              .str.replace(MOD_TAG_REGEX, "", regex=True)
        )
        mask = df["PeptideSequence_clean"].str.match(VALID_SEQ_REGEX)
        df_clean = df.loc[mask].copy()
        after = len(df_clean)

        # Overwrite original sequence column and drop helper
        df_clean["PeptideSequence"] = df_clean["PeptideSequence_clean"]
        df_clean.drop(columns="PeptideSequence_clean", inplace=True)

        total_before += before
        total_after  += after

        # Save cleaned table
        df_clean.to_parquet(out_path, index=False, compression="snappy")
        logging.info(
            "%s: %d → %d sequences (dropped %d)",
            pq.name, before, after, before - after
        )

# ────────────────────────────────────────────────────────────────
# 3 ▸ Summary
# ────────────────────────────────────────────────────────────────
elapsed = datetime.now() - _start_time
logging.info(
    "🏁 Completed in %s | Files: %d processed (%d skipped) | Sequences: %d → %d",
    elapsed, total_files, total_skipped, total_before, total_after
)
