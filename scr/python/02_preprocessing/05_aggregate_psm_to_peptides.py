#!/usr/bin/env python
"""
Aggregate CPTAC TMT10/11 PSM files to peptide level.

A memory-friendly pipeline that:
- Streams large PSM files in chunks
- Tags reference samples based on predefined IDs
- Outputs aggregated results in Parquet format

Steps:
1. Parse CLI arguments and configure logging.
2. Define study directories and reference ID lists.
3. Load and normalize channel-to-sample mapping.
4. Stream each PSM file in chunks:
   a. Clean misaligned data
   b. Pivot to long format and merge metadata
   c. Tag reference samples
   d. Aggregate intensities at peptide level
5. Write per-folder Parquet outputs
6. Log summary statistics.
"""
import re
import argparse
import logging
from pathlib import Path
from collections import defaultdict
from datetime import datetime
from glob import glob

import pandas as pd
from tqdm.auto import tqdm

# ------------------------------------------------------------------
# 0. CLI arguments & logging setup
# ------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description=(
        "Aggregate CPTAC TMT10/11 PSM files to peptide level "
        "(streaming, chunked, Parquet output, tagging references)."
    )
)
parser.add_argument(
    "--chunksize", type=int, default=1_000_000,
    help="Number of rows per chunk (default: 1e6)"
)
parser.add_argument(
    "--outdir", default="data/processed/aggregates_withREF",
    help="Output directory for Parquet files"
)
args = parser.parse_args()

LOG_DIR = Path("reports") / "logs"
LOG_DIR.mkdir(parents=True, exist_ok=True)
_fmt = "%(asctime)s | %(levelname)-7s | %(message)s"
handlers = [
    logging.StreamHandler(),
    logging.FileHandler(LOG_DIR / "01_psm_run.log", mode="w"),
    logging.FileHandler(LOG_DIR / "01_psm_err.log", mode="w")
]
handlers[-1].setLevel(logging.ERROR)
for h in handlers:
    h.setFormatter(logging.Formatter(_fmt, "%Y-%m-%d %H:%M:%S"))
logging.basicConfig(level=logging.INFO, handlers=handlers)

_start_time = datetime.now()

# ------------------------------------------------------------------
# 1. Study directories & reference IDs
# ------------------------------------------------------------------
ROOT      = Path("data/raw/TMT_peptide_spectral_matches")
META_ROOT = Path("data/raw/metadata")
TMT10_STUDIES = {
    "PDC000110":"Prospective_Ovarian_JHU_Proteome",
    "PDC000116":"Prospective_Colon_PNNL_Proteome_Qeplus",
    "PDC000120":"Prospective_Breast_BI_Proteome",
    "PDC000125":"CPTAC_UCEC_Discovery_Study_-_Proteome",
    "PDC000127":"CPTAC_CCRCC_Discovery_Study_-_Proteome",
    "PDC000153":"CPTAC_LUAD_Discovery_Study_-_Proteome",
}
TMT11_STUDIES = {
    "PDC000204":"CPTAC_GBM_Discovery_Study_-_Proteome",
    "PDC000221":"CPTAC_HNSCC_Discovery_Study_-_Proteome",
    "PDC000234":"CPTAC_LSCC_Discovery_Study_-_Proteome",
    "PDC000270":"CPTAC_PDA_Discovery_Study_-_Proteome",
}
REFERENCE_IDS = {
    # PDC study ID → list of substrings marking reference samples
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

# ------------------------------------------------------------------
# 2. Channel-map loader & header fixer
# ------------------------------------------------------------------
def find_channel_map(pdc_id, plex):
    """Locate the Label_to_Sample_Mapping_File for a given study."""
    pattern = META_ROOT/plex/f"*{pdc_id}*Label_to_Sample_Mapping_File*.xlsx"
    hits = [Path(p) for p in glob(str(pattern)) if not Path(p).name.startswith("~$")]
    if not hits:
        raise FileNotFoundError(f"No channel map for {pdc_id} in {plex}")
    hits.sort(key=lambda p: ("_fixed" not in p.name, p.name))
    return hits[0]

def _fix_tmt11_headers(df):
    """Standardize TMT11 channel headers by renaming columns."""
    colmap = {}
    for c in df.columns:
        m = re.match(r"^(TMT(?:10|-11))[- ](\d{3}[NC]?)\s*(.*)$", c)
        if not m: continue
        _, chan, rest = m.groups()
        chan = "126" if chan == "126C" else chan
        colmap[c] = f"TMT11-{chan}" + (f" {rest.strip()}" if rest else "")
    return df.rename(columns=colmap)

def load_channel_map(xls, plex):
    """Read the mapping Excel and pivot to long format with sample info."""
    raw = pd.read_excel(xls)
    raw.columns = [
        re.sub(
            rf"^({plex}-\d+[NC]?)\s*(Specimen\s*Label|SpecimenLabel)",
            r"\1 Aliquot ID", c
        ) for c in raw.columns
    ]
    if plex == "TMT11":
        raw = _fix_tmt11_headers(raw)
    long = (
        raw
        .melt(
            id_vars=[c for c in raw.columns if not c.startswith(f"{plex}-")],
            var_name="raw_col", value_name="value"
        )
        .assign(
            TMT_channel=lambda d: d["raw_col"].str.extract(rf"({plex}-\S+)")[0],
            info_type=lambda d: d["raw_col"]
                .str.replace(rf"{plex}-\S+\s*", "", regex=True)
                .replace("", "Participant ID")
        )
        .pivot_table(
            index=["Folder Name", "TMT_channel"],
            columns="info_type", values="value", aggfunc="first"
        )
        .reset_index()
        .assign(
            Date=lambda d: d["Folder Name"].str.extract(r"(\d{8})")[0],
            Prefix=lambda d: d["Folder Name"].str.extract(r"^(\d{2}CPTAC)")[0]
        )
        .rename_axis(None, axis=1)
    )
    return long

def check_map_duplicates(cmap, pdc_id, plex):
    """Warn if any Folder Name × Channel combinations appear twice."""
    dups = (
        cmap.groupby(["Folder Name", "TMT_channel"])
            .size().reset_index(name="count")
            .query("count > 1")
    )
    if not dups.empty:
        logging.warning("Duplicates in %s %s", plex, pdc_id)

# ------------------------------------------------------------------
# 3. Aggregation helper
# ------------------------------------------------------------------
CHANNELS = {
    "TMT10": ["126","127N","127C","128N","128C","129N","129C","130N","130C","131"],
    "TMT11": ["126","127N","127C","128N","128C","129N","129C","130N","130C","131N","131C"],
}
def channel_regex(plex):
    """Build regex matching TMT channel columns for a given plex."""
    return rf"^{plex}-({'|'.join(CHANNELS[plex])})$"

def aggregate_psm_stream(psm, cmap, plex, pdc_id, chunksize):
    """
    Stream a PSM file in chunks, tag references, and aggregate to peptides.

    Returns a DataFrame grouped by Participant ID and PeptideSequence.
    """
    folder = psm.parent.name
    chan_re = channel_regex(plex)

    # 1) Determine header row
    with psm.open("r") as f:
        header_idx = next(
            (i for i, line in enumerate(f) if line.strip() and not line.startswith("#")),
            None
        )
        if header_idx is None:
            raise RuntimeError(f"No header in {psm}")

    cols = pd.read_table(psm, sep="\t", nrows=0).columns.tolist()
    reader_kwargs = dict(
        sep="\t", header=header_idx, names=cols,
        dtype=str, keep_default_na=False, chunksize=chunksize
    )

    frames = []
    for chunk in pd.read_table(psm, **reader_kwargs):
        # Repair misaligned PDC116/PDC120 sequences
        if chunk["PeptideSequence"].isin({"0","1"}).all():
            rt = chunk["RTAtPrecursorHalfElution"]
            num_rt = rt.str.extract(r"^([+-]?\d+(?:\.\d+)?)")[0]
            chunk["RTAtPrecursorHalfElution"] = pd.to_numeric(num_rt, errors="coerce")
            chunk["PeptideSequence"] = rt

        # Melt to long format and merge channel map
        long = (
            chunk
            .melt(
                id_vars=[c for c in chunk.columns if not re.match(chan_re, c)],
                var_name="TMT_channel", value_name="intensity_raw"
            )
            .loc[lambda d: d["TMT_channel"].str.match(chan_re)]
            .assign(
                **{"Folder Name": folder},
                intensity=lambda d: pd.to_numeric(
                    d["intensity_raw"].str.extract(r"^(\d+(?:\.\d+)?)")[0],
                    errors="coerce"
                )
            )
            .merge(cmap, on=["Folder Name", "TMT_channel"], how="left")
        )

        logging.info(
            f"[{pdc_id}/{folder}] Found Participant IDs: "
            f"{long['Participant ID'].dropna().unique().tolist()}"
        )

        # Tag reference samples by substring match
        ref_ids = REFERENCE_IDS.get(pdc_id, [])
        if ref_ids:
            mask = False
            for rid in ref_ids:
                mask |= long["Participant ID"].str.contains(rid, case=False, na=False)
            long.loc[mask, "Participant ID"] = f"REF_{folder}"
            logging.info(f"Tagged {mask.sum()} rows as REF_{folder}")

        # Aggregate per chunk
        agg = (
            long
            .groupby(["Participant ID", "PeptideSequence"], as_index=False)
            .agg(
                total_intensity=("intensity", "sum"),
                mean_intensity=("intensity", "mean"),
                n_psm=("intensity", "size")
            )
        )
        frames.append(agg)

    # Final aggregation across all chunks
    df = pd.concat(frames, ignore_index=True)
    return (
        df
        .groupby(["Participant ID", "PeptideSequence"], as_index=False)
        .agg(
            total_intensity=("total_intensity", "sum"),
            mean_intensity=("mean_intensity", "mean"),
            n_psm=("n_psm", "sum")
        )
    )

# ------------------------------------------------------------------
# 4. Main processing loop
# ------------------------------------------------------------------
out_root = Path(args.outdir)
out_root.mkdir(exist_ok=True)

total_files = total_rows = 0
for plex, studies in (("TMT10", TMT10_STUDIES), ("TMT11", TMT11_STUDIES)):
    logging.info("▶ Processing %s (%d studies)", plex, len(studies))
    for pdc_id, tail in tqdm(studies.items(), desc=f"{plex} studies"):
        study_root = ROOT / f"{plex}_files/{pdc_id}_{tail}"
        cmap_file  = find_channel_map(pdc_id, plex)
        cmap       = load_channel_map(cmap_file, plex)
        check_map_duplicates(cmap, pdc_id, plex)

        out_dir = out_root / pdc_id
        out_dir.mkdir(exist_ok=True)
        done = {p.stem for p in out_dir.glob("*.parquet")}

        psm_files = sorted(study_root.rglob("*.psm"))
        if not psm_files:
            logging.warning("No `.psm` files found in %s", study_root)
            continue

        buckets = defaultdict(list)
        for psm in tqdm(psm_files, desc=pdc_id, leave=False):
            folder = psm.parent.name
            if folder in done:
                continue
            agg = aggregate_psm_stream(psm, cmap, plex, pdc_id, args.chunksize)
            buckets[folder].append(agg)
            total_rows  += agg.shape[0]
            total_files += 1

        for folder, dfs in buckets.items():
            out_path = out_dir / f"{folder}.parquet"
            (
                pd.concat(dfs, ignore_index=True)
                .groupby(["Participant ID", "PeptideSequence"], as_index=False)
                .agg(
                    total_intensity=("total_intensity", "sum"),
                    mean_intensity=("mean_intensity", "mean"),
                    n_psm=("n_psm", "sum")
                )
                .to_parquet(out_path, index=False, compression="snappy")
            )
        logging.info("✓ %s — wrote %d runs → %s", pdc_id, len(buckets), out_dir)

elapsed = datetime.now() - _start_time
logging.info("🏁 Done in %s | %d files → %d rows", elapsed, total_files, total_rows)
