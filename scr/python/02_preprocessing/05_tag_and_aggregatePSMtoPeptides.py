#!/usr/bin/env python
"""
01_process_tmt_psm_referenceSamples.py  —  memory-friendly CPTAC TMT10/11 pipeline
        (streaming, chunked, Parquet output, tagging references by folder)
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

# ──────────────────────────────────────────────────────────────────────────────
# 0. CLI arguments & logging setup
# ──────────────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(
    description=(
        "Aggregate CPTAC TMT10/11 PSM files to peptide level "
        "(streaming, chunked, Parquet output, tagging references by folder)."
    )
)
parser.add_argument("--chunksize", type=int, default=1_000_000,
                    help="Rows per chunk (default 1e6)")
parser.add_argument("--outdir", default="data/processed/aggregates",
                    help="Folder for Parquet outputs")
args = parser.parse_args()

LOG_DIR = Path("reports") / "logs"
LOG_DIR.mkdir(parents=True, exist_ok=True)
fmt = "%(asctime)s | %(levelname)-7s | %(message)s"
datefmt = "%Y-%m-%d %H:%M:%S"
handlers = [
    logging.StreamHandler(),
    logging.FileHandler(LOG_DIR / "05_psm_run.log", mode="w"),
    logging.FileHandler(LOG_DIR / "05_psm_err.log", mode="w")
]
handlers[-1].setLevel(logging.ERROR)
for h in handlers:
    h.setFormatter(logging.Formatter(fmt, datefmt))
logging.basicConfig(level=logging.INFO, handlers=handlers)

_start_time = datetime.now()

# ──────────────────────────────────────────────────────────────────────────────
# 1. Study directories & reference IDs
# ──────────────────────────────────────────────────────────────────────────────
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

# Keys = PDC study ID; values = list of Participant IDs marking references
REFERENCE_IDS = {
    "PDC000110": ["pooled sample","QC1","QC2","QC3","QC4","QC5","QC6","QC7","QC8","QC9","Ref"],
    "PDC000116": ["ColonRef"],
    "PDC000120": ["Internal Reference - Pooled Sample","RetroIR"],
    "PDC000125": ["Ref", "NX1", "NX2", "NX3", "NX4", "NX5", "NX6", "NX7", "NX8", "NX9"],
    "PDC000127": ["pooled sample","NCI7-1","NCI7-2","NCI7-3","NCI7-4","NCI7-5","QC1","QC2","QC3","QC4","QC5","QC6","QC7","QC8","QC9"],
    "PDC000153": ["Internal Reference - Pooled Sample","Tumor Only IR","Normal Only IR","Taiwanese IR"],
    "PDC000204": ["ref", "GTEX-NPJ7-0011-R10A-SM-HAKXW", "GTEX-P44H-0011-R10A-SM-HAKXX", "GTEX-Q2AG-0011-R10A-SM-HAKXT", "GTEX-QVJO-0011-R10A-SM-HAKXV", 
                  "GTEX-R55F-0011-R10A-SM-HAKXY", "GTEX-RN5K-0011-R10A-SM-HAKXU", "GTEX-RU72-0011-R10A-SM-HAKXS", "GTEX-UTHO-0011-R10A-SM-HAKY2", 
                  "GTEX-WVLH-0011-R10A-SM-HAKXZ", "GTEX-Y8DK-0011-R10A-SM-HAKY1"],
    "PDC000221": ["pooled sample","QC1","QC2","QC3","QC4","QC5","QC6","QC7","QC9","LungTumor1","LungTumor2","LungTumor3"],
    "PDC000234": ["LSCC Global CR","LUAD Global CR (pool 2)","JHU HNSCC CR","LUAD Global CR (pool 1)","LSCC Tumor ONLY CR"],
    "PDC000270": ["pooled sample","QC1","QC2","QC3","QC4","QC5","QC6","KoreanReference1","WU-PDA1","WU-pooled sample", "KoreanReference2", "KoreanReference3"],
}



# ──────────────────────────────────────────────────────────────────────────────
# 2. Channel-map loader & header fixer
# ──────────────────────────────────────────────────────────────────────────────
def find_channel_map(pdc_id, plex):
    pattern = META_ROOT/plex/f"*{pdc_id}*Label_to_Sample_Mapping_File*.xlsx"
    hits = [Path(p) for p in glob(str(pattern)) if not Path(p).name.startswith("~$")]
    if not hits:
        raise FileNotFoundError(f"No channel map for {pdc_id} in {plex}")
    hits.sort(key=lambda p: ("_fixed" not in p.name, p.name))
    return hits[0]

def _fix_tmt11_headers(df):
    colmap = {}
    for c in df.columns:
        m = re.match(r"^(TMT(?:10|-11))[- ](\d{3}[NC]?)\s*(.*)$", c)
        if not m: continue
        _, chan, rest = m.groups()
        chan = "126" if chan=="126C" else chan
        colmap[c] = f"TMT11-{chan}" + (f" {rest.strip()}" if rest else "")
    return df.rename(columns=colmap)

def load_channel_map(xls, plex):
    raw = pd.read_excel(xls)
    raw.columns = [
        re.sub(
            rf"^({plex}-\d+[NC]?)\s*(Specimen\s*Label|SpecimenLabel)",
            r"\1 Aliquot ID", c
        ) for c in raw.columns
    ]
    if plex=="TMT11":
        raw = _fix_tmt11_headers(raw)
    long = (
        raw.melt(
            id_vars=[c for c in raw.columns if not c.startswith(f"{plex}-")],
            var_name="raw_col", value_name="value"
        )
        .assign(
            TMT_channel=lambda d: d["raw_col"].str.extract(rf"({plex}-\S+)")[0],
            info_type=lambda d: (
                d["raw_col"]
                 .str.replace(rf"{plex}-\S+\s*", "", regex=True)
                 .replace("", "Participant ID")
            )
        )
        .pivot_table(
            index=["Folder Name","TMT_channel"],
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
    dups = (
        cmap.groupby(["Folder Name","TMT_channel"])
            .size().reset_index(name="count")
            .query("count>1")
    )
    if not dups.empty:
        logging.warning("Duplicates in %s %s", plex, pdc_id)

# ──────────────────────────────────────────────────────────────────────────────
# 3. Aggregation helper
# ──────────────────────────────────────────────────────────────────────────────
CHANNELS = {
    "TMT10":["126","127N","127C","128N","128C","129N","129C","130N","130C","131"],
    "TMT11":["126","127N","127C","128N","128C","129N","129C","130N","130C","131N","131C"],
}
def channel_regex(plex):
    return rf"^{plex}-({'|'.join(CHANNELS[plex])})$"

def aggregate_psm_stream(psm, cmap, plex, chunksize):
    folder = psm.parent.name
    chan_re = channel_regex(plex)

    # infer header row
    with psm.open("r") as f:
        header_idx = None
        for i,line in enumerate(f):
            if line.strip() and not line.startswith("#"):
                header_idx = i
                break
        if header_idx is None:
            raise RuntimeError(f"No header in {psm}")

    cols = pd.read_table(psm, sep="\t", nrows=0).columns.tolist()
    reader_kwargs = dict(
        sep="\t", header=header_idx, names=cols,
        dtype=str, keep_default_na=False, chunksize=chunksize
    )

    frames = []
    for chunk in pd.read_table(psm, **reader_kwargs):
        # repair mis‐alignment if needed
        if chunk["PeptideSequence"].isin({"0","1"}).all():
            rt = chunk["RTAtPrecursorHalfElution"]
            num = rt.str.extract(r"^([+-]?\d+(?:\.\d+)?)")[0]
            chunk["RTAtPrecursorHalfElution"] = pd.to_numeric(num, errors="coerce")
            chunk["PeptideSequence"] = rt

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
            .merge(cmap, on=["Folder Name","TMT_channel"], how="left")
        )

        # aggregate per chunk
        agg = (
            long
            .groupby(["Participant ID","PeptideSequence"], as_index=False)
            .agg(
                total_intensity=("intensity","sum"),
                mean_intensity=("intensity","mean"),
                n_psm=("intensity","size")
            )
        )
        frames.append(agg)

    df = pd.concat(frames, ignore_index=True)
    return (
        df.groupby(["Participant ID","PeptideSequence"], as_index=False)
          .agg(
              total_intensity=("total_intensity","sum"),
              mean_intensity=("mean_intensity","mean"),
              n_psm=("n_psm","sum")
          )
    )

# ──────────────────────────────────────────────────────────────────────────────
# 4. Main processing loop
# ──────────────────────────────────────────────────────────────────────────────
out_root = Path(args.outdir)
out_root.mkdir(exist_ok=True)

total_files = total_rows = 0
for plex, studies in (("TMT10",TMT10_STUDIES),("TMT11",TMT11_STUDIES)):
    logging.info("▶ %s — %d studies", plex, len(studies))
    for pdc_id, tail in tqdm(studies.items(), desc=f"{plex} studies"):
        study_root = ROOT / f"{plex}_files/{pdc_id}_{tail}"

        # 1) load map
        xls  = find_channel_map(pdc_id, plex)
        cmap = load_channel_map(xls, plex)

        # 2) tag references in the map by folder_name
        refs = REFERENCE_IDS.get(pdc_id, [])
        if refs:
            m = cmap["Participant ID"].isin(refs)
            if m.any():
                cmap.loc[m, "Participant ID"] = (
                    cmap.loc[m, "Folder Name"].apply(lambda f: f"REF_{f}")
                )
                logging.info(
                    f"[{pdc_id}] Tagged {m.sum()} channel-map entries as REF_<folder>"
                )

        check_map_duplicates(cmap, pdc_id, plex)

        out_dir = out_root / pdc_id
        out_dir.mkdir(exist_ok=True)
        done = {p.stem for p in out_dir.glob("*.parquet")}

        psm_files = sorted(study_root.rglob("*.psm"))
        if not psm_files:
            logging.warning("No .psm in %s", study_root)
            continue

        buckets = defaultdict(list)
        for psm in tqdm(psm_files, desc=pdc_id, leave=False):
            folder = psm.parent.name
            if folder in done:
                continue
            agg = aggregate_psm_stream(psm, cmap, plex, args.chunksize)
            buckets[folder].append(agg)
            total_rows  += agg.shape[0]
            total_files += 1

        for folder, dfs in buckets.items():
            outp = out_dir / f"{folder}.parquet"
            pd.concat(dfs, ignore_index=True) \
              .groupby(["Participant ID","PeptideSequence"], as_index=False) \
              .agg(
                  total_intensity=("total_intensity","sum"),
                  mean_intensity=("mean_intensity","mean"),
                  n_psm=("n_psm","sum")
              ) \
              .to_parquet(outp, index=False, compression="snappy")
        logging.info("✓ %s — %d folders → %s", pdc_id, len(buckets), out_dir)

elapsed = datetime.now() - _start_time
logging.info("🏁 Done in %s | %d files → %d rows",
             elapsed, total_files, total_rows)
