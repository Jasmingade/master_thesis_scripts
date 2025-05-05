#!/usr/bin/env python
# ────────────────────────────────────────────────────────────────
#  process_tmt_psm.py  —  memory-friendly CPTAC TMT10/11 pipeline
#          (streaming, chunked, Parquet output)
# ────────────────────────────────────────────────────────────────

import re
import argparse
import logging
from pathlib import Path
from collections import defaultdict
from datetime import datetime
from glob import glob

import pandas as pd
from tqdm.auto import tqdm

# ────────────────────────────────────────────────────────────────
#  0 ▸ CLI & logging
# ────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(
    description="Aggregate CPTAC TMT10/11 PSM files to peptide level "
                "(streaming, chunked, Parquet output)."
)
parser.add_argument("--chunksize", type=int, default=1_000_000,
                    help="Rows per chunk (default 1e6)")
parser.add_argument("--outdir", default="aggregates",
                    help="Folder for Parquet outputs")
args = parser.parse_args()

LOG_DIR = Path("reports") / "logs"
LOG_DIR.mkdir(parents=True, exist_ok=True)

fmt = "%(asctime)s | %(levelname)-7s | %(message)s"
datefmt = "%Y-%m-%d %H:%M:%S"
handlers = [
    logging.StreamHandler(),
    logging.FileHandler(LOG_DIR / "01_psm_run.log", mode="w"),
    logging.FileHandler(LOG_DIR / "01_psm_err.log", mode="w")
]
handlers[-1].setLevel(logging.ERROR)
for h in handlers:
    h.setFormatter(logging.Formatter(fmt, datefmt))
logging.basicConfig(level=logging.INFO, handlers=handlers)

_t0 = datetime.now()

# ────────────────────────────────────────────────────────────────
#  1 ▸ Paths & studies
# ────────────────────────────────────────────────────────────────
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

# ────────────────────────────────────────────────────────────────
#  2 ▸ Channel-map loader
# ────────────────────────────────────────────────────────────────
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

# ────────────────────────────────────────────────────────────────
#  3 ▸ helpers & constants
# ────────────────────────────────────────────────────────────────
CHANNELS = {
    "TMT10":["126","127N","127C","128N","128C","129N","129C","130N","130C","131"],
    "TMT11":["126","127N","127C","128N","128C","129N","129C","130N","130C","131N","131C"],
}
def channel_regex(plex):
    return rf"^{plex}-({'|'.join(CHANNELS[plex])})$"

# ────────────────────────────────────────────────────────────────
#  4 ▸ aggregate_psm_stream (find header as first non-comment)
# ────────────────────────────────────────────────────────────────
def aggregate_psm_stream(psm, cmap, plex, chunksize):
    folder = psm.parent.name
    chan_re = channel_regex(plex)

    # 1) Read header line (first non-empty, non-comment)
    with psm.open("r") as f:
        header_idx = None
        for i, line in enumerate(f):
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            header_line = line
            header_idx = i
            break
        if header_idx is None:
            raise RuntimeError(f"No header found in {psm}")

    cols = header_line.split("\t")

    reader_kwargs = dict(
        sep="\t",
        header=None,
        names=cols,
        skiprows=header_idx+1,
        dtype=str,
        keep_default_na=False,
        chunksize=chunksize,
    )

    frames = []
    for chunk in pd.read_table(psm, **reader_kwargs):
        # ———————— repair mis-aligned PDC116/PDC120 ————————
        if chunk["PeptideSequence"].isin({"0","1"}).all():
            rt_raw = chunk["RTAtPrecursorHalfElution"]
            # extract numeric RT
            numeric_rt = rt_raw.str.extract(r"^([+-]?\d+(?:\.\d+)?)")[0]
            chunk["RTAtPrecursorHalfElution"] = pd.to_numeric(numeric_rt, errors="coerce")
            # move full original (mass+sequence) string into PeptideSequence
            chunk["PeptideSequence"] = rt_raw

        logging.info("[%s] %s – sample PeptideSequence: %s",
                     plex, folder,
                     chunk["PeptideSequence"].dropna().unique()[:5].tolist())

        long = (
            chunk
            .melt(
                id_vars=[c for c in chunk.columns if not re.match(chan_re, c)],
                var_name="TMT_channel",
                value_name="intensity_raw"
            )
            .loc[lambda d: d["TMT_channel"].str.match(chan_re)]
            .assign(**{"Folder Name": folder},
                    intensity=lambda d: pd.to_numeric(
                        d["intensity_raw"].str.extract(r"^(\d+(?:\.\d+)?)")[0],
                        errors="coerce"
                    ))
            .merge(cmap, on=["Folder Name","TMT_channel"], how="left")
            .groupby(["Participant ID","PeptideSequence"], as_index=False)
            .agg(
                total_intensity=("intensity","sum"),
                mean_intensity=("intensity","mean"),
                n_psm=("intensity","size")
            )
        )
        frames.append(long)

    df = pd.concat(frames, ignore_index=True)
    return (
        df.groupby(["Participant ID","PeptideSequence"], as_index=False)
          .agg(
              total_intensity=("total_intensity","sum"),
              mean_intensity=("mean_intensity","mean"),
              n_psm=("n_psm","sum")
          )
    )

# ────────────────────────────────────────────────────────────────
#  5 ▸ main loop + resume
# ────────────────────────────────────────────────────────────────
out_root = Path(args.outdir)
out_root.mkdir(exist_ok=True)

tot_files = tot_rows = 0
for plex, studies in (("TMT10", TMT10_STUDIES), ("TMT11", TMT11_STUDIES)):
    logging.info("▶ %s — %d studies", plex, len(studies))
    for pdc_id, tail in tqdm(studies.items(), desc=f"{plex} studies"):
        study_root = ROOT / f"{plex}_files" / f"{pdc_id}_{tail}"

        xls  = find_channel_map(pdc_id, plex)
        cmap = load_channel_map(xls, plex)
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
            tot_rows += agg.shape[0]
            tot_files += 1

        for folder, dfs in buckets.items():
            outp = out_dir / f"{folder}.parquet"
            pd.concat(dfs, ignore_index=True) \
              .groupby(["Participant ID","PeptideSequence"], as_index=False) \
              .agg(total_intensity=("total_intensity","sum"),
                   mean_intensity=("mean_intensity","mean"),
                   n_psm=("n_psm","sum")) \
              .to_parquet(outp, index=False, compression="snappy")

        logging.info("✓ %s — %d folders → %s", pdc_id, len(buckets), out_dir)

# ────────────────────────────────────────────────────────────────
#  6 ▸ summary
# ────────────────────────────────────────────────────────────────
dt = datetime.now() - _t0
logging.info("🏁 Done in %s | %d files → %d rows",
             dt, tot_files, tot_rows)
logging.info("Parquet tables in %s", out_root.resolve())
