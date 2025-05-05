#!/usr/bin/env python
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

def load_reference_matrix(study_dir: Path):
    """
    Read all *_matrix.parquet under this single study_dir,
    pull out the one REF_ column per run, and return
    a DataFrame (peptide × runs).
    """
    series_list = []
    for pq in sorted(study_dir.rglob("*_matrix.parquet")):
        df = pd.read_parquet(pq)
        # find the one REF column
        ref_cols = [c for c in df.columns if c.startswith("REF_")]
        if not ref_cols:
            continue
        ref = ref_cols[0]
        s = df[ref].copy()
        # make a nice run name: drop the "_matrix" suffix
        run_name = pq.stem.replace("_matrix", "")
        s.name = run_name
        series_list.append(s)

    if not series_list:
        # no REF columns here
        return None

    # align on peptides (union of indices)
    return pd.concat(series_list, axis=1)

# before your loop, build a little lookup of full study names:
STUDY_LABELS = {
    "PDC000110":"Prospective_Ovarian_JHU_Proteome",
    "PDC000116":"Prospective_Colon_PNNL_Proteome_Qeplus",
    "PDC000120":"Prospective_Breast_BI_Proteome",
    "PDC000125":"CPTAC_UCEC_Discovery_Study_-_Proteome",
    "PDC000127":"CPTAC_CCRCC_Discovery_Study_-_Proteome",
    "PDC000153":"CPTAC_LUAD_Discovery_Study_-_Proteome",
    "PDC000204":"CPTAC_GBM_Discovery_Study_-_Proteome",
    "PDC000221":"CPTAC_HNSCC_Discovery_Study_-_Proteome",
    "PDC000234":"CPTAC_LSCC_Discovery_Study_-_Proteome",
    "PDC000270":"CPTAC_PDA_Discovery_Study_-_Proteome",
}

# then inside the loop, after you’ve got plex_dir.name and study_dir.name:
study_id  = study_dir.name
plex      = plex_dir.name
label     = STUDY_LABELS.get(study_id, study_id)

ax.set_title(f"{plex} – {label}\nInternal‐Reference Spearman ρ", fontsize=8)


if __name__=="__main__":
    mat_root     = Path("data/for_analysis/fromPython_data_withREF")
    qc_dir       = Path("reports/qc")
    scatter_dir  = qc_dir / "scatter"
    heatmap_dir  = qc_dir / "heatmaps"

    # make sure our output folders exist
    scatter_dir.mkdir(parents=True, exist_ok=True)
    heatmap_dir.mkdir(parents=True, exist_ok=True)

    for plex_dir in sorted(mat_root.iterdir()):
        if not plex_dir.is_dir(): continue
        for study_dir in sorted(plex_dir.iterdir()):
            if not study_dir.is_dir(): continue

            ref_df = load_reference_matrix(study_dir)
            if ref_df is None or ref_df.shape[1] < 2:
                print(f"⚠️  skipping {plex_dir.name}/{study_dir.name} (need ≥2 REF runs)")
                continue

            # heatmap of Spearman correlations
            corr = ref_df.corr(method="spearman")
            fig, ax = plt.subplots(figsize=(6,6))
            im = ax.imshow(corr, vmin=0.9, vmax=1.0, cmap="viridis")
            cbar = fig.colorbar(im, ax=ax, label="Spearman ρ")
            ax.set_xticks(range(len(corr)))
            ax.set_xticklabels(corr.columns, rotation=90, fontsize=6)
            ax.set_yticks(range(len(corr)))
            ax.set_yticklabels(corr.index, fontsize=6)

            # more descriptive title
            study_id = study_dir.name
            plex     = plex_dir.name
            label    = STUDY_LABELS.get(study_id, study_id)
            ax.set_title(f"{plex} – {label}\nInternal‐Reference Spearman ρ", fontsize=8)

            fig.tight_layout()
            heatmap_path = heatmap_dir / f"{plex}_{study_id}_REF_heatmap.png"
            fig.savefig(heatmap_path)
            plt.close(fig)

            
            # scatter of first two REF runs
            xcol, ycol = ref_df.columns[:2]
            fig, ax = plt.subplots(figsize=(5,5))
            ax.scatter(ref_df[xcol], ref_df[ycol], s=1, alpha=0.3)
            ax.set_xlabel(xcol)
            ax.set_ylabel(ycol)
            ax.set_title(f"REF: {xcol} vs {ycol}")
            fig.tight_layout()
            scatter_path = scatter_dir / f"{plex_dir.name}_{study_dir.name}_REF_scatter.png"
            fig.savefig(scatter_path)
            print("WROTE scatter:", scatter_path)
            plt.close(fig)