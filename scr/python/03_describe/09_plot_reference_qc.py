#!/usr/bin/env python
"""
09_plot_reference_qc.py —  
Comprehensive QC for proteomics internal-reference runs:
- Save Spearman correlation matrices as CSV for downstream inspection
- Clustered Spearman correlation heatmap (prettified, with NaN handling and fallback)
- Pairwise scatter with identity line
- PCA of runs
- Peptide coefficient-of-variation histogram (auto-adjust bins or skip if insufficient data)
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from pathlib import Path
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

# Mapping of CPTAC study IDs to human-readable labels
STUDY_LABELS = {
    "PDC000110": "Prospective_Ovarian_JHU_Proteome",
    "PDC000116": "Prospective_Colon_PNNL_Proteome_Qeplus",
    "PDC000120": "Prospective_Breast_BI_Proteome",
    "PDC000125": "CPTAC_UCEC_Discovery_Study_-_Proteome",
    "PDC000127": "CPTAC_CCRCC_Discovery_Study_-_Proteome",
    "PDC000153": "CPTAC_LUAD_Discovery_Study_-_Proteome",
    "PDC000204": "CPTAC_GBM_Discovery_Study_-_Proteome",
    "PDC000221": "CPTAC_HNSCC_Discovery_Study_-_Proteome",
    "PDC000234": "CPTAC_LSCC_Discovery_Study_-_Proteome",
    "PDC000270": "CPTAC_PDA_Discovery_Study_-_Proteome",
}


def shorten_labels(full_run_names):
    """
    Given a list (or Index) of full run names like:
      "01CPTAC_OVprospective_Proteome_JHU_20161209"
    return just the date substring "20161209".  
    If no date‐style substring is found, fall back to the last underscore portion.
    """
    short = []
    for rn in full_run_names:
        parts = rn.split("_")
        last = parts[-1]
        if len(last) == 8 and last.isdigit():
            short.append(last)
        else:
            short.append(last)
    return short


def plot_pretty_clustermap(
    ref_df: pd.DataFrame,
    out_path: Path,
    title: str,
    rotate_labels: bool = True
):
    """
    ref_df:       DataFrame of shape (n_peptides × n_runs), whose columns are "Run" names.
    out_path:     Path where the .png will be saved.
    title:        Main title for the plot.
    rotate_labels: If True, x tick labels are rotated 45°; otherwise they remain horizontal.
    """
    # 1) Drop any run columns that have zero variance
    stds = ref_df.std(axis=0, skipna=True)
    nonconst = stds[stds > 0].index.tolist()
    if len(nonconst) < 2:
        print(f"⚠️  Skipping clustermap {title}: need ≥2 non-constant runs")
        return
    subset = ref_df[nonconst]

    # 2) Compute Spearman correlation matrix
    corr = subset.corr(method="spearman")
    # Replace any NaNs/infinite with the minimum finite correlation
    if not np.isfinite(corr.values).all():
        minval = corr.values[np.isfinite(corr.values)].min()
        corr = corr.fillna(minval)

    # 3) Shorten each run name to its date (or final underscore block)
    short_names = shorten_labels(corr.columns.tolist())
    corr.index = short_names
    corr.columns = short_names

    # 4) Choose figure size: width & height scale with n_runs
    n = len(short_names)
    fig_size = max(3, n * 0.2)

    # 5) Build a custom colorbar axes to avoid overlap
    fig = plt.figure(figsize=(fig_size + 1, fig_size + 0.6))
    grid = plt.GridSpec(1, 40, left=0.05, right=0.92, top=0.88, bottom=0.10)
    ax_heat = fig.add_subplot(grid[:, :35])   # main heatmap
    ax_cbar = fig.add_subplot(grid[:, 37:40]) # colorbar

    # 6) Hierarchical clustering to reorder rows/columns
    dist = 1.0 - corr.values
    np.fill_diagonal(dist, 0.0)
    condensed = squareform(dist, checks=False)
    linkage_matrix = linkage(condensed, method="average")
    ordered_idx = leaves_list(linkage_matrix)

    # Reorder correlation matrix
    corr_ord = corr.iloc[ordered_idx, :].iloc[:, ordered_idx]

    # 7) Plot the heatmap
    im = ax_heat.imshow(
        corr_ord.values,
        vmin=0.70,       # fixed lower bound for QC consistency
        vmax=1.00,       # fixed upper bound
        cmap="viridis",
        aspect="equal",
        interpolation="nearest"
    )
    tick_labels = [corr_ord.columns[i] for i in range(len(ordered_idx))]

    ax_heat.set_xticks(np.arange(len(tick_labels)))
    ax_heat.set_yticks(np.arange(len(tick_labels)))
    ax_heat.set_xticklabels(
        tick_labels,
        fontsize=6,
        rotation=45 if rotate_labels else 0,
        ha="right"
    )
    ax_heat.set_yticklabels(tick_labels, fontsize=6)
    ax_heat.set_title(title, fontsize=10, pad=8)

    # 8) Add colorbar with fixed tick positions
    cbar = plt.colorbar(im, cax=ax_cbar, orientation="vertical", label="Spearman ρ")
    cbar.set_ticks([0.70, 0.80, 0.90, 1.00])
    cbar.ax.set_yticklabels(["0.70", "0.80", "0.90", "1.00"])
    cbar.ax.yaxis.set_minor_locator(plt.NullLocator())

    # 9) Final layout adjustments
    plt.tight_layout(rect=[0, 0, 1, 1])
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def load_reference_matrix(study_dir: Path, join: str = "inner") -> pd.DataFrame:
    """
    Load all *_matrix.parquet under study_dir,
    extract the REF_ channel per run, and return
    a DataFrame (peptide × runs), aligned by index.

    join: 'inner' for intersection of peptides,
          'outer' for union (with NaNs).
    """
    series_list = []
    for pq in sorted(study_dir.rglob("*_matrix.parquet")):
        df = pd.read_parquet(pq)
        ref_cols = [c for c in df.columns if c.startswith("REF_")]
        if not ref_cols:
            continue
        col = ref_cols[0]
        s = df[col].copy()
        run_name = pq.stem.replace("_matrix", "")
        s.name = run_name
        series_list.append(s)

    if not series_list:
        return None

    return pd.concat(series_list, axis=1, join=join)


def plot_pair_scatter(ref_df: pd.DataFrame, out_path: Path, title: str):
    """Scatter of first two REF runs with identity line & equal aspect"""
    xcol, ycol = ref_df.columns[:2]
    x = ref_df[xcol]
    y = ref_df[ycol]
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(x, y, s=1, alpha=0.3)
    # Identity line
    lims = [
        np.nanmin([x.min(), y.min()]),
        np.nanmax([x.max(), y.max()])
    ]
    ax.plot(lims, lims, "--", linewidth=0.5)
    ax.set_aspect("equal", "box")
    ax.set_xlabel(xcol, fontsize=7)
    ax.set_ylabel(ycol, fontsize=7)
    ax.set_title(title, fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_pca(ref_df: pd.DataFrame, out_path: Path, title: str):
    """PCA of REF runs (2 components)"""
    data = ref_df.fillna(ref_df.mean(axis=0)).T
    pca = PCA(n_components=2)
    coords = pca.fit_transform(data)
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(coords[:, 0], coords[:, 1], s=10)
    for i, label in enumerate(data.index):
        ax.text(coords[i, 0], coords[i, 1], label, fontsize=6)
    ax.set_title(title + " (PCA)", fontsize=8)
    ax.set_xlabel("PC1", fontsize=7)
    ax.set_ylabel("PC2", fontsize=7)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_cv_histogram(ref_df: pd.DataFrame, out_path: Path, title: str):
    """
    Histogram of per-peptide coefficient of variation.
    Automatically reduces the number of bins if “Too many bins” error arises,
    or skips entirely if insufficient non‐zero data.
    """
    cv = ref_df.std(axis=1) / ref_df.mean(axis=1)
    cv_nonan = cv.dropna()
    cv_nonzero = cv_nonan[cv_nonan > 0.0]

    if len(cv_nonzero) < 2:
        print(f"⚠️  Skipping CV histogram {title}: need ≥2 non-zero CV values")
        return

    bins = 50
    while bins > 1:
        try:
            fig, ax = plt.subplots(figsize=(4, 3))
            ax.hist(cv_nonzero, bins=bins, alpha=0.7)
            ax.set_title(title + " (CV distribution)", fontsize=8)
            ax.set_xlabel("CV", fontsize=7)
            ax.set_ylabel("Count", fontsize=7)
            fig.tight_layout()
            fig.savefig(out_path, dpi=150, bbox_inches="tight")
            plt.close(fig)
            return
        except ValueError as ve:
            if "Too many bins" in str(ve):
                bins = bins // 2
            else:
                print(f"⚠️  Skipping CV histogram {title}: {ve}")
                return

    print(f"⚠️  Skipping CV histogram {title}: could not find suitable bin count.")


if __name__ == "__main__":
    # Directories
    mat_root   = Path("data/for_analysis/abundance_matrices")
    qc_dir     = Path("reports/qc")
    dirs = {
        'heatmap': qc_dir / 'heatmaps',
        'scatter': qc_dir / 'scatter',
        'pca':     qc_dir / 'pca',
        'cv':      qc_dir / 'cv',
        'corrcsv': qc_dir / 'correlation_csv'
    }
    # Create output folders
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)

    # Loop over plexes (TMT10, TMT11) and studies
    for plex_dir in sorted(mat_root.iterdir()):
        if not plex_dir.is_dir():
            continue
        for study_dir in sorted(plex_dir.iterdir()):
            if not study_dir.is_dir():
                continue

            study_id  = study_dir.name
            plex_name = plex_dir.name
            label     = STUDY_LABELS.get(study_id, study_id)

            # Load only the REF_ columns for this study (inner join of peptides)
            ref_df = load_reference_matrix(study_dir, join="inner")
            if ref_df is None or ref_df.shape[1] < 2:
                print(f"⚠️  Skipping {plex_name}/{label}: need ≥2 REF runs")
                continue

            print(f"Processing {plex_name}/{label}: "
                  f"{ref_df.shape[0]} peptides × {ref_df.shape[1]} REF runs")
            title_base = f"{study_id}_{plex_name} – {label}"

            # --- Save raw Spearman correlation matrix as CSV ---
            # (Drop zero-variance runs, compute Spearman rho, fill NaNs, then save)
            stds = ref_df.std(axis=0, skipna=True)
            nonconst = stds[stds > 0].index.tolist()
            if len(nonconst) >= 2:
                corr_df = ref_df[nonconst].corr(method="spearman")
                if not np.isfinite(corr_df.values).all():
                    minval = corr_df.values[np.isfinite(corr_df.values)].min()
                    corr_df = corr_df.fillna(minval)
                csv_out = dirs['corrcsv'] / f"{plex_name}_{study_id}_REF_correlation_matrix.csv"
                corr_df.to_csv(csv_out, index=True)
                print(f"  • WROTE correlation CSV: {csv_out}")
            else:
                print(f"⚠️  Skipping CSV for {plex_name}/{label}: need ≥2 non-constant runs")

            # 1) Clustered Spearman correlation heatmap (prettified)
            hpath = dirs['heatmap'] / f"{plex_name}_{study_id}_REF_corr_pretty.png"
            plot_pretty_clustermap(
                ref_df,
                hpath,
                title_base + "\nSpearman ρ (REF runs)"
            )

            # 2) Pairwise scatter (first two runs)
            spath = dirs['scatter'] / f"{plex_name}_{study_id}_REF_scatter.png"
            plot_pair_scatter(ref_df, spath, "REF: first two runs")

            # 3) PCA plot of REF runs
            ppath = dirs['pca'] / f"{plex_name}_{study_id}_REF_pca.png"
            plot_pca(ref_df, ppath, title_base)

            # 4) CV histogram
            cvpath = dirs['cv'] / f"{plex_name}_{study_id}_REF_cv.png"
            plot_cv_histogram(ref_df, cvpath, title_base)

            print(f"  • Figures and CSV written for {plex_name}/{label}\n")
