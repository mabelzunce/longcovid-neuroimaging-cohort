#!/usr/bin/env python3
"""
PET intersubject metabolic covariance analysis for Long COVID vs Controls (CERMEP only).

What this script computes
-------------------------
1) Group-level ROI-by-ROI connectivity (Pearson correlation across subjects)
   - Long COVID (Grupo == 'LC')
   - Controls from CERMEP only (Grupo == 'Control' and SITE == 'CERMEP')

2) Whole-network metrics for each group
   - Mean correlation (off-diagonal)
   - Mean absolute correlation (off-diagonal)
   - Positive edge density
   - Graph metrics on thresholded positive network (top-k density):
     global efficiency, clustering, transitivity, average node strength.

3) Frontal-focused metrics (ROIs starting with 'FL-')
   - Within-frontal mean correlation
   - Frontal-to-nonfrontal mean correlation

4) Permutation-based group comparisons for scalar metrics.

5) Edge-level differential connectivity with Fisher-z difference and permutation p-values,
   plus FDR correction.

Outputs are written to data/results/pet_network_analysis and data/plots/pet_network_analysis.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns


META_COLS = {"ID", "Grupo", "SITE", "edad", "sexo"}


@dataclass
class GroupData:
    name: str
    df: pd.DataFrame
    X: np.ndarray


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction."""
    p = np.asarray(p_values, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]

    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)

    out = np.empty_like(q)
    out[order] = q
    return out


def offdiag_values(mat: np.ndarray) -> np.ndarray:
    m = mat.shape[0]
    tri = np.triu_indices(m, k=1)
    return mat[tri]


def correlation_matrix(X: np.ndarray) -> np.ndarray:
    """ROI-by-ROI Pearson correlation using subjects as observations."""
    if X.ndim != 2:
        raise ValueError("X must be 2D: subjects x ROIs")
    if X.shape[0] < 3:
        raise ValueError("Need at least 3 subjects for stable correlations")
    return np.corrcoef(X, rowvar=False)


def threshold_positive_network(corr: np.ndarray, density: float = 0.2) -> np.ndarray:
    """Keep top positive edges according to requested density."""
    n = corr.shape[0]
    upper = np.triu_indices(n, k=1)
    vals = corr[upper].copy()

    positive_mask = vals > 0
    pos_vals = vals[positive_mask]

    adj = np.zeros_like(corr)
    if pos_vals.size == 0:
        return adj

    n_possible = len(vals)
    n_keep = max(1, int(round(density * n_possible)))
    n_keep = min(n_keep, pos_vals.size)

    thresh = np.partition(pos_vals, -n_keep)[-n_keep]
    keep_upper = vals >= thresh
    keep_upper &= vals > 0

    adj[upper] = vals * keep_upper
    adj = adj + adj.T
    return adj


def graph_metrics_from_corr(corr: np.ndarray, threshold_density: float = 0.2) -> Dict[str, float]:
    """Compute graph metrics from thresholded positive weighted adjacency."""
    adj = threshold_positive_network(corr, density=threshold_density)
    G = nx.from_numpy_array(adj)

    # Remove isolated nodes only for some metrics requiring connected components handling
    weighted_degrees = np.array([d for _, d in G.degree(weight="weight")], dtype=float)

    # Global efficiency in NetworkX is unweighted. We compute on binarized graph.
    G_bin = nx.from_numpy_array((adj > 0).astype(int))

    metrics = {
        "n_nodes": float(G.number_of_nodes()),
        "n_edges": float(G.number_of_edges()),
        "density": float(nx.density(G_bin)),
        "mean_strength": float(np.nanmean(weighted_degrees)),
        "global_efficiency": float(nx.global_efficiency(G_bin)),
        "avg_clustering_weighted": float(nx.average_clustering(G, weight="weight")),
        "transitivity": float(nx.transitivity(G_bin)),
        "mean_corr": float(np.nanmean(offdiag_values(corr))),
        "mean_abs_corr": float(np.nanmean(np.abs(offdiag_values(corr)))),
        "positive_edge_fraction": float(np.mean(offdiag_values(corr) > 0)),
    }
    return metrics


def frontal_metrics(corr: np.ndarray, roi_names: List[str]) -> Dict[str, float]:
    roi_names = list(roi_names)
    frontal_idx = np.array([i for i, r in enumerate(roi_names) if r.startswith("FL-")], dtype=int)
    non_frontal_idx = np.array([i for i, r in enumerate(roi_names) if not r.startswith("FL-")], dtype=int)

    if frontal_idx.size < 2:
        raise ValueError("Not enough frontal ROIs found (prefix 'FL-').")

    ff = corr[np.ix_(frontal_idx, frontal_idx)]
    fn = corr[np.ix_(frontal_idx, non_frontal_idx)] if non_frontal_idx.size > 0 else np.array([np.nan])

    ff_vals = offdiag_values(ff)

    return {
        "n_frontal_rois": float(frontal_idx.size),
        "within_frontal_mean_corr": float(np.nanmean(ff_vals)),
        "within_frontal_mean_abs_corr": float(np.nanmean(np.abs(ff_vals))),
        "frontal_to_nonfrontal_mean_corr": float(np.nanmean(fn)),
        "frontal_to_nonfrontal_mean_abs_corr": float(np.nanmean(np.abs(fn))),
    }


def permutation_test_metric(
    X1: np.ndarray,
    X2: np.ndarray,
    metric_func,
    n_perm: int = 1000,
    seed: int = 123,
) -> Tuple[float, float]:
    """Two-sided permutation test on scalar metric of correlation matrix."""
    rng = np.random.default_rng(seed)

    obs = metric_func(correlation_matrix(X1)) - metric_func(correlation_matrix(X2))

    X_all = np.vstack([X1, X2])
    n1 = X1.shape[0]
    n = X_all.shape[0]

    diffs = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        perm = rng.permutation(n)
        g1 = X_all[perm[:n1], :]
        g2 = X_all[perm[n1:], :]
        diffs[i] = metric_func(correlation_matrix(g1)) - metric_func(correlation_matrix(g2))

    p = (np.sum(np.abs(diffs) >= np.abs(obs)) + 1.0) / (n_perm + 1.0)
    return float(obs), float(p)


def permutation_edge_tests(
    X1: np.ndarray,
    X2: np.ndarray,
    roi_names: List[str],
    n_perm: int = 500,
    seed: int = 123,
) -> pd.DataFrame:
    """Permutation tests for each ROI-ROI edge difference in correlation."""
    rng = np.random.default_rng(seed)

    c1 = correlation_matrix(X1)
    c2 = correlation_matrix(X2)

    n_roi = c1.shape[0]
    tri = np.triu_indices(n_roi, k=1)

    obs_diff = c1[tri] - c2[tri]

    X_all = np.vstack([X1, X2])
    n1 = X1.shape[0]
    n = X_all.shape[0]

    perm_diffs = np.empty((n_perm, obs_diff.size), dtype=float)
    for i in range(n_perm):
        perm = rng.permutation(n)
        g1 = X_all[perm[:n1], :]
        g2 = X_all[perm[n1:], :]
        cp1 = correlation_matrix(g1)
        cp2 = correlation_matrix(g2)
        perm_diffs[i, :] = cp1[tri] - cp2[tri]

    p_vals = (np.sum(np.abs(perm_diffs) >= np.abs(obs_diff), axis=0) + 1.0) / (n_perm + 1.0)
    q_vals = benjamini_hochberg(p_vals)

    rows = []
    for k, (i, j) in enumerate(zip(*tri)):
        rows.append(
            {
                "roi_i": roi_names[i],
                "roi_j": roi_names[j],
                "diff_corr_lc_minus_control": obs_diff[k],
                "p_perm": p_vals[k],
                "q_fdr": q_vals[k],
            }
        )

    out = pd.DataFrame(rows).sort_values("p_perm", ascending=True)
    return out


def save_heatmap(mat: np.ndarray, title: str, out_png: Path, vmin=-1, vmax=1):
    plt.figure(figsize=(10, 8))
    sns.heatmap(mat, cmap="coolwarm", vmin=vmin, vmax=vmax, center=0, square=True, cbar_kws={"label": "r"})
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def build_group_data(df: pd.DataFrame, roi_cols: List[str]) -> Dict[str, GroupData]:
    control = df[(df["Grupo"] == "Control") & (df["SITE"] == "CERMEP")].copy()
    lc = df[df["Grupo"] == "LC"].copy()

    if control.empty:
        raise ValueError("No Control subjects found at SITE=CERMEP")
    if lc.empty:
        raise ValueError("No Long COVID subjects found (Grupo=LC)")

    control_x = control[roi_cols].apply(pd.to_numeric, errors="coerce").dropna(axis=0).to_numpy()
    lc_x = lc[roi_cols].apply(pd.to_numeric, errors="coerce").dropna(axis=0).to_numpy()

    if control_x.shape[0] < 5 or lc_x.shape[0] < 5:
        raise ValueError(
            f"Too few usable subjects after NA removal. Control={control_x.shape[0]}, LC={lc_x.shape[0]}"
        )

    return {
        "Control_CERMEP": GroupData(name="Control_CERMEP", df=control, X=control_x),
        "LongCOVID": GroupData(name="LongCOVID", df=lc, X=lc_x),
    }


def main():
    parser = argparse.ArgumentParser(description="PET intersubject network analysis")
    parser.add_argument(
        "--input-csv",
        default="/home/martin/data/UNSAM/CovidProject/longcovid-neuroimaging-cohort/data/processed_pet/PET_CERMEP_CEUNIM_harmonized.csv",
    )
    parser.add_argument(
        "--out-results",
        default="/home/martin/data/UNSAM/CovidProject/longcovid-neuroimaging-cohort/data/results/pet_network_analysis",
    )
    parser.add_argument(
        "--out-plots",
        default="/home/martin/data/UNSAM/CovidProject/longcovid-neuroimaging-cohort/data/plots/pet_network_analysis",
    )
    parser.add_argument("--n-perm-metrics", type=int, default=1000)
    parser.add_argument("--n-perm-edges", type=int, default=500)
    parser.add_argument("--threshold-density", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=123)
    args = parser.parse_args()

    in_csv = Path(args.input_csv)
    out_results = Path(args.out_results)
    out_plots = Path(args.out_plots)
    out_results.mkdir(parents=True, exist_ok=True)
    out_plots.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(in_csv)

    # ROI columns = numeric PET regions excluding metadata
    roi_cols = [c for c in df.columns if c not in META_COLS]
    # keep only numeric-looking columns
    roi_cols = [c for c in roi_cols if pd.api.types.is_numeric_dtype(df[c]) or pd.to_numeric(df[c], errors="coerce").notna().mean() > 0.9]

    groups = build_group_data(df, roi_cols)

    # Connectivity matrices
    corr_control = correlation_matrix(groups["Control_CERMEP"].X)
    corr_lc = correlation_matrix(groups["LongCOVID"].X)
    diff = corr_lc - corr_control

    pd.DataFrame(corr_control, index=roi_cols, columns=roi_cols).to_csv(out_results / "connectivity_control_CERMEP.csv")
    pd.DataFrame(corr_lc, index=roi_cols, columns=roi_cols).to_csv(out_results / "connectivity_longcovid.csv")
    pd.DataFrame(diff, index=roi_cols, columns=roi_cols).to_csv(out_results / "connectivity_diff_lc_minus_control.csv")

    # Network metrics
    m_control = graph_metrics_from_corr(corr_control, threshold_density=args.threshold_density)
    m_lc = graph_metrics_from_corr(corr_lc, threshold_density=args.threshold_density)

    metrics_df = pd.DataFrame(
        [
            {"group": "Control_CERMEP", **m_control},
            {"group": "LongCOVID", **m_lc},
        ]
    )
    metrics_df.to_csv(out_results / "network_metrics_by_group.csv", index=False)

    # Frontal metrics
    f_control = frontal_metrics(corr_control, roi_cols)
    f_lc = frontal_metrics(corr_lc, roi_cols)
    frontal_df = pd.DataFrame(
        [
            {"group": "Control_CERMEP", **f_control},
            {"group": "LongCOVID", **f_lc},
        ]
    )
    frontal_df.to_csv(out_results / "frontal_metrics_by_group.csv", index=False)

    # Permutation tests for selected scalar metrics
    metric_functions = {
        "mean_abs_corr": lambda c: float(np.nanmean(np.abs(offdiag_values(c)))),
        "global_efficiency": lambda c: graph_metrics_from_corr(c, threshold_density=args.threshold_density)["global_efficiency"],
        "avg_clustering_weighted": lambda c: graph_metrics_from_corr(c, threshold_density=args.threshold_density)["avg_clustering_weighted"],
        "within_frontal_mean_corr": lambda c: frontal_metrics(c, roi_cols)["within_frontal_mean_corr"],
        "frontal_to_nonfrontal_mean_corr": lambda c: frontal_metrics(c, roi_cols)["frontal_to_nonfrontal_mean_corr"],
    }

    perm_rows = []
    X_lc = groups["LongCOVID"].X
    X_ctrl = groups["Control_CERMEP"].X

    for metric_name, func in metric_functions.items():
        obs_diff, p = permutation_test_metric(
            X_lc,
            X_ctrl,
            metric_func=func,
            n_perm=args.n_perm_metrics,
            seed=args.seed,
        )
        perm_rows.append(
            {
                "metric": metric_name,
                "observed_diff_lc_minus_control": obs_diff,
                "p_perm": p,
            }
        )

    perm_df = pd.DataFrame(perm_rows)
    perm_df["q_fdr"] = benjamini_hochberg(perm_df["p_perm"].to_numpy())
    perm_df.to_csv(out_results / "permutation_metric_tests.csv", index=False)

    # Edge-level tests
    edge_df = permutation_edge_tests(
        X_lc,
        X_ctrl,
        roi_names=roi_cols,
        n_perm=args.n_perm_edges,
        seed=args.seed,
    )
    edge_df.to_csv(out_results / "edge_differences_permutation.csv", index=False)
    edge_df.head(30).to_csv(out_results / "edge_differences_top30.csv", index=False)

    # Frontal-only edge summary
    frontal_edges = edge_df[
        edge_df["roi_i"].str.startswith("FL-") | edge_df["roi_j"].str.startswith("FL-")
    ].copy()
    frontal_edges.to_csv(out_results / "edge_differences_frontal_involved.csv", index=False)

    # Plots
    save_heatmap(corr_control, "Control (CERMEP) ROI connectivity", out_plots / "heatmap_connectivity_control_CERMEP.png")
    save_heatmap(corr_lc, "Long COVID ROI connectivity", out_plots / "heatmap_connectivity_longcovid.png")
    save_heatmap(diff, "Connectivity difference (LC - Control)", out_plots / "heatmap_connectivity_difference_lc_minus_control.png")

    # Plot selected scalar metrics by group
    selected = ["mean_abs_corr", "global_efficiency", "avg_clustering_weighted"]
    melt_df = metrics_df[["group"] + selected].melt(id_vars="group", var_name="metric", value_name="value")
    plt.figure(figsize=(9, 4))
    sns.barplot(data=melt_df, x="metric", y="value", hue="group")
    plt.title("Whole-network metrics by group")
    plt.tight_layout()
    plt.savefig(out_plots / "bar_network_metrics_by_group.png", dpi=300)
    plt.close()

    # Frontal metrics plot
    fsel = ["within_frontal_mean_corr", "frontal_to_nonfrontal_mean_corr"]
    fplot = frontal_df[["group"] + fsel].melt(id_vars="group", var_name="metric", value_name="value")
    plt.figure(figsize=(9, 4))
    sns.barplot(data=fplot, x="metric", y="value", hue="group")
    plt.title("Frontal-involved connectivity metrics")
    plt.tight_layout()
    plt.savefig(out_plots / "bar_frontal_metrics_by_group.png", dpi=300)
    plt.close()

    # Save analysis log
    with open(out_results / "analysis_summary.txt", "w", encoding="utf-8") as f:
        f.write("PET intersubject network analysis completed\n")
        f.write(f"Input CSV: {in_csv}\n")
        f.write(f"N controls at CERMEP (raw rows): {groups['Control_CERMEP'].df.shape[0]}\n")
        f.write(f"N long COVID (raw rows): {groups['LongCOVID'].df.shape[0]}\n")
        f.write(f"N controls used after NA filtering: {groups['Control_CERMEP'].X.shape[0]}\n")
        f.write(f"N long COVID used after NA filtering: {groups['LongCOVID'].X.shape[0]}\n")
        f.write(f"N ROIs used: {len(roi_cols)}\n")
        f.write(f"Threshold density for graph metrics: {args.threshold_density}\n")
        f.write(f"Permutation tests (metrics): {args.n_perm_metrics}\n")
        f.write(f"Permutation tests (edges): {args.n_perm_edges}\n")

    print("Done.")
    print(f"Results: {out_results}")
    print(f"Plots:   {out_plots}")


if __name__ == "__main__":
    main()
