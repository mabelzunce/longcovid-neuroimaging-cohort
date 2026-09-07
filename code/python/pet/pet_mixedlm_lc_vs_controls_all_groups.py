from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests


def _project_root() -> Path:
    # File is expected at code/python/pet/
    return Path(__file__).resolve().parents[3]


def _build_long_dataframe(df: pd.DataFrame) -> Tuple[pd.DataFrame, List[str]]:
    required_cols = ["ID", "group", "site", "sex", "age"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Keep only LC and controls from all sites (three cohorts in this dataset).
    keep = (df["group"].eq("LC")) | (df["group"].eq("Controls"))
    sub_df = df.loc[keep].copy()

    excluded = {"ID", "group", "site", "sex", "age", "weight", "analysis_group"}
    candidate_roi_cols = [c for c in sub_df.columns if c not in excluded]

    roi_cols: List[str] = []
    for col in candidate_roi_cols:
        converted = pd.to_numeric(sub_df[col], errors="coerce")
        if converted.notna().sum() > 0:
            sub_df[col] = converted
            roi_cols.append(col)

    if not roi_cols:
        raise ValueError("No numeric ROI columns detected.")

    long_df = sub_df.melt(
        id_vars=["ID", "group", "site", "sex", "age"],
        value_vars=roi_cols,
        var_name="ROI",
        value_name="SUVR",
    )

    long_df = long_df.rename(columns={"ID": "subject"})
    long_df["age"] = pd.to_numeric(long_df["age"], errors="coerce")
    long_df["age_c"] = long_df["age"] - long_df["age"].mean(skipna=True)

    # Ensure stable factor references for interpretation.
    long_df["group"] = pd.Categorical(long_df["group"], categories=["Controls", "LC"])
    long_df["site"] = pd.Categorical(long_df["site"])
    long_df["sex"] = pd.Categorical(long_df["sex"])
    long_df["ROI"] = pd.Categorical(long_df["ROI"], categories=roi_cols)
    long_df["subject"] = long_df["subject"].astype(str)

    long_df = long_df.dropna(subset=["SUVR", "group", "ROI", "age_c", "sex", "site", "subject"]).copy()

    return long_df, roi_cols


def _find_group_lc_term(params_index: pd.Index) -> str:
    matches = [p for p in params_index if p.startswith("C(group") and "[T.LC]" in p and ":" not in p]
    if not matches:
        raise ValueError("Could not find LC group main-effect term in model parameters.")
    return matches[0]


def _find_group_roi_interaction(params_index: pd.Index, roi: str) -> str | None:
    # For non-reference ROIs, statsmodels encodes one interaction coefficient.
    for name in params_index:
        if "C(group" in name and "[T.LC]" in name and ":C(ROI" in name and f"[T.{roi}]" in name:
            return name
        if "C(ROI" in name and f"[T.{roi}]" in name and ":C(group" in name and "[T.LC]" in name:
            return name
    return None


def _roi_lc_contrasts(result, roi_cols: List[str]) -> pd.DataFrame:
    params = result.fe_params
    cov = result.cov_params().loc[params.index, params.index]
    names = list(params.index)
    name_to_idx: Dict[str, int] = {n: i for i, n in enumerate(names)}

    group_term = _find_group_lc_term(params.index)

    records = []
    for roi in roi_cols:
        vec = np.zeros(len(names), dtype=float)
        vec[name_to_idx[group_term]] = 1.0

        inter_term = _find_group_roi_interaction(params.index, roi)
        if inter_term is not None:
            vec[name_to_idx[inter_term]] = 1.0

        est = float(np.dot(vec, params.values))
        var = float(np.dot(vec, np.dot(cov.values, vec)))
        se = float(np.sqrt(var)) if var >= 0 else np.nan
        z = est / se if se and np.isfinite(se) and se > 0 else np.nan
        p = float(2.0 * norm.sf(abs(z))) if np.isfinite(z) else np.nan

        records.append(
            {
                "ROI": roi,
                "estimate_lc_minus_control": est,
                "se": se,
                "z": z,
                "p": p,
                "is_reference_roi": inter_term is None,
            }
        )

    out = pd.DataFrame(records).sort_values("p").reset_index(drop=True)
    if len(out) > 0:
        rej, p_fdr, _, _ = multipletests(out["p"], alpha=0.05, method="fdr_bh")
        out["p_fdr"] = p_fdr
        out["significant_fdr"] = rej
    else:
        out["p_fdr"] = pd.Series(dtype=float)
        out["significant_fdr"] = pd.Series(dtype=bool)

    return out


def main() -> None:
    repo_root = _project_root()

    input_csv = repo_root / "data" / "processed_pet" / "merged_summary_csv" / "summary_suvr_mean_harmonized.csv"
    output_dir = repo_root / "results" / "pet" / "ancova_rois"
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_csv.exists():
        raise FileNotFoundError(f"Input file not found: {input_csv}")

    df = pd.read_csv(input_csv)
    pet_long, roi_cols = _build_long_dataframe(df)

    print("Subjects by site/group:")
    print(df.groupby(["site", "group"]).size())
    print(f"Rows in long dataframe: {len(pet_long):,}")
    print(f"ROIs in model: {len(roi_cols)}")

    roi_ref = roi_cols[0]
    formula = (
        f'SUVR ~ C(group, Treatment(reference="Controls")) * '
        f'C(ROI, Treatment(reference="{roi_ref}")) + '
        'age_c + C(sex) + C(site)'
    )

    # Equivalent structure to: lmer(SUVR ~ group * ROI + age_c + sex + site + (1 | subject), data = pet_long)
    model = smf.mixedlm(formula, data=pet_long, groups=pet_long["subject"])
    result = model.fit(reml=False, method="lbfgs", maxiter=200)

    summary_path = output_dir / "mixedlm_lc_vs_controls_all_groups_summary.txt"
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write(str(result.summary()))

    fixed_effects = (
        pd.DataFrame(
            {
                "term": result.params.index,
                "coef": result.params.values,
                "se": result.bse.values,
                "z": result.tvalues.values,
                "p": result.pvalues.values,
            }
        )
        .sort_values("p")
        .reset_index(drop=True)
    )
    fixed_effects.to_csv(output_dir / "mixedlm_lc_vs_controls_all_groups_fixed_effects.csv", index=False)

    roi_contrasts = _roi_lc_contrasts(result, roi_cols)
    roi_contrasts.to_csv(output_dir / "mixedlm_lc_vs_controls_all_groups_roi_contrasts.csv", index=False)

    sig = roi_contrasts.loc[roi_contrasts["significant_fdr"]].copy()
    sig.to_csv(output_dir / "mixedlm_lc_vs_controls_all_groups_roi_contrasts_significant.csv", index=False)

    print(f"Saved model summary: {summary_path}")
    print(
        "Saved ROI contrasts: "
        f"{output_dir / 'mixedlm_lc_vs_controls_all_groups_roi_contrasts.csv'}"
    )
    print(f"Significant ROIs (FDR < 0.05): {len(sig)}")


if __name__ == "__main__":
    main()