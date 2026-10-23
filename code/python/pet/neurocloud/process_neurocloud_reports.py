# -*- coding: utf-8 -*-
# Reads neurocloud reports and computes within-group z-scores, statistics, and plots.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm

from scipy.stats import ttest_1samp
import warnings, os, math
import numpy as np
import pandas as pd
from scipy.stats import ttest_1samp, wilcoxon

warnings.filterwarnings('ignore')

paths_neurocloud_reports = "/home/martin/data_imaging/CovidProject/FullStudy/PET/InformesNeurocloud/"
#paths_informes_covid =  "/content/drive/MyDrive/NeurocloudPET/InformesNeurocloud"
#paths_informes_fleni = "/content/drive/MyDrive/NeurocloudPET/InformesNeurocloudFleni"
rois_to_lobe_csv = os.path.join('./data/pet/', 'hammers_roi_to_lobe.csv')
output_path = os.path.join('./data/processed_pet/', 'neurocloud/')
plots_path = os.path.join('./results/pet/', 'neurocloud/')
os.makedirs(output_path, exist_ok=True)
os.makedirs(plots_path, exist_ok=True)

# Load rois to lobe mapping
rois_to_lobe = pd.read_csv(rois_to_lobe_csv)

# Set plotting style
sns.set_style('whitegrid')
plt.rcParams['figure.figsize'] = (12, 6)


def load_long_zscores(base_dir,
                      csv_rel_path=("data_es", "es", "rois_analysis.csv"),
                      sep=';'):
    """
    Load ROI z-scores from multiple subjects and return a long-format DataFrame.

    Parameters
    ----------
    base_dir : str
        Base directory containing subject folders.
    csv_rel_path : tuple
        Relative path to the CSV inside each subject folder.
    sep : str
        CSV separator.

    Returns
    -------
    pd.DataFrame
        Long-format dataframe with columns:
        SUBJECT, ROI_NAME, HEMI, ZSCORE
    """

    dfs = []

    for subj in sorted(os.listdir(base_dir)):
        subj_path = os.path.join(base_dir, subj)

        if not os.path.isdir(subj_path):
            continue

        csv_path = os.path.join(subj_path, *csv_rel_path)

        if not os.path.exists(csv_path):
            continue

        df = pd.read_csv(csv_path, sep=sep)

        long_df = pd.DataFrame({
            'SUBJECT': subj,
            'ROI_NAME': df['ROI_NAME'],
            'ZSCORE_LH': df['ZSCORE_LH'],
            'ZSCORE_RH': df['ZSCORE_RH']
        })

        long_df = long_df.melt(
            id_vars=['SUBJECT', 'ROI_NAME'],
            value_vars=['ZSCORE_LH', 'ZSCORE_RH'],
            var_name='HEMI',
            value_name='ZSCORE'
        )

        long_df['HEMI'] = long_df['HEMI'].str.replace('ZSCORE_', '', regex=False)

        dfs.append(long_df)

    if not dfs:
        raise ValueError("No CSV files were found.")

    return pd.concat(dfs, ignore_index=True)

def plot_zscore_hist(ax, data,
                     bins=15,
                     x_range=(-3, 3),
                     title=None):
    data = data.dropna()

    ax.hist(
        data,
        bins=bins,
        range=x_range,
        alpha=0.8
    )

    ax.axvline(0, linestyle='--', linewidth=1)

    if title is not None:
        ax.set_title(title, fontsize=9)

    ax.grid(True, alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

def plot_hist_grid(df_wide,
                   n_cols=6,
                   bins=15,
                   x_range=(-3, 3),
                   figsize_per_subplot=(4, 3)):
    """
    Plot a grid of z-score histograms for each column in df_wide.
    """

    cols = df_wide.columns
    n = len(cols)

    n_rows = math.ceil(n / n_cols)

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(figsize_per_subplot[0] * n_cols,
                 figsize_per_subplot[1] * n_rows),
        sharex=True,
        sharey=True
    )

    axes = axes.flatten()

    for ax, col in zip(axes, cols):
        plot_zscore_hist(
            ax,
            df_wide[col],
            bins=bins,
            x_range=x_range,
            title=col
        )

    # apagar subplots vacíos
    for ax in axes[len(cols):]:
        ax.axis('off')

    fig.supxlabel("Z-score")
    fig.supylabel("Frecuencia")

    plt.tight_layout()
    plt.show()


def compute_roi_stats(df_wide, min_n=10):
    """
    One-sample stats vs 0 for each ROI:
    - Mean
    - t-test (one-sided)
    - Wilcoxon (one-sided)
    - Cohen's d
    """

    rows = []

    for roi in df_wide.columns:
        data = df_wide[roi].dropna()
        n = len(data)

        if n < min_n:
            continue

        mean_z = data.mean()
        std_z = data.std(ddof=1)

        # t-test (one-sided)
        tstat, pval_t = ttest_1samp(data, 0.0, alternative="greater")

        # Wilcoxon (one-sided, robust)
        try:
            wstat, pval_w = wilcoxon(data, alternative="greater")
        except ValueError:
            wstat, pval_w = np.nan, np.nan

        # Cohen's d (one-sample)
        cohen_d = mean_z / std_z if std_z > 0 else np.nan

        rows.append({
            "Region": roi,
            "N": n,
            "MEAN_Z": mean_z,
            "STD_Z": std_z,
            "COHEN_D": cohen_d,
            "TSTAT": tstat,
            "PVAL_T": pval_t,
            "WSTAT": wstat,
            "PVAL_W": pval_w
        })

    return pd.DataFrame(rows)

from statsmodels.stats.multitest import multipletests

def add_fdr(df, alpha=0.05):
    """
    Add FDR-corrected p-values to a ROI-level DataFrame.
    """

    reject, pval_fdr, _, _ = multipletests(
        df["PVAL_TWO"],
        alpha=alpha,
        method="fdr_bh"
    )

    df = df.copy()
    df["PVAL_FDR"] = pval_fdr
    df["FDR_SIG"] = reject

    return df

import numpy as np
import pandas as pd
from scipy.stats import ttest_1samp, wilcoxon

def compute_roi_stats(df_wide, min_n=10):
    """
    One-sample stats vs 0 for each ROI:
    - Mean
    - t-test (one-sided)
    - Wilcoxon (one-sided)
    - Cohen's d
    """

    rows = []

    for roi in df_wide.columns:
        data = df_wide[roi].dropna()
        n = len(data)

        if n < min_n:
            continue

        mean_z = data.mean()
        median_z = data.median()
        std_z = data.std(ddof=1)

        # t-test (one-sided)
        tstat, pval_two = ttest_1samp(data, 0.0)  # two-sided por defecto

        pval_one = pval_two / 2

        # Wilcoxon
        try:
           wstat, pval_w = wilcoxon(data, alternative="two-sided")
        except ValueError:
            wstat, pval_w = np.nan, np.nan

        # Cohen's d (one-sample)
        cohen_d = mean_z / std_z if std_z > 0 else np.nan

        rows.append({
            "Region": roi,
            "N": n,
            "MEAN_Z": mean_z,
            "MEDIAN_Z": median_z,
            "STD_Z": std_z,
            "COHEN_D": cohen_d,
            "TSTAT": tstat,
            "PVAL_T": pval_one,
            "WSTAT": wstat,
            "PVAL_W": pval_w
        })

    return pd.DataFrame(rows)

def classify_effect(df):
    """
    Classify ROIs as HIPER or HIPO based on mean z-score.
    """

    df = df.copy()

    df["Effect"] = np.where(
        df["MEAN_Z"] > 0,
        "HIPER",
        "HIPO"
    )

    return df

def add_percentage_extreme_column(
    df_all_covid,
    df_wide,
    threshold=1.96,
    new_col="PERC_1_96"
):
    """
    Adds a column with the percentage of subjects showing extreme z-scores:
    - HIPER  -> z >  threshold
    - HIPO   -> z < -threshold

    Parameters
    ----------
    df_all_covid : DataFrame
        ROI-level stats, must contain columns ['Region', 'Effect']
    df_wide : DataFrame
        Wide dataframe with z-scores (rows = subjects, columns = ROIs)
    threshold : float
        Z-score threshold (default = 1.96)
    new_col : str
        Name of the new column to add

    Returns
    -------
    DataFrame
        df_all_covid with the new percentage column added
    """

    percentages = []

    for _, row in df_all_covid.iterrows():
        roi = row["Region"]
        effect = row["Effect"]

        if roi not in df_wide.columns:
            percentages.append(np.nan)
            continue

        data = df_wide[roi].dropna()

        if len(data) == 0:
            percentages.append(np.nan)
            continue

        if effect == "HIPER":
            perc = (data > threshold).mean() * 100
        elif effect == "HIPO":
            perc = (data < -threshold).mean() * 100
        else:
            perc = np.nan

        percentages.append(perc)

    df_all_covid[new_col] = percentages

    return df_all_covid



def volcano_plot(
    df,
    effect_col="MEDIAN_Z",
    pval_col="PVAL_FDR",
    alpha=0.05,
    sesoi=0.5,
    X_label = "Median z-score",
    title="Volcano plot (Wilcoxon)",
    figsize=(8, 6),

):
    """
    Volcano plot for ROI-level statistics.
    """

    df = df.copy()

    # Evitar problemas con log(0)
    df[pval_col] = df[pval_col].replace(0, np.nan)

    x = df[effect_col]
    y = -np.log10(df[pval_col])

    # Criterio combinado
    # significant = (df[pval_col] < alpha) & (np.abs(df[effect_col]) >= sesoi)
    significant = (df[pval_col] < alpha)


    plt.figure(figsize=figsize)

    # Puntos no relevantes
    plt.scatter(
        x[~significant],
        y[~significant],
        color="lightgray",
        alpha=0.7,
        label="Not selected"
    )

    # Puntos relevantes
    plt.scatter(
        x[significant],
        y[significant],
        color="crimson",
        alpha=0.9,
        label="Significant"
    )

    # Líneas de corte
    plt.axhline(-np.log10(alpha), color="black", linestyle="--", linewidth=1)
    plt.axvline(+sesoi, color="black", linestyle="--", linewidth=1)
    plt.axvline(-sesoi, color="black", linestyle="--", linewidth=1)

    plt.xlabel(X_label)
    plt.ylabel("-log10(p-value)")
    plt.title(title)
    plt.legend(frameon=False)

    plt.tight_layout()
    plt.show()



# main

df_long_covid = load_long_zscores(paths_neurocloud_reports)

# Pivot to wide format: one row per SUBJECT, one column per ROI_NAME and hemisphere
# (e.g. "Hipocampo_LH", "Hipocampo_RH")
df_wide_covid = df_long_covid.pivot_table(
    index="SUBJECT",
    columns=["ROI_NAME", "HEMI"],
    values="ZSCORE",
    aggfunc="first"
)
df_wide_covid.columns = [f"{roi}_{hemi}" for roi, hemi in df_wide_covid.columns]
df_wide_covid = df_wide_covid.reset_index()

# Save the raw data for neurocloud long covid z-scores:
df_wide_covid.to_csv(os.path.join(output_path, "neurocloud_long_covid_zscores.csv"), index=False)

# Do the ttest and wilcoxon test for each ROI, and compute Cohen's d
df_long_covid = compute_roi_stats(df_wide_covid.set_index("SUBJECT"))

df_long_covid = classify_effect(df_long_covid)

# Sort using p-value and Cohen's d
df_long_covid = df_long_covid.sort_values(by=["PVAL_T", "COHEN_D"], ascending=[True, False])

df_long_covid.to_csv(os.path.join(output_path, "ttest_neurocloud_long_covid.csv"), index=False)


# Change the names to the standard hammers names in the nerucloud output
neurocloud_to_hammers_df = pd.read_csv(os.path.join('./data/pet/', 'neurocloud_to_hammers_rois_match.csv'))
rename_dict = dict(zip(neurocloud_to_hammers_df["neurocloud_roi"], neurocloud_to_hammers_df["hammers_roi"]))
df_long_covid["Region"] = df_long_covid["Region"].map(rename_dict).fillna(df_long_covid["Region"])

df_long_covid.to_csv(
    os.path.join(output_path, "ttest_neurocloud_long_covid_hammers_names.csv"),
    index=False
)

volcano_plot(
    df=df_long_covid,
    effect_col="COHEN_D",
    pval_col="PVAL_W",
    alpha=0.05,
    sesoi=0.5,
    title="Volcano plot – Long COVID (Wilcoxon)"
)

volcano_plot(
    df=df_long_covid,
    effect_col="COHEN_D",
    pval_col="PVAL_T",
    alpha=0.05,
    sesoi=0.5,
    X_label = "Cohen D",
    title="Volcano plot – Long COVID (Test T)"
)