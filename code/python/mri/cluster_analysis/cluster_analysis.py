from statsmodels.stats.multicomp import pairwise_tukeyhsd
import statsmodels.api as sm
from statsmodels.formula.api import ols
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from scipy.stats import chi2_contingency

import math, os
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def plot_all_variables(df, vars_to_plot, results_anova=None, save_path=None, title=None):
    """
    Generate a grid of boxplots (one subplot per variable).
    Layout adjusts dynamically based on the number of variables.
    """
    n_vars = len(vars_to_plot)

    # Determine grid size: try to keep it balanced
    ncols = 3 if n_vars > 3 else n_vars
    nrows = math.ceil(n_vars / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols*5, nrows*4), squeeze=False)
    if title:
        fig.suptitle(title, fontsize=16, weight="bold")

    for i, var in enumerate(vars_to_plot):
        row, col = divmod(i, ncols)
        ax = axes[row][col]

        sns.boxplot(x="Cluster", y=var, data=df, palette="Set2",
                    order=["1","2","CONTROL"], ax=ax)
        sns.swarmplot(x="Cluster", y=var, data=df, color="black", alpha=0.6,
                      size=3, order=["1","2","CONTROL"], ax=ax)

        ax.set_title(var, fontsize=12, weight="bold")
        ax.set_xlabel("")
        ax.set_ylabel("")

        if results_anova is not None:
            pval = results_anova.loc[results_anova["Variable"]==var, "p_value_anova"].values
            if len(pval) > 0:
                ax.text(
                    0.05, 0.95,
                    f"ANOVA p = {pval[0]:.3f}",
                    transform=ax.transAxes,
                    ha="left", va="top", fontsize=9,
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", lw=0.5)
                )

    # Remove unused subplots
    for j in range(i+1, nrows*ncols):
        row, col = divmod(j, ncols)
        axes[row][col].axis("off")

    plt.tight_layout(rect=[0, 0, 1, 0.97])

    if save_path:
        plt.savefig(save_path, dpi=300)


def cluster_data(
    data: pd.DataFrame,
    group_filter: str = None,
    free_group: str = None,
    exclude_vars = ["BDNF", "M6a", "NFL", "TotalGM_CBF", "DepresionEscalaUnificada", "AnsiedadEscalaUnificada"],
    centers: int = 2,
    nstart: int = 25,
    random_state: int = 123
):
    """
    Perform K-means clustering on a dataset of subjects, with flexibility for:
      1. Clustering the entire dataset.
      2. Clustering only one group (e.g., "CONTROL" or "COVID").
      3. Clustering everything, but forcing one group to remain fixed (`free_group`).

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame containing:
        - "ID": subject identifier.
        - "Grupo": group label ("CONTROL" / "COVID").
        - Numeric variables for clustering.

    group_filter : str, optional
        If provided ("CONTROL" or "COVID"), clustering is performed only on that group.
        If None, clustering includes all subjects.

    free_group : str, optional
        If provided ("CONTROL" or "COVID"), clustering is run on the whole dataset,
        but for subjects in this group the cluster assignment is overridden to be
        exactly `free_group`.

    exclude_vars : list, optional
        Variables to be excluded from clustering (default = protein markers + CBF).

    centers : int, optional
        Number of clusters (default = 2).

    nstart : int, optional
        Number of random initializations for K-means (default = 25).

    random_state : int, optional
        Random seed for reproducibility (default = 123).

    Returns
    -------
    dict with:
        - "merged_data": pd.DataFrame with all subjects and a new "Cluster" column.
                         Subjects with missing values remain with Cluster = NaN.
                         If `free_group` is used, that group is forced to its label.
        - "km_result": sklearn KMeans object (centers, inertia, labels, etc.).
        - "clustering_data_clean": DataFrame with the variables actually used
                                   for clustering (after exclusions and NA removal).
    """
    # 1. Filter group if requested
    if group_filter is not None:
        data_subset = data[data["Grupo"] == group_filter].copy()
    else:
        data_subset = data.copy()

    # 2. Prepare clustering matrix (set ID as index)
    clustering_data = data_subset.drop(columns=["Grupo"]).set_index("ID")

    # Keep only numeric columns
    numeric_cols = clustering_data.select_dtypes(include=["number"]).columns
    clustering_data_numeric = clustering_data[numeric_cols]

    # Scale numeric features
    clustering_data_scaled = StandardScaler().fit_transform(clustering_data_numeric)
    clustering_data_scaled = pd.DataFrame(
        clustering_data_scaled,
        columns=numeric_cols,
        index=clustering_data.index  # IDs preserved
    )

    # 3. Drop excluded variables
    clustering_data_clean = clustering_data_scaled.drop(
        columns=[c for c in exclude_vars if c in clustering_data_scaled.columns],
        errors="ignore"
    ).dropna(axis=0)

    # 4. Run KMeans
    clean_ids = clustering_data_clean.index
    km = KMeans(n_clusters=centers, n_init=nstart, random_state=random_state)
    km.fit(clustering_data_clean)

    # 5. Merge cluster labels back
    cluster_df = pd.DataFrame({"ID": clean_ids, "Cluster": (km.labels_ + 1).astype(str)})
    merged_data = pd.merge(data_subset, cluster_df, on="ID", how="left")

    # 6. Apply free_group override
    if free_group is not None:
        merged_data.loc[merged_data["Grupo"] == free_group, "Cluster"] = free_group

    return {
        "merged_data": merged_data,
        "km_result": km,
        "clustering_data_clean": clustering_data_clean
    }


def run_anova(df, var):
    """
    Corre un ANOVA de una variable continua (var) respecto a Cluster.
    Maneja los NaN SOLO de esa variable.
    """
    df = df.copy()
    df["Cluster"] = df["Cluster"].astype(str)

    # Usar solo filas válidas para la variable de interés
    df_var = df.dropna(subset=[var])

    # Ajustar modelo
    model = ols(f"{var} ~ C(Cluster)", data=df_var).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)

    print(f"\n=== ANOVA {var} ~ Cluster ===")
    print(anova_table)

    return anova_table, model


def run_chisq_gender(df):
    df = df.copy()
    df["Cluster"] = df["Cluster"].astype(str)

    contingency = pd.crosstab(df["Cluster"], df["Genero"])
    chi2, p, dof, expected = chi2_contingency(contingency)

    print("\n=== Chi-cuadrado Género ~ Cluster ===")
    print("Tabla de contingencia:")
    print(contingency, "\n")
    print(f"Chi2 = {chi2:.3f}")
    print(f"p-value = {p:.4f}")
    print(f"Grados de libertad = {dof}\n")
    print("Valores esperados (bajo independencia):")
    print(pd.DataFrame(expected,
                       index=contingency.index,
                       columns=contingency.columns))

    return contingency, chi2, p, dof, expected


def plot_age_gender_bmi(df, save_path=None):
    df = df.copy()
    df["Cluster"] = df["Cluster"].astype(str)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # --- Boxplot Edad ---
    sns.boxplot(x="Cluster", y="Edad", data=df,
                palette="Set2", order=["1", "2", "CONTROL"], ax=axes[0])
    sns.swarmplot(x="Cluster", y="Edad", data=df,
                  color="black", alpha=0.6, size=3,
                  order=["1", "2", "CONTROL"], ax=axes[0])
    axes[0].set_title("Edad por Cluster", fontsize=14, weight="bold")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("Edad")

    # --- Barplot % Género ---
    sns.barplot(x="Cluster", y="Genero_Perc", hue="Genero", data=df,
                palette="Set2", order=["1", "2", "CONTROL"], ax=axes[1])
    axes[1].set_title("% Género por Cluster", fontsize=14, weight="bold")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("%")

    # --- Boxplot BMI ---
    sns.boxplot(x="Cluster", y="BMI", data=df,
                palette="Set2", order=["1", "2", "CONTROL"], ax=axes[2])
    sns.swarmplot(x="Cluster", y="BMI", data=df,
                  color="black", alpha=0.6, size=3,
                  order=["1", "2", "CONTROL"], ax=axes[2])
    axes[2].set_title("BMI por Cluster", fontsize=14, weight="bold")
    axes[2].set_xlabel("")
    axes[2].set_ylabel("BMI")

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"✅ Gráfico guardado en {save_path}")





# === Paths ===
path_data_analysis = os.path.abspath("../../../../data/")
path_results_data_analysis = os.path.join(path_data_analysis, "results")
path_plots_data_analysis = os.path.join(path_data_analysis, "plots/")
path_cuestionarios = os.path.join(path_data_analysis, "study_summary.xlsx")

if not os.path.exists(path_results_data_analysis):
    os.makedirs(path_results_data_analysis)
if not os.path.exists(path_plots_data_analysis):
    os.makedirs(path_plots_data_analysis)

# Input
matrix_path = os.path.join(path_data_analysis, "kmeans_analysis/matrix_kmeans.csv")

# Outputs
cluster_matrix = os.path.join(path_results_data_analysis, "cluster_data.csv")
anova_csv = os.path.join(path_results_data_analysis, "anova_results.csv")
tukey_csv = os.path.join(path_results_data_analysis, "tukey_results.csv")


# === Read input matrix ===
matrix_kmeans = pd.read_csv(matrix_path)


# === Define variable categories ===
volumes_vars = [
    "frontalpole_vol_rh","superiorfrontal_vol_rh","rostralmiddlefrontal_vol_rh",
    "parsopercularis_vol_rh","fusiform_vol_rh","inferiortemporal_vol_rh",
    "lingual_vol_lh","Left_Cerebellum_Cortex","Right_Cerebellum_Cortex"
]

thickness_vars = [
    "fusiform_thick_rh","postcentral_thick_rh","lingual_thick_lh",
    "postcentral_thick_lh","supramarginal_thick_lh"
]

wmh_vars = ["WMHTotal","WMHPeriventricular","WMHDeep"]
asl_vars = ["TotalGM_sCOV","TotalGM_CBF"]
scales_proteins_vars = ["DepresionEscalaUnificada","AnsiedadEscalaUnificada","NFL","BDNF","M6a"]
tests_vars = ["MOCA_perc","TMT_A_perc"]


# === Run clustering ===
res_all = cluster_data(
    data=matrix_kmeans,
    free_group=None,   # cluster all subjects
    centers=2,
)

# Extract merged DataFrame
final_matrix_all = res_all["merged_data"]

# Reorder columns so ID, Grupo, Cluster come first
cols = ["ID", "Grupo", "Cluster"] + [c for c in final_matrix_all.columns if c not in ["ID", "Grupo", "Cluster"]]
final_matrix_all = final_matrix_all[cols]


# 1. Directly extract controls
final_matrix_control = matrix_kmeans[matrix_kmeans["Grupo"] == "CONTROL"].copy()
final_matrix_control["Cluster"] = "CONTROL"

# 2. Cluster only COVID subjects
final_matrix_covid = cluster_data(
    matrix_kmeans,
    group_filter="COVID",
    centers=2
)["merged_data"]

# 3. Concatenate both
df_clusters = pd.concat([final_matrix_control, final_matrix_covid], ignore_index=True)

# 4. Clean cluster column: remove NaN and convert to string
df_clusters = df_clusters[df_clusters["Cluster"].notna()].copy()
df_clusters["Cluster"] = df_clusters["Cluster"].astype(str)

# Optionally flip labels
df_clusters["Cluster"] = df_clusters["Cluster"].replace({"1": "2", "2": "1"})

# Select numeric variables
numeric_vars = df_clusters.select_dtypes(include="number").columns.difference(["ID"])
numeric_vars = [
    c for c in df_clusters.columns
    if c not in ["ID", "Grupo", "Cluster"]
    and pd.api.types.is_numeric_dtype(df_clusters[c])
]

# Save CSV with ID, Group, and Cluster
clusters_ids = df_clusters[["ID", "Grupo", "Cluster"]]
clusters_ids.to_csv(os.path.join(path_results_data_analysis, "clusters_ids.csv"), index=False)

# === Run ANOVA and Tukey ===
results_anova = []
results_tukey = {}

for var in numeric_vars:
    mask = df_clusters[var].notna() & df_clusters["Cluster"].notna()

    # ANOVA
    model = ols(f"Q('{var}') ~ C(Cluster)", data=df_clusters[mask]).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)
    pval_anova = anova_table["PR(>F)"][0]
    results_anova.append({"Variable": var, "p_value_anova": pval_anova})

    # Tukey only if ANOVA < 0.1
    if pval_anova < 0.1:
        tukey = pairwise_tukeyhsd(
            endog=df_clusters.loc[mask, var].astype(float),
            groups=df_clusters.loc[mask, "Cluster"].astype(str),
            alpha=0.05
        )
        results_tukey[var] = pd.DataFrame(
            data=tukey.summary().data[1:],
            columns=tukey.summary().data[0]
        )

results_anova = pd.DataFrame(results_anova)
results_anova.to_csv(anova_csv, index=False)

# Concatenate Tukey results into one DataFrame
tukey_all = []
for var, df_tukey in results_tukey.items():
    df_temp = df_tukey.copy()
    df_temp["Variable"] = var
    tukey_all.append(df_temp)
tukey_all = pd.concat(tukey_all, ignore_index=True)
tukey_all.to_csv(tukey_csv, index=False)


# === Save plots by category ===
plot_all_variables(df_clusters, volumes_vars, results_anova,
                   save_path=os.path.join(path_plots_data_analysis, "volumes.png"),
                   title="Volumes")

plot_all_variables(df_clusters, thickness_vars, results_anova,
                   save_path=os.path.join(path_plots_data_analysis, "thickness.png"),
                   title="Cortical Thickness")

plot_all_variables(df_clusters, wmh_vars, results_anova,
                   save_path=os.path.join(path_plots_data_analysis, "wmh.png"),
                   title="White Matter Hyperintensities (WMH)")

plot_all_variables(df_clusters, asl_vars, results_anova,
                   save_path=os.path.join(path_plots_data_analysis, "asl.png"),
                   title="ASL Measures: sCOV and CBF")

plot_all_variables(df_clusters, scales_proteins_vars, results_anova,
                   save_path=os.path.join(path_plots_data_analysis, "scales_proteins.png"),
                   title="Clinical Scales and Proteins")

plot_all_variables(df_clusters, tests_vars, results_anova,
                   save_path=os.path.join(path_plots_data_analysis, "tests.png"),
                   title="Cognitive Tests")

# Save final clustered matrix
df_clusters.to_csv(cluster_matrix, index=False)

## Demographics

# === Merge con cuestionarios ===
df_cuestionarios = pd.read_csv(path_cuestionarios)

df_merged = pd.merge(df_cuestionarios, df_clusters, on="ID", how="inner")

# Calcular % de género por cluster
gender_perc = (
    df_merged.groupby(["Cluster", "Genero"])
             .size()
             .reset_index(name="count")
)
gender_perc["Genero_Perc"] = gender_perc.groupby("Cluster")["count"].transform(
    lambda x: 100 * x / x.sum()
)

df_merged = df_merged.merge(
    gender_perc[["Cluster", "Genero", "Genero_Perc"]],
    on=["Cluster", "Genero"], how="left"
)

# ANOVAs
anova_edad, model_edad = run_anova(df_merged, "Edad")
anova_bmi, model_bmi   = run_anova(df_merged, "BMI")

# Chi-square para Género
contingency, chi2, p, dof, expected = run_chisq_gender(df_merged)

# Caso especial: setear BMI en NaN para CP0171
df_merged.loc[df_merged["ID"] == "CP0171", "BMI"] = np.nan

# Plot Edad + Género + BMI
plot_age_gender_bmi(
    df_merged,
    save_path=os.path.join(path_plots_data_analysis, "edad_genero_bmi.png")
)









