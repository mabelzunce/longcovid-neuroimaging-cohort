import pandas as pd
import numpy as np
import re

# Paths
path_data_analysis = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Freesurfer/DataAnalysis/"
path_results_data_analysis = path_data_analysis + "Results/"
path_plots_data_analysis = path_data_analysis + "Plots/"

path_cuestionarios = path_data_analysis + "ResumenTotal_2_06.csv"
path_sienax2 = path_data_analysis + "SienaxResults2.csv"
path_brain_vol = path_data_analysis + "brain_volumes.csv"
path_segmentation = path_data_analysis + "segmentation.csv"
path_parcellation = path_data_analysis + "parcellation.csv"
path_parcellation_DKT = path_data_analysis + "parcellation_DKT.csv"
path_parcellation_thick = path_data_analysis + "parcellation_thick.csv"
path_parcellation_thick_dkt = path_data_analysis + "parcellation_thick_DKT.csv"

path_bianca_file = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/WMH/ProcessedWMH/TotalPvAndDeepValuesBiancaThr0_7.csv"
path_socioeconomics_cardiovascular_index = path_data_analysis + "ASL/SocialAndCardiovascularIndex.csv"
path_gray_matter_perfusion = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/ASLProcessedCOVID/derivatives/ExploreASL/Population/Stats/mean_qCBF_StandardSpace_TotalGM_n=203_15-Jan-2025_PVC2.tsv"
path_scov_without_pvc = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/ASLProcessedCOVID/derivatives/ExploreASL/Population/Stats/CoV_qCBF_StandardSpace_TotalGM_n=203_15-Jan-2025_PVC0.tsv"

# IDs a excluir
values_to_exclude_from_the_study = ['CP0011', 'CP0015', 'CP0106', 'CP0087', 'CP0144', 'CP0193']

values_to_exclude_mri = ['CP0035', 'CP0192', 'CP0196']
# Concatenate lists
values_to_exclude = list(set(values_to_exclude_from_the_study) | set(values_to_exclude_mri))

# === Leer Cuestionarios ===
cuestionarios_excel = pd.read_csv(path_cuestionarios)

# Excluir filas
cuestionarios_excel = cuestionarios_excel[~cuestionarios_excel["ID"].isin(values_to_exclude)]

# === Leer otros archivos ===
sienax2 = pd.read_csv(path_sienax2, decimal=",")
brain_vol = pd.read_csv(path_brain_vol)
segmentation = pd.read_csv(path_segmentation)
parcellation = pd.read_csv(path_parcellation)
parcellation_DKT = pd.read_csv(path_parcellation_DKT)
parcellation_thick = pd.read_csv(path_parcellation_thick)
parcellation_thick_DKT = pd.read_csv(path_parcellation_thick_dkt)
socioeconomics_cardiovascular_index = pd.read_csv(path_socioeconomics_cardiovascular_index)


# === Leer WMH ===
bianca = pd.read_csv(path_bianca_file, names=["ID","Total_WMH","PV_WMH","Deep_WMH"])


# === Lista de IDs con artefactos vasculares ===
vascular_artifacts = [
    'CP0031','CP0044','CP0053','CP0054','CP0061','CP0096','CP0107','CP0108','CP0118',
    'CP0123','CP0142','CP0147','CP0154','CP0167','CP0176','CP0178','CP0180','CP0183',
    'CP0188','CP0193','CP0199','CP0205','CP0213','CP0229','CP0238','CP0245','CP0247','CP0238'
]

# === CBF (PVC2) ===

perfusion_csv = pd.read_csv(path_gray_matter_perfusion, sep="\t")

# Eliminar primera fila como en R [-1, ]
perfusion_csv = perfusion_csv.iloc[1:, :]

# Limpiar participant_id
perfusion_csv["participant_id"] = (
    perfusion_csv["participant_id"]
    .str.replace("^sub-", "", regex=True)
    .str.replace("_1$", "", regex=True)
)

# Reemplazar "n/a" por NaN
perfusion_csv = perfusion_csv.replace("n/a", np.nan)

# Convertir TotalGM_B a numérico
perfusion_csv["TotalGM_B"] = pd.to_numeric(perfusion_csv["TotalGM_B"], errors="coerce")

# Setear NA en artefactos
perfusion_csv.loc[perfusion_csv["participant_id"].isin(vascular_artifacts), "TotalGM_B"] = np.nan


# === sCOV (PVC0) ===

perfusion_scov_without_pvc = pd.read_csv(path_scov_without_pvc, sep="\t")
perfusion_scov_without_pvc = perfusion_scov_without_pvc.iloc[1:, :]

# Limpiar participant_id
perfusion_scov_without_pvc["participant_id"] = (
    perfusion_scov_without_pvc["participant_id"]
    .str.replace("^sub-", "", regex=True)
    .str.replace("_1$", "", regex=True)
)

# Reemplazar "n/a" por NaN
perfusion_scov_without_pvc = perfusion_scov_without_pvc.replace("n/a", np.nan)

# Convertir TotalGM_B a numérico
perfusion_scov_without_pvc["TotalGM_B"] = pd.to_numeric(perfusion_scov_without_pvc["TotalGM_B"], errors="coerce")


# === Regiones Freesurfer de interés ===
regions_segmentation = ["Left-Cerebellum-Cortex", "Right-Cerebellum-Cortex"]

right_regions_parcellation = [
    "frontalpole", "superiorfrontal", "rostralmiddlefrontal",
    "parsopercularis", "fusiform", "inferiortemporal"
]
left_regions_parcellation = ["lingual"]

right_regions_parcellation_thick = ["fusiform", "postcentral"]
left_regions_parcellation_thick = ["lingual", "postcentral", "supramarginal"]

# === Variables de cuestionarios ===
proteins_columns_select = ["M6a", "BDNF", "NFL"]
test_select = ["MOCA_perc", "TMT-A_perc", "DepresionEscalaUnificada", "AnsiedadEscalaUnificada"]
WMHvolumes = ["WMHTotal", "WMHPeriventricular", "WMHDeep"]
#
# --- 1. Normalización de columnas numéricas ---
for col in ["M6a", "BDNF", "NFL", "MOCA_perc", "IndiceSocioeconomico"]:
    cuestionarios_excel[col] = (
        cuestionarios_excel[col]
        .astype(str)
        .str.replace(",", ".", regex=False)
        .astype(float)
    )

cols_matrix = ['ID'] + ['Grupo'] + proteins_columns_select + test_select + WMHvolumes
matrix_kmeans = cuestionarios_excel[cols_matrix].copy()

# --- 4. Merge con sCOV ---
scov = perfusion_scov_without_pvc[['participant_id', 'TotalGM_B']].rename(columns={
    'participant_id': 'ID',
    'TotalGM_B': 'TotalGM_sCOV'
})
matrix_kmeans = pd.merge(matrix_kmeans, scov, on="ID", how="left")

# --- 5. Merge con CBF ---
cbf = perfusion_csv[['participant_id', 'TotalGM_B']].rename(columns={
    'participant_id': 'ID',
    'TotalGM_B': 'TotalGM_CBF'
})
matrix_kmeans = pd.merge(matrix_kmeans, cbf, on="ID", how="left")

# --- 6. Merge con Cardiovascular Index ---
matrix_kmeans = pd.merge(
    matrix_kmeans,
    socioeconomics_cardiovascular_index[['ID', 'CardiovascularRiskIndex']],
    on="ID",
    how="left"
)

# --- 7. Preparar datos FreeSurfer ---
# Right Hemisphere Volume
right_hemisphere_vol_data = (
    parcellation_DKT.query("Hemisphere == 'rh'")
    .loc[:, ['subject'] + right_regions_parcellation]
    .rename(columns={col: f"{col}_vol_rh" for col in right_regions_parcellation})
)

# Left Hemisphere Volume
left_hemisphere_vol_data = (
    parcellation_DKT.query("Hemisphere == 'lh'")
    .loc[:, ['subject'] + left_regions_parcellation]
    .rename(columns={col: f"{col}_vol_lh" for col in left_regions_parcellation})
)

# Right Hemisphere Thickness
right_hemisphere_thick_data = (
    parcellation_thick_DKT.query("Hemisphere == 'rh'")
    .loc[:, ['subject'] + right_regions_parcellation_thick]
    .rename(columns={col: f"{col}_thick_rh" for col in right_regions_parcellation_thick})
)

# Left Hemisphere Thickness
left_hemisphere_thick_data = (
    parcellation_thick_DKT.query("Hemisphere == 'lh'")
    .loc[:, ['subject'] + left_regions_parcellation_thick]
    .rename(columns={col: f"{col}_thick_lh" for col in left_regions_parcellation_thick})
)

# Cerebellum
cerebellum_data = segmentation[['subject'] + regions_segmentation]

# --- 8. Unir FreeSurfer ---
all_freesurfer = (
    right_hemisphere_vol_data
    .merge(left_hemisphere_vol_data, on="subject", how="outer")
    .merge(right_hemisphere_thick_data, on="subject", how="outer")
    .merge(left_hemisphere_thick_data, on="subject", how="outer")
    .merge(cerebellum_data, on="subject", how="outer")
)
# --- 9. Merge con matrix_kmeans ---
matrix_kmeans = pd.merge(matrix_kmeans, all_freesurfer, left_on="ID", right_on="subject", how="left")

# Eliminar la columna 'subject' porque ya es redundante con 'ID'
matrix_kmeans = matrix_kmeans.drop(columns=["subject"])

matrix_kmeans = matrix_kmeans.rename(
    columns=lambda x: re.sub(r'[^0-9a-zA-Z]+', '_', x)
)

# WMH
# CP0209
matrix_kmeans.loc[
    matrix_kmeans["ID"] == "CP0209",
    ["WMHTotal", "WMHPeriventricular", "WMHDeep"]
] = [83, 26, 57]

# CP0219
matrix_kmeans.loc[
    matrix_kmeans["ID"] == "CP0219",
    ["WMHTotal", "WMHPeriventricular", "WMHDeep"]
] = [708, 489, 219]

# Forzar a NaN en M6a para IDs específicos
ids_to_nan = ["CP0223","CP0224","CP0229","CP0230","CP0235","CP0156"]
ids_to_nan_BDNF = ["CP0002", "CP0223"]

matrix_kmeans.loc[matrix_kmeans["ID"].isin(ids_to_nan), "M6a"] = np.nan
matrix_kmeans.loc[matrix_kmeans["ID"].isin(ids_to_nan_BDNF), "BDNF"] = np.nan


# IDs para forzar a NaN en TotalGM_sCOV
ids_to_nan_scov = ['CP0009','CP0093','CP0101','CP0140','CP0105','CP0117','CP0216','CP0227','CP0062']

matrix_kmeans.loc[matrix_kmeans["ID"].isin(ids_to_nan_scov), "TotalGM_sCOV"] = np.nan


# Definir path donde guardar
output_path_matrix = f"{path_results_data_analysis}/matrix_kmeans.csv"
matrix_kmeans.to_csv(output_path_matrix, index=False)
