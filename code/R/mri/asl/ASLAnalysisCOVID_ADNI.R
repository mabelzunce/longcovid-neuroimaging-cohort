# Load necessary libraries
library(tidyverse)
library(readxl)
library(car)
library(gridExtra)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(effects)
library(ggplot2)
library(patchwork)

source("functionsANCOVA.R")

# Data Analysis 
path_data_analysis <- "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Freesurfer/DataAnalysis/ASL/"
path_results_data_analysis <- paste(path_data_analysis, "Results/", sep="")
path_plots_data_analysis <-  paste(path_data_analysis, "Plots/", sep="")

# COVID
path_asl_analysis_covid <- "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/ASLProcessedCOVID/derivatives/ExploreASL/Population/Stats/"
path_cuestionarios <- paste(path_data_analysis, "RespuestasCuestionariosCovid_15_01.csv", sep="")
raw_perfussion_data_covid <-"mean_qCBF_StandardSpace_TotalGM_n=203_15-Jan-2025_PVC2.tsv"
path_dkt_perfusion_covid <-paste(path_asl_analysis_covid, raw_perfussion_data_covid , sep="")

# ADNI
path_asl_analysis_adni <- "/media/sol/Expansion/ADNI/ADNI_ExploreASL_3D/derivatives/ExploreASL/Population/Stats/"
path_info_images_adni <- "/media/sol/Expansion/ADNI/all_processed_images_info.csv"
raw_perfussion_data_adni <-"mean_qCBF_StandardSpace_TotalGM_n=198_28-Jan-2025_PVC2.tsv"
path_dkt_perfusion_adni <-paste(path_asl_analysis_adni, raw_perfussion_data_adni , sep="")


# COVID
# Values to exclude 
values_to_exclude_covid <- c('CP0011', 'CP0015','CP0009', 'CP0035', 'CP0093', 'CP0101', 'CP0140',
                       'CP0035', 'CP0087', 'CP0105', 'CP0117', 'CP0144', 'CP0192','CP0216', 'CP0227', 'CP0062')

vascular_artifacts_covid <- c('CP0031','CP0044', 'CP0053', 'CP0054', 'CP0061', 'CP0096','CP0107','CP0108','CP0118',
                        'CP0123', 'CP0142', 'CP0147', 'CP0154','CP0167', 'CP0176','CP0178',
                        'CP0180', 'CP0183', 'CP0188','CP0193','CP0199', 'CP0205', 'CP0213',
                        'CP0229', 'CP0238','CP0245', 'CP0247', 'CP0238')

values_to_exclude_covid <- paste0("sub-", values_to_exclude_covid, "_1")
vascular_artifacts_covid <- paste0("sub-", vascular_artifacts_covid, "_1")

# Read files
cuestionarios_excel <- read.csv(path_cuestionarios)
cuestionarios_excel$ID <- paste0("sub-", cuestionarios_excel$ID,"_1")
perfusion_csv_covid <- read.delim(path_dkt_perfusion_covid)

# Exclude columns with +10 NaN in Perfusion
perfusion_csv_covid[perfusion_csv_covid == "n/a"] <- NA
cols_with_nan_perfusion <- colSums(is.na(perfusion_csv_covid)) > 10
columns_to_exclude_perfusion <- names(perfusion_csv_covid)[cols_with_nan_perfusion]
perfusion_csv_covid <- perfusion_csv_covid[, !cols_with_nan_perfusion]
colnames(perfusion_csv_covid) <- make.names(colnames(perfusion_csv_covid), unique = TRUE)

# Columnas con regiones de la segmentación 
column_names_perfusion <- names(perfusion_csv_covid)
columns_regions_perfusion <- column_names_perfusion[6:length(column_names_perfusion)]
columns_regions_perfusion <- columns_regions_perfusion[!grepl("^Wm.", columns_regions_perfusion)]

# Unir con variables necesarias de cuestionarios
cols_to_merge_covid <- c("Edad", "Grupo", "ID", "Genero", "BMI", "FechaCognitivo")
cols_to_merge_covid <- intersect(cols_to_merge_covid, colnames(cuestionarios_excel))  

merged_df_perfusion_covid <- merge(
  cuestionarios_excel[, cols_to_merge_covid], 
  perfusion_csv_covid, 
  by.x = "ID", 
  by.y = "participant_id", 
  all.y = TRUE
) %>% 
  slice(-1) %>%  # Eliminar la primera fila
  mutate(across(all_of(columns_regions_perfusion), ~ suppressWarnings(as.numeric(.))))  # Regiones %>%

merged_df_perfusion_with_ata_covid <- subset(merged_df_perfusion_covid, !(ID %in% values_to_exclude_covid))
merged_df_perfusion_covid <- subset(merged_df_perfusion_with_ata_covid, !(ID %in% vascular_artifacts_covid)) 


# ADNI Data 
# Info images 
info_images_adni <-read.csv(path_info_images_adni)

# Unir con variables necesarias de la información adicional
cols_to_merge_adni <- c("Age", "Group", "Sex", "Subject", "Site")
cols_to_merge_adni <- intersect(cols_to_merge_adni, colnames(info_images_adni))  

# Limpiar tsv perfusión
perfusion_csv_adni <- read.delim(path_dkt_perfusion_adni)
# perfusion_csv_adni <- perfusion_csv_adni[, !cols_with_nan_perfusion]
colnames(perfusion_csv_adni) <- make.names(colnames(perfusion_csv_adni), unique = TRUE)

# Values to exclude 
# ASL images
asl_images_path <- "/media/sol/Expansion/ADNI/ADNI_ExploreASL_3D/derivatives/ExploreASL/Population/ASLCheck/"

artifact_images_path <- paste(asl_images_path, "3_ArtifactContrast", sep="")
vascular_images_path <- paste(asl_images_path, "2_VascularContrast", sep="")
processed_images_path <- paste(asl_images_path, "1_CBFContrast", sep="")


artifact_images <- list.files(path = artifact_images_path, pattern = "\\.jpg$", full.names = FALSE)
vascular_images <- list.files(path = vascular_images_path, pattern = "\\.jpg$", full.names = FALSE)
processed_images <- list.files(path = processed_images_path, pattern = "\\.jpg$", full.names = FALSE)

patron <- "(sub[^_]*_1)"

artifact_ids <- str_extract(artifact_images, patron)
vascular_ids <- str_extract(vascular_images, patron)
processed_ids <- str_extract(processed_images, patron)

merged_df_perfusion_adni <- merge(
  info_images_adni[, cols_to_merge_adni], 
  perfusion_csv_adni, 
  by.x = "Subject", 
  by.y = "participant_id", 
  all.y = TRUE
) %>% 
  slice(-1)  # Eliminar la primera fila


merged_df_perfusion_adni_vascular <- merged_df_perfusion_adni %>%
  filter(Subject %in% vascular_ids)

merged_df_perfusion_adni_artifact <- merged_df_perfusion_adni %>%
  filter(Subject %in% artifact_ids)

merged_df_perfusion_adni_processed <- merged_df_perfusion_adni %>%
  filter(Subject %in% processed_ids)

count_vascular <- merged_df_perfusion_adni_vascular %>%
  group_by(Group) %>%
  summarise(Conteo = n())

count_artifact <- merged_df_perfusion_adni_artifact %>%
  group_by(Group) %>%
  summarise(Conteo = n())

count_processed <- merged_df_perfusion_adni_processed %>%
  group_by(Group) %>%
  summarise(Conteo = n())
# 
# print("Vasculares")
# print(count_vascular)
# 
# print("Artefactos")
# print(count_artifact)
# 
# print("Procesadas")
# print(count_processed)
# 
merged_df_perfusion_adni_processed <- merged_df_perfusion_adni_processed %>%
  mutate(
    TotalGM_L = as.numeric(TotalGM_L)
  )

merged_df_perfusion_adni_processed <- merged_df_perfusion_adni_processed %>%
  mutate(Group = ifelse(Group == "CN", "CN", "CI")) %>%
  mutate(Group = factor(Group, levels = c("CN", "CI")))


boxplots_paper <- list()
#
#
boxplots_paper[[1]] <- ggplot(merged_df_perfusion_covid, aes(x = Grupo, y =  TotalGM_L , fill = Grupo)) +
  geom_boxplot(alpha = 0.5) +  # Add transparency to boxplot
  geom_jitter(color = "black", size = 2, alpha = 0.3, width = 0.2) +  # Add jittered points
  labs(title = "PL_supramarginal_gyrus_l",
       x = "Group",
       y = "CBF") +
  theme_minimal() +
  guides(fill = FALSE) +# Remove legend
  theme(
    plot.title = element_text(size = 10, hjust=0.5),  # Cambia el tamaño del título
  ) +
  coord_cartesian(ylim = c(30, 50))  # Establecer escala fija en el eje Y

boxplots_paper[[2]] <- ggplot(merged_df_perfusion_adni_processed, aes(x = Group, y =  TotalGM_L , fill = Group)) +
  geom_boxplot(alpha = 0.5) +  # Add transparency to boxplot
  geom_jitter(color = "black", size = 2, alpha = 0.3, width = 0.2) +  # Add jittered points
  labs(title = "PL_supramarginal_gyrus_l",
       x = "Group",
       y = "CBF") +
  theme_minimal() +
  guides(fill = FALSE) +# Remove legend
  theme(
    plot.title = element_text(size = 10, hjust=0.5)  # Cambia el tamaño del título
  ) +
  coord_cartesian(ylim = c(30, 50))  # Establecer escala fija en el eje Y

# Combinar los dos boxplots en una sola figura
combined_boxplots <- boxplots_paper[[1]] + boxplots_paper[[2]] +
  plot_layout(ncol = 2) +  # Organizar en 2 columnas
  plot_annotation(title = "PL_supramarginal_gyrus_L",
                  subtitle = "COVID Vs ADNI",
                  theme = theme(plot.title = element_text(size = 14, hjust = 0.5),
                                plot.subtitle = element_text(size = 12, hjust = 0.5)))

print(combined_boxplots)
