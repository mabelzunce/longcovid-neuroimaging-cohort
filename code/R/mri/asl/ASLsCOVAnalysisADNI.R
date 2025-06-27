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


# ADNI
path_asl_analysis_adni <- "/media/sol/Expansion/ADNI/ADNI_ExploreASL_3D/derivatives/ExploreASL/Population/Stats/"
path_info_images_adni <- "/media/sol/Expansion/ADNI/all_processed_images_info.csv"
raw_perfussion_data_adni <-"CoV_qCBF_StandardSpace_TotalGM_n=198_28-Jan-2025_PVC0.tsv"
path_dkt_perfusion_adni <-paste(path_asl_analysis_adni, raw_perfussion_data_adni , sep="")


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

# Linear Regression
# Ensure Genero is a factor
merged_df_perfusion_adni$Sex <- as.factor(merged_df_perfusion_adni$Sex)

filtered_df_adni <- merged_df_perfusion_adni %>% filter(!is.na(TotalGM_B))
# 
# # Define dependent variabl
filtered_df_adni$TotalGM_B <- as.numeric(filtered_df_adni$TotalGM_B)
filtered_df_adni$TotalGM_B <- filtered_df_adni$TotalGM_B * 100

filtered_df_adni <- filtered_df_adni %>% filter(!is.na(Group))


# Plot 3: Scatter plot with regression line (Age vs TotalGM_B, colored by Grupo)
p3 <- ggplot(filtered_df_adni, aes(x = Age, y = TotalGM_B, color = Group)) +
  geom_point(alpha = 0.8, size = 3) +  # Adjust transparency and size
  scale_color_viridis_d() +
  scale_shape_manual(values = c(15, 17)) +  # Square (15) and Triangle (17)
  labs(title = "Age vs sCOVgm (without PVC)",
       x = "Age", y = "sCOV(%)", shape = "Vascular Artifact", color = "Grupo") +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"),
        legend.position = "top")  # Move legend to the top for better readability

ggsave(filename = paste0(path_plots_data_analysis, "scatter_Age_vs_sCOV_without_PVC_adni.png"), 
       plot = p3, width = 8, height = 6, dpi = 300, bg = "white")

filtered_df_adni_2 <- filtered_df_adni %>% filter(TotalGM_B != 0)

# Plot 3: Scatter plot with regression line (Age vs TotalGM_B, colored by Grupo)
p5 <- ggplot(filtered_df_adni_2, aes(x = Age, y = TotalGM_B, color = Group)) +
  geom_point(alpha = 0.8, size = 3) +  # Adjust transparency and size
  scale_color_viridis_d() +
  scale_shape_manual(values = c(15, 17)) +  # Square (15) and Triangle (17)
  labs(title = "Age vs sCOVgm (without PVC)",
       x = "Age", y = "sCOV(%)", shape = "Vascular Artifact", color = "Grupo") +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"),
        legend.position = "top")  # Move legend to the top for better readability

ggsave(filename = paste0(path_plots_data_analysis, "scatter_Age_vs_sCOV_without_PVC_adni_2.png"), 
       plot = p5, width = 8, height = 6, dpi = 300, bg = "white")

filtered_df_adni_3 <- filtered_df_adni_2 %>% filter(Subject != "sub-941S4365_1")

# Plot 3: Scatter plot with regression line (Age vs TotalGM_B, colored by Grupo)
p6 <- ggplot(filtered_df_adni_3, aes(x = Age, y = TotalGM_B, color = Group)) +
  geom_point(alpha = 0.8, size = 3) +  # Adjust transparency and size
  scale_color_viridis_d() +
  scale_shape_manual(values = c(15, 17)) +  # Square (15) and Triangle (17)
  labs(title = "Age vs sCOVgm (without PVC)",
       x = "Age", y = "sCOV(%)", shape = "Vascular Artifact", color = "Grupo") +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"),
        legend.position = "top")  # Move legend to the top for better readability

ggsave(filename = paste0(path_plots_data_analysis, "scatter_Age_vs_sCOV_without_PVC_adni_3.png"), 
       plot = p6, width = 8, height = 6, dpi = 300, bg = "white")


# Crear una nueva variable de grupo basada en los artefactos
filtered_df_adni_4 <- filtered_df_adni_3 %>%
  mutate(
    Artifact_Group = case_when(
      Subject %in% artifact_ids ~ "Artifact Images",
      Subject %in% vascular_ids ~ "Vascular Artifact",
      TRUE ~ "No Artifacts"
    )
  )

# Verificar la asignación de los grupos
table(filtered_df_adni_4$Artifact_Group)

# Scatter plot con los nuevos grupos de artefactos
p7 <- ggplot(filtered_df_adni_4, aes(x = Age, y = TotalGM_B, color = Artifact_Group)) +
  geom_point(alpha = 0.8, size = 3) +  # Ajustar transparencia y tamaño
  scale_color_viridis_d() +
  scale_shape_manual(values = c(16, 17, 18)) +  # Círculo, Triángulo, Diamante
  labs(title = "Age vs sCOVgm (without PVC) by Artifact Group",
       x = "Age", y = "sCOV(%)", shape = "Artifact Group", color = "Artifact Group") +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"),
        legend.position = "top")  # Mover la leyenda arriba

# Guardar la imagen
ggsave(filename = paste0(path_plots_data_analysis, "scatter_Age_vs_sCOV_without_PVC_adni_artifact_groups.png"), 
       plot = p7, width = 8, height = 6, dpi = 300, bg = "white")
