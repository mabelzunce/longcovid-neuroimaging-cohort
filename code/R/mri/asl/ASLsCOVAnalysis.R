# Load necessary libraries
library(mgcv)
library(tidyverse)
library(readxl)
library(car)
library(gridExtra)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(effects)
library(randomForest)

setwd("/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Freesurfer/DataAnalysis/CodeR")
source("functionsASL.R") # Para que funcione debo estar en el mismo directorio

# Define the paths
path_data_analysis <- "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Freesurfer/DataAnalysis/ASL/"
#path_data_analysis <- "/home/solcat/work/data_analysis_paper/"
# path_data_analysis <- "C:/Users/Tomas/Downloads/ASLSolCOVID/"

path_asl_analysis <- "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/ASLProcessedCOVID/derivatives/ExploreASL/Population/Stats/"


path_results_data_analysis <- paste(path_data_analysis, "Results/", sep="")
path_plots_data_analysis <-  paste(path_data_analysis, "Plots/", sep="")
path_cuestionarios <- paste(path_data_analysis, "RespuestasCuestionariosCovid_15_01.csv", sep="")
path_cuestionarios <- paste(path_data_analysis, "ResumenCuestionarios_30_01.csv", sep="")
path_socialIndex <- paste(path_data_analysis, "SocialAndCardiovascularIndex.csv", sep="")


mean_perfussion_data <-"median_qCBF_StandardSpace_TotalGM_n=203_15-Jan-2025_PVC0.tsv"
cov_perfussion_data_without_PVC <- "CoV_qCBF_StandardSpace_TotalGM_n=203_15-Jan-2025_PVC0.tsv"
cov_perfussion_data_with_PVC <- "CoV_qCBF_StandardSpace_TotalGM_n=203_15-Jan-2025_PVC2.tsv"


path_dkt_perfusion_mean <-paste(path_asl_analysis, mean_perfussion_data, sep="")
path_dkt_perfusion_cov_without_pvc <-paste(path_asl_analysis, cov_perfussion_data_without_PVC, sep="")
path_dkt_perfusion_cov_with_pvc <-paste(path_asl_analysis, cov_perfussion_data_with_PVC, sep="")


# Tests para evaluar correlaciones 

tests_variables <- c("EQ.VAS", "FAS", "PSQI",
                     "TMT.A", "TMT.A_perc",
                     "TMT.B", "TMT.B_perc",
                     "WMS.R.DIR", "WMS.R.DIR_perc",
                     "WMS.R.INV", "WMS.R.INV_perc",
                     "STROOP_P", "STROOP_P_perc",
                     "STROOP_C", "STROOP_C_perc",
                     "STROOP_P.C","STROOP_P.C_perc",
                     "STROOP_P.C_INTERF_perc", "STROOP_P.C_INTERF",
                     "MOCA")

columnas_antecedentes <- c("PresionAlta", "Diabetes", "Asma", "ColesterolAlto", "InfartoCardiaco", "AnginaPecho",
                           "EmboliaTrombos", "MedicacionPresionArterial", "Aspirinas", "Antiagregantes")

escalas_depresion_ansiedad <- c("DepresionCategoriaUnificada", "AnsiedadEscalaUnificada")


# Exclude rows
values_to_exclude <- c('CP0011', 'CP0015','CP0009', 'CP0035', 'CP0093', 'CP0101', 'CP0140',
                       'CP0035', 'CP0087', 'CP0105', 'CP0117', 'CP0144', 'CP0192', 'CP0193','CP0216', 'CP0062', 'CP0227')

vascular_artifacts <- c('CP0031','CP0044', 'CP0053', 'CP0054', 'CP0061', 'CP0096','CP0107','CP0108','CP0118',
                        'CP0123', 'CP0142', 'CP0147', 'CP0154','CP0167', 'CP0176','CP0178',
                        'CP0180', 'CP0183', 'CP0188','CP0193','CP0199', 'CP0205', 'CP0213',
                        'CP0229', 'CP0238','CP0245', 'CP0247', 'CP0238')

values_to_exclude <- paste0("sub-", values_to_exclude, "_1")
vascular_artifacts <- paste0("sub-", vascular_artifacts, "_1")


# Read files
cuestionarios_excel <- read.csv(path_cuestionarios)
cuestionarios_excel$ID <- paste0("sub-", cuestionarios_excel$ID,"_1")

perfusion_mean_csv <- read.delim(path_dkt_perfusion_mean)
perfusion_cov_csv_without_pvc <- read.delim(path_dkt_perfusion_cov_without_pvc)
perfusion_cov_csv_with_pvc <- read.delim(path_dkt_perfusion_cov_with_pvc)

social_index <- read.csv(path_socialIndex)
social_index$ID <- paste0("sub-", social_index$ID,"_1")

# Exclude columns with +10 NaN in Perfusion (mean)
perfusion_mean_csv[perfusion_mean_csv == "n/a"] <- NA
cols_with_nan_perfusion <- colSums(is.na(perfusion_mean_csv)) > 10
columns_to_exclude_perfusion <- names(perfusion_mean_csv)[cols_with_nan_perfusion]
perfusion_mean_csv <- perfusion_mean_csv[, !cols_with_nan_perfusion]
colnames(perfusion_mean_csv) <- make.names(colnames(perfusion_mean_csv), unique = TRUE)

# sCOV
columns_to_exclude_perfusion <- names(perfusion_mean_csv)[cols_with_nan_perfusion]

# Columnas con regiones de la segmentación 
column_names_perfusion <- names(perfusion_mean_csv)
columns_regions_perfusion <- column_names_perfusion[6:length(column_names_perfusion)]
columns_regions_perfusion <- columns_regions_perfusion[!grepl("^Wm.", columns_regions_perfusion)]


# Save filtrado perfusion 
perfusion_cov_csv_without_pvc_cp <- subset(perfusion_cov_csv_without_pvc, !(participant_id %in% values_to_exclude))
perfusion_cov_csv_without_pvc_cp$participant_id <- gsub("^sub-|_1$", "", perfusion_cov_csv_without_pvc_cp$participant_id)
write.csv(perfusion_cov_csv_without_pvc_cp, paste0(path_results_data_analysis, 'perfusion_sCOV_PVC0.csv'), row.names = FALSE)


# Unir con variables necesarias de cuestionarios
cols_to_merge <- c("Edad", "Grupo", "ID", "Genero", "BMI", tests_variables, columnas_antecedentes)
cols_to_merge <- intersect(cols_to_merge, colnames(cuestionarios_excel)) 

merged_df_perfusion_mean <- merge(
  cuestionarios_excel[, cols_to_merge], 
  perfusion_mean_csv, 
  by.x = "ID", 
  by.y = "participant_id", 
  all.y = TRUE
) %>% 
  slice(-1) %>%  # Eliminar la primera fila
  mutate(across(all_of(tests_variables), ~ suppressWarnings(as.numeric(.)))) %>%  # Variables tests
  mutate(across(all_of(columns_regions_perfusion), ~ suppressWarnings(as.numeric(.))))  # Regiones %>%

merged_df_perfusion_cov_without_pvc <- merge(
  cuestionarios_excel[, cols_to_merge], 
  perfusion_cov_csv_without_pvc, 
  by.x = "ID", 
  by.y = "participant_id", 
  all.y = TRUE
) %>% 
  slice(-1) %>%  # Eliminar la primera fila
  mutate(across(all_of(tests_variables), ~ suppressWarnings(as.numeric(.)))) %>%  # Variables tests
  mutate(across(all_of(columns_regions_perfusion), ~ suppressWarnings(as.numeric(.))))  # Regiones %>%

merged_df_perfusion_cov_with_pvc <- merge(
  cuestionarios_excel[, cols_to_merge], 
  perfusion_cov_csv_with_pvc, 
  by.x = "ID", 
  by.y = "participant_id", 
  all.y = TRUE
) %>% 
  slice(-1) %>%  # Eliminar la primera fila
  mutate(across(all_of(tests_variables), ~ suppressWarnings(as.numeric(.)))) %>%  # Variables tests
  mutate(across(all_of(columns_regions_perfusion), ~ suppressWarnings(as.numeric(.))))  # Regiones %>%

merged_df_perfusion_cov_without_pvc <- merge(
  social_index, 
  merged_df_perfusion_cov_without_pvc, 
  by.x = "ID", 
  by.y = "ID", 
  all.y = TRUE
) %>% 
  slice(-1)

# Mean
merged_df_perfusion_mean_with_ata <- subset(merged_df_perfusion_mean, !(ID %in% values_to_exclude))
merged_df_perfusion_mean <- subset(merged_df_perfusion_mean_with_ata, !(ID %in% vascular_artifacts)) 

# sCoV without PVC
merged_df_perfusion_cov_with_ata_without_pvc <- subset(merged_df_perfusion_cov_without_pvc, !(ID %in% values_to_exclude))
merged_df_perfusion_cov_without_pvc <- subset(merged_df_perfusion_cov_with_ata_without_pvc, !(ID %in% vascular_artifacts)) 

# sCoV with PVC
merged_df_perfusion_cov_with_ata_with_pvc <- subset(merged_df_perfusion_cov_with_pvc, !(ID %in% values_to_exclude))
merged_df_perfusion_cov_with_pvc <- subset(merged_df_perfusion_cov_with_ata_with_pvc, !(ID %in% vascular_artifacts)) 

# Linear Regression
# Ensure Genero is a factor
merged_df_perfusion_cov_with_ata_without_pvc$Genero <- as.factor(merged_df_perfusion_cov_with_ata_without_pvc$Genero)

# Define dependent variable
merged_df_perfusion_cov_with_ata_without_pvc$TotalGM_B <- merged_df_perfusion_cov_with_ata_without_pvc$TotalGM_B * 100
merged_df_perfusion_cov_with_ata_with_pvc$TotalGM_B <- merged_df_perfusion_cov_with_ata_with_pvc$TotalGM_B *100 

merged_df_perfusion_cov_with_ata_without_pvc$TotalGM_R <- merged_df_perfusion_cov_with_ata_without_pvc$TotalGM_R * 100
merged_df_perfusion_cov_with_ata_with_pvc$TotalGM_R <- merged_df_perfusion_cov_with_ata_with_pvc$TotalGM_R *100 

merged_df_perfusion_cov_with_ata_without_pvc$TotalGM_L <- merged_df_perfusion_cov_with_ata_without_pvc$TotalGM_L * 100
merged_df_perfusion_cov_with_ata_with_pvc$TotalGM_l <- merged_df_perfusion_cov_with_ata_with_pvc$TotalGM_L *100 

# Define a new column to classify vascular artifacts
merged_df_perfusion_cov_with_ata_without_pvc <- merged_df_perfusion_cov_with_ata_without_pvc %>%
  mutate(VascularArtifact = ifelse(ID %in% vascular_artifacts, "With Vascular Artifact", "Without Vascular Artifact"))

merged_df_perfusion_cov_with_ata_with_pvc <- merged_df_perfusion_cov_with_ata_with_pvc %>%
  mutate(VascularArtifact = ifelse(ID %in% vascular_artifacts, "With Vascular Artifact", "Without Vascular Artifact"))


# Convert the new variable to a factor for ggplot aesthetics
merged_df_perfusion_cov_with_ata_without_pvc$VascularArtifact <- factor(merged_df_perfusion_cov_with_ata_without_pvc$VascularArtifact, levels = c("Without Vascular Artifact", "With Vascular Artifact"))
merged_df_perfusion_cov_with_ata_with_pvc$VascularArtifact <- factor(merged_df_perfusion_cov_with_ata_with_pvc$VascularArtifact, levels = c("Without Vascular Artifact", "With Vascular Artifact"))

# Linear regression 
model_TotalGM_B_group <- lm(TotalGM_B ~ Edad + Genero + Grupo, data = merged_df_perfusion_cov_with_ata_without_pvc)
summary(model_TotalGM_B_group)

model_TotalGM_R_group <- lm(TotalGM_R ~ Edad + Genero + Grupo, data = merged_df_perfusion_cov_with_ata_without_pvc)
summary(model_TotalGM_R_group)

model_TotalGM_L_group <- lm(TotalGM_L ~ Edad + Genero + Grupo, data = merged_df_perfusion_cov_with_ata_without_pvc)
summary(model_TotalGM_L_group)


plot_scov_b <- ggplot(merged_df_perfusion_cov_with_ata_without_pvc, aes(x = Grupo, y = TotalGM_B, fill = Grupo)) +
  geom_boxplot(alpha = 0.5) +  # Add transparency to boxplot
  geom_jitter(color = "black", size = 2, alpha = 0.3, width = 0.2) +  # Add jittered points
  labs(title = "Spatial Coefficient of Variation (sCOV)of Total Gray Matter", x = "Group", y = "sCOV (%)")+
  theme_minimal() +
  guides(fill = FALSE)# +   Remove legend
  # theme(
  #   plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),  # Increase title size
  #   axis.text = element_text(size = 14),  # Increase axis text size
  #   axis.title = element_text(size = 16)  # Increase axis title size
  # )
ggsave(filename = paste0(path_plots_data_analysis, "boxplot_TotalGM_B_without_pvc_by_group.png"), plot = plot_scov_b, width = 3, height = 6, dpi = 300, bg = "white")


plot_scov_l <- ggplot(merged_df_perfusion_cov_with_ata_without_pvc, aes(x = Grupo, y = TotalGM_L, fill = Grupo)) +
  geom_boxplot(alpha = 0.5) +  # Add transparency to boxplot
  geom_jitter(color = "black", size = 2, alpha = 0.3, width = 0.2) +  # Add jittered points
  labs(title = "Left Grey Matter sCOV", x = "Group", y = "sCOV (%)")+
  theme_minimal() +
  guides(fill = FALSE)# +   Remove legend
# theme(
#   plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),  # Increase title size
#   axis.text = element_text(size = 14),  # Increase axis text size
#   axis.title = element_text(size = 16)  # Increase axis title size
# )
ggsave(filename = paste0(path_plots_data_analysis, "boxplot_TotalGM_L_without_pvc_by_group.png"), plot = plot_scov_l, width = 3, height = 6, dpi = 300, bg = "white")


plot_scov_r <- ggplot(merged_df_perfusion_cov_with_ata_without_pvc, aes(x = Grupo, y = TotalGM_R, fill = Grupo)) +
  geom_boxplot(alpha = 0.5) +  # Add transparency to boxplot
  geom_jitter(color = "black", size = 2, alpha = 0.3, width = 0.2) +  # Add jittered points
  labs(title = "Right Grey Matter sCOV", x = "Group", y = "sCOV (%)")+
  theme_minimal() +
  guides(fill = FALSE)# +   Remove legend
# theme(
#   plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),  # Increase title size
#   axis.text = element_text(size = 14),  # Increase axis text size
#   axis.title = element_text(size = 16)  # Increase axis title size
# )
ggsave(filename = paste0(path_plots_data_analysis, "boxplot_TotalGM_R_without_pvc_by_group.png"), plot = plot_scov_r, width = 3, height = 6, dpi = 300, bg = "white")


plot_scov_regression_b <- ggplot(merged_df_perfusion_cov_with_ata_without_pvc, aes(x = Edad, y = TotalGM_B, color = Grupo)) +
  geom_point(alpha = 0.8, size = 3) +  # Adjust transparency and size of points
  geom_smooth(method = "lm", se = TRUE, aes(group = Grupo)) +  # Add a regression line for each group
  scale_color_viridis_d() +  # Color scale for points and lines
  scale_shape_manual(values = c(15, 17)) +  # Define shapes, if needed
  labs(title = "Age vs sCOV",
       x = "Age", y = "sCOV(%)", shape = "Vascular Artifact", color = "Grupo") +
  theme_minimal(base_size = 14) +  # Minimal theme with base size for text
  theme(panel.background = element_rect(fill = "white"), 
        plot.background = element_rect(fill = "white"),
        legend.position = "top")  # Adjust legend position for better readability

# Save the plot to a file
ggsave(filename = paste0(path_plots_data_analysis, "scatter_Age_vs_sCOV_without_PVC_b.png"), 
       plot = p3_new, width = 8, height = 6, dpi = 300, bg = "white")


plot_scov_regression_l <- ggplot(merged_df_perfusion_cov_with_ata_without_pvc, aes(x = Edad, y = TotalGM_L, color = Grupo)) +
  geom_point(alpha = 0.8, size = 3) +  # Adjust transparency and size of points
  geom_smooth(method = "lm", se = TRUE, aes(group = Grupo)) +  # Add a regression line for each group
  scale_color_viridis_d() +  # Color scale for points and lines
  scale_shape_manual(values = c(15, 17)) +  # Define shapes, if needed
  labs(title = "Age vs sCOV",
       x = "Age", y = "sCOV(%)", shape = "Vascular Artifact", color = "Grupo") +
  theme_minimal(base_size = 14) +  # Minimal theme with base size for text
  theme(panel.background = element_rect(fill = "white"), 
        plot.background = element_rect(fill = "white"),
        legend.position = "top")  # Adjust legend position for better readability

# Save the plot to a file
ggsave(filename = paste0(path_plots_data_analysis, "scatter_Age_vs_sCOV_without_PVC_L.png"), 
       plot = plot_scov_regression_l, width = 8, height = 6, dpi = 300, bg = "white")



plot_scov_regression_r <- ggplot(merged_df_perfusion_cov_with_ata_without_pvc, aes(x = Edad, y = TotalGM_R, color = Grupo)) +
  geom_point(alpha = 0.8, size = 3) +  # Adjust transparency and size of points
  geom_smooth(method = "lm", se = TRUE, aes(group = Grupo)) +  # Add a regression line for each group
  scale_color_viridis_d() +  # Color scale for points and lines
  scale_shape_manual(values = c(15, 17)) +  # Define shapes, if needed
  labs(title = "Age vs sCOV",
       x = "Age", y = "sCOV(%)", shape = "Vascular Artifact", color = "Grupo") +
  theme_minimal(base_size = 14) +  # Minimal theme with base size for text
  theme(panel.background = element_rect(fill = "white"), 
        plot.background = element_rect(fill = "white"),
        legend.position = "top")  # Adjust legend position for better readability

# Save the plot to a file
ggsave(filename = paste0(path_plots_data_analysis, "scatter_Age_vs_sCOV_without_PVC_R.png"), 
       plot = plot_scov_regression_l, width = 8, height = 6, dpi = 300, bg = "white")


##PLOTS FINALES SCOV
plot_gm_r <- ggplot(merged_df_perfusion_cov_with_ata_without_pvc, aes(x = Grupo, y = get("TotalGM_L"), fill = Grupo)) +
  geom_boxplot(alpha = 0.5) +  # Add transparency to boxplot
  geom_jitter(color = "black", size = 2, alpha = 0.3, width = 0.2) +  # Add jittered points
  labs(title = "Left Grey Matter sCOV", x = "Group", y = "sCOV (%)")+
  theme_minimal() +
  guides(fill = FALSE) +  # Remove legend
  theme(
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),  # Increase title size
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )+
  ylim(30, 90)  # Set y-axis limits


# Save the plot
plot_filename <- paste(path_plots_data_analysis, "Total_GM_L_all_sCOV.jpg", sep="")
ggsave(plot_filename, plot = plot_gm_r , device = "jpeg", bg = "white", width = 4, height = 7)


plot_gm_r <- ggplot(merged_df_perfusion_cov_with_ata_without_pvc, aes(x = Grupo, y = get("TotalGM_R"), fill = Grupo)) +
  geom_boxplot(alpha = 0.5) +  # Add transparency to boxplot
  geom_jitter(color = "black", size = 2, alpha = 0.3, width = 0.2) +  # Add jittered points
  labs(title = "Right Grey Matter sCOV", x = "Group", y = "sCOV (%)") +
  theme_minimal() +
  guides(fill = FALSE) +  # Remove legend
  theme(
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),  # Increase title size
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )+
  ylim(30, 90)  # Set y-axis limits

# Save the plot
plot_filename <- paste(path_plots_data_analysis, "Total_GM_R_all_SCOV.jpg", sep="")
ggsave(plot_filename, plot = plot_gm_r , device = "jpeg", bg = "white", width = 4, height = 7)


# Correlacion 
correlation_data <- merged_df_perfusion_cov_with_ata_without_pvc %>%
  dplyr::select(TotalGM_B, dplyr::all_of(tests_variables))%>%
  mutate(across(everything(), as.numeric))  # Asegurar que todas las variables sean numéricas

# Calcular la matriz de correlación de Pearson
cor_matrix <- cor(correlation_data, use = "complete.obs", method = "pearson")

other_variables_to_correlated <- c("Social.Index",  "CardiovascularRiskIndex" ,"WMH_vol", "WMH_count", "Edad")

correlation_data_2 <- merged_df_perfusion_cov_with_ata_without_pvc %>%
  dplyr::select(TotalGM_B, dplyr::all_of(other_variables_to_correlated))%>%
  mutate(across(everything(), as.numeric))  # Asegurar que todas las variables sean numéricas

# Calcular la matriz de correlación de Pearson
cor_matrix_2 <- cor(correlation_data_2, use = "complete.obs", method = "pearson")

# # Add ADNI images
# 
# # ADNI
# path_asl_analysis_adni <- "/media/sol/Expansion/ADNI/ADNI_ExploreASL_3D/derivatives/ExploreASL/Population/Stats/"
# path_info_images_adni <- "/media/sol/Expansion/ADNI/all_processed_images_info.csv"
# raw_perfussion_data_adni <-"mean_qCBF_StandardSpace_TotalGM_n=198_28-Jan-2025_PVC2.tsv"
# path_dkt_perfusion_adni <-paste(path_asl_analysis_adni, raw_perfussion_data_adni , sep="")
# 
# 
# 
# # ADNI
# path_asl_analysis_adni <- "/media/sol/Expansion/ADNI/ADNI_ExploreASL_3D/derivatives/ExploreASL/Population/Stats/"
# path_info_images_adni <- "/media/sol/Expansion/ADNI/all_processed_images_info.csv"
# raw_perfussion_data_adni <-"CoV_qCBF_StandardSpace_TotalGM_n=198_28-Jan-2025_PVC0.tsv"
# path_dkt_perfusion_adni <-paste(path_asl_analysis_adni, raw_perfussion_data_adni , sep="")
# 
# 
# # ADNI Data 
# # Info images 
# info_images_adni <-read.csv(path_info_images_adni)
# 
# # Unir con variables necesarias de la información adicional
# cols_to_merge_adni <- c("Age", "Group", "Sex", "Subject", "Site")
# cols_to_merge_adni <- intersect(cols_to_merge_adni, colnames(info_images_adni))  
# 
# # Limpiar tsv perfusión
# perfusion_csv_adni <- read.delim(path_dkt_perfusion_adni)
# # perfusion_csv_adni <- perfusion_csv_adni[, !cols_with_nan_perfusion]
# colnames(perfusion_csv_adni) <- make.names(colnames(perfusion_csv_adni), unique = TRUE)
# 
# # Values to exclude 
# # ASL images
# asl_images_path <- "/media/sol/Expansion/ADNI/ADNI_ExploreASL_3D/derivatives/ExploreASL/Population/ASLCheck/"
# 
# artifact_images_path <- paste(asl_images_path, "3_ArtifactContrast", sep="")
# vascular_images_path <- paste(asl_images_path, "2_VascularContrast", sep="")
# processed_images_path <- paste(asl_images_path, "1_CBFContrast", sep="")
# 
# 
# artifact_images <- list.files(path = artifact_images_path, pattern = "\\.jpg$", full.names = FALSE)
# vascular_images <- list.files(path = vascular_images_path, pattern = "\\.jpg$", full.names = FALSE)
# processed_images <- list.files(path = processed_images_path, pattern = "\\.jpg$", full.names = FALSE)
# 
# patron <- "(sub[^_]*_1)"
# 
# artifact_ids <- str_extract(artifact_images, patron)
# vascular_ids <- str_extract(vascular_images, patron)
# processed_ids <- str_extract(processed_images, patron)
# 
# merged_df_perfusion_adni <- merge(
#   info_images_adni[, cols_to_merge_adni], 
#   perfusion_csv_adni, 
#   by.x = "Subject", 
#   by.y = "participant_id", 
#   all.y = TRUE
# ) %>% 
#   slice(-1)  # Eliminar la primera fila
# library(dplyr)
# 
# 
# # Lista de columnas a convertir a numérico
# region_columns <- c("Age", "GM_vol", "WM_vol", "CSF_vol", 
#                     "GM_ICVRatio", "GMWM_ICVRatio", "MeanMotion", 
#                     "TotalGM_B", "TotalGM_L", "TotalGM_R")
# 
# # Convertir solo las columnas de regiones a numérico
# merged_df_perfusion_adni <- merged_df_perfusion_adni %>%
#   rename(ID = Subject) %>%  # Si necesitas renombrar la columna "Subject" a "ID"
#   mutate(across(all_of(region_columns), ~ suppressWarnings(as.numeric(.))))
# 
# merged_df_perfusion_adni$TotalGM_B <- merged_df_perfusion_adni$TotalGM_B * 100
# 
# 
# merged_combined_df <- bind_rows(merged_df_perfusion_adni, merged_df_perfusion_cov_with_ata_without_pvc)
# 
# filtered_merged_combined_df <- merged_combined_df %>% 
#     filter(TotalGM_B != 0) %>%
#     filter(ID != "sub-941S4365_1")
# 
# 
# library(ggplot2)
# library(RColorBrewer)
# 
# # Scatter plot con una paleta más distintiva
# p6 <- ggplot(filtered_merged_combined_df, aes(x = Age, y = TotalGM_B, color = Group)) +
#   geom_point(alpha = 0.8, size = 3) +  
#   scale_color_brewer(palette = "Set1") +  # Paleta con colores más contrastantes
#   scale_shape_manual(values = c(15, 17)) +  
#   labs(title = "Age vs sCOVgm (without PVC)",
#        x = "Age", y = "sCOV(%)", shape = "Vascular Artifact", color = "Grupo") +
#   theme_minimal(base_size = 14) +
#   theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"),
#         legend.position = "top")  
# 
# ggsave(filename = paste0(path_plots_data_analysis, "scatter_Age_vs_sCOV_without_PVC_adni_and_covid.png"), 
#        plot = p6, width = 8, height = 6, dpi = 300, bg = "white")
# 
# 
# variables_to_correlate <- c("WMH_vol", "WMH_count", "Age")
# 
# correlation_data_2 <- filtered_merged_combined_df %>%
#   dplyr::select(TotalGM_B, dplyr::all_of(variables_to_correlate))%>%
#   mutate(across(everything(), as.numeric))  # Asegurar que todas las variables sean numéricas
# 
# # Calcular la matriz de correlación de Pearson
# cor_matrix_2 <- cor(correlation_data_2, use = "complete.obs", method = "pearson")
# 
# # Filtrar solo los participantes de 60 años o más
# filtered_merged_combined_df_old <- filtered_merged_combined_df %>%
#   filter(Age >= 65)
# 
# correlation_data_old <- filtered_merged_combined_df_old %>%
#   dplyr::select(TotalGM_B, dplyr::all_of(variables_to_correlate))%>%
#   mutate(across(everything(), as.numeric))  # Asegurar que todas las variables sean numéricas
# 
# # Calcular la matriz de correlación de Pearson
# cor_matrix_old <- cor(correlation_data_old, use = "complete.obs", method = "pearson")
# 
# p6 <- ggplot(filtered_merged_combined_df_old, aes(x = Age, y = TotalGM_B, color = Group)) +
#   geom_point(alpha = 0.8, size = 3) +
#   scale_color_brewer(palette = "Set1") +  # Paleta con colores más contrastantes
#   scale_shape_manual(values = c(15, 17)) +
#   labs(title = "Age vs sCOVgm (without PVC)",
#        x = "Age", y = "sCOV(%)", shape = "Vascular Artifact", color = "Grupo") +
#   theme_minimal(base_size = 14) +
#   theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"),
#         legend.position = "top")
# 
# ggsave(filename = paste0(path_plots_data_analysis, "scatter_Age_vs_sCOV_without_PVC_adni_and_covid_old.png"),
#        plot = p6, width = 8, height = 6, dpi = 300, bg = "white")
# 
# 
