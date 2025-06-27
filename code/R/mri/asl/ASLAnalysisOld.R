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
library(dplyr)


setwd("/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Freesurfer/DataAnalysis/CodeR")
source("functionsASL.R") # Para que funcione debo estar en el mismo directorio

# Define the paths
path_data_analysis <- "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Freesurfer/DataAnalysis/ASL/"
path_asl_analysis <- "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/Nifti/Final_data_csv/DataISMRM/"

path_results_data_analysis <- paste(path_data_analysis, "Results/", sep="")
path_plots_data_analysis <-  paste(path_data_analysis, "Plots/", sep="")
path_cuestionarios <- paste("/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Freesurfer/DataAnalysis/ResumenTotal_04_04.csv", sep="")

# path_cuestionarios <- paste(path_data_analysis, "RespuestasCuestionariosCovid_24_10.csv", sep="")
# path_voluntarios <- '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/VoluntariosProyectoCovidProlongado.csv'
# 
# raw_perfussion_data <-"mean_qCBF_StandardSpace_Desikan_Killiany_MNI_SPM12_PVC2.csv"
raw_perfussion_data <-"mean_qCBF_StandardSpace_Hammers_PVC2.csv"
raw_gray_matter_data <-"mean_qCBF_StandardSpace_TotalGM_2024_PVC2.csv"
path_scov_ismrm_data <- paste(path_asl_analysis, "CoV_qCBF_StandardSpace_TotalGM_ISRM.csv",sep="")

path_dkt_perfusion <-paste(path_asl_analysis, raw_perfussion_data , sep="")
path_gray_matter_perfusion <-paste(path_asl_analysis, raw_gray_matter_data , sep="")

tests_variables <- c("EQ.VAS", "FAS", "PSQI",
                     "TMT.A", "TMT.A_perc",
                     "TMT.B", "TMT.B_perc",
                     "WMS.R.DIR", "WMS.R.DIR_perc",
                     "WMS.R.INV", "WMS.R.INV_perc",
                     "STROOP_P", "STROOP_P_perc",
                     "STROOP_C", "STROOP_C_perc",
                     "STROOP_P.C","STROOP_P.C_perc",
                     "STROOP_P.C_INTERF_perc", "STROOP_P.C_INTERF",
                     "MOCA", "MOCA_perc")

tests_variables_perc_and_scales <- tests_variables[grepl("_perc$", tests_variables) | 
                                                     tests_variables %in% c("EQ.VAS", "FAS", "PSQI")]


# Exclude rows
# values_to_exclude <- c('CP0011', 'CP0062', 'CP0105', 'CP0216', 
#                        'CP0075', 'CP0076', 'CP0078', 
#                        'CP0079', 'CP0140', 'CP0144','CP0196')


values_to_exclude <- c('CP0011', 'CP0015','CP0009', 'CP0035', 'CP0093', 'CP0101', 'CP0140','CP0106',
                       'CP0087', 'CP0105', 'CP0117', 'CP0144', 'CP0192','CP0216', 'CP0193','CP0227', 'CP0062',
                       'CP0075', 'CP0076', 'CP0078', 'CP0079', 'CP0196')

vascular_artifacts <- c('CP0031','CP0044', 'CP0053', 'CP0054', 'CP0061', 'CP0096','CP0107','CP0108','CP0118',
                        'CP0123', 'CP0142', 'CP0147', 'CP0154','CP0167', 'CP0176','CP0178',
                        'CP0180', 'CP0183', 'CP0188','CP0193','CP0199', 'CP0205', 'CP0213',
                        'CP0229', 'CP0238','CP0245', 'CP0247')

# regions_volume <- c('Ctx.lingual_L', 'Ctx.inferiorparietal_L', 'Ctx.inferiortemporal_R')
# regions_thick <- c('Ctx.lingual_L', 'Ctx.postcentral_B',  'Ctx.superiorparietal_L', 'Ctx.precuneus_L', 
#                    'Ctx.supramarginal_L', 'Ctx.bankssts_L', 'Ctx.precuneus_R', 'Ctx.superiortemporal_R'
#                    )


# Read files
cuestionarios_excel <- read.csv(path_cuestionarios)
perfusion_csv <- read.csv(path_dkt_perfusion)
gray_matter_csv <- read.csv(path_gray_matter_perfusion)
perfusion_csv <- read.csv(path_dkt_perfusion)
perfusion_scov_csv <- read.csv(path_scov_ismrm_data)

# Exclude columns with +10 NaN in Perfusion
perfusion_csv[perfusion_csv == "n/a"] <- NA
cols_with_nan_perfusion <- colSums(is.na(perfusion_csv)) > 10
columns_to_exclude_perfusion <- names(perfusion_csv)[cols_with_nan_perfusion]
perfusion_csv <- perfusion_csv[, !cols_with_nan_perfusion]
colnames(perfusion_csv) <- make.names(colnames(perfusion_csv), unique = TRUE)

# Columnas con regiones de la segmentación 
column_names_perfusion <- names(perfusion_csv)
columns_regions_perfusion <- column_names_perfusion[3:length(column_names_perfusion) - 1]
columns_regions_perfusion <- columns_regions_perfusion[!grepl("^Wm.", columns_regions_perfusion)]
columns_regions_perfusion <- columns_regions_perfusion[9:length(columns_regions_perfusion)]

# Unir con variables necesarias de cuestionarios
cols_to_merge <- c("Edad", "Grupo", "ID", "Genero", "BMI",tests_variables)
cols_to_merge <- intersect(cols_to_merge, colnames(cuestionarios_excel))  

merged_df_perfusion <- merge(
  cuestionarios_excel[, cols_to_merge], 
  perfusion_csv, 
  by.x = "ID", 
  by.y = "ID", 
  all.y = TRUE
) 

cols_to_merge_gm <- c('ID','TotalGM_B', 'TotalGM_L', 'TotalGM_R')

merged_df_perfusion <- merge(
  merged_df_perfusion, 
  gray_matter_csv[, cols_to_merge_gm],
  by.x = "ID", 
  by.y = "ID", 
  all.y = TRUE
) 

# scov
perfusion_scov_csv <- perfusion_scov_csv %>%
  rename("TotalGM_B_scov" = "TotalGM_B",
         "TotalGM_L_scov" = "TotalGM_L", 
         "TotalGM_R_scov" = "TotalGM_R")


merged_df_perfusion <- merge(
  merged_df_perfusion, 
  perfusion_scov_csv, 
  by.x = "ID", 
  by.y = "ID", 
  all.x = TRUE
) 

merged_df_perfusion_with_ata <- subset(merged_df_perfusion, !(ID %in% values_to_exclude))
merged_df_perfusion <- subset(merged_df_perfusion_with_ata, !(ID %in% vascular_artifacts)) 

colnames(merged_df_perfusion_with_ata) <- make.names(colnames(merged_df_perfusion_with_ata), unique = TRUE)
colnames(merged_df_perfusion) <- make.names(colnames(merged_df_perfusion), unique = TRUE)


# merged_df_perfusion_with_ata <- subset(merged_df_perfusion_with_ata, select = -TL_hippocampus_R)
# merged_df_perfusion_with_ata <- subset(merged_df_perfusion_with_ata, select = -TL_superior_temporal_gyrus_anterior_part_R )

# multiple_ancova_correcting_scov(merged_df_perfusion_with_ata, columns_regions_perfusion)
multiple_ancova_old(merged_df_perfusion_with_ata, columns_regions_perfusion)
multiple_ancova_old(merged_df_perfusion_with_ata, c('TotalGM_L_scov', 'TotalGM_R_scov'))

multiple_ancova_old(merged_df_perfusion_with_ata, columns_regions_perfusion)
# multiple_ancova_old(merged_df_perfusion_with_ata, c('TotalGM_B', 'TotalGM_L', 'TotalGM_R'))


# Reemplazar las columnas de perfusión con sus valores normalizados
# normalized_df_perfusion <- merged_df_perfusion %>%
#   mutate(across(all_of(columns_regions_perfusion), ~ . / TotalGM_B))
# 
# normalized_df_perfusion_with_ata <- merged_df_perfusion_with_ata %>%
#   mutate(across(all_of(columns_regions_perfusion), ~ . / TotalGM_B))
# 
# multiple_ancova_old(normalized_df_perfusion_with_ata, columns_regions_perfusion)
# 
# outliers_covid_df <- detect_outliers_covid(merged_df_perfusion, columns_regions_perfusion)
# 
# # Create a boxplot for each region
# ggplot(outliers_covid_df, aes(x = Region, y = Value)) +
#   geom_boxplot(outlier.colour = "red", outlier.shape = 16, fill = "lightblue") +
#   labs(title = "Lower Outliers in COVID Group",
#        x = "Brain Region",
#        y = "Perfusion Value") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for readability

regions <- c("TL_superior_temporal_gyrus_anterior_part_R", "PL_superior_parietal_gyrus_R",
             "FL_anterior_orbital_gyrus_L", "FL_middle_frontal_gyrus_L",
             "TL_superior_temporal_gyrus_middle_part_L", "PL_postcentral_gyrus_R")

# # Loop through each region and create a boxplot
for (region in regions) {
  plot <- ggplot(merged_df_perfusion, aes(x = Grupo, y = get(region), fill = Grupo)) +
    geom_boxplot(alpha = 0.5) +  # Add transparency to boxplot
    geom_jitter(color = "black", size = 2, alpha = 0.3, width = 0.2) +  # Add jittered points
    # geom_text(aes(label = ID), position = position_jitter(width = 0.2, height = 0),
    #           vjust = -1, size = 2.5, alpha = 0.8) +
    labs(title = region, x = "Group", y = "CBF (mL/100g/min)") +
    theme_minimal() +
    guides(fill = FALSE) +  # Remove legend
    theme(
      plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),  # Increase title size
      axis.text = element_text(size = 14),  # Increase axis text size
      axis.title = element_text(size = 16)  # Increase axis title size
    )

  # Save the plot
  plot_filename <- paste(path_plots_data_analysis, region, "_AAIC_NEW.jpg", sep="")
  ggsave(plot_filename, plot = plot, device = "jpeg", bg = "white", width = 4, height = 7)
}

# Total Grey Matter
plot_gm_l <- ggplot(merged_df_perfusion, aes(x = Grupo, y = get("TotalGM_L"), fill = Grupo)) +
  geom_boxplot(alpha = 0.5) +  # Add transparency to boxplot
  geom_jitter(color = "black", size = 2, alpha = 0.3, width = 0.2) +  # Add jittered points
  labs(title = "Left Grey Matter", x = "Group", y = "CBF (mL/100g/min)")+
  theme_minimal() +
  ylim(23, 35) +
  guides(fill = FALSE) +  # Remove legend
  theme(
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),  # Increase title size
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )

# Save the plot
plot_filename <- paste(path_plots_data_analysis, "Total_GM_L_best.jpg", sep="")
ggsave(plot_filename, plot = plot_gm_l , device = "jpeg", bg = "white", width = 4, height = 7)


plot_gm_r <- ggplot(merged_df_perfusion, aes(x = Grupo, y = get("TotalGM_R"), fill = Grupo)) +
  geom_boxplot(alpha = 0.5) +  # Add transparency to boxplot
  geom_jitter(color = "black", size = 2, alpha = 0.3, width = 0.2) +  # Add jittered points
  ylim(23, 35) +
  labs(title = "Right Grey Matter", x = "Group", y = "CBF (mL/100g/min)") +
  theme_minimal() +
  guides(fill = FALSE) +  # Remove legend
  theme(
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),  # Increase title size
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )

# Save the plot
plot_filename <- paste(path_plots_data_analysis, "Total_GM_R_best.jpg", sep="")
ggsave(plot_filename, plot = plot_gm_r , device = "jpeg", bg = "white", width = 4, height = 7)

#scovmultiple_ancova_old(merged_df_perfusion_scov_with_ata, c('TotalGM_B', 'TotalGM_L', 'TotalGM_R'))

merged_df_perfusion_with_ata$TotalGM_R_scov <- merged_df_perfusion_with_ata$TotalGM_R_scov * 100
merged_df_perfusion_with_ata$TotalGM_L_scov <- merged_df_perfusion_with_ata$TotalGM_L_scov * 100


plot_gm_l <- ggplot(merged_df_perfusion_with_ata, aes(x = Grupo, y = get("TotalGM_L_scov"), fill = Grupo)) +
  geom_boxplot(alpha = 0.5) +  # Add transparency to boxplot
  geom_jitter(color = "black", size = 2, alpha = 0.3, width = 0.2) +  # Add jittered points
  labs(title = "Left Grey Matter", x = "Group", y = "sCOV (%)")+
  theme_minimal() +
  ylim(30, 100) +
  guides(fill = FALSE) +  # Remove legend
  theme(
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),  # Increase title size
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )

# Save the plot
plot_filename <- paste(path_plots_data_analysis, "Total_GM_L_ISMRM.jpg", sep="")
ggsave(plot_filename, plot = plot_gm_l , device = "jpeg", bg = "white", width = 4, height = 7)


plot_gm_r <- ggplot(merged_df_perfusion_with_ata, aes(x = Grupo, y = get("TotalGM_R_scov"), fill = Grupo)) +
  geom_boxplot(alpha = 0.5) +  # Add transparency to boxplot
  geom_jitter(color = "black", size = 2, alpha = 0.3, width = 0.2) +  # Add jittered points
  ylim(30, 100) +
  labs(title = "Right Grey Matter", x = "Group", y = "sCOV (%)") +
  theme_minimal() +
  guides(fill = FALSE) +  # Remove legend
  theme(
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),  # Increase title size
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )

# Save the plot
plot_filename <- paste(path_plots_data_analysis, "Total_GM_R_ISMRM.jpg", sep="")
ggsave(plot_filename, plot = plot_gm_r , device = "jpeg", bg = "white", width = 4, height = 7)

# Linear regression
merged_df_perfusion_with_ata$Genero <- as.factor(merged_df_perfusion_with_ata$Genero)

model_TotalGM_R_group <- lm(TotalGM_R_scov ~ Edad + Genero + Grupo , data = merged_df_perfusion_with_ata)
summary(model_TotalGM_R_group)

model_TotalGM_L_group <- lm(TotalGM_L_scov ~ Edad + Genero + Grupo , data = merged_df_perfusion_with_ata)
summary(model_TotalGM_L_group)


# Dispersion Right 
plot_dispersion_R <- ggplot(merged_df_perfusion_with_ata, aes(x = Edad, y = TotalGM_R_scov, color = Grupo)) +
  geom_point(alpha = 0.8, size = 3) +  # Adjust transparency and size
  scale_color_viridis_d() 
  labs(title = "Age vs sCOV",
       x = "Age", y = "sCOV(%)", color = "Grupo") +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"),
        legend.position = "top")  # Move legend to the top for better readability

ggsave(filename = paste0(path_plots_data_analysis, "scatter_Age_vs_sCOV_R_ismrm.png"), 
       plot = plot_dispersion_R, width = 8, height = 6, dpi = 300, bg = "white")



plot_dispersion_L <- ggplot(merged_df_perfusion_with_ata, aes(x = Edad, y = TotalGM_L_scov, color = Grupo)) +
  geom_point(alpha = 0.8, size = 3) +  # Adjust transparency and size
  scale_color_viridis_d() +
  labs(title = "Age vs sCOV",
       x = "Age", y = "sCOV(%)", color = "Grupo") +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"),
        legend.position = "top")  # Move legend to the top for better readability

ggsave(filename = paste0(path_plots_data_analysis, "scatter_Age_vs_sCOV_L_ismrm.png"), 
       plot = plot_dispersion_L, width = 8, height = 6, dpi = 300, bg = "white")

# Correlation between tests and sCOV
char_vars <- c("WMS.R.DIR_perc", "WMS.R.INV_perc", 
               "STROOP_P_perc", "STROOP_C_perc", 
               "STROOP_P.C_perc", "STROOP_P.C_INTERF_perc")

for (var in char_vars) {
  merged_df_perfusion_with_ata[[var]] <- gsub("<", "", merged_df_perfusion_with_ata[[var]])
  merged_df_perfusion_with_ata[[var]] <- as.numeric(merged_df_perfusion_with_ata[[var]])
}

cor_scov <- cor(merged_df_perfusion_with_ata[, c(tests_variables_perc_and_scales, "TotalGM_L_scov", "TotalGM_R_scov" )], use="pairwise.complete.obs")

# Definir la ruta de guardado del gráfico
output_path_cor_SCOV_covid <- paste0(path_plots_data_analysis, "correlation_plot_cor_SCOV_ismrm.png")

# Guardar gráfico unificado con colores pastel
png(output_path_cor_SCOV_covid, width=800, height=800)

corrplot(cor_scov, method="circle", type="lower", order="original", 
         col=colorRampPalette(c("lightblue", "white", "lightcoral"))(200),
         tl.col="black", tl.srt=45, addCoef.col="black")

title("Correlation Matrix: Standardised Scales (COVID)", cex.main=1.5)

dev.off()

