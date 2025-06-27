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

source("functionsANCOVA.R")
source("functionsASL.R") # Para que funcione debo estar en el mismo directorio

# Define the paths
path_data_analysis <- "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Freesurfer/DataAnalysis/ASL/"
path_asl_analysis <- "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/ASLProcessedCOVID/derivatives/ExploreASL/Population/Stats/"

path_results_data_analysis <- paste(path_data_analysis, "Results/", sep="")
path_plots_data_analysis <-  paste(path_data_analysis, "Plots/", sep="")
path_cuestionarios <- paste(path_data_analysis, "RespuestasCuestionariosCovid_15_01.csv", sep="")

cov_perfussion_data_without_PVC <- "CoV_qCBF_StandardSpace_TotalGM_n=203_15-Jan-2025_PVC0.tsv"


raw_perfussion_data <-"mean_qCBF_StandardSpace_Hammers_n=203_15-Jan-2025_PVC2.tsv"
raw_gray_matter_data <-"mean_qCBF_StandardSpace_TotalGM_n=203_15-Jan-2025_PVC2.tsv"

# raw_perfussion_data <-"mean_qCBF_StandardSpace_TotalGM_n=203_15-Jan-2025_PVC2.tsv"

path_dkt_perfusion <-paste(path_asl_analysis, raw_perfussion_data , sep="")
path_gray_matter_perfusion <-paste(path_asl_analysis, raw_gray_matter_data , sep="")
path_scov_perfusion <-paste(path_asl_analysis, cov_perfussion_data_without_PVC , sep="")


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


# Exclude rows
values_to_exclude <- c('CP0011', 'CP0015','CP0009', 'CP0035', 'CP0093', 'CP0101', 'CP0140',
                       'CP0035', 'CP0087', 'CP0105', 'CP0117', 'CP0144', 'CP0192','CP0216', 'CP0227', 'CP0062')

vascular_artifacts <- c('CP0031','CP0044', 'CP0053', 'CP0054', 'CP0061', 'CP0096','CP0107','CP0108','CP0118',
                        'CP0123', 'CP0142', 'CP0147', 'CP0154','CP0167', 'CP0176','CP0178',
                         'CP0180', 'CP0183', 'CP0188','CP0193','CP0199', 'CP0205', 'CP0213',
                        'CP0229', 'CP0238','CP0245', 'CP0247', 'CP0238')

values_to_exclude <- paste0("sub-", values_to_exclude, "_1")
vascular_artifacts <- paste0("sub-", vascular_artifacts, "_1")

# Read files
cuestionarios_excel <- read.csv(path_cuestionarios)
cuestionarios_excel$ID <- paste0("sub-", cuestionarios_excel$ID,"_1")
perfusion_csv <- read.delim(path_dkt_perfusion)
gray_matter_csv <- read.delim(path_gray_matter_perfusion)

# Exclude columns with +10 NaN in Perfusion
perfusion_csv[perfusion_csv == "n/a"] <- NA
cols_with_nan_perfusion <- colSums(is.na(perfusion_csv)) > 10
columns_to_exclude_perfusion <- names(perfusion_csv)[cols_with_nan_perfusion]
perfusion_csv <- perfusion_csv[, !cols_with_nan_perfusion]
colnames(perfusion_csv) <- make.names(colnames(perfusion_csv), unique = TRUE)


# Columnas con regiones de la segmentación 
column_names_perfusion <- names(perfusion_csv)
columns_regions_perfusion <- column_names_perfusion[6:length(column_names_perfusion)]
columns_regions_perfusion <- columns_regions_perfusion[!grepl("^Wm.", columns_regions_perfusion)]

# Save filtrado perfusion 
perfusion_csv_with_ata_cp <- subset(perfusion_csv, !(participant_id %in% values_to_exclude))
perfusion_csv_cp <- subset(perfusion_csv_with_ata_cp, !(participant_id %in% vascular_artifacts))

perfusion_csv_with_ata_cp$participant_id <- gsub("^sub-|_1$", "", perfusion_csv_with_ata_cp$participant_id)
perfusion_csv_cp$participant_id <- gsub("^sub-|_1$", "", perfusion_csv_cp$participant_id)

write.csv(perfusion_csv_with_ata_cp, paste0(path_results_data_analysis, 'perfusion_csv_with_ATA.csv'))
write.csv(perfusion_csv_cp, paste0(path_results_data_analysis, 'perfusion_csv_without_ATA.csv'))


# Unir con variables necesarias de cuestionarios
cols_to_merge <- c("Edad", "Grupo", "ID", "Genero", "BMI", "FechaCognitivo", tests_variables)
cols_to_merge <- intersect(cols_to_merge, colnames(cuestionarios_excel))  

merged_df_perfusion <- merge(
  cuestionarios_excel[, cols_to_merge], 
  perfusion_csv, 
  by.x = "ID", 
  by.y = "participant_id", 
  all.y = TRUE
) %>% 
  slice(-1) %>%  # Eliminar la primera fila
  mutate(across(all_of(tests_variables), ~ suppressWarnings(as.numeric(.)))) %>%  # Variables tests
  mutate(across(all_of(columns_regions_perfusion), ~ suppressWarnings(as.numeric(.))))  # Regiones %>%

cols_to_merge_gm <- c('participant_id','TotalGM_B', 'TotalGM_L', 'TotalGM_R')

merged_df_perfusion <- merge(
  merged_df_perfusion, 
  gray_matter_csv[,cols_to_merge_gm], 
  by.x = "ID", 
  by.y = "participant_id", 
  all.y = TRUE
) 


merged_df_perfusion_with_ata <- subset(merged_df_perfusion, !(ID %in% values_to_exclude))
merged_df_perfusion <- subset(merged_df_perfusion_with_ata, !(ID %in% vascular_artifacts)) 
merged_df_perfusion <- merged_df_perfusion[-1, ]



df_covid <- merged_df_perfusion %>% filter(Grupo == "COVID")

table(merged_df_perfusion$Grupo)
table(merged_df_perfusion_with_ata$Grupo)

colnames(merged_df_perfusion) <- make.names(colnames(merged_df_perfusion), unique = TRUE)

# merged_df_perfusion_sin_24 <- merged_df_perfusion_with_ata %>%
#   filter(!grepl("24$", FechaCognitivo))

# gam_analysis(merged_df_perfusion, columns_regions_perfusion)
# multiple_ancova(merged_df_perfusion, columns_regions_perfusion)
# kruskal_wallis_analysis(merged_df_perfusion, tests_variables)

# gam_analysis(merged_df_perfusion_with_ata, columns_regions_perfusion)
# multiple_ancova(merged_df_perfusion_with_ata, columns_regions_perfusion)
# kruskal_wallis_analysis(merged_df_perfusion_with_ata, tests_variables)
multiple_ancova_correcting_scov(merged_df_perfusion_with_ata, columns_regions_perfusion)


# Partial Correlation
# partial_correlation_tests(merged_df_perfusion, columns_regions_perfusion, tests_variables)
# partial_correlation_tests(merged_df_perfusion_with_ata, columns_regions_perfusion, tests_variables)
# partial_correlation_tests(df_covid, columns_regions_perfusion, tests_variables)

  # boxplots_paper <- list()
  # #
  # #
  # boxplots_paper[[1]] <- ggplot(merged_df_perfusion_sin_24, aes(x = Grupo, y =  PL_supramarginal_gyrus_B , fill = Grupo)) +
  #   geom_boxplot(alpha = 0.5) +  # Add transparency to boxplot
  #   geom_jitter(color = "black", size = 2, alpha = 0.3, width = 0.2) +  # Add jittered points
  #   labs(title = "PL_supramarginal_gyrus_B",
  #        x = "Group",
  #        y = "sCoV") +
  #   theme_minimal() +
  #   guides(fill = FALSE) +# Remove legend
  #   theme(
  #     plot.title = element_text(size = 10, hjust=0.5),  # Cambia el tamaño del título
  #   )
  # 
  # 
  # boxplots_paper[[2]] <- ggplot(merged_df_perfusion_sin_24, aes(x = Grupo, y =  PL_supramarginal_gyrus_B, fill = Grupo)) +
  #   geom_boxplot(alpha = 0.5) +  # Add transparency to boxplot
  #   geom_jitter(color = "black", size = 2, alpha = 0.3, width = 0.2) +  # Add jittered points
  #   geom_text(aes(label = ID), position = position_jitter(width = 0.2, height = 0),
  #                vjust = -1, size = 2.5, alpha = 0.8) +
  #   labs(title = "PL_supramarginal_gyrus_B",
  #        x = "Group",
  #        y = "sCoV") +
  #   theme_minimal() +
  #   guides(fill = FALSE) +# Remove legend
  #   theme(
  #     plot.title = element_text(size = 10, hjust=0.5),  # Cambia el tamaño del título
  #   )
  # ggsave('postcentral_supramarginal_new.jpg', bg="white")

  # png(paste(path_plots_data_analysis, "boxplot_gm_cbf_scov.png",sep=""), width = 12, height = 6, units = "in", res = 300)
  # grid_arrange <- do.call(gridExtra::grid.arrange, c(boxplots_paper, ncol = 2))
  # dev.off()

merged_df_perfusion$TotalGM_L <- as.numeric(as.character(merged_df_perfusion$TotalGM_L))

plot_gm_l <- ggplot(merged_df_perfusion, aes(x = Grupo, y = get("TotalGM_L"), fill = Grupo)) +
  geom_boxplot(alpha = 0.5) +  # Add transparency to boxplot
  geom_jitter(color = "black", size = 2, alpha = 0.3, width = 0.2) +  # Add jittered points
  labs(title = "Left Grey Matter", x = "Group", y = "CBF (mL/100g/min)")+
  theme_minimal() +
  guides(fill = FALSE) +  # Remove legend
  theme(
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),  # Increase title size
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )

# Save the plot
plot_filename <- paste(path_plots_data_analysis, "Total_GM_L_all.jpg", sep="")
ggsave(plot_filename, plot = plot_gm_l , device = "jpeg", bg = "white", width = 4, height = 7)

merged_df_perfusion$TotalGM_R <- as.numeric(as.character(merged_df_perfusion$TotalGM_R))

plot_gm_r <- ggplot(merged_df_perfusion, aes(x = Grupo, y = get("TotalGM_R"), fill = Grupo)) +
  geom_boxplot(alpha = 0.5) +  # Add transparency to boxplot
  geom_jitter(color = "black", size = 2, alpha = 0.3, width = 0.2) +  # Add jittered points
  labs(title = "Right Grey Matter", x = "Group", y = "CBF (mL/100g/min)") +
  theme_minimal() +
  guides(fill = FALSE) +  # Remove legend
  theme(
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),  # Increase title size
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )

# Save the plot
plot_filename <- paste(path_plots_data_analysis, "Total_GM_R_all.jpg", sep="")
ggsave(plot_filename, plot = plot_gm_r , device = "jpeg", bg = "white", width = 4, height = 7)
