# Load necessary libraries
library(tidyverse)
library(readxl)
library(car)
library(gridExtra)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(effects)
library(stats)
library(corrplot) 
library(patchwork)


setwd("/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Freesurfer/DataAnalysis/CodeR")
source("functionsANCOVA.R")

# Define the paths
path_data_analysis <- "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Freesurfer/DataAnalysis/"
#path_data_analysis <- "/home/solcat/work/data_analysis_paper/"
path_results_data_analysis <- paste(path_data_analysis, "Results/", sep="")
path_plots_data_analysis <-  paste(path_data_analysis, "Plots/", sep="")

path_cuestionarios <- paste(path_data_analysis, "ResumenRespuestas_14_04.csv", sep="")
path_sienax2 <-paste(path_data_analysis, "SienaxResults2.csv", sep="")
path_brain_vol <- paste(path_data_analysis, "brain_volumes.csv", sep="")
path_segmentation <- paste(path_data_analysis, "segmentation.csv", sep="")
path_parcellation <- paste(path_data_analysis, "parcellation.csv", sep="")
path_parcellation_DKT <- paste(path_data_analysis, "parcellation_DKT.csv", sep="")
path_parcellation_thick <- paste(path_data_analysis, "parcellation_thick.csv", sep="")
path_parcellation_thick_dkt <- paste(path_data_analysis, "parcellation_thick_DKT.csv", sep="")

path_bianca_file <- "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/WMH/ProcessedWMH/TotalPvAndDeepValuesBiancaThr0_7.csv"
path_socioeconomics_cardiovascular_index <- "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Freesurfer/DataAnalysis/ASL/SocialAndCardiovascularIndex.csv"

# Exclude rows
values_to_exclude <- c('CP0011', 'CP0015', 'CP0035', 'CP0144', 'CP0106', 'CP0192', 'CP0193', 'CP0087','CP0196')

tests_variables_perc <- c("EQ.VAS", "FAS", "PSQI",
                     "TMT.A_perc",
                     "TMT.B_perc",
                     "WMS.R.DIR_perc",
                     "WMS.R.INV_perc",
                     "STROOP_P_perc",
                     "STROOP_C_perc",
                     "STROOP_P.C_perc",
                     "STROOP_P.C_INTERF_perc",
                     "MOCA")



cuestionarios_excel <- read.csv(path_cuestionarios)
cuestionarios_excel <- subset(cuestionarios_excel, !(ID %in% values_to_exclude))





# Read files (Fresesurfer)
sienax2_csv <-read.csv(path_sienax2,dec = ",")
brain_vol_csv <- read.csv(path_brain_vol)
segmentation_csv <- read.csv(path_segmentation)
parcellation_csv <- read.csv(path_parcellation)
parcellation_DKT_csv <- read.csv(path_parcellation_DKT)
parcellation_thick_csv <- read.csv(path_parcellation_thick)
parcellation_thick_DKT_csv <- read.csv(path_parcellation_thick_dkt)

# Read WMH
bianca_csv <-read.csv(path_bianca_file)
colnames(bianca_csv) <- c("ID", "Total_WMH", "PV_WMH", "Deep_WMH")

socioeconomics_cardiovascular_index_csv <- read.csv(path_socioeconomics_cardiovascular_index)

# Read CBF

vascular_artifacts <- c('CP0031','CP0044', 'CP0053', 'CP0054', 'CP0061', 'CP0096','CP0107','CP0108','CP0118',
                        'CP0123', 'CP0142', 'CP0147', 'CP0154','CP0167', 'CP0176','CP0178',
                        'CP0180', 'CP0183', 'CP0188','CP0193','CP0199', 'CP0205', 'CP0213',
                        'CP0229', 'CP0238','CP0245', 'CP0247', 'CP0238')


path_gray_matter_perfusion <- "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/ASLProcessedCOVID/derivatives/ExploreASL/Population/Stats/mean_qCBF_StandardSpace_TotalGM_n=203_15-Jan-2025_PVC2.tsv"
perfusion_csv <- read.delim(path_gray_matter_perfusion)
perfusion_csv <- perfusion_csv[-1, ]
perfusion_csv$participant_id <- gsub("^sub-|_1$", "", perfusion_csv$participant_id)
perfusion_csv[perfusion_csv == "n/a"] <- NA
perfusion_csv$TotalGM_B <- as.numeric(perfusion_csv$TotalGM_B)

# NA in Artifacts
perfusion_csv$"TotalGM_B"[perfusion_csv$participant_id %in% vascular_artifacts] <- NA

# Read sCOV
# Definir el path al archivo de perfusión de sCOV sin PVC
path_scov_without_pvc <- "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/ASLProcessedCOVID/derivatives/ExploreASL/Population/Stats/CoV_qCBF_StandardSpace_TotalGM_n=203_15-Jan-2025_PVC0.tsv"

perfusion_scov_without_pvc <- read.delim(path_scov_without_pvc)
perfusion_scov_without_pvc <- perfusion_scov_without_pvc[-1, ]
perfusion_scov_without_pvc$participant_id <- gsub("^sub-|_1$", "", perfusion_scov_without_pvc$participant_id)
perfusion_scov_without_pvc[perfusion_scov_without_pvc == "n/a"] <- NA

perfusion_scov_without_pvc$TotalGM_B <- as.numeric(perfusion_scov_without_pvc$TotalGM_B)


# Signficative Freesurfer
regions_segmentation <- c('Left.Cerebellum.Cortex', 'Right.Cerebellum.Cortex')
right_regions_parcellation <- c('frontalpole', 'superiorfrontal',
                                'rostralmiddlefrontal', 'parsopercularis', 'fusiform', 'inferiortemporal' )
left_regions_parcellation <- c('lingual')

right_regions_parcellation_thick <- c('fusiform', 'postcentral')
left_regions_parcellation_thick <- c('lingual', 'postcentral', 'supramarginal')

# Cuestionarios (proteinas + tmt a + moca)
proteins_columns_select <- c('M6a','BDNF', 'NFL')
test_select <- c('MOCA_perc', 'TMT.A_perc', 'FAS', 'DepresionEscalaUnificada', 'AnsiedadEscalaUnificada')

# # Filtro los del grupo COVID
# cuestionarios_excel <- cuestionarios_excel %>%
#   filter(Grupo == "COVID")


cuestionarios_excel$M6a <- as.numeric(gsub(",", ".", as.character(cuestionarios_excel$M6a)))
cuestionarios_excel$BDNF <- as.numeric(gsub(",", ".", as.character(cuestionarios_excel$BDNF)))
cuestionarios_excel$NFL <- as.numeric(gsub(",", ".", as.character(cuestionarios_excel$NFL)))
cuestionarios_excel$MOCA_perc <- as.numeric(gsub(",", ".", as.character(cuestionarios_excel$MOCA_perc)))


matrix_kmeans <- cuestionarios_excel[,c('ID', proteins_columns_select, test_select,'Grupo')]

# Merge WMH
matrix_kmeans <- merge(matrix_kmeans, bianca_csv, by = "ID", all.x = TRUE)

# Merge sCOV
matrix_kmeans <- merge(matrix_kmeans, perfusion_scov_without_pvc[,c('participant_id','TotalGM_B')], by.x = "ID", 
                       by.y = "participant_id", 
                       all.x = TRUE)

names(matrix_kmeans)[names(matrix_kmeans) == 'TotalGM_B'] <- 'TotalGM_sCOV'

# Merge CBF
matrix_kmeans <- merge(matrix_kmeans, perfusion_csv[,c('participant_id','TotalGM_B')], by.x = "ID",
                       by.y = "participant_id",
                       all.x = TRUE)

names(matrix_kmeans)[names(matrix_kmeans) == 'TotalGM_B'] <- 'TotalGM_CBF'

# Merge Cardiovascular INDEX
matrix_kmeans <- merge(matrix_kmeans, socioeconomics_cardiovascular_index_csv[, c('ID','CardiovascularRiskIndex')], by = "ID", all.x = TRUE)

# Freesurfer
# Right Hemisphere Volume
right_hemisphere_vol_data <- parcellation_DKT_csv %>%
  filter(Hemisphere == 'rh') %>%
  dplyr::select(subject, all_of(right_regions_parcellation)) %>%
  rename_with(~ paste0(., "_vol_rh"), all_of(right_regions_parcellation))

# Left Hemisphere Volume
left_hemisphere_vol_data <- parcellation_DKT_csv %>%
  filter(Hemisphere == 'lh') %>%
  dplyr::select(subject, all_of(left_regions_parcellation)) %>%
  rename_with(~ paste0(., "_vol_lh"), all_of(left_regions_parcellation))

# Right Hemisphere thick
right_hemisphere_thick_data <- parcellation_thick_DKT_csv %>%
  filter(Hemisphere == 'rh') %>%
  dplyr::select(subject, all_of(right_regions_parcellation_thick)) %>%
  rename_with(~ paste0(., "_thick_rh"), all_of(right_regions_parcellation_thick))

# Left Hemisphere thick
left_hemisphere_thick_data <- parcellation_thick_DKT_csv %>%
  filter(Hemisphere == 'lh') %>%
  dplyr::select(subject, all_of(left_regions_parcellation_thick)) %>%
  rename_with(~ paste0(., "_thick_lh"), all_of(left_regions_parcellation_thick))

# Cerebellum
cerebellum_data <- segmentation_csv %>%
  dplyr::select(subject, all_of(regions_segmentation)) 

# Combinar los datos del hemisferio derecho e izquierdo
all_freesurfer <- right_hemisphere_vol_data %>%
  full_join(left_hemisphere_vol_data, by = "subject") %>%
  full_join(right_hemisphere_thick_data, by = "subject") %>%
  full_join(left_hemisphere_thick_data, by = "subject") %>%
  full_join(cerebellum_data, by = "subject")

# Merge Freesurfer
matrix_kmeans <- merge(matrix_kmeans, all_freesurfer, by.x = "ID", by.y = "subject", all.x = TRUE)

# 
# # Correlate FS data 
# cor_matrix_freesurfer <- cor(all_freesurfer[ , colnames(all_freesurfer) != "subject"], use = "pairwise.complete.obs")
# 
# output_path_cor_matrix_fs <- paste0(path_plots_data_analysis, "freesurfer_correlation_tests.png")
# 
# # Guardar gráfico unificado con colores pastel
# png(output_path_cor_matrix_fs, width=800, height=800)
# 
# corrplot(
#   cor_matrix_freesurfer , method="circle", type="lower", order="original", 
#   col=colorRampPalette(c("lightblue", "white", "lightcoral"))(200),
#   tl.col="black", tl.srt=45, addCoef.col="black")
# dev.off()

# Sintomas cognitivos
# Definir las columnas cognitivas
# cols_cognitivas <- c("Fatiga", "Confusion", "Comunicacion", "Memoria", "Atencion")
# 
# # Crear nueva matriz binarizada
# matriz_sintomas_cognitivos <- cuestionarios_excel[, c("ID", cols_cognitivas)]
# matriz_sintomas_cognitivos[cols_cognitivas] <- lapply(matriz_sintomas_cognitivos[cols_cognitivas],
#                                                       function(x) ifelse(is.na(x), 0, ifelse(x == 1, 1, 0)))
# 
# # Agregar índice a la nueva matriz
# matriz_sintomas_cognitivos$Indice_Sintomas_Cognitivos <- rowSums(matriz_sintomas_cognitivos[, cols_cognitivas], na.rm = TRUE)
# 
# # Unir
# matrix_kmeans <- merge(matrix_kmeans, matriz_sintomas_cognitivos, by.x = "ID", by.y = "ID", all.x = TRUE)

bianca_names <- setdiff(colnames(bianca_csv), "ID")
freesurfer_names <- setdiff(colnames(all_freesurfer), "subject")

exclude_vars <- c(bianca_names, freesurfer_names, "CardiovascularRiskIndex", 
                  "MOCA_perc", "TMT.A_perc", "FAS", 
                  "TotalGM_sCOV")


# # 1. Clusterizar BDNF + M6A
final_matrix_all_BDNF_M6A <- cluster_data_flexible_v2(data = matrix_kmeans,
                                             group_filter = NULL,    # Si es NULL, usa todo
                                             value_to_exclude = c(exclude_vars, "TotalGM_CBF"),
                                             proteins = c("NFL"),
                                             centers = 2, 
                                             nstart = 25) 


final_matrix_all_BDNF_M6A_sin_escalas <- cluster_data_flexible_v2(data = matrix_kmeans,
                                                                  group_filter = NULL,    # Si es NULL, usa todo
                                                                  value_to_exclude = c(exclude_vars, "TotalGM_CBF", "DepresionEscalaUnificada", "AnsiedadEscalaUnificada"),
                                                                  proteins = c("NFL"),
                                                                  centers = 2, 
                                                                  nstart = 25) 

# 
# variables_to_plot <- setdiff(names(final_matrix_all_BDNF_M6A), c("ID", "Grupo", "Cluster", proteins, colnames(bianca_csv)))
# variables_to_plot_2 <- setdiff(names(final_matrix_all_BDNF_M6A), c("ID", "Grupo", "Cluster")) 


matrix_kmeans_no_na_BDNF_M6A <- matrix_kmeans %>%
  dplyr::filter(!is.na(final_matrix_all_BDNF_M6A$merged_data$Cluster))

matrix_kmeans_no_na_BDNF_M6A_sin_escalas <- matrix_kmeans %>%
  dplyr::filter(!is.na(final_matrix_all_BDNF_M6A_sin_escalas$merged_data$Cluster))

final_matrix_covid_plot_BDNF_M6A <- cluster_data_flexible_v2(data = matrix_kmeans_no_na_BDNF_M6A,
                                                    group_filter = "COVID",    # Si es NULL, usa todo
                                                    value_to_exclude = c(exclude_vars, "TotalGM_CBF"),
                                                    proteins = c("NFL"),
                                                    centers = 2, 
                                                    nstart = 25) 


final_matrix_covid_plot_BDNF_M6A_sin_escalas <- cluster_data_flexible_v2(data = matrix_kmeans_no_na_BDNF_M6A_sin_escalas,
                                                             group_filter = "COVID",    # Si es NULL, usa todo
                                                             value_to_exclude = c(exclude_vars, "TotalGM_CBF", "DepresionEscalaUnificada", "AnsiedadEscalaUnificada" ),
                                                             proteins = c("NFL"),
                                                             centers = 2, 
                                                             nstart = 25) 



final_matrix_control_plot_BDNF_M6A <- cluster_data_flexible_v2(data = matrix_kmeans_no_na_BDNF_M6A,
                                                      group_filter = "CONTROL",    # Si es NULL, usa todo
                                                      value_to_exclude = c(exclude_vars, "TotalGM_CBF"),
                                                      proteins = c("NFL"),
                                                      centers = 2, 
                                                      nstart = 25) 


# BDNF + M6A COVID 
n_total_covid <- nrow(final_matrix_covid_plot_BDNF_M6A$clustering_data_clean)
cat("Total de sujetos analizados (COVID):", n_total_covid, "\n")

n_total_control <- nrow(final_matrix_control_plot_BDNF_M6A$clustering_data_clean)
cat("Total de sujetos analizados (Control):", n_total_control, "\n")

n_total_all <- nrow(final_matrix_all_BDNF_M6A$clustering_data_clean)
cat("Total de sujetos analizados (All):", n_total_all, "\n")

n_total_all_sin_escalas <- nrow(final_matrix_all_BDNF_M6A_sin_escalas$clustering_data_clean)
cat("Total de sujetos analizados sin escalas (All):", n_total_all_sin_escalas, "\n")


# Find ID de diferencia 
ids_validos_con_escalas <- final_matrix_all_BDNF_M6A$merged_data %>%
  filter(!is.na(Cluster)) %>%
  pull(ID)  # reemplazá 'subject_id' con el nombre real de tu columna ID

# Para matriz sin escalas
ids_validos_sin_escalas <- final_matrix_all_BDNF_M6A_sin_escalas$merged_data %>%
  filter(!is.na(Cluster)) %>%
  pull(ID)

setdiff(ids_validos_sin_escalas, ids_validos_con_escalas)


final_matrix_covid_plot_BDNF_M6A$merged_data %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(n = n())

final_matrix_covid_plot_BDNF_M6A_sin_escalas$merged_data %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(n = n())


crear_panel_covid(
  data = final_matrix_covid_plot_BDNF_M6A$merged_data,
  filename_prefix = "panel_COVID_BDNF_M6A",
  panel_title = "COVID - M6A + BDNF - Clinical and Protein Markers"
)

# Sin Escalas
crear_panel_covid(
  data = final_matrix_covid_plot_BDNF_M6A_sin_escalas$merged_data,
  filename_prefix = "panel_COVID_BDNF_M6A_sin_escalas",
  panel_title = "COVID - M6A + BDNF - Protein Markers"
)

# BDNF + M6A All
n_total_covid <- nrow(final_matrix_all_BDNF_M6A$clustering_data_clean)
cat("Total de sujetos analizados:", n_total_covid, "\n")

final_matrix_all_BDNF_M6A$merged_data <- final_matrix_all_BDNF_M6A$merged_data %>%
  dplyr::filter(!is.na(Cluster))

final_matrix_all_BDNF_M6A_sin_escalas$merged_data <- final_matrix_all_BDNF_M6A_sin_escalas$merged_data %>%
  dplyr::filter(!is.na(Cluster))

final_matrix_all_BDNF_M6A$merged_data %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(n = n())


final_matrix_all_BDNF_M6A_sin_escalas$merged_data %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(n = n())

crear_panel_covid(
  data = final_matrix_all_BDNF_M6A$merged_data,
  filename_prefix = "panel_ALL_BDNF_M6A",
  panel_title = "ALL - M6A + BDNF - Clinical and Protein Markers"
)


crear_panel_covid(
  data = final_matrix_all_BDNF_M6A_sin_escalas$merged_data,
  filename_prefix = "panel_ALL_BDNF_M6A_sin_escalas",
  panel_title = "ALL - M6A + BDNF - Protein Markers"
)



# # 1. Clusterizar BDNF
final_matrix_all_BDNF <- cluster_data_flexible_v2(data = matrix_kmeans,
                                                      group_filter = NULL,    # Si es NULL, usa todo
                                                      value_to_exclude = c(exclude_vars, "TotalGM_CBF"),
                                                      proteins = c("NFL", "M6a"),
                                                      centers = 2, 
                                                      nstart = 25) 

matrix_kmeans_no_na_BDNF <- matrix_kmeans %>%
  dplyr::filter(!is.na(final_matrix_all_BDNF$merged_data$Cluster))

final_matrix_covid_plot_BDNF <- cluster_data_flexible_v2(data = matrix_kmeans_no_na_BDNF,
                                                             group_filter = "COVID",    # Si es NULL, usa todo
                                                             value_to_exclude = c(exclude_vars, "TotalGM_CBF"),
                                                             proteins = c("NFL", "M6a"),
                                                             centers = 2, 
                                                             nstart = 25) 

final_matrix_control_plot_BDNF <- cluster_data_flexible_v2(data = matrix_kmeans_no_na_BDNF,
                                                               group_filter = "CONTROL",    # Si es NULL, usa todo
                                                               value_to_exclude = c(exclude_vars, "TotalGM_CBF"),
                                                               proteins = c("NFL","M6a"),
                                                               centers = 2, 
                                                               nstart = 25) 

# BDNF 
n_total_covid <- nrow(final_matrix_covid_plot_BDNF$clustering_data_clean)
cat("Total de sujetos analizados (COVID):", n_total_covid, "\n")

n_total_control <- nrow(final_matrix_control_plot_BDNF$clustering_data_clean)
cat("Total de sujetos analizados (Control):", n_total_control, "\n")

n_total_all <- nrow(final_matrix_all_BDNF$clustering_data_clean)
cat("Total de sujetos analizados (All):", n_total_all, "\n")

final_matrix_covid_plot_BDNF$merged_data %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(n = n())

crear_panel_covid(
  data = final_matrix_covid_plot_BDNF$merged_data,
  filename_prefix = "panel_COVID_BDNF",
  panel_title = "COVID - BDNF - Clinical and Protein Markers"
)

# BDNF - All
n_total_all <- nrow(final_matrix_all_BDNF$clustering_data_clean)
cat("Total de sujetos analizados:", n_total_all, "\n")

final_matrix_all_BDNF$merged_data <- final_matrix_all_BDNF$merged_data %>%
  dplyr::filter(!is.na(Cluster))

final_matrix_all_BDNF$merged_data %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(n = n())


crear_panel_covid(
  data = final_matrix_all_BDNF$merged_data,
  filename_prefix = "panel_ALL_BDNF",
  panel_title = "ALL - BDNF - Clinical and Protein Markers"
)




# # 1. Clusterizar M6a
final_matrix_all_M6a <- cluster_data_flexible_v2(data = matrix_kmeans,
                                                  group_filter = NULL,    # Si es NULL, usa todo
                                                  value_to_exclude = c(exclude_vars, "TotalGM_CBF"),
                                                  proteins = c("NFL", "BDNF"),
                                                  centers = 2, 
                                                  nstart = 25) 

matrix_kmeans_no_na_M6a <- matrix_kmeans %>%
  dplyr::filter(!is.na(final_matrix_all_M6a$merged_data$Cluster))

final_matrix_covid_plot_M6a <- cluster_data_flexible_v2(data = matrix_kmeans_no_na_M6a,
                                                         group_filter = "COVID",    # Si es NULL, usa todo
                                                         value_to_exclude = c(exclude_vars, "TotalGM_CBF"),
                                                         proteins = c("NFL", "BDNF"),
                                                         centers = 2, 
                                                         nstart = 25) 

final_matrix_control_plot_M6a <- cluster_data_flexible_v2(data = matrix_kmeans_no_na_BDNF,
                                                           group_filter = "CONTROL",    # Si es NULL, usa todo
                                                           value_to_exclude = c(exclude_vars, "TotalGM_CBF"),
                                                           proteins = c("NFL","BDNF"),
                                                           centers = 2, 
                                                           nstart = 25) 

# M6a 
n_total_covid <- nrow(final_matrix_covid_plot_M6a$clustering_data_clean)
cat("Total de sujetos analizados (COVID):", n_total_covid, "\n")

n_total_control <- nrow(final_matrix_control_plot_M6a$clustering_data_clean)
cat("Total de sujetos analizados (CONTROL):", n_total_control, "\n")

n_total_all <- nrow(final_matrix_all_M6a$clustering_data_clean)
cat("Total de sujetos analizados (All):", n_total_all, "\n")

final_matrix_covid_plot_M6a$merged_data %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(n = n())

crear_panel_covid(
  data = final_matrix_covid_plot_M6a$merged_data,
  filename_prefix = "panel_COVID_M6a",
  panel_title = "COVID - M6a - Clinical and Protein Markers"
)

# M6a - All
n_total_all <- nrow(final_matrix_all_M6a$clustering_data_clean)
cat("Total de sujetos analizados:", n_total_all, "\n")

final_matrix_all_M6a$merged_data <- final_matrix_all_M6a$merged_data %>%
  dplyr::filter(!is.na(Cluster))

final_matrix_all_M6a$merged_data %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(n = n())


crear_panel_covid(
  data = final_matrix_all_M6a$merged_data,
  filename_prefix = "panel_ALL_M6a",
  panel_title = "ALL - M6a - Clinical and Protein Markers"
)




# # 1. Clusterizar NFL
final_matrix_all_NFL <- cluster_data_flexible_v2(data = matrix_kmeans,
                                                 group_filter = NULL,    # Si es NULL, usa todo
                                                 value_to_exclude = c(exclude_vars, "TotalGM_CBF"),
                                                 proteins = c("M6a", "BDNF"),
                                                 centers = 2, 
                                                 nstart = 25) 

matrix_kmeans_no_na_NFL <- matrix_kmeans %>%
  dplyr::filter(!is.na(final_matrix_all_NFL$merged_data$Cluster))

final_matrix_covid_plot_NFL <- cluster_data_flexible_v2(data = matrix_kmeans_no_na_NFL,
                                                        group_filter = "COVID",    # Si es NULL, usa todo
                                                        value_to_exclude = c(exclude_vars, "TotalGM_CBF"),
                                                        proteins = c("M6a", "BDNF"),
                                                        centers = 2, 
                                                        nstart = 25) 

final_matrix_control_plot_NFL <- cluster_data_flexible_v2(data = matrix_kmeans_no_na_NFL,
                                                          group_filter = "CONTROL",    # Si es NULL, usa todo
                                                          value_to_exclude = c(exclude_vars, "TotalGM_CBF"),
                                                          proteins = c("M6a","BDNF"),
                                                          centers = 2, 
                                                          nstart = 25) 

# NFL 
n_total_covid <- nrow(final_matrix_covid_plot_NFL$clustering_data_clean)
cat("Total de sujetos analizados (COVID):", n_total_covid, "\n")

n_total_control <- nrow(final_matrix_control_plot_NFL$clustering_data_clean)
cat("Total de sujetos analizados (Control):", n_total_control, "\n")

n_total_all <- nrow(final_matrix_all_NFL$clustering_data_clean)
cat("Total de sujetos analizados (All):", n_total_all, "\n")


final_matrix_covid_plot_NFL$merged_data %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(n = n())

crear_panel_covid(
  data = final_matrix_covid_plot_NFL$merged_data,
  filename_prefix = "panel_COVID_NFL",
  panel_title = "COVID - NFL - Clinical and Protein Markers"
)

# NFL - All
n_total_all <- nrow(final_matrix_all_NFL$clustering_data_clean)
cat("Total de sujetos analizados:", n_total_all, "\n")

final_matrix_all_NFL$merged_data <- final_matrix_all_NFL$merged_data %>%
  dplyr::filter(!is.na(Cluster))

final_matrix_all_NFL$merged_data %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(n = n())


crear_panel_covid(
  data = final_matrix_all_NFL$merged_data,
  filename_prefix = "panel_ALL_NFL",
  panel_title = "ALL - NFL - Clinical and Protein Markers"
)


## QUITANDO OUTLIERS CP0223 Y CP0224

matrix_kmeans_new <- matrix_kmeans[!(matrix_kmeans$ID %in% c("CP0223", "CP0224")), ]


# # 1. Clusterizar BDNF + M6A
final_matrix_all_BDNF_M6A <- cluster_data_flexible_v2(data = matrix_kmeans_new,
                                             group_filter = NULL,    # Si es NULL, usa todo
                                             value_to_exclude = c(exclude_vars, "TotalGM_CBF"),
                                             proteins = c("NFL"),
                                             centers = 2, 
                                             nstart = 25) 

final_matrix_all_BDNF_M6A_sin_escalas <- cluster_data_flexible_v2(data = matrix_kmeans_new,
                                                      group_filter = NULL,    # Si es NULL, usa todo
                                                      value_to_exclude = c(exclude_vars, "TotalGM_CBF", "AnsiedadEscalaUnificada", "DepresionEscalaUnificada"),
                                                      proteins = c("NFL"),
                                                      centers = 2, 
                                                      nstart = 25) 


# 
# variables_to_plot <- setdiff(names(final_matrix_all_BDNF_M6A), c("ID", "Grupo", "Cluster", proteins, colnames(bianca_csv)))
# variables_to_plot_2 <- setdiff(names(final_matrix_all_BDNF_M6A), c("ID", "Grupo", "Cluster")) 


matrix_kmeans_no_na_BDNF_M6A <- matrix_kmeans_new %>%
  dplyr::filter(!is.na(final_matrix_all_BDNF_M6A$merged_data$Cluster))


matrix_kmeans_no_na_BDNF_M6A_sin_escalas <- matrix_kmeans_new %>%
  dplyr::filter(!is.na(final_matrix_all_BDNF_M6A_sin_escalas$merged_data$Cluster))

final_matrix_covid_plot_BDNF_M6A <- cluster_data_flexible_v2(data = matrix_kmeans_no_na_BDNF_M6A,
                                                    group_filter = "COVID",    # Si es NULL, usa todo
                                                    value_to_exclude = c(exclude_vars, "TotalGM_CBF"),
                                                    proteins = c("NFL"),
                                                    centers = 2, 
                                                    nstart = 25)

final_matrix_covid_plot_BDNF_M6A_sin_escalas <- cluster_data_flexible_v2(data = matrix_kmeans_no_na_BDNF_M6A_sin_escalas,
                                                             group_filter = "COVID",    # Si es NULL, usa todo
                                                             value_to_exclude = c(exclude_vars, "TotalGM_CBF", "AnsiedadEscalaUnificada", "DepresionEscalaUnificada"),
                                                             proteins = c("NFL"),
                                                             centers = 2, 
                                                             nstart = 25)



final_matrix_control_plot_BDNF_M6A <- cluster_data_flexible_v2(data = matrix_kmeans_no_na_BDNF_M6A,
                                                      group_filter = "CONTROL",    # Si es NULL, usa todo
                                                      value_to_exclude = c(exclude_vars, "TotalGM_CBF"),
                                                      proteins = c("NFL"),
                                                      centers = 2, 
                                                      nstart = 25) 


# BDNF + M6A COVID 
n_total_covid <- nrow(final_matrix_covid_plot_BDNF_M6A$clustering_data_clean)
cat("Total de sujetos analizados (COVID):", n_total_covid, "\n")

n_total_control <- nrow(final_matrix_control_plot_BDNF_M6A$clustering_data_clean)
cat("Total de sujetos analizados (Control):", n_total_control, "\n")

n_total_all <- nrow(final_matrix_all_BDNF_M6A$clustering_data_clean)
cat("Total de sujetos analizados (All):", n_total_all, "\n")



final_matrix_covid_plot_BDNF_M6A$merged_data %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(n = n())

final_matrix_covid_plot_BDNF_M6A_sin_escalas$merged_data %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(n = n())

crear_panel_covid(
  data = final_matrix_covid_plot_BDNF_M6A$merged_data,
  filename_prefix = "panel_COVID_BDNF_M6A_without_outliers",
  panel_title = "COVID - M6A + BDNF - without outliers - Clinical and Protein Markers"
)

crear_panel_covid(
  data = final_matrix_covid_plot_BDNF_M6A_sin_escalas$merged_data,
  filename_prefix = "panel_COVID_BDNF_M6A_without_outliers_sin_escalas",
  panel_title = "COVID - M6A + BDNF - without outliers - Protein Markers"
)

crear_panel_cognitivas(data = final_matrix_covid_plot_BDNF_M6A$merged_data,
                       filename_prefix = "panel_COVID_BDNF_M6A_without_outliers_cognitivas",
                       panel_title = "COVID - M6A + BDNF - without outliers - Cognitive Tests")

crear_panel_cognitivas(data = final_matrix_covid_plot_BDNF_M6A_sin_escalas$merged_data,
                       filename_prefix = "panel_COVID_BDNF_M6A_without_outliers_cognitivas_sin_escalas",
                       panel_title = "COVID - M6A + BDNF - without outliers - Cognitive Tests")

crear_panel_WMH(data = final_matrix_covid_plot_BDNF_M6A$merged_data,
                filename_prefix = "panel_COVID_BDNF_M6A_without_outliers_wmh",
                panel_title = "COVID - M6A + BDNF - without outliers - WMH")

crear_panel_WMH(data = final_matrix_covid_plot_BDNF_M6A_sin_escalas$merged_data,
                filename_prefix = "panel_COVID_BDNF_M6A_without_outliers_wmh_sin_escalas",
                panel_title = "COVID - M6A + BDNF - without outliers - WMH")


crear_panel_ASL(data = final_matrix_covid_plot_BDNF_M6A$merged_data,
                filename_prefix = "panel_COVID_BDNF_M6A_without_outliers_ASL",
                panel_title = "COVID - M6A + BDNF - without outliers - ASL")

crear_panel_ASL(data = final_matrix_covid_plot_BDNF_M6A_sin_escalas$merged_data,
                filename_prefix = "panel_COVID_BDNF_M6A_without_outliers_ASL_sin_escalas",
                panel_title = "COVID - M6A + BDNF - without outliers - ASL")



# BDNF + M6A All
n_total_covid <- nrow(final_matrix_all_BDNF_M6A$clustering_data_clean)
cat("Total de sujetos analizados:", n_total_covid, "\n")


final_matrix_all_BDNF_M6A$merged_data <- final_matrix_all_BDNF_M6A$merged_data %>%
  dplyr::filter(!is.na(Cluster))


final_matrix_all_BDNF_M6A_sin_escalas$merged_data <- final_matrix_all_BDNF_M6A_sin_escalas$merged_data %>%
  dplyr::filter(!is.na(Cluster))

final_matrix_all_BDNF_M6A$merged_data %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(n = n())

final_matrix_all_BDNF_M6A_sin_escalas$merged_data %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(n = n())


crear_panel_covid(
  data = final_matrix_all_BDNF_M6A$merged_data,
  filename_prefix = "panel_ALL_BDNF_M6A_without_outliers",
  panel_title = "ALL - M6A + BDNF - without outliers - Clinical and Protein Markers"
)


crear_panel_covid(
  data = final_matrix_all_BDNF_M6A_sin_escalas$merged_data,
  filename_prefix = "panel_ALL_BDNF_M6A_without_outliers_sin_escalas",
  panel_title = "ALL - M6A + BDNF - without outliers - Protein Markers"
)


crear_panel_cognitivas(data = final_matrix_all_BDNF_M6A$merged_data,
                       filename_prefix = "panel_ALL_BDNF_M6A_without_outliers_cognitivas",
                       panel_title = "ALL - M6A + BDNF - without outliers - Cognitive Tests")


crear_panel_cognitivas(data = final_matrix_all_BDNF_M6A_sin_escalas$merged_data,
                       filename_prefix = "panel_ALL_BDNF_M6A_without_outliers_cognitivas_sin_escalas",
                       panel_title = "ALL - M6A + BDNF - without outliers - Cognitive Tests")


crear_panel_WMH(data = final_matrix_all_BDNF_M6A$merged_data,
                filename_prefix = "panel_ALL_BDNF_M6A_without_outliers_wmh",
                panel_title = "All - M6A + BDNF - without outliers - WMH")

crear_panel_WMH(data = final_matrix_all_BDNF_M6A_sin_escalas$merged_data,
                filename_prefix = "panel_ALL_BDNF_M6A_without_outliers_wmh_sin_escalas",
                panel_title = "All - M6A + BDNF - without outliers - WMH")


crear_panel_ASL(data = final_matrix_all_BDNF_M6A$merged_data,
                filename_prefix = "panel_ALL_BDNF_M6A_without_outliers_ASL",
                panel_title = "All - M6A + BDNF - without outliers - ASL")

crear_panel_ASL(data = final_matrix_all_BDNF_M6A_sin_escalas$merged_data,
                filename_prefix = "panel_ALL_BDNF_M6A_without_outliers_ASL_sin_escalas",
                panel_title = "All - M6A + BDNF - without outliers - ASL")


# # 1. Clusterizar BDNF
final_matrix_all_BDNF <- cluster_data_flexible_v2(data = matrix_kmeans_new,
                                                      group_filter = NULL,    # Si es NULL, usa todo
                                                      value_to_exclude = c(exclude_vars, "TotalGM_CBF"),
                                                      proteins = c("NFL", "M6a"),
                                                      centers = 2, 
                                                      nstart = 25) 

matrix_kmeans_no_na_BDNF <- matrix_kmeans_new %>%
  dplyr::filter(!is.na(final_matrix_all_BDNF$merged_data$Cluster))

final_matrix_covid_plot_BDNF <- cluster_data_flexible_v2(data = matrix_kmeans_no_na_BDNF,
                                                             group_filter = "COVID",    # Si es NULL, usa todo
                                                             value_to_exclude = c(exclude_vars, "TotalGM_CBF"),
                                                             proteins = c("NFL", "M6a"),
                                                             centers = 2, 
                                                             nstart = 25) 

final_matrix_control_plot_BDNF <- cluster_data_flexible_v2(data = matrix_kmeans_no_na_BDNF,
                                                               group_filter = "CONTROL",    # Si es NULL, usa todo
                                                               value_to_exclude = c(exclude_vars, "TotalGM_CBF"),
                                                               proteins = c("NFL","M6a"),
                                                               centers = 2, 
                                                               nstart = 25) 

# BDNF 
n_total_covid <- nrow(final_matrix_covid_plot_BDNF$clustering_data_clean)
cat("Total de sujetos analizados (COVID):", n_total_covid, "\n")

n_total_control <- nrow(final_matrix_control_plot_BDNF$clustering_data_clean)
cat("Total de sujetos analizados (Control):", n_total_control, "\n")

n_total_all <- nrow(final_matrix_all_BDNF$clustering_data_clean)
cat("Total de sujetos analizados (All):", n_total_all, "\n")



final_matrix_covid_plot_BDNF$merged_data %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(n = n())

crear_panel_covid(
  data = final_matrix_covid_plot_BDNF$merged_data,
  filename_prefix = "panel_COVID_BDNF_without_outliers",
  panel_title = "COVID - BDNF - without outliers - Clinical and Protein Markers"
)

# BDNF - All
n_total_all <- nrow(final_matrix_all_BDNF$clustering_data_clean)
cat("Total de sujetos analizados:", n_total_all, "\n")

final_matrix_all_BDNF$merged_data <- final_matrix_all_BDNF$merged_data %>%
  dplyr::filter(!is.na(Cluster))

final_matrix_all_BDNF$merged_data %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(n = n())


crear_panel_covid(
  data = final_matrix_all_BDNF$merged_data,
  filename_prefix = "panel_ALL_BDNF_without_outliers",
  panel_title = "ALL - BDNF - Without outliers - Clinical and Protein Markers"
)




# # 1. Clusterizar M6a
# Eliminar sujetos problemáticos
matrix_kmeans_filtered <- matrix_kmeans_new %>%
  dplyr::filter(!(ID %in% c("CP0223", "CP0224")))

# 1. Clusterizar M6a (ALL)
final_matrix_all_M6a <- cluster_data_flexible_v2(
  data = matrix_kmeans_filtered,
  group_filter = NULL,    
  value_to_exclude = c(exclude_vars, "TotalGM_CBF"),
  proteins = c("NFL", "BDNF"),
  centers = 2, 
  nstart = 25
)

# Quitar sujetos sin cluster
matrix_kmeans_no_na_M6a <- matrix_kmeans_filtered %>%
  dplyr::filter(!is.na(final_matrix_all_M6a$merged_data$Cluster))

# 2. Clusterizar M6a (COVID)
final_matrix_covid_plot_M6a <- cluster_data_flexible_v2(
  data = matrix_kmeans_no_na_M6a,
  group_filter = "COVID",
  value_to_exclude = c(exclude_vars, "TotalGM_CBF"),
  proteins = c("NFL", "BDNF"),
  centers = 2, 
  nstart = 25
)

# 3. Clusterizar M6a (CONTROL)
final_matrix_control_plot_M6a <- cluster_data_flexible_v2(
  data = matrix_kmeans_no_na_M6a,
  group_filter = "CONTROL",
  value_to_exclude = c(exclude_vars, "TotalGM_CBF"),
  proteins = c("NFL", "BDNF"),
  centers = 2, 
  nstart = 25
)


# COVID
n_total_covid <- nrow(final_matrix_covid_plot_M6a$clustering_data_clean)
cat("Total de sujetos analizados (COVID):", n_total_covid, "\n")

final_matrix_covid_plot_M6a$merged_data %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(n = n())

crear_panel_covid(
  data = final_matrix_covid_plot_M6a$merged_data,
  filename_prefix = "panel_COVID_M6a_without_outliers",
  panel_title = "COVID - M6a - without outliers - Clinical and Protein Markers"
)

# CONTROL
n_total_control <- nrow(final_matrix_control_plot_M6a$clustering_data_clean)
cat("Total de sujetos analizados (CONTROL):", n_total_control, "\n")

final_matrix_control_plot_M6a$merged_data %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(n = n())

# ALL
n_total_all <- nrow(final_matrix_all_M6a$clustering_data_clean)
cat("Total de sujetos analizados (ALL):", n_total_all, "\n")

final_matrix_all_M6a$merged_data <- final_matrix_all_M6a$merged_data %>%
  dplyr::filter(!is.na(Cluster))

final_matrix_all_M6a$merged_data %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(n = n())

crear_panel_covid(
  data = final_matrix_all_M6a$merged_data,
  filename_prefix = "panel_ALL_M6a_without_outliers",
  panel_title = "ALL - M6a - without outliers - Clinical and Protein Markers"
)




# # 1. Clusterizar NFL
final_matrix_all_NFL <- cluster_data_flexible_v2(data = matrix_kmeans_new,
                                                 group_filter = NULL,    # Si es NULL, usa todo
                                                 value_to_exclude = c(exclude_vars, "TotalGM_CBF"),
                                                 proteins = c("M6a", "BDNF"),
                                                 centers = 2, 
                                                 nstart = 25) 

matrix_kmeans_no_na_NFL <- matrix_kmeans_new %>%
  dplyr::filter(!is.na(final_matrix_all_NFL$merged_data$Cluster))

final_matrix_covid_plot_NFL <- cluster_data_flexible_v2(data = matrix_kmeans_no_na_NFL,
                                                        group_filter = "COVID",    # Si es NULL, usa todo
                                                        value_to_exclude = c(exclude_vars, "TotalGM_CBF"),
                                                        proteins = c("M6a", "BDNF"),
                                                        centers = 2, 
                                                        nstart = 25) 

final_matrix_control_plot_NFL <- cluster_data_flexible_v2(data = matrix_kmeans_no_na_NFL,
                                                          group_filter = "CONTROL",    # Si es NULL, usa todo
                                                          value_to_exclude = c(exclude_vars, "TotalGM_CBF"),
                                                          proteins = c("M6a","BDNF"),
                                                          centers = 2, 
                                                          nstart = 25) 

# NFL 
n_total_covid <- nrow(final_matrix_covid_plot_NFL$clustering_data_clean)
cat("Total de sujetos analizados (COVID):", n_total_covid, "\n")

n_total_control <- nrow(final_matrix_control_plot_NFL$clustering_data_clean)
cat("Total de sujetos analizados (Control):", n_total_control, "\n")

n_total_all <- nrow(final_matrix_all_NFL$clustering_data_clean)
cat("Total de sujetos analizados (All:", n_total_all, "\n")

final_matrix_covid_plot_NFL$merged_data %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(n = n())

crear_panel_covid(
  data = final_matrix_covid_plot_NFL$merged_data,
  filename_prefix = "panel_COVID_NFL_without_outliers",
  panel_title = "COVID - NFL - Without outliers -  Clinical and Protein Markers"
)

# NFL - All
n_total_all <- nrow(final_matrix_all_NFL$clustering_data_clean)
cat("Total de sujetos analizados:", n_total_all, "\n")

final_matrix_all_NFL$merged_data <- final_matrix_all_NFL$merged_data %>%
  dplyr::filter(!is.na(Cluster))

final_matrix_all_NFL$merged_data %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(n = n())


crear_panel_covid(
  data = final_matrix_all_NFL$merged_data,
  filename_prefix = "panel_ALL_NFL_without_outliers",
  panel_title = "ALL - NFL - Without Outliers - Clinical and Protein Markers"
)

# 
# # Kmeans using proteins and tests


matrix_kmeans <- cuestionarios_excel[,c('ID', proteins_columns_select, test_select,'Grupo')]


# final_matrix_all_M6a <- cluster_data_flexible_v2(
#   data = matrix_kmeans_filtered,
#   group_filter = NULL,    
#   value_to_exclude = c(exclude_vars, "TotalGM_CBF"),
#   proteins = c("NFL", "BDNF"),
#   centers = 2, 
#   nstart = 25
# )


