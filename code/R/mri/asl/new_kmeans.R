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
path_plots_data_analysis <-  paste(path_data_analysis, "Plots/ClustersTests/", sep="")

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
proteins_columns_select_1 <- c('M6a','BDNF', 'NFL')

test_select <- c('MOCA_perc', 'TMT.A_perc', 'FAS', 'DepresionEscalaUnificada', 'AnsiedadEscalaUnificada')

# # Filtro los del grupo COVID
# cuestionarios_excel <- cuestionarios_excel %>%
#   filter(Grupo == "COVID")


cuestionarios_excel$M6a <- as.numeric(gsub(",", ".", as.character(cuestionarios_excel$M6a)))
cuestionarios_excel$BDNF <- as.numeric(gsub(",", ".", as.character(cuestionarios_excel$BDNF)))
cuestionarios_excel$NFL <- as.numeric(gsub(",", ".", as.character(cuestionarios_excel$NFL)))
cuestionarios_excel$MOCA_perc <- as.numeric(gsub(",", ".", as.character(cuestionarios_excel$MOCA_perc)))

matrix_kmeans <- cuestionarios_excel[,c('ID', proteins_columns_select_1,'Grupo')]

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
tests_variables_perc <- c("TMT.A_perc")

for (test in tests_variables_perc) {
  print(paste("Iterando con:", test))
  
  # Copiar base del data frame original sin la columna del test aún
  matrix_kmeans <- cuestionarios_excel[, c("ID", proteins_columns_select_1, "Grupo")]
  
    # Agregar la columna del test actual
  matrix_kmeans[[test]] <- cuestionarios_excel[[test]]
  
  # Convertir el test a numérico antes de agregar
  test_values <- as.numeric(cuestionarios_excel[[test]])
  matrix_kmeans[[test]] <- test_values
  
  matrix_kmeans_new <- matrix_kmeans[!(matrix_kmeans$ID %in% c("CP0223","CP0224", "CP0229", "CP0230", "CP0235", "CP0156")), ]
  
  
  
  # 
  # Ejecutar clustering solo con M6a y el test
  resultado <- cluster_data_flexible_v2(
    data = matrix_kmeans_new,
    group_filter = "COVID",
    proteins = NULL,
    value_to_exclude = c("BDNF",  "NFL"),
    centers = 2,
    nstart = 25
  )
  #
  
  
  # Filtrar filas con cluster no NA
  merged_data_filtrado <- resultado$merged_data %>%
    dplyr::filter(!is.na(Cluster))
  
  
  cat("\nResumen de clusters para", test, ":\n")
  merged_data_filtrado %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    print()
  
  # Crear título y prefijo dinámico
  prefijo <- paste0("M6a_", test, "_without outliers")
  titulo  <- paste("M6a -", test, "without outliers")
  
  tests_a_agregar <- setdiff(test_select, colnames(merged_data_filtrado))
  
  # Crear un mini-dataframe con solo ID y los tests faltantes
  tests_df <- cuestionarios_excel[, c("ID", tests_a_agregar)]
  
  # Hacer join por ID
  merged_data_filtrado <- merged_data_filtrado %>%
    left_join(tests_df, by = "ID")
  
  # Crear panel
  crear_panel_covid(
    data = merged_data_filtrado,
    filename_prefix = prefijo,
    panel_title = titulo
  )
  
  # Crear boxplot individual
  boxplot_test <- crear_boxplot(
    data = merged_data_filtrado,
    variable = test,
    title = paste("Boxplot de", test, "por clusters")
  )
  
  # Guardar boxplot individual
  ggsave(
    filename = paste0(path_plots_data_analysis, "boxplot_M6a_", test, ".png"),
    plot = boxplot_test,
    width = 8,
    height = 6
  )
  
  colores_clusters <- c("1" = "#440154", "CONTROL" = "#21908C", "2" = "#FDE725")
  
  
  # Crear gráfico de dispersión M6a vs test
  scatter_plot <- ggplot(merged_data_filtrado, aes_string(x = "M6a", y = test, color = "Cluster")) +
    geom_point(size = 3, alpha = 0.7) +
    labs(
      title = paste("Dispersión M6a vs", test),
      x = "M6a",
      y = test
    ) +
    theme_minimal(base_size = 14) +
    scale_color_manual(values = colores_clusters) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold")
    )
  # Filtrar solo Cluster 1 y 2
  data_filtrada <- merged_data_filtrado %>%
    dplyr::filter(Cluster %in% c(1, 2))
  
  # Crear gráfico
  scatter_plot_sin_control <- ggplot(data_filtrada, aes_string(x = "M6a", y = test, color = "Cluster")) +
    geom_point(size = 3, alpha = 0.7) +
    labs(
      title = paste("Dispersión M6a vs", test, "– Clusters 1 y 2"),
      x = "M6a",
      y = test
    ) +
    theme_bw(base_size = 14) +
    scale_color_manual(values = colores_clusters) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold")
    )
  

  
  # Crear gráfico combinado
  grafico_combinado <- scatter_plot + scatter_plot_sin_control +
    plot_layout(ncol = 2) +
    plot_annotation(title = paste("Comparación: M6a vs", test),
                    theme = theme(plot.title = element_text(face = "bold", size = 16)))
  
  # Guardar gráfico combinado
  ggsave(
    filename = paste0(path_plots_data_analysis, "comparacion_M6a_vs_", test, ".png"),
    plot = grafico_combinado,
    width = 16,  # doble de ancho
    height = 6,
    bg = "white"
  )
  
  
  # Asegurarse que Cluster es factor
  merged_data_filtrado$Cluster <- as.factor(merged_data_filtrado$Cluster)
  
  # ANOVA
  anova_result <- aov(M6a ~ Cluster, data = merged_data_filtrado)
  
  # Mostrar resumen
  cat("\nANOVA para M6a entre Cluster 1, 2 y Control (", test, "):\n")
  print(summary(anova_result))
  
  cat("\nTukey post hoc para M6a:\n")
  print(TukeyHSD(anova_result))
  
}

