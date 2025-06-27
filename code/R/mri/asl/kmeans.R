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

path_cuestionarios <- paste(path_data_analysis, "ResumenTotal_04_04.csv", sep="")
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
test_select <- c('MOCA_perc', 'TMT.A_perc', 'FAS','DepresionEscalaUnificada', 'AnsiedadEscalaUnificada')

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

# Correlate FS data 
cor_matrix_freesurfer <- cor(all_freesurfer[ , colnames(all_freesurfer) != "subject"], use = "pairwise.complete.obs")

output_path_cor_matrix_fs <- paste0(path_plots_data_analysis, "freesurfer_correlation_tests.png")

# Guardar gráfico unificado con colores pastel
png(output_path_cor_matrix_fs, width=800, height=800)

corrplot(
  cor_matrix_freesurfer , method="circle", type="lower", order="original", 
         col=colorRampPalette(c("lightblue", "white", "lightcoral"))(200),
         tl.col="black", tl.srt=45, addCoef.col="black")
dev.off()

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

# Ejemplos de uso:

# # 1. Clusterizar todo el dataset (incluyendo controles) usando 2 clusters
final_matrix_all <- cluster_data(data = matrix_kmeans, free_group = NULL, centers = 2, cognitive_vars=NULL)
final_matrix_all <- final_matrix_all %>% select(ID, Grupo, Cluster, everything())

# 2. Clusterizar todo el dataset pero asignar "CONTROL" a los sujetos cuyo Grupo es CONTROL.
# final_matrix_covid <- cluster_data(data = matrix_kmeans, free_group = "CONTROL", cognitive_vars=NULL, centers = 2)
# final_matrix_covid <- final_matrix_covid %>% select(ID, Grupo, Cluster, everything())
# 
# final_matrix_control <- cluster_data(data = matrix_kmeans, free_group = "COVID", cognitive_vars=NULL, centers = 2)
# final_matrix_control <- final_matrix_control %>% select(ID, Grupo, Cluster, everything())

# final_matrix_control <- cluster_data_by_group(data = matrix_kmeans,
#                                               group_filter = "CONTROL",
#                                               centers = 2)
# 
# final_matrix_covid <- cluster_data_by_group(data = matrix_kmeans,
#                                             group_filter = "COVID",
#                                             centers = 2)
# 
# 
# levels(final_matrix_control$Cluster) <- paste0("C", levels(final_matrix_control$Cluster))
# levels(final_matrix_covid$Cluster) <- paste0("V", levels(final_matrix_covid$Cluster))
# 
# final_matrix_both <- bind_rows(final_matrix_control, final_matrix_covid)
# final_matrix_both <- final_matrix_both %>% select(ID, Grupo, Cluster, everything())
# 
# final_matrix_both <- final_matrix_both[order(final_matrix_both$ID), ]
# rownames(final_matrix_both) <- NULL
# 
# final_matrix_both <- final_matrix_both %>% filter(!is.na(Cluster))
# 
# # Seleccionar las columnas a graficar (excluyendo ID, Grupo y Cluster)

proteins = c("BDNF", "M6a", "NFL", "TotalGM_CBF")

variables_to_plot <- setdiff(names(final_matrix_all), c("ID", "Grupo", "Cluster",proteins))
variables_to_plot_2 <- setdiff(names(final_matrix_all), c("ID", "Grupo", "Cluster")) 

# Resta entre centroides para ir probando 
final_matrix_control <- cluster_data_by_group(data = matrix_kmeans,
                                              group_filter = "CONTROL",
                                              centers = 2)

final_matrix_covid <- cluster_data_by_group(data = matrix_kmeans,
                                            group_filter = "COVID",
                                            centers = 2)
centroids_control <- final_matrix_control$km_result$centers
diff_centroids_control <- centroids_control[1,] - centroids_control[2,]
abs_diff_centroids_control <- sort(abs(diff_centroids_control),decreasing=TRUE)


print("CONTROL:")
print(abs_diff_centroids_control)

centroids_covid <- final_matrix_covid$km_result$centers
diff_centroids_covid <- centroids_covid[1,] - centroids_covid[2,]
abs_diff_centroids_covid <- sort(abs(diff_centroids_covid),decreasing=TRUE)

print("COVID:")
print(abs_diff_centroids_covid)

# Inicializar un data frame para almacenar los valores de withinss
withinss_df <- data.frame(Variable = character(),
                          Group = character(),
                          Cluster = character(),
                          Withinss = numeric(),
                          stringsAsFactors = FALSE)

# Extract ID of cluster 1 CONTROL
cluster1_ids_control <- final_matrix_control[["merged_data"]] %>% 
  filter(Cluster == "1") %>% 
  pull(ID)


cluster1_ids_covid <- final_matrix_covid[["merged_data"]] %>% 
  filter(Cluster == "1") %>% 
  pull(ID)


# Iterar sobre todos los variables y calcular la varianza
for (var in variables_to_plot) {
  
  final_matrix_control_temp <- cluster_data_by_group(data = matrix_kmeans,
                                                group_filter = "CONTROL",
                                                value_to_exclude = var,
                                                centers = 2)
  
  final_matrix_covid_temp <- cluster_data_by_group(data = matrix_kmeans,
                                              group_filter = "COVID",
                                              value_to_exclude = var,
                                              centers = 2)
  
  # Extract ID of cluster 1 CONTROL
  cluster1_ids_control_tmp <- final_matrix_control_temp[["merged_data"]] %>% 
    filter(Cluster == "1") %>% 
    pull(ID)
  
  
  cluster1_ids_covid_tmp <- final_matrix_covid_temp[["merged_data"]] %>% 
    filter(Cluster == "1") %>% 
    pull(ID)
  
  

  
  # Calcular el porcentaje de coincidencia para cada grupo
  percentage_control <- (length(intersect(cluster1_ids_control, cluster1_ids_control_tmp)) / length(cluster1_ids_control)) * 100
  percentage_covid   <- (length(intersect(cluster1_ids_covid, cluster1_ids_covid_tmp)) / length(cluster1_ids_covid)) * 100
  
  cat("Variable:", var, "\n")
  cat("CONTROL: ", round(percentage_control, 2), "%\n")
  cat("COVID: ", round(percentage_covid, 2), "%\n")
  
  # Si el porcentaje en CONTROL es menor al 50%, se intercambian los nombres de cluster 
  # y se reordena el vector withinss (intercambiando los dos elementos)
  if (percentage_control < 50) {
    final_matrix_control_temp[["merged_data"]] <- final_matrix_control_temp[["merged_data"]] %>%
      mutate(Cluster = ifelse(Cluster == "1", "2", ifelse(Cluster == "2", "1", Cluster)))
    
    # Reordenar el vector withinss (intercambiando el primer y el segundo valor)
    final_matrix_control_temp[["km_result"]]$withinss <- final_matrix_control_temp[["km_result"]]$withinss[c(2, 1)]
    
    cat("Se intercambiaron los nombres de cluster y se reordenaron los valores del vector dentro de km_result para CONTROL\n")
  }
  
  # Si el porcentaje en COVID es menor al 50%, se intercambian los nombres de cluster 
  # y se reordena el vector withinss para COVID
  if (percentage_covid < 50) {
    final_matrix_covid_temp[["merged_data"]] <- final_matrix_covid_temp[["merged_data"]] %>%
      mutate(Cluster = ifelse(Cluster == "1", "2", ifelse(Cluster == "2", "1", Cluster)))
    
    # Reordenar el vector withinss (intercambiando el primer y el segundo valor)
    final_matrix_covid_temp[["km_result"]]$withinss <- final_matrix_covid_temp[["km_result"]]$withinss[c(2, 1)]
    
    cat("Se intercambiaron los nombres de cluster y se reordenaron los valores del vector dentro de km_result para COVID\n")
  }
  
  km_control <- final_matrix_control_temp$km_result
  km_covid   <- final_matrix_covid_temp$km_result
  
    
  # Extraer la suma de cuadrados intra-cluster (withinss)
  w_control <- km_control$withinss   # Vector de dos valores para CONTROL
  w_covid   <- km_covid$withinss     # Vector de dos valores para COVID
 
  
  # Crear un data frame temporal con los resultados para la variable 'var'
  temp_df <- data.frame(
    Variable = rep(var, times = 4),
    Group = rep(c("CONTROL", "COVID"), each = 2),
    Cluster = rep(c("1", "2"), times = 2),
    Withinss = c(w_control, w_covid),
    stringsAsFactors = FALSE
  )
  
  # Unir al data frame final
  withinss_df <- rbind(withinss_df, temp_df)
  
    
}



# Agregar variabilidad total para cada cluster
final_matrix_control <- cluster_data_by_group(data = matrix_kmeans,
                                              group_filter = "CONTROL",
                                              centers = 2)

final_matrix_covid <- cluster_data_by_group(data = matrix_kmeans,
                                            group_filter = "COVID",
                                            centers = 2)

km_control_final <- final_matrix_control$km_result
km_covid_final   <- final_matrix_covid$km_result

# Extraer la suma de cuadrados intra-cluster (withinss)
w_control_final <- km_control_final$withinss   # Vector de dos valores para CONTROL
w_covid_final   <- km_covid_final$withinss     # Vector de dos valores para COVID


total_withinss_df <- data.frame(
  Variable = rep("Total", times = 4),
  Group = rep(c("CONTROL", "COVID"), each = 2),
  Cluster = rep(c("1", "2"), times = 2),
  Withinss = c(w_control_final, w_covid),
  stringsAsFactors = FALSE
)

withinss_df <- rbind(withinss_df, total_withinss_df)


# Filtrar para cada grupo:
withinss_control_df <- withinss_df %>% filter(Group == "CONTROL")
withinss_covid_df   <- withinss_df %>% filter(Group == "COVID")

# Calcular límites comunes en el eje Y (por ejemplo, iniciando en 0)
y_min <- 0
y_max_control <- max(withinss_control_df$Withinss, na.rm = TRUE)
y_max_covid <- max(withinss_covid_df$Withinss, na.rm = TRUE)

# Extraer el valor de referencia para la línea horizontal
withinss_total_control_1 <- withinss_control_df %>% 
  filter(Cluster == 1 & Variable == "Total") %>% 
  pull(Withinss)

withinss_total_control_2 <- withinss_control_df %>% 
  filter(Cluster == 2 & Variable == "Total") %>% 
  pull(Withinss)


# Gráfico para el grupo CONTROL
p_control <- ggplot(withinss_control_df, aes(x = Variable, y = Withinss, fill = factor(Cluster))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  scale_y_continuous(limits = c(y_min, y_max_control)) +
  geom_hline(yintercept = withinss_total_control_1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = withinss_total_control_2, linetype = "dashed", color = "black") +
  labs(title = "Within-cluster Sum of Squares - CONTROL",
       x = "Variable",
       y = "Withinss") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Extraer el valor de referencia para la línea horizontal
withinss_total_covid_1 <- withinss_covid_df %>% 
  filter(Cluster == 1 & Variable == "Total") %>% 
  pull(Withinss)

withinss_total_covid_2 <- withinss_covid_df %>% 
  filter(Cluster == 2 & Variable == "Total") %>% 
  pull(Withinss)


# Gráfico para el grupo COVID
p_covid <- ggplot(withinss_covid_df, aes(x = Variable, y = Withinss, fill = factor(Cluster))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = withinss_total_covid_1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = withinss_total_covid_2, linetype = "dashed", color = "black") +
  scale_y_continuous(limits = c(y_min, y_max_covid)) +
  labs(title = "Within-cluster Sum of Squares - COVID",
       x = "Variable",
       y = "Withinss") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Crear el gráfico combinado
combined_plot <- grid.arrange(p_control, p_covid, ncol = 1)

# Guardarlo con ggsave
ggsave(
  filename = paste0(path_plots_data_analysis, "withinss_groups.png"),
  plot = combined_plot,
  width = 10, height = 12, dpi = 300
)

## New Analysis: excluding values
values_to_exclude_control <- c('Deep_WMH', 'FAS', 'frontalpole_vol_rh', 
                               'parsopercularis_vol_rh', 'postcentral_thick_rh', 
                               'postcentral_thick_lh', 'PV_WMH', 
                               'rostralmiddlefrontal_vol_rh', 'supramarginal_thick_lh', 'TMT.A_perc', 'Total_WMH')


values_to_exclude_covid <- c('CardiovascularRiskIndex','Deep_WMH', 'FAS', 'frontalpole_vol_rh', 
                               'fusiform_thick_rh', 'fusiform_vol_rh', 
                               'postcentral_thick_lh', 'postcentral_thick_rh', 
                               'supramarginal_thick_lh', 'TMT.A_perc', 'TotalGM_sCOV')


final_matrix_control_excluding_values <- cluster_data_by_group(data = matrix_kmeans,
                                                   group_filter = "CONTROL",
                                                   value_to_exclude = values_to_exclude_control,
                                                   centers = 2)

final_matrix_covid_excluding_values <- cluster_data_by_group(data = matrix_kmeans,
                                                 group_filter = "COVID",
                                                 value_to_exclude = values_to_exclude_covid,
                                                 centers = 2)

withinss_df_excluding <- data.frame(Variable = character(),
                          Group = character(),
                          Cluster = character(),
                          Withinss = numeric(),
                          stringsAsFactors = FALSE)

# Extract ID of cluster 1 CONTROL
cluster1_ids_control <- final_matrix_control_excluding_values[["merged_data"]] %>% 
  filter(Cluster == "1") %>% 
  pull(ID)


cluster1_ids_covid <- final_matrix_covid_excluding_values[["merged_data"]] %>% 
  filter(Cluster == "1") %>% 
  pull(ID)

colnames_control <- setdiff(colnames(final_matrix_control_excluding_values[["merged_data"]]), c("ID", "Cluster", "Grupo", values_to_exclude_control, proteins))
colnames_covid <- setdiff(colnames(final_matrix_covid_excluding_values[["merged_data"]]), c("ID", "Cluster", "Grupo", values_to_exclude_covid, proteins))

# Numero de controles en merged data 

n_controles_con_cluster <- sum(
  !is.na(final_matrix_control_excluding_values[["merged_data"]]$Cluster)
)

n_covid_con_cluster <- sum(
  !is.na(final_matrix_covid_excluding_values[["merged_data"]]$Cluster)
)

# Iterar sobre todos los variables y calcular la varianza
for (var in colnames_control) {
  
  final_matrix_control_temp <- cluster_data_by_group(data = matrix_kmeans,
                                                     group_filter = "CONTROL",
                                                     value_to_exclude = c(var, values_to_exclude_control),
                                                     centers = 2)
  
  # Extract ID of cluster 1 CONTROL
  cluster1_ids_control_tmp <- final_matrix_control_temp[["merged_data"]] %>% 
    filter(Cluster == "1") %>% 
    pull(ID)

  
  # Calcular el porcentaje de coincidencia para cada grupo
  percentage_control <- (length(intersect(cluster1_ids_control, cluster1_ids_control_tmp)) / length(cluster1_ids_control)) * 100

  # Si el porcentaje en CONTROL es menor al 50%, se intercambian los nombres de cluster 
  # y se reordena el vector withinss (intercambiando los dos elementos)
  if (percentage_control < 50) {
    final_matrix_control_temp[["merged_data"]] <- final_matrix_control_temp[["merged_data"]] %>%
      mutate(Cluster = ifelse(Cluster == "1", "2", ifelse(Cluster == "2", "1", Cluster)))
    
    # Reordenar el vector withinss (intercambiando el primer y el segundo valor)
    final_matrix_control_temp[["km_result"]]$withinss <- final_matrix_control_temp[["km_result"]]$withinss[c(2, 1)]
    
    cat("Se intercambiaron los nombres de cluster y se reordenaron los valores del vector dentro de km_result para CONTROL\n")
  }
  

  km_control <- final_matrix_control_temp$km_result
  w_control <- (km_control$withinss)/n_controles_con_cluster  # Vector de dos valores para CONTROL

  
  # Crear un data frame temporal con los resultados para la variable 'var'
  temp_df <- data.frame(
    Variable = rep(var, times = 2),
    Group = rep(c("CONTROL"), each = 2),
    Cluster = c("1", "2"),
    Withinss = c(w_control),
    stringsAsFactors = FALSE
  )
  
  # Unir al data frame final
  withinss_df_excluding <- rbind(withinss_df_excluding, temp_df)
  
  
}

# COVID
# Iterar sobre todos los variables y calcular la varianza
for (var in colnames_covid) {
  
  final_matrix_covid_temp <- cluster_data_by_group(data = matrix_kmeans,
                                                     group_filter = "COVID",
                                                     value_to_exclude = c(var,values_to_exclude_covid),
                                                     centers = 2)
  
  # Extract ID of cluster 1 CONTROL
  cluster1_ids_covid_tmp <- final_matrix_covid_temp[["merged_data"]] %>% 
    filter(Cluster == "1") %>% 
    pull(ID)
  
  
  # Calcular el porcentaje de coincidencia para cada grupo
  
  # Si el porcentaje en CONTROL es menor al 50%, se intercambian los nombres de cluster 
  # y se reordena el vector withinss (intercambiando los dos elementos)
  if (percentage_covid < 50) {
  final_matrix_covid_temp[["merged_data"]] <- final_matrix_covid_temp[["merged_data"]] %>%
    mutate(Cluster = ifelse(Cluster == "1", "2", ifelse(Cluster == "2", "1", Cluster)))
  
  # Reordenar el vector withinss (intercambiando el primer y el segundo valor)
  final_matrix_covid_temp[["km_result"]]$withinss <- final_matrix_covid_temp[["km_result"]]$withinss[c(2, 1)]


}

  
  km_covid <- final_matrix_covid_temp$km_result
  w_covid <- (km_covid$withinss)/n_covid_con_cluster   # Vector de dos valores para CONTROL
  
  
  # Crear un data frame temporal con los resultados para la variable 'var'
  temp_df <- data.frame(
    Variable = rep(var, times = 2),
    Group = rep(c("COVID"), each = 2),
    Cluster = c("1", "2"),
    Withinss = c(w_covid),
    stringsAsFactors = FALSE
  )
  
  # Unir al data frame final
  withinss_df_excluding <- rbind(withinss_df_excluding, temp_df)
  
  
}

# Agregar variabilidad total para cada cluster
final_matrix_control_excluding_values <- cluster_data_by_group(data = matrix_kmeans,
                                              group_filter = "CONTROL",
                                              value_to_exclude = values_to_exclude_control,
                                              centers = 2)

final_matrix_covid_excluding_values <- cluster_data_by_group(data = matrix_kmeans,
                                            group_filter = "COVID",
                                            value_to_exclude =  values_to_exclude_covid, 
                                            centers = 2)

km_control_final <- final_matrix_control_excluding_values$km_result
km_covid_final   <- final_matrix_covid_excluding_values$km_result

# Extraer la suma de cuadrados intra-cluster (withinss)
w_control_final <- (km_control_final$withinss)/n_controles_con_cluster   # Vector de dos valores para CONTROL
w_covid_final   <- (km_covid_final$withinss)/n_covid_con_cluster     # Vector de dos valores para COVID


total_withinss_df <- data.frame(
  Variable = rep("Total", times = 4),
  Group = rep(c("CONTROL", "COVID"), each = 2),
  Cluster = rep(c("1", "2"), times = 2),
  Withinss = c(w_control_final, w_covid_final),
  stringsAsFactors = FALSE
)

withinss_df_excluding <- rbind(withinss_df_excluding, total_withinss_df)


# Filtrar para cada grupo:
withinss_control_excluding_df <- withinss_df_excluding %>% filter(Group == "CONTROL")

withinss_covid_excluding_df   <- withinss_df_excluding %>% filter(Group == "COVID")

# Calcular límites comunes en el eje Y (por ejemplo, iniciando en 0)
y_min <- 0
y_max_control <- max(withinss_control_excluding_df$Withinss, na.rm = TRUE)
y_max_covid <- max(withinss_covid_excluding_df$Withinss, na.rm = TRUE)

# Extraer el valor de referencia para la línea horizontal
withinss_total_control_1_excluding <- withinss_control_excluding_df %>% 
  filter(Cluster == 1 & Variable == "Total") %>% 
  pull(Withinss)

withinss_total_control_2_excluding <- withinss_control_excluding_df %>% 
  filter(Cluster == 2 & Variable == "Total") %>% 
  pull(Withinss)


# Gráfico para el grupo CONTROL
p_control <- ggplot(withinss_control_excluding_df, aes(x = Variable, y = Withinss, fill = factor(Cluster))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  scale_y_continuous(limits = c(y_min, y_max_control)) +
  geom_hline(yintercept = withinss_total_control_1_excluding, linetype = "dashed", color = "black") +
  geom_hline(yintercept = withinss_total_control_2_excluding, linetype = "dashed", color = "black") +
  labs(title = "Within-cluster Sum of Squares - CONTROL",
       x = "Variable",
       y = "Withinss") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Extraer el valor de referencia para la línea horizontal
withinss_total_covid_1_excluding <- withinss_covid_excluding_df %>% 
  filter(Cluster == 1 & Variable == "Total") %>% 
  pull(Withinss)

withinss_total_covid_2_excluding <- withinss_covid_excluding_df %>% 
  filter(Cluster == 2 & Variable == "Total") %>% 
  pull(Withinss)


# Gráfico para el grupo COVID
p_covid <- ggplot(withinss_covid_excluding_df, aes(x = Variable, y = Withinss, fill = factor(Cluster))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = withinss_total_covid_1_excluding, linetype = "dashed", color = "black") +
  geom_hline(yintercept = withinss_total_covid_2_excluding, linetype = "dashed", color = "black") +
  scale_y_continuous(limits = c(y_min, y_max_covid)) +
  labs(title = "Within-cluster Sum of Squares - COVID",
       x = "Variable",
       y = "Withinss") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Crear el gráfico combinado
combined_plot <- grid.arrange(p_control, p_covid, ncol = 1)

# Guardarlo con ggsave
ggsave(
  filename = paste0(path_plots_data_analysis, "withinss_groups_excluding.png"),
  plot = combined_plot,
  width = 10, height = 12, dpi = 300)



# Reemplazar la columna Cluster por "CONTROL" en la parte de control
final_matrix_control_excluding_values_new <- final_matrix_control_excluding_values[["merged_data"]]
final_matrix_control_excluding_values_new$Cluster <- "CONTROL"


# Unir ambas matrices COVID
final_matrix_covid_final <- dplyr::bind_rows(
  final_matrix_control_excluding_values_new,
  final_matrix_covid_excluding_values[["merged_data"]]
)

# Reemplazar la columna Cluster por "COVID" en la parte de control
final_matrix_covid_excluding_values_new <- final_matrix_covid_excluding_values[["merged_data"]]
final_matrix_covid_excluding_values_new$Cluster <- "COVID"


# Unir ambas matrices CONTROL
final_matrix_control_final <- dplyr::bind_rows(
  final_matrix_covid_excluding_values_new,
  final_matrix_control_excluding_values[["merged_data"]]
)


final_matrix_covid_final <- final_matrix_covid_final %>%
  dplyr::filter(
    !is.na(Cluster) 
  )


final_matrix_control_final <- final_matrix_control_final %>%
  dplyr::filter(
    !is.na(Cluster) 
  )

titulos_proteinas <- c("M6a", "BDNF", "NFL")
clinicas <- c("DepresionEscalaUnificada", "AnsiedadEscalaUnificada")
titulos_clinicas <- c("Depresión", "Ansiedad")


# Función auxiliar para generar y guardar paneles
for (i in seq_along(titulos_proteinas)) {
  
  prot <- titulos_proteinas[i]
  titulo_prot <- titulos_proteinas[i]
  
  # Variables a graficar en el panel
  vars <- c(clinicas, titulo_prot)
  titulos <- c(titulos_clinicas, titulo_prot)
  
  ## === COVID ===
  plots_covid <- list()
  for (j in seq_along(vars)) {
    if (vars[j] == "M6a") {
      ylim <- 25
    } else if (vars[j] == "BDNF") {
      ylim <- 0.3
    } else {
      ylim <- NULL
    }
    
    plots_covid[[j]] <- crear_boxplot(
      data = final_matrix_covid_final,
      variable = vars[j],
      title = titulos[j],
      ylim_sup = ylim, 
      mostrar_ids = FALSE
    )
  }
  
  panel_covid <- (plots_covid[[1]] + plots_covid[[2]] + plots_covid[[3]]) +
    plot_annotation(
      title = paste("COVID -", titulo_prot),
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
  
  ggsave(
    filename = paste0(path_plots_data_analysis, "panel_", prot, "_COVID_new_lim.png"),
    plot = panel_covid,
    width = 10, height = 5, dpi = 300
  )
  
  ## === CONTROL ===
  plots_control <- list()
  for (j in seq_along(vars)) {
    
    if (vars[j] == "M6a") {
      ylim <- 25
    } else if (vars[j] == "BDNF") {
      ylim <- 0.3
    } else {
      ylim <- NULL
    }
    
    
    plots_control[[j]] <- crear_boxplot(
      data = final_matrix_control_final,
      variable = vars[j],
      title = titulos[j],
      ylim_sup = ylim,
      mostrar_ids = FALSE
    )
  }
  
  panel_control <- (plots_control[[1]] + plots_control[[2]] + plots_control[[3]]) +
    plot_annotation(
      title = paste("CONTROL -", titulo_prot),
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
  
  ggsave(
    filename = paste0(path_plots_data_analysis, "panel_", prot, "_CONTROL_new_lim.png"),
    plot = panel_control,
    width = 10, height = 5, dpi = 300
  )
}


