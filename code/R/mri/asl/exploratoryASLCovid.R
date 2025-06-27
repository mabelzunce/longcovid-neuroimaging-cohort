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

source("functionsASL.R") # Para que funcione debo estar en el mismo directorio

# Define the paths
path_data_analysis <- "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Freesurfer/DataAnalysis/ASL/"
path_asl_analysis <- "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/ASLProcessedCOVID/derivatives/ExploreASL/Population/Stats/"

path_results_data_analysis <- paste(path_data_analysis, "Results/", sep="")
path_plots_data_analysis <-  paste(path_data_analysis, "Plots/", sep="")
path_cuestionarios <- paste(path_data_analysis, "RespuestasCuestionarios_27_01.csv", sep="")

raw_perfussion_data <-"mean_qCBF_StandardSpace_Hammers_n=203_15-Jan-2025_PVC2.tsv"
# raw_perfussion_data <-"mean_qCBF_StandardSpace_TotalGM_n=203_15-Jan-2025_PVC2.tsv"

path_dkt_perfusion <-paste(path_asl_analysis, raw_perfussion_data , sep="")

# Exclude rows
values_to_exclude <- c('CP0011', 'CP0015','CP0009', 'CP0035', 'CP0093', 'CP0101', 'CP0140',
                       'CP0035', 'CP0087', 'CP0105', 'CP0117', 'CP0144', 'CP0192','CP0216', 'CP0227', 'CP0062', 'CP0232','CP0226')

# vascular_artifacts <- c('CP0031','CP0044', 'CP0053', 'CP0054', 'CP0061', 'CP0096','CP0107','CP0108','CP0118',
#                         'CP0123', 'CP0142', 'CP0147', 'CP0154','CP0167', 'CP0176','CP0178',
#                         'CP0180', 'CP0183', 'CP0188','CP0193','CP0199', 'CP0205', 'CP0213',
#                         'CP0229', 'CP0238','CP0245', 'CP0247', 'CP0238')

values_to_exclude <- paste0("sub-", values_to_exclude, "_1")
# vascular_artifacts <- paste0("sub-", vascular_artifacts, "_1")


# Read files
cuestionarios_excel <- read.csv(path_cuestionarios)
cuestionarios_excel$ID <- paste0("sub-", cuestionarios_excel$ID,"_1")
perfusion_csv <- read.delim(path_dkt_perfusion)

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

# Convertir a Nan lo que no lo es para limpiar fecha de infección
cuestionarios_excel <- cuestionarios_excel %>%
  mutate(
    `Días.entre.primera.infección.y.fecha.de.cuestionario` = ifelse(
      `Días.entre.primera.infección.y.fecha.de.cuestionario` == "#VALUE!",
      NaN,
      as.numeric(`Días.entre.primera.infección.y.fecha.de.cuestionario`)
    ),
    `Días.entre.última.infección.y.fecha.de.cuestionario` = ifelse(
      `Días.entre.última.infección.y.fecha.de.cuestionario` == "#VALUE!",
      NaN,
      as.numeric(`Días.entre.última.infección.y.fecha.de.cuestionario`)
    )
  )

# Generar una columna con un valor único
cuestionarios_excel <- cuestionarios_excel %>%
  mutate(
    Dias_entre_infeccion_y_cuestionario = coalesce(
      `Días.entre.última.infección.y.fecha.de.cuestionario`,
      `Días.entre.primera.infección.y.fecha.de.cuestionario`
    )
  )

# Unir con variables necesarias de cuestionarios
cols_to_merge <- c("Grupo", "ID", "X1.1", "X2.2", "BMI", "Dias_entre_infeccion_y_cuestionario", "Días.entre.primera.infección.y.fecha.de.cuestionario", "Fecha.test.cognitivo")
cols_to_merge <- intersect(cols_to_merge, colnames(cuestionarios_excel))  

merged_df_perfusion <- merge(
  cuestionarios_excel[, cols_to_merge], 
  perfusion_csv, 
  by.x = "ID", 
  by.y = "participant_id", 
  all.y = TRUE
) %>% 
  slice(-1) %>%  # Eliminar la primera fila
  mutate(across(all_of(columns_regions_perfusion), ~ suppressWarnings(as.numeric(.))))  # Regiones %>%

merged_df_perfusion_with_ata <- subset(merged_df_perfusion, !(ID %in% values_to_exclude))
merged_df_perfusion <- subset(merged_df_perfusion_with_ata, !(ID %in% vascular_artifacts)) 

# NaN de valores que no son logicos
ids_a_nan <- c("sub-CP0226_1", "sub-CP0232_1", "sub-CP0182_1")
merged_df_perfusion$Días.entre.primera.infección.y.fecha.de.cuestionario[merged_df_perfusion$ID %in% ids_a_nan] <- NaN

# Filtrar por grupo covid 
df_covid <- merged_df_perfusion %>% filter(Grupo == "COVID")
df_control <- merged_df_perfusion %>% filter(Grupo == "CONTROL")



# Hacer gráfico de dispersión 
library(ggpmisc)

scatter_plot <- ggplot(
  df_covid, 
  aes(
    x = `Dias_entre_infeccion_y_cuestionario`, 
    y = FL_middle_frontal_gyrus_L,
    color = Grupo
  )
) +
  geom_point(size = 3, alpha = 0.8) +
  # Etiquetas sobre cada punto
  geom_text(
    aes(label = ID),
    position = position_jitter(width = 0.2, height = 0),
    vjust = -1, size = 2.5, alpha = 0.8
  ) +
  # Línea de tendencia (regresión lineal)
  geom_smooth(
    method = "lm", 
    formula = y ~ x, 
    se = TRUE,          # TRUE para mostrar banda de error
    color = "black"     # Cambia el color si deseas
  ) +
  # Mostrar la ecuación y el R^2
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., sep = " ~~~ ")),
    formula = y ~ x,
    parse = TRUE,        # Interpretar como expresión matemática
    label.x.npc = "right", # Ajusta posición de la etiqueta en X
    label.y.npc = "top",   # Ajusta posición de la etiqueta en Y
    color = "black"
  ) +
  labs(
    title = "Relación entre días post-infección y FL_middle_frontal_gyrus_L",
    x = "Días entre infección y cuestionario",
    y = "FL_middle_frontal_gyrus_L"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

# Mostrar en el visor
print(scatter_plot)


df_covid_box <- df_covid %>%
  mutate(
    Year = if_else(
      grepl("2024", Fecha.test.cognitivo),
      "2024",
      "2023"
    )
  )

df_control_box <- df_control %>%
  mutate(
    Year = if_else(
      grepl("2024", Fecha.test.cognitivo),
      "2024",
      "2023"
    )
  )


# Crear un dataframe combinado con una columna "Cohorte" que indique el grupo
df_combined <- bind_rows(
  df_covid_box  %>% mutate(Cohorte = "COVID"),
  df_control_box %>% mutate(Cohorte = "CONTROL")
)

bx <- ggplot(df_combined, aes(x = Year, y = FL_middle_frontal_gyrus_L, fill = Year)) +
  geom_boxplot(alpha = 0.5) +
  geom_jitter(color = "black", size = 2, alpha = 0.3, width = 0.2) +
  facet_wrap(~ Cohorte, scales = "fixed") + 
  # scales="fixed" es la opción por defecto, 
  # que asegura la misma escala en el eje Y para ambos paneles.
  labs(
    title = "Comparación COVID vs CONTROL por Año de Test Cognitivo",
    x = "Año (Fecha de Test Cognitivo)",
    y = "PL_supramarginal_gyrus_B"
  ) +
  theme_minimal() +
  guides(fill = FALSE)  # Quita la leyenda si no es necesaria

print(bx)
