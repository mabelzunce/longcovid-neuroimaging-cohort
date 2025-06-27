library(shiny)
library(bslib)
library(dplyr)
library(ggplot2)
library(ggExtra)
library(readr)

# Define paths
path_data_analysis <- "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Freesurfer/DataAnalysis/ASL/"
path_asl_analysis <- "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/ASLProcessedCOVID/derivatives/ExploreASL/Population/Stats/"

path_cuestionarios <- paste(path_data_analysis, "RespuestasCuestionariosCovid_15_01.csv", sep="")
path_dkt_perfusion <- paste(path_asl_analysis, "mean_qCBF_StandardSpace_Desikan_Killiany_MNI_SPM12_n=203_15-Jan-2025_PVC2.tsv", sep="")

tests_variables <- c("EQ.VAS", "FAS", "PSQI",
                     "TMT.A", "TMT.A_perc", "TMT.B", "TMT.B_perc",
                     "WMS.R.DIR", "WMS.R.DIR_perc", "WMS.R.INV", "WMS.R.INV_perc",
                     "STROOP_P", "STROOP_P_perc", "STROOP_C", "STROOP_C_perc",
                     "STROOP_P.C","STROOP_P.C_perc",
                     "STROOP_P.C_INTERF_perc", "STROOP_P.C_INTERF",
                     "MOCA", "MOCA_perc", "HADSAnsiedad", "HADSDepresion",
                     "PHQ.9", "GAD7")

# Exclude rows
values_to_exclude <- c('CP0011', 'CP0015','CP0009', 'CP0035', 'CP0087', 'CP0105', 'CP0117', 'CP0144', 'CP0192','CP0216', 'CP0062')
vascular_artifacts <- c('CP0031','CP0044', 'CP0053', 'CP0054', 'CP0061','CP0107','CP0108',
                        'CP0123', 'CP0142', 'CP0147', 'CP0154', 'CP0105','CP0176',
                        'CP0178', 'CP0180', 'CP0183', 'CP0188', 'CP0205', 'CP0213')

values_to_exclude <- paste0("sub-", values_to_exclude, "_1")
vascular_artifacts <- paste0("sub-", vascular_artifacts, "_1")

# Leer archivos
cuestionarios_excel <- read.csv(path_cuestionarios)
cuestionarios_excel$ID <- paste0("sub-", cuestionarios_excel$ID,"_1")
perfusion_csv <- read_tsv(path_dkt_perfusion, show_col_types = FALSE)

# Excluir columnas con +10 NaN en Perfusion
perfusion_csv[perfusion_csv == "n/a"] <- NA
cols_with_nan_perfusion <- colSums(is.na(perfusion_csv)) > 10
perfusion_csv <- perfusion_csv[, !cols_with_nan_perfusion]

# Columnas con regiones de la segmentación 
columns_regions_perfusion <- setdiff(names(perfusion_csv)[6:length(names(perfusion_csv))], 
                                     grep("^Wm.", names(perfusion_csv), value = TRUE))

# Fusionar datos
merged_df_perfusion <- merge(
  cuestionarios_excel, perfusion_csv, 
  by.x = "ID", by.y = "participant_id", all.y = TRUE
) %>%
  slice(-1) %>% 
  mutate(across(everything(), ~ na_if(trimws(as.character(.)), "NaN"))) %>%
  drop_na(tests_variables)

# Filtrar sujetos a excluir
merged_df_perfusion <- subset(merged_df_perfusion, !(ID %in% values_to_exclude))
merged_df_perfusion <- subset(merged_df_perfusion, !(ID %in% vascular_artifacts))

# Variables numéricas
numeric_vars <- c(tests_variables, columns_regions_perfusion, "Edad")
numeric_vars <- numeric_vars[numeric_vars %in% names(merged_df_perfusion)]

numeric_df <- merged_df_perfusion %>%
  select(all_of(numeric_vars))  # Selecciona solo las columnas de numeric_vars

# UI
ui <- page_sidebar(
  sidebar = sidebar(
    varSelectInput("boxvar", "Select Variable", numeric_df, selected = names(numeric_df)[1]),  # Selector de variable
    checkboxGroupInput(
      "Group", "Filter by Group",
      choices = unique(merged_df_perfusion$Grupo), 
      selected = unique(merged_df_perfusion$Grupo)
    ),
    hr()
  ),
  plotOutput("boxplot")
)

# Server
server <- function(input, output, session) {
  subsetted <- reactive({
    req(input$Group)
    merged_df_perfusion |> filter(Grupo %in% input$Group)
  })
  
  output$boxplot <- renderPlot({
    req(input$boxvar)
    
    ggplot(subsetted(), aes(x = factor(Grupo), y = !!input$boxvar, fill = Grupo)) +
      geom_boxplot(alpha = 0.5, width = 0.3, color = "black", outlier.shape = NA) +  # ✅ Más ancho, sin outliers
      geom_jitter(position = position_jitter(width = 0.2, height = 0), 
                  color = "black", size = 2, alpha = 0.6) +  # ✅ Puntos más separados
      scale_fill_manual(values = c("CONTROL" = "#F4A582", "COVID" = "#92C5DE")) +  # ✅ Colores personalizados
      labs(title = paste("Distribución de", input$boxvar, "por grupo"),
           x = "Grupo",
           y = "CBF") +
      theme_minimal() +
      guides(fill = FALSE) +  # ✅ Oculta la leyenda de color
      theme(
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),  # ✅ Título más visible
        axis.title = element_text(size = 11, face = "bold"),  # ✅ Ejes en negrita
        axis.text = element_text(size = 10)  # ✅ Tamaño del texto de los ejes
      )
  }, res = 100)
  
  
}

shinyApp(ui, server)
shinyApp(ui, server)