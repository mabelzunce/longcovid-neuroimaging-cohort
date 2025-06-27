library(shiny)
library(bslib)
library(dplyr)
library(ggplot2)
library(ggExtra)
library(readr)

# Define paths
path_data_analysis <- "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Freesurfer/DataAnalysis/ASL/"
path_asl_analysis <- "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/ASLProcessedCOVID/derivatives/ExploreASL/Population/Stats/"

path_results_data_analysis <- paste(path_data_analysis, "Results/", sep="")
path_plots_data_analysis <-  paste(path_data_analysis, "Plots/", sep="")
path_cuestionarios <- paste(path_data_analysis, "RespuestasCuestionariosCovid_15_01.csv", sep="")

path_dkt_perfusion <-paste(path_asl_analysis, "mean_qCBF_StandardSpace_Desikan_Killiany_MNI_SPM12_n=203_15-Jan-2025_PVC2.tsv", sep="")

tests_variables <- c("EQ.VAS", "FAS", "PSQI",
                     "TMT.A", "TMT.A_perc",
                     "TMT.B", "TMT.B_perc",
                     "WMS.R.DIR", "WMS.R.DIR_perc",
                     "WMS.R.INV", "WMS.R.INV_perc",
                     "STROOP_P", "STROOP_P_perc",
                     "STROOP_C", "STROOP_C_perc",
                     "STROOP_P.C","STROOP_P.C_perc",
                     "STROOP_P.C_INTERF_perc", "STROOP_P.C_INTERF",
                     "MOCA", "MOCA_perc",
                     "HADSAnsiedad", "HADSDepresion",
                     "PHQ.9", "GAD7",
                     "M6a",	"BDNF",	"NFL")

# Exclude rows
values_to_exclude <- c('CP0011', 'CP0015','CP0009', 'CP0035', 'CP0087', 'CP0105', 'CP0117', 'CP0144', 'CP0192','CP0216', 'CP0062')

vascular_artifacts <- c('CP0031','CP0044', 'CP0053', 'CP0054', 'CP0061','CP0107','CP0108',
                        'CP0123', 'CP0142', 'CP0147', 'CP0154', 'CP0105','CP0176',
                        'CP0178', 'CP0180', 'CP0183', 'CP0188', 'CP0205', 'CP0213')

values_to_exclude <- paste0("sub-", values_to_exclude, "_1")
vascular_artifacts <- paste0("sub-", vascular_artifacts, "_1")

# Read files
cuestionarios_excel <- read.csv(path_cuestionarios)
cuestionarios_excel$ID <- paste0("sub-", cuestionarios_excel$ID,"_1")
perfusion_csv <- read_tsv(path_dkt_perfusion)

# Exclude columns with +10 NaN in Perfusion
perfusion_csv[perfusion_csv == "n/a"] <- NA
cols_with_nan_perfusion <- colSums(is.na(perfusion_csv)) > 10
columns_to_exclude_perfusion <- names(perfusion_csv)[cols_with_nan_perfusion]
perfusion_csv <- perfusion_csv[, !cols_with_nan_perfusion]

# Columnas con regiones de la segmentación 
column_names_perfusion <- names(perfusion_csv)
columns_regions_perfusion <- column_names_perfusion[6:length(column_names_perfusion)]
columns_regions_perfusion <- columns_regions_perfusion[!grepl("^Wm.", columns_regions_perfusion)]

# Unir con variables necesarias de cuestionarios
cols_to_merge <- c("Edad", "Grupo", "ID", "Genero", "BMI", tests_variables)
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
    varSelectInput("xvar", "X variable", numeric_df, selected = names(numeric_df)[1]),
    varSelectInput("yvar", "Y variable", numeric_df, selected = names(numeric_df)[2]),
    checkboxGroupInput(
      "Group", "Filter by Group",
      choices = unique(merged_df_perfusion$Grupo), 
      selected = unique(merged_df_perfusion$Grupo)
    ),
    hr(),
    checkboxInput("by_group", "Show Groups", TRUE),
    checkboxInput("show_labels", "Show IDs", FALSE)  # ✅ Nuevo checkbox para mostrar ID
  ),
  plotOutput("scatter")
)

# Server
server <- function(input, output, session) {
  subsetted <- reactive({
    req(input$Group)
    merged_df_perfusion |> filter(Grupo %in% input$Group)
  })
  
  output$scatter <- renderPlot({
    req(input$xvar, input$yvar)
    
    p <- ggplot(subsetted(), aes(x = !!input$xvar, y = !!input$yvar)) +
      geom_point(aes(color = if (input$by_group) Grupo else NULL)) +  
      theme_minimal() +
      theme(legend.position = "bottom")
    
    # ✅ Agregar etiquetas con ID si la opción está activada
    if (input$show_labels) {
      p <- p + geom_text(aes(label = ID), vjust = -1, size = 3)
    }
    
    p
  }, res = 100)
}

shinyApp(ui, server)

