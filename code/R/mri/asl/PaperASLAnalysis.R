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

library(Hmisc)
library(corrplot)

setwd("/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Freesurfer/DataAnalysis/CodeR")

source("functionsANCOVA.R")
source("functionsASL.R") 

# ------------------ PRE PROCESSING  ------------------

# Define Variables
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

tests_variables_perc <- c("TMT.A_perc", "TMT.B_perc",
                          "WMS.R.DIR_perc", "WMS.R.INV_perc",
                          "STROOP_P_perc", "STROOP_C_perc",
                          "STROOP_P.C_perc", "STROOP_P.C_INTERF_perc",
                          "MOCA_perc")


columnas_antecedentes <- c("PresionAlta", "Diabetes", "Asma", "ColesterolAlto", "InfartoCardiaco", "AnginaPecho",
                           "EmboliaTrombos", "MedicacionPresionArterial", "Aspirinas", "Antiagregantes")

escalas_depresion_ansiedad <- c("DepresionCategoriaUnificada", "AnsiedadEscalaUnificada")
wmh_variables <-  c("WMHTotal", "WMHPeriventricular", "WMHDeep")

# Subjects CLassification
# Exclude rows
values_to_exclude <- c('CP0011', 'CP0015','CP0009', 'CP0035', 'CP0093', 'CP0101', 'CP0140', "CP0106", "CP0193",
                       'CP0035', 'CP0087', 'CP0105', 'CP0117', 'CP0144', 'CP0192','CP0216', 'CP0227', 'CP0062')

vascular_artifacts <- c('CP0031','CP0044', 'CP0053', 'CP0054', 'CP0061', 'CP0096','CP0107','CP0108','CP0118',
                        'CP0123', 'CP0142', 'CP0147', 'CP0154','CP0167', 'CP0176','CP0178',
                        'CP0180', 'CP0183', 'CP0188','CP0193','CP0199', 'CP0205', 'CP0213',
                        'CP0229', 'CP0238','CP0245', 'CP0247')

values_to_exclude <- paste0("sub-", values_to_exclude, "_1")
vascular_artifacts <- paste0("sub-", vascular_artifacts, "_1")


# sCOV
cov_perfussion_data_total_gm_without_PVC <- "CoV_qCBF_StandardSpace_TotalGM_n=203_15-Jan-2025_PVC0.tsv"
cov_perfussion_data_hammers_without_PVC <- "CoV_qCBF_StandardSpace_Hammers_n=203_15-Jan-2025_PVC0.tsv"
cov_perfussion_data_lobes_without_PVC <- "CoV_qCBF_StandardSpace_MNI_Structural_n=203_15-Jan-2025_PVC0.tsv"

# CBF
raw_perfussion_data_Hammers <-"mean_qCBF_StandardSpace_Hammers_n=203_15-Jan-2025_PVC2.tsv"
raw_gray_matter_data_total <-"mean_qCBF_StandardSpace_TotalGM_n=203_15-Jan-2025_PVC2.tsv"
raw_perfussion_data_lobes <-"mean_qCBF_StandardSpace_MNI_Structural_n=203_15-Jan-2025_PVC2.tsv"

# Define the paths
path_data_analysis <- "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Freesurfer/DataAnalysis/ASL/"
path_asl_analysis <- "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/ASLProcessedCOVID/derivatives/ExploreASL/Population/Stats/"
path_bianca_file <- "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/WMH/ProcessedWMH/TotalPvAndDeepValuesBiancaThr0_7.csv"

path_results_data_analysis <- paste(path_data_analysis, "Results/", sep="")
path_plots_data_analysis <-  paste(path_data_analysis, "Plots/Paper/", sep="")
path_cuestionarios <- paste(path_data_analysis, "ResumenTotal_2_06.csv", sep="")
path_socialAndCardiovascularIndex <- paste(path_data_analysis, "SocialAndCardiovascularIndex.csv", sep="")


path_perfussion_Hammers <-paste(path_asl_analysis, raw_perfussion_data_Hammers , sep="")
path_perfussion_data_lobes <-paste(path_asl_analysis, raw_perfussion_data_lobes , sep="")
path_gray_matter_perfusion <-paste(path_asl_analysis, raw_gray_matter_data_total , sep="")

path_scov_perfusion_total_gm <-paste(path_asl_analysis, cov_perfussion_data_total_gm_without_PVC , sep="")
path_scov_perfusion_Hammers <-paste(path_asl_analysis, cov_perfussion_data_hammers_without_PVC , sep="")
path_scov_perfussion_data_lobes <-paste(path_asl_analysis, cov_perfussion_data_lobes_without_PVC , sep="")

# Read Files
cuestionarios_excel <- read.csv(path_cuestionarios)
social_index <- read.csv(path_socialAndCardiovascularIndex)

# CBF
cbf_Hammers_csv <- read.delim(path_perfussion_Hammers)
cbf_gray_matter_csv <- read.delim(path_gray_matter_perfusion)
cbf_lobes_csv <- read.delim(path_perfussion_data_lobes)

# sCOV 
scov_Hammers_csv <- read.delim(path_scov_perfusion_Hammers)
scov_gray_matter_csv <- read.delim(path_scov_perfusion_total_gm)
scov_lobes_csv <- read.delim(path_scov_perfussion_data_lobes)

# Exclude columns with +10 NaN in Perfusion
cbf_Hammers_csv[cbf_Hammers_csv == "n/a"] <- NA
cols_with_nan_perfusion <- colSums(is.na(cbf_Hammers_csv)) > 10
columns_to_exclude_perfusion <- names(cbf_Hammers_csv)[cols_with_nan_perfusion]

cbf_Hammers_csv <- cbf_Hammers_csv[, !cols_with_nan_perfusion]
scov_Hammers_csv <- scov_Hammers_csv[, !cols_with_nan_perfusion]

colnames(cbf_Hammers_csv) <- make.names(colnames(cbf_Hammers_csv), unique = TRUE)
colnames(scov_Hammers_csv) <- make.names(colnames(scov_Hammers_csv), unique = TRUE)

cols_to_remove_B <- grep("_B$", colnames(cbf_Hammers_csv), value = TRUE)
cols_to_remove_B_lobes <- grep("_B$", colnames(scov_lobes_csv), value = TRUE)



# Limpiar Hammers de regiones subcorticales, y valores bilaterales
cols_to_remove <- c("cerebellum_L", "cerebellum_R",
  "caudate_nucleus_L", "caudate_nucleus_R",
  "putamen_L", "putamen_R",
  "thalamus_L", "thalamus_R",
  "TL_anterior_temporal_lobe_lateral_part_L"
)

cols_to_remove_lobes <- c("Cerebellum_L", "Cerebellum_R", "Caudate_L", "Caudate_R", 
                          "Putamen_R", "Putamen_L", "Thalamus_L", "Thalamus_R")

cbf_Hammers_csv <- cbf_Hammers_csv[, !(names(cbf_Hammers_csv) %in% cols_to_remove)]
cbf_Hammers_csv <- cbf_Hammers_csv[, !(colnames(cbf_Hammers_csv) %in% cols_to_remove_B)]

scov_Hammers_csv <- scov_Hammers_csv[, !(names(scov_Hammers_csv) %in% cols_to_remove)]
scov_Hammers_csv <- scov_Hammers_csv[, !(colnames(scov_Hammers_csv) %in% cols_to_remove_B)]

cbf_lobes_csv <- cbf_lobes_csv[, !(names(cbf_lobes_csv) %in% cols_to_remove_lobes)]
cbf_lobes_csv <- cbf_lobes_csv[, !(colnames(cbf_lobes_csv) %in% cols_to_remove_B_lobes)]

scov_lobes_csv <- scov_lobes_csv[, !(names(scov_lobes_csv) %in% cols_to_remove_lobes)]
scov_lobes_csv <- scov_lobes_csv[, !(colnames(scov_lobes_csv) %in% cols_to_remove_B_lobes)]
                               
# n_col_Hammers <- length(cbf_Hammers_csv[14:length(cbf_Hammers_csv)]
# 


# Classify between CBF and Scov (Estoy renombrando todas las columnas, para despues darme cuenta sin CBF o sCOV)
colnames(cbf_Hammers_csv)[14:ncol(cbf_Hammers_csv)] <- paste0(
  colnames(cbf_Hammers_csv)[14:ncol(cbf_Hammers_csv)],
  "_cbf"
)

colnames(cbf_gray_matter_csv)[14:ncol(cbf_gray_matter_csv)] <- paste0(
  colnames(cbf_gray_matter_csv)[14:ncol(cbf_gray_matter_csv)],
  "_cbf"
)

colnames(cbf_lobes_csv)[14:ncol(cbf_lobes_csv)] <- paste0(
  colnames(cbf_lobes_csv)[14:ncol(cbf_lobes_csv)],
  "_cbf"
)

# sCOV
colnames(scov_Hammers_csv)[14:ncol(scov_Hammers_csv)] <- paste0(
  colnames(scov_Hammers_csv)[14:ncol(scov_Hammers_csv)],
  "_scov"
)

colnames(scov_gray_matter_csv)[14:ncol(scov_gray_matter_csv)] <- paste0(
  colnames(scov_gray_matter_csv)[14:ncol(scov_gray_matter_csv)],
  "_scov"
)

colnames(scov_lobes_csv)[14:ncol(scov_lobes_csv)] <- paste0(
  colnames(scov_lobes_csv)[14:ncol(scov_lobes_csv)],
  "_scov"
)

# Merge Cuestionarios Cardiovascular Index
merged_cuestionarios <- merge(
  cuestionarios_excel,
  social_index[, c("ID", "CardiovascularRiskIndex")],  # Solo ID y la columna deseada
  by = "ID",
  all.x = TRUE
)

# Merge CBF
merged_cuestionarios$ID <- paste0("sub-", merged_cuestionarios$ID,"_1")

colnames_cbf_gray_matter <- colnames(cbf_gray_matter_csv)[14:length(colnames(cbf_gray_matter_csv))]
colnames_cbf_Hammers <-  colnames(cbf_Hammers_csv)[14:length(colnames(cbf_Hammers_csv))]
colnames_cbf_lobes <- colnames(cbf_lobes_csv)[14:length(colnames(cbf_lobes_csv))]

# Merge scov
colnames_scov_gray_matter <- colnames(scov_gray_matter_csv)[14:length(colnames(scov_gray_matter_csv))]
colnames_scov_Hammers <-  colnames(scov_Hammers_csv)[14:length(colnames(scov_Hammers_csv))]
colnames_scov_lobes <- colnames(scov_lobes_csv)[14:length(colnames(scov_lobes_csv))]


all_cbf_scov_columns <- c(
  colnames_cbf_gray_matter,
  colnames_cbf_Hammers,
  colnames_cbf_lobes,
  colnames_scov_gray_matter,
  colnames_scov_Hammers,
  colnames_scov_lobes
)

merged_asl <- reduce(
  list(
    cbf_gray_matter_csv[, c("participant_id", colnames_cbf_gray_matter)],
    cbf_Hammers_csv[, c("participant_id", colnames_cbf_Hammers)],
    cbf_lobes_csv[, c("participant_id", colnames_cbf_lobes)],
    scov_gray_matter_csv[, c("participant_id", colnames_scov_gray_matter)],
    scov_Hammers_csv[, c("participant_id", colnames_scov_Hammers)],
    scov_lobes_csv[, c("participant_id", colnames_scov_lobes)]
  ),
  ~ merge(.x, .y, by = "participant_id", all = TRUE)
)


# Merge con cuestionarios - ASL
merged_cuestionarios <- merge(
  merged_cuestionarios,
  merged_asl[, c("participant_id", all_cbf_scov_columns)],
  by.x = "ID",
  by.y = "participant_id",
  all.y = TRUE
)

merged_cuestionarios[merged_cuestionarios == "n/a" | merged_cuestionarios == "NaN"] <- NA

merged_cuestionarios_with_ata <- subset(merged_cuestionarios, !(ID %in% values_to_exclude))
merged_cuestionarios_without_ata <- subset(merged_cuestionarios_with_ata, !(ID %in% vascular_artifacts)) 

merged_cuestionarios_with_ata[all_cbf_scov_columns] <- lapply(merged_cuestionarios_with_ata[all_cbf_scov_columns], function(x) as.numeric(as.character(x)))
merged_cuestionarios_without_ata[all_cbf_scov_columns] <- lapply(merged_cuestionarios_without_ata[all_cbf_scov_columns], function(x) as.numeric(as.character(x)))

cbf_Hammers_csv_clean <- subset(cbf_Hammers_csv, !(participant_id %in% values_to_exclude))
cbf_Hammers_csv_clean <- cbf_Hammers_csv_clean[-1, ]

scov_Hammers_csv_clean <- subset(scov_Hammers_csv, !(participant_id %in% values_to_exclude))
scov_Hammers_csv_clean <- scov_Hammers_csv_clean[-1, ]

cbf_Hammers_csv_clean <- merge(
  cbf_Hammers_csv_clean,
  merged_cuestionarios[, c("ID", "Grupo", "Edad", "Genero")],
  by.x = "participant_id",
  by.y = "ID",
  all.x = TRUE
)

scov_Hammers_csv_clean <- merge(
  scov_Hammers_csv_clean,
  merged_cuestionarios[, c("ID", "Grupo", "Edad", "Genero")],
  by.x = "participant_id",
  by.y = "ID",
  all.x = TRUE
)

cbf_lobes_csv_clean <- subset(cbf_lobes_csv, !(participant_id %in% values_to_exclude))
cbf_lobes_csv_clean <- cbf_lobes_csv_clean[-1, ]

scov_lobes_csv_clean <- subset(scov_lobes_csv, !(participant_id %in% values_to_exclude))
scov_lobes_csv_clean <- scov_lobes_csv_clean[-1, ]



cbf_lobes_csv_clean <- merge(
  cbf_lobes_csv_clean,
  merged_cuestionarios[, c("ID", "Grupo", "Edad", "Genero")],
  by.x = "participant_id",
  by.y = "ID",
  all.x = TRUE
)

scov_lobes_csv_clean <- merge(
  scov_lobes_csv_clean,
  merged_cuestionarios[, c("ID", "Grupo", "Edad", "Genero")],
  by.x = "participant_id",
  by.y = "ID",
  all.x = TRUE
)

# Problematics rows in WMS and Stroop

# Vector con las columnas problemáticas
cols_to_fix <- c("WMS.R.DIR_perc", "WMS.R.INV_perc", "STROOP_P_perc",
                 "STROOP_C_perc", "STROOP_P.C_perc", "STROOP_P.C_INTERF_perc")

# Reemplazar "<5" por "5", ">95" por 95  y convertir a numérico
merged_cuestionarios_with_ata[cols_to_fix] <- lapply(
  merged_cuestionarios_with_ata[cols_to_fix], function(x) {
    x <- as.character(x)
    x[x == "<5"] <- "5"
    x[x == "<2"] <- "2"
    x[x == ">3"] <- "3"
    x[x == "<3"] <- "3"
    x[x == ">95"] <- "95"
    x[x == "<95"] <- "95"
    as.numeric(x)
  }
)

sapply(merged_cuestionarios_with_ata[cols_to_fix], function(x) sum(is.na(x)))
# ------------------  CBF Analysis  ------------------

# multiple_ancova_correcting_scov(merged_cuestionarios_with_ata, colnames_cbf_Hammers)
# multiple_ancova_correcting_scov(merged_cuestionarios_with_ata, colnames_cbf_gray_matter)

model_TotalGM_B_group <- lm(TotalGM_B_scov ~ Edad + Genero + Grupo, data = merged_cuestionarios_with_ata)
print(summary(model_TotalGM_B_group))

# ------------------  sCOV Analysis  ------------------

# ------------------ Regional sCOV --------------------
scov_Hammers_csv_clean$regional_sCOV <- apply(
  cbf_Hammers_csv_clean[, colnames_cbf_Hammers],
  1,
  function(x) sd(as.numeric(x), na.rm = TRUE) / mean(as.numeric(x), na.rm = TRUE)
)

# kruskal.test(regional_sCOV ~ Grupo, data = cbf_Hammers_csv_clean)
# 
# for (region in colnames_scov_Hammers) {
#   formula <- as.formula(paste(region, "~ Grupo"))
#   test_result <- kruskal.test(formula, data = scov_Hammers_csv_clean)
#   
#   if (test_result$p.value < 0.05) {
#     cat("\n", region, "- p-value:", test_result$p.value, "\n")
#   }
# }
# 
# for (region in colnames_scov_Hammers) {
#   formula <- as.formula(paste(region, "~ Grupo"))
#   test_result <- kruskal.test(formula, data = scov_Hammers_csv_clean)
#   
#   if (test_result$p.value < 0.05) {
#     cat("\n", region, "- p-value:", test_result$p.value, "\n")
#   }
# }
# colnames_scov_Hammers <- setdiff(colnames_scov_Hammers, c("TL_superior_temporal_gyrus_middle_part_R_scov", "FL_anterior_orbital_gyrus_L_scov", "FL_anterior_orbital_gyrus_R_scov", "PL_postcentral_gyrus_L_scov"))

scov_Hammers_csv_clean[colnames_scov_Hammers] <- lapply(
  scov_Hammers_csv_clean[colnames_scov_Hammers],
  function(x) as.numeric(as.character(x))
)

scov_lobes_csv_clean[colnames_scov_lobes] <- lapply(
  scov_lobes_csv_clean[colnames_scov_lobes],
  function(x) as.numeric(as.character(x))
)


scov_lobes_csv_clean[, colnames_scov_lobes] <- scov_lobes_csv_clean[, colnames_scov_lobes] * 100

multiple_ancova(scov_Hammers_csv_clean, c(colnames_scov_Hammers, "regional_sCOV"))
multiple_ancova(scov_lobes_csv_clean, c(colnames_scov_lobes))


# 
# 
# for (region in colnames_scov_lobes) {
#   formula <- as.formula(paste(region, "~ Grupo"))
#   test_result <- kruskal.test(formula, data = scov_lobes_csv_clean)
#   
#   cat("\n", region, "- p-value:", test_result$p.value, "\n")
# }

# ------------------ DESCRIPTIVE STATISTICS ------------------
# Indice de riesgo cardiovascular
vars_cardiovascular <- c("PresionAlta", "Diabetes", "ColesterolAlto", "InfartoCardiaco", "AnginaPecho", "EmboliaTrombosis", "MedicacionPresionArterial")
merged_cuestionarios_with_ata$cardiovascular_risk_index <- rowSums(
  merged_cuestionarios_with_ata[vars_cardiovascular] == "SI", 
  na.rm = TRUE
) # Este no lo estoy usando, es para chequear

# Numero por grupo
results_groups <- table(merged_cuestionarios_with_ata$Grupo)
print(results_groups)

# Variables a analizar
vars_numericas <- c("Edad", "BMI")
vars_tests <- c(tests_variables_perc)
vars_wmh <- c("WMHTotal", "WMHPeriventricular","WMHDeep")
vars_qol <- c("FAS","PSQI", "EQ.VAS")
vars_categoricas <- c("Genero")

# Calcular medias y desvíos por grupo para las variables numericas 
stats_numericas <- merged_cuestionarios_with_ata %>%
  group_by(Grupo) %>%
  summarise(across(all_of(vars_numericas),
                   list(mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE)),
                   .names = "{.col}_{.fn}"))

# # Calcular mediana e IQR por grupo para tests cognitivos
# Paso 1: Calcular mediana, Q1 y Q3
iqr_table <- merged_cuestionarios_with_ata %>%
  group_by(Grupo) %>%
  summarise(across(all_of(vars_tests), list(
    med = ~median(., na.rm = TRUE),
    q1  = ~quantile(., 0.25, na.rm = TRUE),
    q3  = ~quantile(., 0.75, na.rm = TRUE)
  ), .names = "{.col}_{.fn}")) %>%
  ungroup()

# # Variables categóricas
# Paso 1: Ver qué variables categóricas tienen NA
na_cp_table <- merged_cuestionarios_with_ata %>%
  dplyr::select(ID, all_of(vars_cardiovascular)) %>%
  pivot_longer(-ID, names_to = "Variable", values_to = "Valor") %>%
  filter(is.na(Valor)) %>%
  group_by(Variable) %>%
  summarise(NA_CP = paste(ID, collapse = ", "), .groups = "drop")

print("Variables con NA:")
print(na_cp_table)

# Paso 2: Calcular N y porcentaje excluyendo NA
cat_table_clean <- merged_cuestionarios_with_ata %>%
  dplyr::select(ID, Grupo, all_of(vars_cardiovascular)) %>%
  pivot_longer(-c(ID, Grupo), names_to = "Variable", values_to = "Valor") %>%
  filter(!is.na(Valor)) %>%
  group_by(Variable, Grupo, Valor) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(Variable, Grupo) %>%
  mutate(Perc = round(100 * N / sum(N), 1)) %>%
  ungroup() %>%
  mutate(Resumen = paste0(N, " (", Perc, "%)")) %>%
  left_join(na_cp_table, by = "Variable")

# Resultado
print(cat_table_clean)



# COMPARACIONES ENTRE GRUPOS
# Kruskal Wallis variables numéricas

# Aplicar Kruskal-Wallisy redondear resultados a dos decimales
kruskal_results <- lapply(c(vars_numericas,vars_tests) , function(var) {
  formula <- as.formula(paste(var, "~ Grupo"))
  test <- kruskal.test(formula, data = merged_cuestionarios_with_ata)
  data.frame(
    Variable = var,
    Statistic = round(test$statistic, 2),
    p_value = round(test$p.value, 2)
  )
})

# Combinar en un único data.frame
kruskal_results_df <- do.call(rbind, kruskal_results)

# Ver resultados
print(kruskal_results_df)


# CHI square para variables categoricas

# Crear tabla vacía para resultados
chi_results <- data.frame(
  Variable = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop para cada variable
for (var in c(vars_categoricas,vars_cardiovascular)) {
  # Crear tabla de contingencia quitando NA
  tab <- table(merged_cuestionarios_with_ata$Grupo, 
               merged_cuestionarios_with_ata[[var]], useNA = "no")
  
  # Verificamos que haya más de una categoría para evitar error
  if (ncol(tab) > 1 && nrow(tab) > 1) {
    test <- chisq.test(tab)
    chi_results <- rbind(chi_results, 
                         data.frame(Variable = var, 
                                    p_value = round(test$p.value, 4)))
  } else {
    chi_results <- rbind(chi_results, 
                         data.frame(Variable = var, 
                                    p_value = NA))
  }
}

# Mostrar resultados
print(chi_results)

# ------------------ Analysis WMH ------------------
merged_cuestionarios_with_ata[vars_wmh] <- lapply(merged_cuestionarios_with_ata[vars_wmh], function(x) as.numeric(as.character(x)))


stats_numericas_wmh <- merged_cuestionarios_with_ata %>%
  group_by(Grupo) %>%
  summarise(across(all_of(vars_wmh),
                   list(mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE)),
                   .names = "{.col}_{.fn}"))
print(stats_numericas_wmh)

kruskal_results_wmh <- lapply(c(vars_wmh) , function(var) {
  formula <- as.formula(paste(var, "~ Grupo"))
  test <- kruskal.test(formula, data = merged_cuestionarios_with_ata)
  data.frame(
    Variable = var,
    Statistic = round(test$statistic, 2),
    p_value = round(test$p.value, 2)
  )
})

# Combinar en un único data.frame
kruskal_results_df_wmh <- do.call(rbind, kruskal_results_wmh)

# Ver resultados
print(kruskal_results_df_wmh)

# ------------------ Correlaciones  ------------------

# Paso 1: Unir variables de interés
vars_correlacion <- c(vars_wmh, vars_tests, vars_qol, "TotalGM_B_scov")

# Subset de datos
cor_data <- merged_cuestionarios_with_ata[, vars_correlacion]

# Paso 2: Calcular matriz de correlación y p-valores con Pearson
rcorr_result <- rcorr(as.matrix(cor_data), type = "pearson")
cor_matrix <- rcorr_result$r
p_matrix <- rcorr_result$P

# Paso 3: Guardar el gráfico como PNG
png(paste0(path_plots_data_analysis, "correlation_matrix_pearson.png"),
    width = 1200, height = 1000, res = 150)

# Paso 4: Graficar corrplot completo (sin filtrar por significancia)
corrplot(cor_matrix,
         method = "color",
         type = "upper",
         order = "hclust",
         tl.col = "black",
         tl.cex = 0.8,
         number.cex = 0.7,
         addCoef.col = "black",
         col = colorRampPalette(c("blue", "white", "red"))(200),
         title = "Pearson Correlation Matrix",
         mar = c(0,0,2,0)
)

dev.off()

# ------------------ Gráficos  ------------------

# Scatter sCOV vs Age
merged_cuestionarios_with_ata$ATA <- ifelse(merged_cuestionarios_with_ata$ID %in% vascular_artifacts, "ATA+", "ATA-")

merged_cuestionarios_with_ata$TotalGM_B_scov <- merged_cuestionarios_with_ata$TotalGM_B_scov * 100

plot_scov_regression_b <- ggplot(merged_cuestionarios_with_ata, aes(x = Edad, y = TotalGM_B_scov, color = Grupo, shape = ATA)) +
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
ggsave(filename = paste0(path_plots_data_analysis, "scatter_Age_vs_sCOV_paper.png"), 
       plot = plot_scov_regression_b, width = 8, height = 6, dpi = 300, bg = "white")

# Con labels de ID 

# Crear el plot con etiquetas de ID
plot_scov_regression_b_labels <- ggplot(merged_cuestionarios_with_ata, 
                                 aes(x = Edad, y = TotalGM_B_scov, color = Grupo)) +
  geom_point(alpha = 0.8, size = 3) +
  geom_smooth(method = "lm", se = TRUE, aes(group = Grupo)) +
  geom_text(aes(label = ID), size = 3, vjust = -0.8) +  # <- etiquetas sin repel
  scale_color_viridis_d() +
  labs(title = "Age vs sCOV",
       x = "Age", y = "sCOV(%)", shape = "Vascular Artifact", color = "Grupo") +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_rect(fill = "white"), 
        plot.background = element_rect(fill = "white"),
        legend.position = "top")

ggsave(filename = paste0(path_plots_data_analysis, "scatter_Age_vs_sCOV_paper_labels.png"), 
       plot = plot_scov_regression_b_labels, width = 8, height = 6, dpi = 300, bg = "white")

# Bar Chart scov lobes 

library(ggsignif)

scov_lobes_tidy <- scov_lobes_csv_clean %>%
  pivot_longer(
    cols = all_of(colnames_scov_lobes),
    names_to = "lobe",
    values_to = "scov"
  )

scov_lobes_tidy <- scov_lobes_tidy %>%
  mutate(
    lobe = str_replace_all(lobe, "_L_scov", " L"),
    lobe = str_replace_all(lobe, "_R_scov", " R"),
    lobe = str_replace_all(lobe, "_scov", "")  # por si quedan nombres sin hemisferio
  )

# Crear el data.frame de significancia
signif_df <- data.frame(
  lobe = c("Frontal L", "Insula L", "Occipital R", "Parietal L", "Parietal R", "Temporal L"),
  pval = c(0.02, 0.001421, 0.0106777, 0.0445590, 0.030872, 0.0209279)
)

# Etiqueta de significancia
signif_df$label <- ifelse(signif_df$pval < 0.01, "**", "*")

# Asignar posiciones X manualmente (asumiendo orden en gráfico)
# Usamos índice base 1 con barras agrupadas → 0.35 de separación aprox.
signif_df$x <- match(signif_df$lobe, levels(factor(scov_lobes_tidy$lobe)))
signif_df$xmin <- signif_df$x - 0.175
signif_df$xmax <- signif_df$x + 0.175

# Agregar altura de línea (ejemplo base, después se puede ajustar dinámicamente)
# Se puede sobrescribir con media + sd + margen
signif_df$y_pos <- c(63, 67, 69, 71, 73, 75)

# El gráfico
ggplot(scov_lobes_tidy, aes(x = lobe, y = scov, fill = Grupo)) +
  
  # Barras de media
  stat_summary(fun = mean, geom = "bar",
               position = position_dodge(width = 0.7),
               width = 0.7,
               color = NA, alpha = 0.85) +
  
  # Barras de error
  stat_summary(fun.data = mean_sdl,
               fun.args = list(mult = 1),
               geom = "errorbar",
               position = position_dodge(width = 0.7),
               width = 0.15,
               size = 0.6,
               color = "black") +
  
  # Eje Y desde 0, sin margen inferior
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  # Paleta de colores
  scale_fill_manual(values = c("CONTROL" = "#77AADD", "COVID" = "#FF7744")) +
  
  # Etiquetas
  labs(
    title = "Regional variability (sCOV) by brain lobe",
    x = "Brain lobe",
    y = "sCOV (%)",
    fill = "Group"
  ) +
  
  # Asteriscos de significancia
  geom_segment(data = signif_df,
               aes(x = xmin, xend = xmax,
                   y = y_pos, yend = y_pos),
               inherit.aes = FALSE,
               color = "black", size = 0.6) +
  
  geom_text(data = signif_df,
            aes(x = (xmin + xmax)/2,
                y = y_pos + 1.5,
                label = label),
            inherit.aes = FALSE,
            size = 5, vjust = 0) +
  
  # Tema
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title.x = element_text(face = "bold", size = 10, vjust = -0.2),
    axis.title.y = element_text(face = "bold", size = 10, vjust = 1.5),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    legend.position = "top",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    axis.line = element_line(color = "black"),
    panel.border = element_blank()
  )


# plot_gm_l <- ggplot(merged_df_perfusion_with_ata, aes(x = Grupo, y = get("TotalGM_L_scov"), fill = Grupo)) +
#   geom_boxplot(alpha = 0.5) +  # Add transparency to boxplot
#   geom_jitter(color = "black", size = 2, alpha = 0.3, width = 0.2) +  # Add jittered points
#   labs(title = "Left Grey Matter", x = "Group", y = "sCOV (%)")+
#   theme_minimal() +
#   ylim(30, 100) +
#   guides(fill = FALSE) +  # Remove legend
#   theme(
#     plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),  # Increase title size
#     axis.text = element_text(size = 14),  # Increase axis text size
#     axis.title = element_text(size = 16)  # Increase axis title size
#   )
# 
# # Save the plot
# plot_filename <- paste(path_plots_data_analysis, "Total_GM_L_ISMRM.jpg", sep="")
# ggsave(plot_filename, plot = plot_gm_l , device = "jpeg", bg = "white", width = 4, height = 7)
# 
# 
# plot_gm_r <- ggplot(merged_df_perfusion_with_ata, aes(x = Grupo, y = get("TotalGM_R_scov"), fill = Grupo)) +
#   geom_boxplot(alpha = 0.5) +  # Add transparency to boxplot
#   geom_jitter(color = "black", size = 2, alpha = 0.3, width = 0.2) +  # Add jittered points
#   ylim(30, 100) +
#   labs(title = "Right Grey Matter", x = "Group", y = "sCOV (%)") +
#   theme_minimal() +
#   guides(fill = FALSE) +  # Remove legend
#   theme(
#     plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),  # Increase title size
#     axis.text = element_text(size = 14),  # Increase axis text size
#     axis.title = element_text(size = 16)  # Increase axis title size
#   )
# 
# # Save the plot
# plot_filename <- paste(path_plots_data_analysis, "Total_GM_R_ISMRM.jpg", sep="")
# ggsave(plot_filename, plot = plot_gm_r , device = "jpeg", bg = "white", width = 4, height = 7)

# Tests

# WMH 