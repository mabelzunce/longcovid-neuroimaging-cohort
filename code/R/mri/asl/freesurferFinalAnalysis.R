# Load necessary libraries
library(tidyverse)
library(readxl)
library(car)
library(gridExtra)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(effects)

source("functionsANCOVA.R")

# Define the paths
path_data_analysis <- "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Freesurfer/DataAnalysis/"
#path_data_analysis <- "/home/solcat/work/data_analysis_paper/"
path_results_data_analysis <- paste(path_data_analysis, "Results/", sep="")
path_plots_data_analysis <-  paste(path_data_analysis, "Plots/", sep="")

path_cuestionarios <- paste(path_data_analysis, "ResumenCuestionarios_30_01.csv", sep="")
path_sienax2 <-paste(path_data_analysis, "SienaxResults2.csv", sep="")
path_brain_vol <- paste(path_data_analysis, "brain_volumes.csv", sep="")
path_segmentation <- paste(path_data_analysis, "segmentation.csv", sep="")
path_parcellation <- paste(path_data_analysis, "parcellation.csv", sep="")
path_parcellation_DKT <- paste(path_data_analysis, "parcellation_DKT.csv", sep="")
path_parcellation_thick <- paste(path_data_analysis, "parcellation_thick.csv", sep="")
path_parcellation_thick_dkt <- paste(path_data_analysis, "parcellation_thick_DKT.csv", sep="")

# Read the Excel file
# cuestionarios_excel <- read_excel(path_cuestionarios)

# Read files
cuestionarios_excel <- read.csv(path_cuestionarios)
sienax2_csv <-read.csv(path_sienax2,dec = ",")
brain_vol_csv <- read.csv(path_brain_vol)
segmentation_csv <- read.csv(path_segmentation)
parcellation_csv <- read.csv(path_parcellation)
parcellation_DKT_csv <- read.csv(path_parcellation_DKT)
parcellation_thick_csv <- read.csv(path_parcellation_thick)
parcellation_thick_DKT_csv <- read.csv(path_parcellation_thick_dkt)

# Exclude rows
values_to_exclude <- c('CP0011', 'CP0015', 'CP0035', 'CP0144', 'CP0106', 'CP0192', 'CP0193', 'CP0087','CP0196')

# Rename columns 
# colnames(cuestionarios_excel)[colnames(cuestionarios_excel) == "FAS-puntuación final"] <- "FAS_puntuacion_final"
# colnames(cuestionarios_excel)[colnames(cuestionarios_excel) == "IPAQ-puntuacion"] <- "IPAQ_puntuacion"
# colnames(cuestionarios_excel)[colnames(cuestionarios_excel) == "PSQI-puntuacion"] <- "PSQI_puntuacion"
# colnames(cuestionarios_excel)[colnames(cuestionarios_excel) == "TMTA-valor de referencia"] <- "TMTA_valor_de_referencia"
# colnames(cuestionarios_excel)[colnames(cuestionarios_excel) == "TMTB-valor de referencia"] <- "TMTB_valor_de_referencia"
# colnames(cuestionarios_excel)[colnames(cuestionarios_excel) == "WMS-R-directos valor de referencia"] <- "WMS_R_directos_valor_de_referencia"
# colnames(cuestionarios_excel)[colnames(cuestionarios_excel) == "WMS-R-inversos valor de referencia"] <- "WMS_R_inversos_valor_de_referencia"
# colnames(cuestionarios_excel)[colnames(cuestionarios_excel) == "STROOP-P/C valor de referencia"] <- "STROOP_P_C_valor_de_referencia"
# colnames(cuestionarios_excel)[colnames(cuestionarios_excel) == "STROOP-interf valor de referencia"] <- "STROOP_interf_valor_de_referencia"
# colnames(cuestionarios_excel)[colnames(cuestionarios_excel) == "MOCA valor de referencia"] <- "MOCA_valor_de_referencia"



# Relevant columns 
# related_symptoms <- colnames(cuestionarios_excel)[26:37]
# standardized_questionnaire <- colnames(cuestionarios_excel)[c(63, 76, 89, 122)]
# cognitive_test <- colnames(cuestionarios_excel)[c(137, 140, 143,  146, 155, 158, 161)]

# Sienax
# Reemplazar las comas por puntos en las columnas numéricas
# sienax2_csv[] <- lapply(sienax2_csv, function(x) {
#   if (is.character(x)) {
#     # Reemplazar las comas por puntos y convertir a numérico
#     as.numeric(gsub(",", ".", x))
#   } else {
#     x
#   }
# })


merged_df_sienax2 <- merge(cuestionarios_excel[, c("Edad", "Grupo", "ID","Genero")], sienax2_csv, by.x="ID", by.y="SubjectNames", all.y=TRUE)
filtered_sienax2 <- subset(merged_df_sienax2, !(ID %in% values_to_exclude))
column_names_sienax2 <- names(sienax2_csv)
column_names_sienax2 <- column_names_sienax2[2:(length(column_names_sienax2))]


ggplot(filtered_sienax2, aes(x = Grupo , y = GreyMatterVolume, fill = Grupo)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2) +
  geom_text(aes(label = ID), vjust = -0.5, size = 3, check_overlap = TRUE) +
  theme(legend.position="top")

# # perform ancova for each dependent variable - Sienax2
for(col_name in column_names_sienax2) {
  print(col_name)
  formula <- as.formula(paste(col_name, "~ Grupo + Edad"))

  # perform ancova with group_variable as a factor and age as covariate
  ancova_result <- lm(formula, data = filtered_sienax2)
  model <- Anova(ancova_result, type="III")

  # Extract p-value from the ANOVA results
  p_value <- model$"Pr(>F)"[2]

  if (p_value < 0.05){
    print(model)
  }

}



# BrainVol and Sienax

# Merge both data frames
merged_df_brainvol <- merge(cuestionarios_excel, brain_vol_csv, by.x="ID", by.y="subject", all.y=TRUE)
filtered_brainvol <- subset(merged_df_brainvol, !(ID %in% values_to_exclude))

# Columns Brainvol -> delete subject, data and Units
column_names_brainvol <- names(brain_vol_csv)
column_names_brainvol <- column_names_brainvol[2:(length(column_names_brainvol)-2)]


# perform ancova for each dependent variable - Brainvol
for(col_name in column_names_brainvol) {
  formula <- as.formula(paste(col_name, "~ Grupo + Edad"))
  
  # perform ancova with group_variable as a factor and age as covariate
  ancova_result <- lm(formula, data = filtered_brainvol)
  model <- Anova(ancova_result, type="III")
    # Extract p-value from the ANOVA results
  
  p_value <- model$"Pr(>F)"[2]

  if (p_value < 0.05){
    print(model)
  }
  
}


# Segmentation
# Exclude Subject
filtered_df_segmentation <- subset(segmentation_csv, !(subject %in% values_to_exclude))

# Etiv from brainvol
df_segmentation <- cbind(filtered_df_segmentation, filtered_brainvol$Estimated.Total.Intracranial.Volume)
colnames(df_segmentation)[ncol(df_segmentation)] <-"Estimated.Total.Intracranial.Volume"

# Edad, Grupo from cuestionarios_excel
merged_df_segmentation <- merge(df_segmentation, cuestionarios_excel[, c("ID", "Edad", "Grupo")], by.x = "subject", by.y = "ID", all.x = TRUE)

column_names_segmentation <- names(merged_df_segmentation)
column_names_segmentation <- column_names_segmentation[2:(length(column_names_segmentation) - 6 )]

right_regions <- column_names_segmentation[grepl("^Right(?!.*hypointensities)", column_names_segmentation,ignore.case = TRUE, perl = TRUE)]
left_regions <- column_names_segmentation[grepl("^Left(?!.*hypointensities)", column_names_segmentation,ignore.case = TRUE, perl = TRUE)]

only_subcortical_regions <- c(right_regions, left_regions)

significant_ROIs_not_age <- list()
significant_ROIs_not_etiv <- list()
significant_ROIs_not_etiv_not_age <- list()

# Create an empty data frame to store significant regions and p-values
p_values_segmentation_df_right <- data.frame(region = character(),
                          mean_control = numeric(), sd_control = numeric(), 
                          mean_covid = numeric(), sd_covid = numeric(), p_value = numeric(), ancova=character() , stringsAsFactors = FALSE)

p_values_segmentation_df_leftt <- data.frame(region = character(),
                                             mean_control = numeric(), sd_control = numeric(), 
                                             mean_covid = numeric(), sd_covid = numeric(), p_value = numeric(), ancova=character() , stringsAsFactors = FALSE)


left_analysis_results <- perform_analysis(left_regions, merged_df_segmentation, p_values_segmentation_df_leftt, "left")
right_analysis_results <- perform_analysis(right_regions, merged_df_segmentation, p_values_segmentation_df_right, "right")

p_values_df_segmentation_left <- left_analysis_results$df
p_values_df_segmentation_right <- right_analysis_results$df

results_segmentation_left <- paste(path_results_data_analysis, "left_segmentation_ancova.csv", sep="")
results_segmentation_right <- paste(path_results_data_analysis, "right_segmentation_ancova.csv", sep="")

p_values_segmentation_df <- rbind(p_values_df_segmentation_left, p_values_df_segmentation_right)

write.csv(p_values_df_segmentation_left , results_segmentation_left, row.names = FALSE)
write.csv(p_values_df_segmentation_right , results_segmentation_right, row.names = FALSE)



# Parcellation DKT

filtered_df_parcellation_DKT <- subset(parcellation_DKT_csv, !(subject %in% values_to_exclude))

# Etiv from brainvol
merge_df_parcellation_DKT <- cbind(filtered_df_parcellation_DKT, filtered_brainvol$Estimated.Total.Intracranial.Volume)
colnames(merge_df_parcellation_DKT)[ncol(merge_df_parcellation_DKT)] <-"Estimated.Total.Intracranial.Volume"

# Edad, Grupo from cuestionarios_excel
merge_df_parcellation_DKT  <- merge(merge_df_parcellation_DKT , cuestionarios_excel[, c("Edad", "Grupo", "ID")], by.x = "subject", by.y = "ID", all.x = TRUE)

# Mean both hemispheres
mean_df_parcellation_DKT <- filtered_df_parcellation_DKT %>%
  group_by(subject) %>%
  summarise_all(mean)

mean_df_parcellation_DKT <- mean_df_parcellation_DKT[1:35]

# Etiv from brainvol
mean_df_parcellation_DKT <- cbind(mean_df_parcellation_DKT, filtered_brainvol$Estimated.Total.Intracranial.Volume)
colnames(mean_df_parcellation_DKT)[ncol(mean_df_parcellation_DKT)] <-"Estimated.Total.Intracranial.Volume"

# Edad, Grupo from cuestionarios_excel
mean_df_parcellation_DKT  <- merge(mean_df_parcellation_DKT , cuestionarios_excel[, c("Edad", "Grupo", "ID")], by.x = "subject", by.y = "ID", all.x = TRUE)

column_names_parcellation_DKT  <- names(mean_df_parcellation_DKT)
column_names_parcellation_DKT <- column_names_parcellation_DKT [2:(length(column_names_parcellation_DKT) - 4)]

# Perform ancova for each dependent ROI - parcellation - DKT 
significant_ROIs_not_age <- list()

# Only Right structures 
right_df_parcellation_dkt <- merge_df_parcellation_DKT %>%
  filter(Hemisphere == "rh")

# Only Left structures 
left_df_parcellation_dkt <- merge_df_parcellation_DKT %>%
  filter(Hemisphere == "lh")


# Create an empty data frame to store significant regions and p-values
p_values_df_parcellation_dkt <- data.frame(region = character(),
                                           mean_control = numeric(), sd_control = numeric(), 
                                           mean_covid = numeric(), sd_covid = numeric(), f_value = numeric(), p_value = numeric(), ancova=character() , stringsAsFactors = FALSE)


p_values_df_parcellation_dkt_left <- data.frame(region = character(),
                                                mean_control = numeric(), sd_control = numeric(), 
                                                mean_covid = numeric(), sd_covid = numeric(), f_value = numeric(), p_value = numeric(), ancova=character() , stringsAsFactors = FALSE)



p_values_df_parcellation_dkt_right <- data.frame(region = character(),
                                                 mean_control = numeric(), sd_control = numeric(), 
                                                 mean_covid = numeric(), sd_covid = numeric(), f_value = numeric(), p_value = numeric(), ancova=character() , stringsAsFactors = FALSE)


# Perform analysis on mean both hemispheres
mean_analysis_results <- perform_analysis(column_names_parcellation_DKT, mean_df_parcellation_DKT, p_values_df_parcellation_dkt, "mean")

# Left
left_analysis_results <- perform_analysis(column_names_parcellation_DKT, left_df_parcellation_dkt, p_values_df_parcellation_dkt_left, "left")

# Right
right_analysis_results <- perform_analysis(column_names_parcellation_DKT, right_df_parcellation_dkt, p_values_df_parcellation_dkt_right, "right")

# Extract results
p_values_df_parcellation_dkt <- mean_analysis_results$df
p_values_df_parcellation_dkt_left <- left_analysis_results$df
p_values_df_parcellation_dkt_right <- right_analysis_results$df

# Merge all volumes 

p_values_df_parcellation_dkt_left$region <- paste0(p_values_df_parcellation_dkt_left$region, "_left")
p_values_df_parcellation_dkt_right$region <- paste0(p_values_df_parcellation_dkt_right$region, "_right")

p_values_volumes_df <- rbind(p_values_df_parcellation_dkt_left, p_values_df_parcellation_dkt_right)
p_values_volumes_df <- rbind(p_values_volumes_df, p_values_segmentation_df)



results_volumes <- paste(path_results_data_analysis, "results_volume.csv", sep="")
write.csv(p_values_volumes_df, results_volumes, row.names = FALSE)


# Cortical Thickness DKT
filtered_df_parcellation_thick_dkt <- subset(parcellation_thick_DKT_csv, !(subject %in% values_to_exclude))

# Mean both hemispheres
mean_df_parcellation_thick_dkt <- filtered_df_parcellation_thick_dkt  %>%
  group_by(subject) %>%
  summarise_all(mean)

# Edad, Grupo from cuestionarios_excel
mean_df_parcellation_thick_dkt  <- merge(mean_df_parcellation_thick_dkt , cuestionarios_excel[, c("Edad", "Grupo", "ID")], by.x = "subject", by.y = "ID", all.x = TRUE)

column_names_parcellation_thick_dkt  <- names(mean_df_parcellation_thick_dkt)
column_names_parcellation_thick_dkt <- column_names_parcellation_thick_dkt[2:(length(column_names_parcellation_thick_dkt) - 5)]

# Perform t-test and save in a csv

# Create an empty data frame to store significant regions and p-values
p_values_df_parcellation_thickness_dkt <- data.frame(region = character(),
                                                     mean_control = numeric(), sd_control = numeric(), 
                                                     mean_covid = numeric(), sd_covid = numeric(), p_value = numeric(), t_statistics = numeric(), df = numeric(), 
                                                     stringsAsFactors = FALSE)

significant_ROIs_cort_thick_dkt <- list()


for(col_name in column_names_parcellation_thick_dkt) {
  
  values_covid <- mean_df_parcellation_thick_dkt[mean_df_parcellation_thick_dkt$Grupo == "COVID", col_name]
  values_control <- mean_df_parcellation_thick_dkt[mean_df_parcellation_thick_dkt$Grupo == "CONTROL", col_name]
  
  # Calculate mean and standard deviation for both groups
  mean_covid <- mean(values_covid , na.rm = TRUE)
  mean_control <- mean(values_control, na.rm = TRUE)
  sd_covid <- sd(values_covid, na.rm = TRUE)
  sd_control <- sd(values_control, na.rm = TRUE)
  
  # Test T 
  result_t_thick_dkt <- t.test(values_covid, values_control)
  p_value <- result_t_thick_dkt$p.value
  degrees_freedom <- result_t_thick_dkt$parameter
  t_stat <- result_t_thick_dkt$statistic
  
  p_values_df_parcellation_thickness_dkt <- bind_rows(p_values_df_parcellation_thickness_dkt, data.frame(region = col_name, mean_covid = mean_covid, sd_covid = sd_covid, 
                                                                                                         mean_control = mean_control, sd_control = sd_control, 
                                                                                                         p_value = p_value, t_statistics = t_stat, df = degrees_freedom))
  
  if (p_value < 0.05) {
    significant_ROIs_cort_thick_dkt <- append(significant_ROIs_cort_thick_dkt, col_name)
  }
}



# Only Left Structures
left_df_parcellation_thick_dkt <- filtered_df_parcellation_thick_dkt %>%
  filter(Hemisphere == "lh")

# Right Structures
right_df_parcellation_thick_dkt <- filtered_df_parcellation_thick_dkt %>%
  filter(Hemisphere == "rh")

left_df_parcellation_thick_dkt  <- merge(left_df_parcellation_thick_dkt, cuestionarios_excel[, c("Edad", "Grupo", "ID")], by.x = "subject", by.y = "ID", all.x = TRUE)
right_df_parcellation_thick_dkt  <- merge(right_df_parcellation_thick_dkt, cuestionarios_excel[, c("Edad", "Grupo", "ID")], by.x = "subject", by.y = "ID", all.x = TRUE)

# Create an empty data frame to store significant regions and p-values
p_values_df_parcellation_thickness_dkt_left <- data.frame(region = character(),
                                                          mean_control = numeric(), sd_control = numeric(), 
                                                          mean_covid = numeric(), sd_covid = numeric(), t_statistics = numeric(), df = numeric(), 
                                                          p_value = numeric(), stringsAsFactors = FALSE)

p_values_df_parcellation_thickness_dkt_right <- data.frame(region = character(),
                                                           mean_control = numeric(), sd_control = numeric(), 
                                                           mean_covid = numeric(), sd_covid = numeric(), t_statistics = numeric(), df = numeric(), 
                                                           p_value = numeric(), stringsAsFactors = FALSE)

significant_ROIs_cort_thick_dkt_left <- list()
significant_ROIs_cort_thick_dkt_right <- list()

for(col_name in column_names_parcellation_thick_dkt) {
  values_covid_left <- left_df_parcellation_thick_dkt[left_df_parcellation_thick_dkt$Grupo == "COVID", col_name]
  values_control_left <- left_df_parcellation_thick_dkt[left_df_parcellation_thick_dkt$Grupo == "CONTROL", col_name]
  
  values_covid_right <- right_df_parcellation_thick_dkt[right_df_parcellation_thick_dkt$Grupo == "COVID", col_name]
  values_control_right <- right_df_parcellation_thick_dkt[right_df_parcellation_thick_dkt$Grupo == "CONTROL", col_name]
  
  # Calculate mean and standard deviation for both groups
  mean_covid_left <- mean(values_covid_left , na.rm = TRUE)
  mean_control_left <- mean(values_control_left, na.rm = TRUE)
  sd_covid_left <- sd(values_covid_left, na.rm = TRUE)
  sd_control_left <- sd(values_control_left, na.rm = TRUE)
  
  mean_covid_right <- mean(values_covid_right , na.rm = TRUE)
  mean_control_right <- mean(values_control_right, na.rm = TRUE)
  sd_covid_right <- sd(values_covid_right, na.rm = TRUE)
  sd_control_right <- sd(values_control_right, na.rm = TRUE)
  
  # Test T - Left
  result_t_thick_left <- t.test(values_covid_left, values_control_left)
  p_value_left <- result_t_thick_left$p.value
  degrees_freedom_left <- result_t_thick_left$parameter
  t_stat_left <- result_t_thick_left$statistic
  
  # Test T - Right
  result_t_thick_right <- t.test(values_covid_right, values_control_right)
  p_value_right <- result_t_thick_right$p.value
  degrees_freedom_right <- result_t_thick_right$parameter
  t_stat_right <- result_t_thick_right$statistic
  
  p_values_df_parcellation_thickness_dkt_left <- bind_rows(p_values_df_parcellation_thickness_dkt_left, data.frame(region = col_name, mean_covid = mean_covid_left, sd_covid = sd_covid_left, 
                                                                                                                   mean_control = mean_control_left, sd_control = sd_control_left, 
                                                                                                                   t_statistics = t_stat_left, 
                                                                                                                   df = degrees_freedom_left, 
                                                                                                                   p_value = p_value_left) )
  
  p_values_df_parcellation_thickness_dkt_right <- bind_rows( p_values_df_parcellation_thickness_dkt_right, data.frame(region = col_name, mean_covid = mean_covid_right, sd_covid = sd_covid_right, 
                                                                                                                      mean_control = mean_control_right, sd_control = sd_control_right, 
                                                                                                                      t_statistics = t_stat_right, 
                                                                                                                      df = degrees_freedom_right, 
                                                                                                                      p_value = p_value_right) )
  
  if (p_value_left < 0.05) {
    significant_ROIs_cort_thick_dkt_left <- append(significant_ROIs_cort_thick_dkt_left, col_name)
  }
  
  if (p_value_right < 0.05) {
    significant_ROIs_cort_thick_dkt_right <- append(significant_ROIs_cort_thick_dkt_right, col_name)
  }
}



results_parcellation_cortical_thickness_dkt <- paste(path_results_data_analysis, "mean_parcellation_cort_thickness_dkt.csv", sep="")
write.csv(p_values_df_parcellation_thickness_dkt, results_parcellation_cortical_thickness_dkt, row.names = FALSE)

results_parcellation_cortical_thickness_left_dkt <- paste(path_results_data_analysis, "left_parcellation_cort_thickness_dkt.csv", sep="")
write.csv(p_values_df_parcellation_thickness_dkt_left, results_parcellation_cortical_thickness_left_dkt, row.names = FALSE)

results_parcellation_cortical_thickness_right_dkt <- paste(path_results_data_analysis, "right_parcellation_cort_thickness_dkt.csv", sep="")
write.csv(p_values_df_parcellation_thickness_dkt_right, results_parcellation_cortical_thickness_right_dkt, row.names = FALSE)


## Boxplot
boxplots_paper <- list()

boxplots_paper[[1]] <- ggplot(merged_df_segmentation, aes(x = Grupo, y =  Left.Cerebellum.Cortex, fill = Grupo)) +
  geom_boxplot(alpha = 0.5) +  # Add transparency to boxplot
  geom_jitter(color = "black", size = 2, alpha = 0.3, width = 0.2) +  # Add jittered points
  labs(title = "Left Cerebellum Cortex",
       x = "Group",
       y = "Volume (cm³)") +
  theme_minimal() +
  guides(fill = FALSE) +# Remove legend
  theme(
    plot.title = element_text(size = 10, hjust=0.5),  # Cambia el tamaño del título
  )




boxplots_paper[[2]] <- ggplot(merged_df_segmentation, aes(x = Grupo, y =  Right.Cerebellum.Cortex, fill = Grupo)) +
  geom_boxplot(alpha = 0.5) +  # Add transparency to boxplot
  geom_jitter(color = "black", size = 2, alpha = 0.3, width = 0.2) +  # Add jittered points
  labs(title = "Right Cerebellum Cortex",
       x = "Group",
       y = "Volume (cm³)") +
  theme_minimal() +
  guides(fill = FALSE) +# Remove legend
  theme(
    plot.title = element_text(size = 10, hjust=0.5),  # Cambia el tamaño del título
  )

grid_arrange <- do.call(gridExtra::grid.arrange, c(boxplots_paper, ncol = 2))

ggsave(paste(path_plots_data_analysis, 'final_results_cerebellum.png'), grid_arrange, width = 6, height = 6)



# Cortical thickness

boxplots_paper <- list()

boxplots_paper[[1]] <- ggplot(left_df_parcellation_dkt, aes(x = Grupo, y = lingual, fill = Grupo)) +
  geom_boxplot(alpha = 0.5) +  # Add transparency to boxplot
  geom_jitter(color = "black", size = 2, alpha = 0.3, width = 0.2) +  # Add jittered points
  labs(title = "Left Lingual Gyrus",
       x = "Group",
       y = "Volume (cm³)") +
  theme_minimal() +
  guides(fill = FALSE) +# Remove legend
  theme(
    plot.title = element_text(size = 10, hjust=0.5),  # Cambia el tamaño del título
  )



boxplots_paper[[2]] <- ggplot(left_df_parcellation_dkt, aes(x = Grupo, y =  postcentral, fill = Grupo)) +
  geom_boxplot(alpha = 0.5) +  # Add transparency to boxplot
  geom_jitter(color = "black", size = 2, alpha = 0.3, width = 0.2) +  # Add jittered points
  labs(title = "Left Postcentral Gyrus",
       x = "Group",
       y = "Volume (cm³)") +
  theme_minimal() +
  guides(fill = FALSE) +# Remove legend
  theme(
    plot.title = element_text(size = 10, hjust=0.5),  # Cambia el tamaño del título
  )


boxplots_paper[[3]] <- ggplot(left_df_parcellation_dkt, aes(x = Grupo, y =  supramarginal, fill = Grupo)) +
  geom_boxplot(alpha = 0.5) +  # Add transparency to boxplot
  geom_jitter(color = "black", size = 2, alpha = 0.3, width = 0.2) +  # Add jittered points
  labs(title = "Left Supramarginal\nGyrus",
       x = "Group",
       y = "Volume (cm³)") +
  theme_minimal() +
  guides(fill = FALSE) +# Remove legend
  theme(
    plot.title = element_text(size = 10, hjust=0.5),  # Cambia el tamaño del título
  )


grid_arrange <- do.call(gridExtra::grid.arrange, c(boxplots_paper, ncol = 3))

ggsave(paste(path_plots_data_analysis, 'final_results_left_thickness.png'), grid_arrange, width = 5, height = 5)




