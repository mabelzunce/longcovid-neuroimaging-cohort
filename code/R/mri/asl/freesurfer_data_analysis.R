# Load necessary libraries
library(tidyverse)
library(readxl)
library(car)
library(gridExtra)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(effects)


# Define the paths
path_data_analysis <- "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Freesurfer/DataAnalysis"
#path_data_analysis <- "/home/solcat/work/data_analysis_paper/"
path_results_data_analysis <- paste(path_data_analysis, "Results/", sep="")
path_plots_data_analysis <-  paste(path_data_analysis, "Plots/", sep="")


path_cuestionarios <- paste(path_data_analysis, "RespuestasCuestionarioEvaluacionCognitivaResonancia.xlsx", sep="")
#path_sienax2 <-paste(path_data_analysis, "SienaxResults2.csv", sep="")
path_brain_vol <- paste(path_data_analysis, "brain_volumes.csv", sep="")
path_segmentation <- paste(path_data_analysis, "segmentation.csv", sep="")

path_parcellation <- paste(path_data_analysis, "parcellation.csv", sep="")
path_parcellation_DKT <- paste(path_data_analysis, "parcellation_DKT.csv", sep="")
path_parcellation_thick <- paste(path_data_analysis, "parcellation_thick.csv", sep="")
path_parcellation_thick_dkt <- paste(path_data_analysis, "parcellation_thick_DKT.csv", sep="")

# Read the Excel file
cuestionarios_excel <- read_excel(path_cuestionarios)

# Read CSV files
#sienax2_csv <-read.csv(path_sienax2)
brain_vol_csv <- read.csv(path_brain_vol)
segmentation_csv <- read.csv(path_segmentation)
parcellation_csv <- read.csv(path_parcellation)
parcellation_DKT_csv <- read.csv(path_parcellation_DKT)
parcellation_thick_csv <- read.csv(path_parcellation_thick)
parcellation_thick_DKT_csv <- read.csv(path_parcellation_thick_dkt)

# Exclude rows
values_to_exclude <- c('CP0011', 'CP0015', 'CP0035', 'CP0144', 'CP0106')



# Rename columns 
colnames(cuestionarios_excel)[colnames(cuestionarios_excel) == "FAS-puntuación final"] <- "FAS_puntuacion_final"
colnames(cuestionarios_excel)[colnames(cuestionarios_excel) == "IPAQ-puntuacion"] <- "IPAQ_puntuacion"
colnames(cuestionarios_excel)[colnames(cuestionarios_excel) == "PSQI-puntuacion"] <- "PSQI_puntuacion"
colnames(cuestionarios_excel)[colnames(cuestionarios_excel) == "TMTA-valor de referencia"] <- "TMTA_valor_de_referencia"
colnames(cuestionarios_excel)[colnames(cuestionarios_excel) == "TMTB-valor de referencia"] <- "TMTB_valor_de_referencia"
colnames(cuestionarios_excel)[colnames(cuestionarios_excel) == "WMS-R-directos valor de referencia"] <- "WMS_R_directos_valor_de_referencia"
colnames(cuestionarios_excel)[colnames(cuestionarios_excel) == "WMS-R-inversos valor de referencia"] <- "WMS_R_inversos_valor_de_referencia"
colnames(cuestionarios_excel)[colnames(cuestionarios_excel) == "STROOP-P/C valor de referencia"] <- "STROOP_P_C_valor_de_referencia"
colnames(cuestionarios_excel)[colnames(cuestionarios_excel) == "STROOP-interf valor de referencia"] <- "STROOP_interf_valor_de_referencia"
colnames(cuestionarios_excel)[colnames(cuestionarios_excel) == "MOCA valor de referencia"] <- "MOCA_valor_de_referencia"



# Relevant columns 
related_symptoms <- colnames(cuestionarios_excel)[26:37]
standardized_questionnaire <- colnames(cuestionarios_excel)[c(63, 76, 89, 122)]
cognitive_test <- colnames(cuestionarios_excel)[c(137, 140, 143,  146, 155, 158, 161)]

# BrainVol and Sienax

# Merge both data frames
merged_df_sienax2 <- merge(cuestionarios_excel, sienax2_csv, by.x="ID", by.y="SubjectNames", all.y=TRUE)
merged_df_brainvol <- merge(cuestionarios_excel, brain_vol_csv, by.x="ID", by.y="subject", all.y=TRUE)


filtered_sienax2 <- subset(merged_df_sienax2, !(ID %in% values_to_exclude))
filtered_brainvol <- subset(merged_df_brainvol, !(ID %in% values_to_exclude))

# Columns Brainvol -> delete subject, data and Units
column_names_brainvol <- names(brain_vol_csv)
column_names_sienax2 <- names(sienax2_csv)

column_names_brainvol <- column_names_brainvol[2:(length(column_names_brainvol)-2)]
column_names_sienax2 <- column_names_sienax2[2:(length(column_names_sienax2))]



ggplot(filtered_sienax2, aes(x = Grupo , y = GreyMatterVolume, fill = Grupo)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2) +
  theme(legend.position="top")

# Independence Assumption (The covariate and the treatment are independent)
anova_result <- aov(Edad ~ Grupo, data = filtered_sienax2)
summary(anova_result)

ancova_result <- aov(GreyMatterVolume ~ Grupo , data = filtered_sienax2)
Anova(ancova_result, type="III")


# perform ancova for each dependent variable - Sienax2
for(col_name in column_names_sienax2) {
  formula <- as.formula(paste(col_name, "~ Grupo + Edad"))
  
  # perform ancova with group_variable as a factor and age as covariate
  ancova_result <- aov(formula, data = filtered_sienax2)
  model <- Anova(ancova_result, type="III")
  
  # Extract p-value from the ANOVA results
  p_value <- model$"Pr(>F)"[2]
  
  print(model)
  # if (p_value < 0.05){
  #   print(model)
  # }
  
}

# perform ancova for each dependent variable - Brainvol
for(col_name in column_names_brainvol) {
  formula <- as.formula(paste(col_name, "~ Grupo + Edad + TMTA_valor_de_referencia"))
  
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

# Independence Assumption (The covariate and the treatment are independent)
anova_result <- aov(Estimated.Total.Intracranial.Volume ~ Grupo, data =  merged_df_segmentation)
summary(anova_result)


column_names_segmentation <- names(merged_df_segmentation)
column_names_segmentation <- column_names_segmentation[2:(length(column_names_segmentation) - 6 )]

right_regions <- column_names_segmentation[grepl("^Right(?!.*hypointensities)", column_names_segmentation,ignore.case = TRUE, perl = TRUE)]
left_regions <- column_names_segmentation[grepl("^Left(?!.*hypointensities)", column_names_segmentation,ignore.case = TRUE, perl = TRUE)]

only_subcortical_regions <- c(right_regions, left_regions)

significant_ROIs_not_age <- list()
significant_ROIs_not_etiv <- list()
significant_ROIs_not_etiv_not_age <- list()

# Create an empty data frame to store significant regions and p-values
p_values_segmentation_df <- data.frame(region = character(),
                          mean_control = numeric(), sd_control = numeric(), 
                          mean_covid = numeric(), sd_covid = numeric(), p_value = numeric(), ancova=character() , stringsAsFactors = FALSE)


# perform ancova for each dependent variable - Segmentation including edad and etiv
for(col_name in only_subcortical_regions) {
  formula <- as.formula(paste(col_name, "~ Grupo + Edad + Estimated.Total.Intracranial.Volume"))
  
  # perform ancova with group_variable as a factor and age as covariate
  ancova_result <- lm(formula, data = merged_df_segmentation)
  model <- Anova(ancova_result, type="III")
  
  # Extract p-value from the ANOVA results
  p_value_etiv <- model$"Pr(>F)"[4]
  p_value_age <- model$"Pr(>F)"[3]
  p_value_group <- model$"Pr(>F)"[2]
  
  
  if (p_value_etiv > 0.05){
    
    if (p_value_age > 0.05){
      significant_ROIs_not_etiv_not_age <- append(significant_ROIs_not_etiv_not_age, col_name)
    } else {
      significant_ROIs_not_etiv <- append(significant_ROIs_not_etiv, col_name)
    }
  } 
  
  if (p_value_age > 0.05) {
    if (p_value_etiv < 0.05){
      significant_ROIs_not_age <- append(significant_ROIs_not_age, col_name)
    }
    
  }
  
  if (p_value_age < 0.05 & p_value_etiv < 0.05){
    
    # Calculate mean and standard deviation for both groups
    mean_covid <- mean(merged_df_segmentation[merged_df_segmentation$Grupo == "COVID", col_name], na.rm = TRUE)
    mean_control <- mean(merged_df_segmentation[merged_df_segmentation$Grupo == "CONTROL", col_name], na.rm = TRUE)
    sd_covid <- sd(merged_df_segmentation[merged_df_segmentation$Grupo == "COVID", col_name], na.rm = TRUE)
    sd_control <- sd(merged_df_segmentation[merged_df_segmentation$Grupo == "CONTROL", col_name], na.rm = TRUE)
    
    # Append to a dataframe
    p_values_segmentation_df <- bind_rows(p_values_segmentation_df, data.frame(region = col_name, mean_covid = mean_covid, sd_covid = sd_covid, 
                                                                               mean_control = mean_control, sd_control = sd_control, 
                                                                               p_value = p_value_group,  ancova="Age, eTIV") )
    
  }
    
}


# perform ancova for each dependent variable - Segmentation including edad
for(col_name in significant_ROIs_not_etiv) {
  formula <- as.formula(paste(col_name, "~ Grupo + Edad"))
  
  # perform ancova with group_variable as a factor and age as covariate
  ancova_result <- lm(formula, data = merged_df_segmentation)
  model <- Anova(ancova_result, type="III")
  
  p_value_group <- model$"Pr(>F)"[2]
  
  # Calculate mean and standard deviation for both groups
  mean_covid <- mean(merged_df_segmentation[merged_df_segmentation$Grupo == "COVID", col_name], na.rm = TRUE)
  mean_control <- mean(merged_df_segmentation[merged_df_segmentation$Grupo == "CONTROL", col_name], na.rm = TRUE)
  sd_covid <- sd(merged_df_segmentation[merged_df_segmentation$Grupo == "COVID", col_name], na.rm = TRUE)
  sd_control <- sd(merged_df_segmentation[merged_df_segmentation$Grupo == "CONTROL", col_name], na.rm = TRUE)
  
  # Append to a dataframe
  p_values_segmentation_df <- bind_rows(p_values_segmentation_df, data.frame(region = col_name, mean_covid = mean_covid, sd_covid = sd_covid, 
                                                                             mean_control = mean_control, sd_control = sd_control, 
                                                                             p_value = p_value_group,  ancova="Age") )
  print(model)
  if (p_value_group < 0.05){
    print(model)
  }
  
}

# perform ancova for each dependent variable - Segmentation including etiv
for(col_name in significant_ROIs_not_age) {
  formula <- as.formula(paste(col_name, "~ Grupo + Estimated.Total.Intracranial.Volume"))
  
  # perform ancova with group_variable as a factor and age as covariate
  ancova_result <- lm(formula, data = merged_df_segmentation)
  model <- Anova(ancova_result, type="III")
  
  p_value_group <- model$"Pr(>F)"[2]
  
  # Calculate mean and standard deviation for both groups
  mean_covid <- mean(merged_df_segmentation[merged_df_segmentation$Grupo == "COVID", col_name], na.rm = TRUE)
  mean_control <- mean(merged_df_segmentation[merged_df_segmentation$Grupo == "CONTROL", col_name], na.rm = TRUE)
  sd_covid <- sd(merged_df_segmentation[merged_df_segmentation$Grupo == "COVID", col_name], na.rm = TRUE)
  sd_control <- sd(merged_df_segmentation[merged_df_segmentation$Grupo == "CONTROL", col_name], na.rm = TRUE)
  
  # Append to a dataframe
  p_values_segmentation_df <- bind_rows(p_values_segmentation_df, data.frame(region = col_name, mean_covid = mean_covid, sd_covid = sd_covid, 
                                                                             mean_control = mean_control, sd_control = sd_control, 
                                                                             p_value = p_value_group,  ancova="eTIV") )
  
  if (p_value_group < 0.05){
    print(model)
  }
  
}

# Right and Left Cerebellum Cortex 
# Right
ancova_result <- lm(Right.Cerebellum.Cortex ~ Grupo + Estimated.Total.Intracranial.Volume , data = merged_df_segmentation)
model <- Anova(ancova_result, type="III")

# Assuming ancova.example is your ANCOVA model
adj.means.ex <- effect("Grupo", ancova_result)

# Create a bar plot of the adjusted means
plot(adj.means.ex, type = "bar")


results_segmentation <- paste(path_results_data_analysis, "results_segmentation.csv", sep="")

write.csv(p_values_segmentation_df, results_segmentation, row.names = FALSE)


# Parcellation
# Exclude Subject
filtered_df_parcellation <- subset(parcellation_csv, !(subject %in% values_to_exclude))

# Etiv from brainvol
merged_df_parcellation <- cbind(filtered_df_parcellation, filtered_brainvol$Estimated.Total.Intracranial.Volume)
colnames(merged_df_parcellation)[ncol(merged_df_parcellation)] <-"Estimated.Total.Intracranial.Volume"

# Edad, Grupo, TMT from Cuestionarios_excel
merged_df_parcellation  <- merge(merged_df_parcellation , cuestionarios_excel[, c("Edad", "Grupo", "ID", "TMTA_valor_de_referencia")], by.x = "subject", by.y = "ID", all.x = TRUE)

# Mean both hemispheres
mean_df_parcellation <- merged_df_parcellation %>%
  group_by(subject) %>%
  summarise_all(mean)


mean_df_parcellation <- mean_df_parcellation[1:75]
mean_df_parcellation <- merge(mean_df_parcellation, cuestionarios_excel[, c("ID", "Edad", "Grupo", "TMTA_valor_de_referencia")], by.x = "subject", by.y = "ID", all.x = TRUE)

# Etiv from brainvol
mean_df_parcellation <- cbind(mean_df_parcellation, filtered_brainvol$Estimated.Total.Intracranial.Volume)
colnames(mean_df_parcellation)[ncol(mean_df_parcellation)] <-"Estimated.Total.Intracranial.Volume"

column_names_parcellation <- names(mean_df_parcellation)
column_names_parcellation<- column_names_parcellation[2:(length(column_names_parcellation) - 4 )]





# Perform ancova for each dependent ROI - parcellation
significant_ROIs_not_age <- list()

# Create an empty data frame to store significant regions and p-values
p_values_df <- data.frame(region = character(),
                          mean_control = numeric(), sd_control = numeric(), 
                          mean_covid = numeric(), sd_covid = numeric(), p_value = numeric(), ancova=character() , stringsAsFactors = FALSE)

for(col_name in column_names_parcellation) {
  formula <- as.formula(paste(col_name, "~ Grupo + Edad + Estimated.Total.Intracranial.Volume"))
  
  # perform ancova with group_variable as a factor and age as covariate
  ancova_result <- lm(formula, data = mean_df_parcellation)
  model <- Anova(ancova_result, type="III")
  # Extract p-value from the ANOVA results
  
  p_value_age <- model$"Pr(>F)"[3]
  p_value_group <- model$"Pr(>F)"[2]
  
  if (p_value_age < 0.05){
    
    # Calculate mean and standard deviation for both groups
    mean_covid <- mean(mean_df_parcellation[mean_df_parcellation$Grupo == "COVID", col_name], na.rm = TRUE)
    mean_control <- mean(mean_df_parcellation[mean_df_parcellation$Grupo == "CONTROL", col_name], na.rm = TRUE)
    sd_covid <- sd(mean_df_parcellation[mean_df_parcellation$Grupo == "COVID", col_name], na.rm = TRUE)
    sd_control <- sd(mean_df_parcellation[mean_df_parcellation$Grupo == "CONTROL", col_name], na.rm = TRUE)
    
    # Append to a dataframe
    p_values_df <- bind_rows(p_values_df, data.frame(region = col_name, mean_covid = mean_covid, sd_covid = sd_covid, 
                                                     mean_control = mean_control, sd_control = sd_control, 
                                                     p_value = p_value_group,  ancova="Age, eTIV") )
    
    } else {
    significant_ROIs_not_age <- append(significant_ROIs_not_age, col_name)
    
    
  }
}

for(col_name in significant_ROIs_not_age) {
  formula <- as.formula(paste(col_name, "~ Grupo + Estimated.Total.Intracranial.Volume"))
  
  # perform ancova with group_variable as a factor and age as covariate
  ancova_result <- lm(formula, data = mean_df_parcellation)
  model <- Anova(ancova_result, type="III")
  
  # Extract p-value from the ANOVA results
  p_value_group <- model$"Pr(>F)"[2]
  
  # Calculate mean and standard deviation for both groups
  mean_covid <- mean(mean_df_parcellation[mean_df_parcellation$Grupo == "COVID", col_name], na.rm = TRUE)
  mean_control <- mean(mean_df_parcellation[mean_df_parcellation$Grupo == "CONTROL", col_name], na.rm = TRUE)
  sd_covid <- sd(mean_df_parcellation[mean_df_parcellation$Grupo == "COVID", col_name], na.rm = TRUE)
  sd_control <- sd(mean_df_parcellation[mean_df_parcellation$Grupo == "CONTROL", col_name], na.rm = TRUE)
  
  # Append to a dataframe
  p_values_df <- bind_rows(p_values_df, data.frame(region = col_name, mean_covid = mean_covid, sd_covid = sd_covid, 
                                                   mean_control = mean_control, sd_control = sd_control, 
                                                   p_value = p_value_group,  ancova="eTIV") )
  
}


results_parcellation <- paste(path_results_data_analysis, "parcellation_ancova.csv", sep="")

write.csv(p_values_df, results_parcellation, row.names = FALSE)



# Parcellation DKT
filtered_df_parcellation_DKT <- subset(parcellation_DKT_csv, !(subject %in% values_to_exclude))

# Etiv from brainvol
merge_df_parcellation_DKT <- cbind(filtered_df_parcellation_DKT, filtered_brainvol$Estimated.Total.Intracranial.Volume)
colnames(merge_df_parcellation_DKT)[ncol(merge_df_parcellation_DKT)] <-"Estimated.Total.Intracranial.Volume"

# Edad, Grupo from cuestionarios_excel
merge_df_parcellation_DKT  <- merge(merge_df_parcellation_DKT , cuestionarios_excel[, c("Edad", "Grupo", "ID", "TMTA_valor_de_referencia")], by.x = "subject", by.y = "ID", all.x = TRUE)
# Mean both hemispheres
mean_df_parcellation_DKT <- filtered_df_parcellation_DKT %>%
  group_by(subject) %>%
  summarise_all(mean)

mean_df_parcellation_DKT <- mean_df_parcellation_DKT[1:35]

# Etiv from brainvol
mean_df_parcellation_DKT <- cbind(mean_df_parcellation_DKT, filtered_brainvol$Estimated.Total.Intracranial.Volume)
colnames(mean_df_parcellation_DKT)[ncol(mean_df_parcellation_DKT)] <-"Estimated.Total.Intracranial.Volume"

# Edad, Grupo from cuestionarios_excel
mean_df_parcellation_DKT  <- merge(mean_df_parcellation_DKT , cuestionarios_excel[, c("Edad", "Grupo", "ID", "TMTA_valor_de_referencia")], by.x = "subject", by.y = "ID", all.x = TRUE)

column_names_parcellation_DKT  <- names(mean_df_parcellation_DKT)
column_names_parcellation_DKT <- column_names_parcellation_DKT [2:(length(column_names_parcellation_DKT) - 4)]




# Perform ancova for each dependent ROI - parcellation - DKT 
significant_ROIs_not_age <- list()

# Create an empty data frame to store significant regions and p-values
p_values_df_parcellation_dkt <- data.frame(region = character(),
                          mean_control = numeric(), sd_control = numeric(), 
                          mean_covid = numeric(), sd_covid = numeric(), p_value = numeric(), ancova=character() , stringsAsFactors = FALSE)

for(col_name in column_names_parcellation_DKT) {
  formula <- as.formula(paste(col_name, "~ Grupo + Edad + Estimated.Total.Intracranial.Volume"))
  
  # perform ancova with group_variable as a factor and age as covariate
  ancova_result <- lm(formula, data = mean_df_parcellation_DKT)
  model <- Anova(ancova_result, type="III")
  # Extract p-value from the ANOVA results
  
  p_value_age <- model$"Pr(>F)"[3]
  p_value_group <- model$"Pr(>F)"[2]
  
  if (p_value_age < 0.05){
    
    # Calculate mean and standard deviation for both groups
    mean_covid <- mean(mean_df_parcellation_DKT[mean_df_parcellation_DKT$Grupo == "COVID", col_name], na.rm = TRUE)
    mean_control <- mean(mean_df_parcellation_DKT[mean_df_parcellation_DKT$Grupo == "CONTROL", col_name], na.rm = TRUE)
    sd_covid <- sd(mean_df_parcellation_DKT[mean_df_parcellation_DKT$Grupo == "COVID", col_name], na.rm = TRUE)
    sd_control <- sd(mean_df_parcellation_DKT[mean_df_parcellation_DKT$Grupo == "CONTROL", col_name], na.rm = TRUE)
    
    # Append to a dataframe
    p_values_df_parcellation_dkt <- bind_rows(p_values_df_parcellation_dkt, data.frame(region = col_name, mean_covid = mean_covid, sd_covid = sd_covid, 
                                                     mean_control = mean_control, sd_control = sd_control, 
                                                     p_value = p_value_group,  ancova="Age, eTIV") )
    
  } else {
    significant_ROIs_not_age <- append(significant_ROIs_not_age, col_name)
    
    
  }
}

for(col_name in significant_ROIs_not_age) {
  formula <- as.formula(paste(col_name, "~ Grupo + Estimated.Total.Intracranial.Volume"))
  
  # perform ancova with group_variable as a factor and age as covariate
  ancova_result <- lm(formula, data = mean_df_parcellation_DKT)
  model <- Anova(ancova_result, type="III")
  
  # Extract p-value from the ANOVA results
  p_value_group <- model$"Pr(>F)"[2]
  
  # Calculate mean and standard deviation for both groups
  mean_covid <- mean(mean_df_parcellation_DKT[mean_df_parcellation_DKT$Grupo == "COVID", col_name], na.rm = TRUE)
  mean_control <- mean(mean_df_parcellation_DKT[mean_df_parcellation_DKT$Grupo == "CONTROL", col_name], na.rm = TRUE)
  sd_covid <- sd(mean_df_parcellation_DKT[mean_df_parcellation_DKT$Grupo == "COVID", col_name], na.rm = TRUE)
  sd_control <- sd(mean_df_parcellation_DKT[mean_df_parcellation_DKT$Grupo == "CONTROL", col_name], na.rm = TRUE)
  
  # Append to a dataframe
  p_values_df_parcellation_dkt <- bind_rows(p_values_df_parcellation_dkt, data.frame(region = col_name, mean_covid = mean_covid, sd_covid = sd_covid, 
                                                   mean_control = mean_control, sd_control = sd_control, 
                                                   p_value = p_value_group,  ancova="eTIV") )
  
}


results_parcellation_dkt <- paste(path_results_data_analysis, "parcellation_ancova_DKT.csv", sep="")

write.csv(p_values_df_parcellation_dkt, results_parcellation_dkt, row.names = FALSE)



# Cortical Thickness
filtered_df_parcellation_thick <- subset(parcellation_thick_csv, !(subject %in% values_to_exclude))

# Mean both hemispheres
mean_df_parcellation_thick <- filtered_df_parcellation_thick  %>%
  group_by(subject) %>%
  summarise_all(mean)

# Edad, Grupo from cuestionarios_excel
mean_df_parcellation_thick  <- merge(mean_df_parcellation_thick , cuestionarios_excel[, c("Edad", "Grupo", "ID")], by.x = "subject", by.y = "ID", all.x = TRUE)

column_names_parcellation_thick  <- names(mean_df_parcellation_thick)
column_names_parcellation_thick <- column_names_parcellation_thick[2:(length(column_names_parcellation_thick) - 5)]

# Perform t-test and save in a csv

# Create an empty data frame to store significant regions and p-values
p_values_df_parcellation_thickness<- data.frame(region = character(),
                                                     mean_control = numeric(), sd_control = numeric(), 
                                                     mean_covid = numeric(), sd_covid = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

significant_ROIs_cort_thick <- list()

for(col_name in column_names_parcellation_thick) {
 
  values_covid <- mean_df_parcellation_thick[mean_df_parcellation_thick$Grupo == "COVID", col_name]
  values_control <- mean_df_parcellation_thick[mean_df_parcellation_thick$Grupo == "CONTROL", col_name]
  
  # Calculate mean and standard deviation for both groups
  mean_covid <- mean(values_covid , na.rm = TRUE)
  mean_control <- mean(values_control, na.rm = TRUE)
  sd_covid <- sd(values_covid, na.rm = TRUE)
  sd_control <- sd(values_control, na.rm = TRUE)
  
  # Test T 
  result_t_thick <- t.test(values_covid, values_control)
  p_value <- result_t_thick$p.value
  
  p_values_df_parcellation_thickness <- bind_rows(p_values_df_parcellation_thickness, data.frame(region = col_name, mean_covid = mean_covid, sd_covid = sd_covid, 
                                                                                                         mean_control = mean_control, sd_control = sd_control, 
                                                                                                         p_value = p_value) )
  
  if (p_value < 0.05) {
    significant_ROIs_cort_thick <- append(significant_ROIs_cort_thick, col_name)
  }
}

results_parcellation_cortical_thickness <- paste(path_results_data_analysis, "parcellation_cort_thickness.csv", sep="")
write.csv(p_values_df_parcellation_thickness, results_parcellation_cortical_thickness, row.names = FALSE)

# Boxplots
boxplots <- list()


for (roi in significant_ROIs_cort_thick) {
  
  values_covid <- mean_df_parcellation_thick[mean_df_parcellation_thick$Grupo == "COVID", roi]
  values_control <- mean_df_parcellation_thick[mean_df_parcellation_thick$Grupo == "CONTROL", roi]
  
  
  # Combine data into a single dataframe
  df <- data.frame(Group = c(rep("COVID", length(values_covid)), rep("Control", length(values_control))),
                   Thickness = c(values_covid, values_control))
  
  # Create boxplot with jittered points
  boxplot <- ggplot(df, aes(x = Group, y = Thickness, fill = Group)) +
    geom_boxplot(alpha = 0.5) +  # Add transparency to boxplot
    geom_jitter(color = "black", size = 2, alpha = 0.3, width = 0.2) +  # Add jittered points
    labs(title = roi,
         x = "Group",
         y = "Cortical Thickness (mm)") +
    theme_minimal() +
    guides(fill = FALSE)  # Remove legend
  
  boxplots[[roi]] <- boxplot
  
}

grid_arrange1 <- do.call(gridExtra::grid.arrange, c(boxplots[1:8], ncol = 4))
grid_arrange2 <- do.call(gridExtra::grid.arrange, c(boxplots[9:12], ncol = 4))
ggsave(paste(path_plots_data_analysis, "boxplots_cort_thickness1.png", sep=""),grid_arrange1, width = 8, height = 6, dpi = 300)
ggsave(paste(path_plots_data_analysis, "boxplots_cort_thickness2.png", sep=""),grid_arrange2, width = 8, height = 6, dpi = 300)


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

# Check Normality 
# Function to create QQ plot
create_qqplot <- function(data, col_name) {
  ggplot(data, aes(sample = !!sym(col_name))) +
    geom_qq() +
    geom_qq_line() +
    labs(title = paste("QQ Plot for", col_name),
         x = "Theoretical Quantiles",
         y = "Sample Quantiles") +
    theme_minimal()
}

# List to store QQ plots
qq_plots1 <- list()
qq_plots2 <- list()


# Create QQ plot for each column
for(col_name in column_names_parcellation_thick_dkt[1:15]) {
  qq_plots1[[col_name]] <- create_qqplot(mean_df_parcellation_thick_dkt, col_name)
}

# Create QQ plot for each column
for(col_name in column_names_parcellation_thick_dkt[16:length(column_names_parcellation_thick_dkt)]) {
  qq_plots2[[col_name]] <- create_qqplot(mean_df_parcellation_thick_dkt, col_name)
}

# Arrange plots in a grid
grid_arrange1 = do.call(gridExtra::grid.arrange, c(qq_plots1, ncol = 5))
grid_arrange2 = do.call(gridExtra::grid.arrange, c(qq_plots2, ncol = 5))

ggsave(paste(path_plots_data_analysis, "qq_plot1.png", sep=""), grid_arrange1, width = 16, height = 10, dpi = 300)
ggsave(paste(path_plots_data_analysis, "qq_plot2.png", sep=""), grid_arrange2, width = 16, height = 10, dpi = 300)

# Perform t-test and save in a csv

# Create an empty data frame to store significant regions and p-values
p_values_df_parcellation_thickness_dkt <- data.frame(region = character(),
                                           mean_control = numeric(), sd_control = numeric(), 
                                           mean_covid = numeric(), sd_covid = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

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
  
  p_values_df_parcellation_thickness_dkt <- bind_rows(p_values_df_parcellation_thickness_dkt, data.frame(region = col_name, mean_covid = mean_covid, sd_covid = sd_covid, 
                                                                                     mean_control = mean_control, sd_control = sd_control, 
                                                                                     p_value = p_value) )
  
  if (p_value < 0.05) {
    significant_ROIs_cort_thick_dkt <- append(significant_ROIs_cort_thick_dkt, col_name)
  }
}

results_parcellation_cortical_thickness_dkt <- paste(path_results_data_analysis, "parcellation_cort_thickness_dkt.csv", sep="")
write.csv(p_values_df_parcellation_thickness_dkt, results_parcellation_cortical_thickness_dkt, row.names = FALSE)

# Boxplots
boxplots <- list()


for (roi in significant_ROIs_cort_thick_dkt) {

  values_covid <- mean_df_parcellation_thick_dkt[mean_df_parcellation_thick_dkt$Grupo == "COVID", roi]
  values_control <- mean_df_parcellation_thick_dkt[mean_df_parcellation_thick_dkt$Grupo == "CONTROL", roi]
  

  # Combine data into a single dataframe
  df <- data.frame(Group = c(rep("COVID", length(values_covid)), rep("Control", length(values_control))),
                   Thickness = c(values_covid, values_control))
  S
  # Create boxplot with jittered points
  boxplot <- ggplot(df, aes(x = Group, y = Thickness, fill = Group)) +
    geom_boxplot(alpha = 0.5) +  # Add transparency to boxplot
    geom_jitter(color = "black", size = 2, alpha = 0.3, width = 0.2) +  # Add jittered points
    labs(title = roi,
         x = "Group",
         y = "Cortical Thickness (mm)") +
    theme_minimal() +
    guides(fill = FALSE)  # Remove legend
  
  boxplots[[roi]] <- boxplot
  
}


grid_arrange <- do.call(gridExtra::grid.arrange, c(boxplots, ncol = 4))
ggsave(paste(path_plots_data_analysis, "boxplots_cort_thickness_dkt.png", sep=""),grid_arrange, width = 8, height = 6, dpi = 300)


## Correlations between cortical thickness and questionaries and test
cognitive_test_id <- append(cognitive_test, "ID")
standardized_questionnaire_id <- append(standardized_questionnaire, "ID")

mean_df_parcellation_thick  <- merge(mean_df_parcellation_thick , cuestionarios_excel[, cognitive_test_id], by.x = "subject", by.y = "ID", all.x = TRUE)
mean_df_parcellation_thick  <- merge(mean_df_parcellation_thick , cuestionarios_excel[, standardized_questionnaire_id], by.x = "subject", by.y = "ID", all.x = TRUE)

mean_df_parcellation_thick_dkt  <- merge(mean_df_parcellation_thick_dkt , cuestionarios_excel[, cognitive_test_id], by.x = "subject", by.y = "ID", all.x = TRUE)
mean_df_parcellation_thick_dkt  <- merge(mean_df_parcellation_thick_dkt , cuestionarios_excel[, standardized_questionnaire_id], by.x = "subject", by.y = "ID", all.x = TRUE)



for (roi in significant_ROIs_cort_thick){
  for (cogn_test in cognitive_test){
    formula <- as.formula(paste(roi, "~", cogn_test))
    
    
    # perform ancova with group_variable as a factor and age as covariate
    ancova_result <- lm(formula, data = mean_df_parcellation_thick)
    model <- Anova(ancova_result, type="III")
    
    p_value <- model$"Pr(>F)"[2]
    if (p_value < 0.05){
      print(model)
    }  
  }
  
  for (quest in standardized_questionnaire){
    formula <- as.formula(paste(roi, "~", quest))
    
    # perform ancova with group_variable as a factor and age as covariate
    ancova_result <- lm(formula, data = mean_df_parcellation_thick)
    model <- Anova(ancova_result, type="III")
    
    p_value <- model$"Pr(>F)"[2]
    if (p_value < 0.05){
      print(model)
    }
  }
}



for (roi in significant_ROIs_cort_thick_dkt){
  for (cogn_test in cognitive_test){
    formula <- as.formula(paste(roi, "~", cogn_test))
    
    # perform ancova with group_variable as a factor and age as covariate
    ancova_result <- lm(formula, data = mean_df_parcellation_thick_dkt)
    model <- Anova(ancova_result, type="III")
    
    p_value <- model$"Pr(>F)"[2]
    if (p_value < 0.05){
      print(model)
    }
  }
  
  for (quest in standardized_questionnaire){
    formula <- as.formula(paste(roi, "~", quest))
    
    # perform ancova with group_variable as a factor and age as covariate
    ancova_result <- lm(formula, data = mean_df_parcellation_thick_dkt)
    model <- Anova(ancova_result, type="III")
    
    p_value <- model$"Pr(>F)"[2]
    if (p_value < 0.05){
      print(model)
    }
  }
}

## Exploring correlations between cerebellum and questionaries and test

merged_df_segmentation_tests  <- merge(merged_df_segmentation , cuestionarios_excel[, cognitive_test_id], by.x = "subject", by.y = "ID", all.x = TRUE)
merged_df_segmentation_tests  <- merge(merged_df_segmentation_tests , cuestionarios_excel[, standardized_questionnaire_id], by.x = "subject", by.y = "ID", all.x = TRUE)

for (cogn_test in cognitive_test){
  formula <- as.formula(paste(  "Right.Cerebellum.Cortex ~", cogn_test))
  
  # perform ancova with group_variable as a factor and age as covariate
  ancova_result <- lm(formula, data = merged_df_segmentation_tests)
  model <- Anova(ancova_result, type="III")
  
  p_value <- model$"Pr(>F)"[2]
  print(summary(ancova_result))
  if (p_value < 0.05){
    print(model)
  }
}

for (quest in standardized_questionnaire){
  formula <- as.formula(paste("Right.Cerebellum.Cortex ~", quest))
  
  # perform ancova with group_variable as a factor and age as covariate
  ancova_result <- lm(formula, data = merged_df_segmentation_tests)
  model <- Anova(ancova_result, type="III")
  
  p_value <- model$"Pr(>F)"[2]
  print(model)
}


for (cogn_test in cognitive_test){
  formula <- as.formula(paste(  "Left.Cerebellum.Cortex ~", cogn_test))
  
  # perform ancova with group_variable as a factor and age as covariate
  ancova_result <- lm(formula, data = merged_df_segmentation_tests)
  model <- Anova(ancova_result, type="III")
  
  p_value <- model$"Pr(>F)"[2]
  if (p_value < 0.05){
    print(model)
  }
}

for (quest in standardized_questionnaire){
  formula <- as.formula(paste("Left.Cerebellum.Cortex ~", quest))
  
  # perform ancova with group_variable as a factor and age as covariate
  ancova_result <- lm(formula, data = merged_df_segmentation_tests)
  model <- Anova(ancova_result, type="III")
  
  p_value <- model$"Pr(>F)"[2]
  
  if (p_value < 0.05){
    print(model)
  }
  
}
