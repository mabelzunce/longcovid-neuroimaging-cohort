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
path_ismrm <- "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/Nifti/Final_data_csv/DataISMRM/"

path_asl_analysis <- "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/ASLProcessedCOVID/derivatives/ExploreASL/Population/Stats/"

path_results_data_analysis <- paste(path_data_analysis, "Results/", sep="")
path_plots_data_analysis <-  paste(path_data_analysis, "Plots/", sep="")

cov_perfussion_data_without_PVC <- "CoV_qCBF_StandardSpace_TotalGM_n=203_15-Jan-2025_PVC0.tsv"
cov_perfussion_data_with_PVC <- "CoV_qCBF_StandardSpace_TotalGM_n=203_15-Jan-2025_PVC2.tsv"

path_scov_general <-paste(path_asl_analysis, cov_perfussion_data_without_PVC , sep="")

raw_gray_matter_data_ismrm <-"mean_qCBF_StandardSpace_TotalGM_2024_PVC2.csv"
path_gray_matter_perfusion_ismrm <-paste(path_ismrm, raw_gray_matter_data_ismrm , sep="")

# ISMRM subjects
gray_matter_csv_ismrm <- read.csv(path_gray_matter_perfusion_ismrm)
subjects_ismrm <- gray_matter_csv_ismrm[["ID"]]
subjects_ismrm <- paste0("sub-", subjects_ismrm, "_1")

# sCOV subjects
scov_data <- read.delim(path_scov_general, sep="\t")
scov_only_ismrm <- scov_data[scov_data$participant_id %in% subjects_ismrm, ]


# All Data

# Renombrar la columna
colnames(scov_only_ismrm)[colnames(scov_only_ismrm) == "participant_id"] <- "ID"

# Modificar todos los valores de la nueva columna
# Eliminar prefijo y sufijo
scov_only_ismrm$ID <- gsub("^sub-|_1$", "", scov_only_ismrm$ID)

path_scov_ismrm_data <- paste(path_ismrm, "CoV_qCBF_StandardSpace_TotalGM_ISRM.csv",sep="")
write.csv(scov_only_ismrm, file = path_scov_ismrm_data, row.names = FALSE)

                                                                                                                                                                                                  








