# Define a function for ANOVA
perform_anova <- function(formula, data) {
  ancova_result <- lm(formula, data = data)
  model <- Anova(ancova_result, type = "III")
  return(model)
}

# Define a function for extracting mean and standard deviation
get_mean_sd <- function(data, group, col_name) {
  mean_val <- mean(data[data$Grupo == group, col_name], na.rm = TRUE)
  sd_val <- sd(data[data$Grupo == group, col_name], na.rm = TRUE)
  return(list(mean = mean_val, sd = sd_val))
}

# Define a function for extracting p-value from the ANOVA results
get_p_value <- function(model) {
  return(model$"Pr(>F)"[2])
}

# Define a function to handle appending results to a dataframe
append_to_dataframe <- function(df, col_name, mean_covid, sd_covid, mean_control, sd_control, p_value, f_value, ancova) {
  df <- bind_rows(df, data.frame(region = col_name, mean_covid = mean_covid, sd_covid = sd_covid, 
                                 mean_control = mean_control, sd_control = sd_control, 
                                 p_value = p_value, f_value = f_value, ancova = ancova))
  return(df)
}

# Perform ANOVA or ANCOVA for each column
perform_analysis <- function(column_names, data, df, hemisphere) {
  for (col_name in column_names) {
    formula_full <- as.formula(paste(col_name, "~ Grupo + Edad + Estimated.Total.Intracranial.Volume"))
    formula_etiv <- as.formula(paste(col_name, "~ Grupo + Estimated.Total.Intracranial.Volume"))
    formula_age <- as.formula(paste(col_name, "~ Grupo + Edad"))
    formula_group <- as.formula(paste(col_name, "~ Grupo"))
    
    model_full <- perform_anova(formula_full, data)
    p_value_group_full <- model_full$"Pr(>F)"[2]
    p_value_age_full <- model_full$"Pr(>F)"[3]
    p_value_etiv_full <- model_full$"Pr(>F)"[4]
    
    # Extract p-values and F-values
    f_value_group_full <- model_full$"F value"[2]
    
    
    mean_covid <- get_mean_sd(data, "COVID", col_name)$mean
    sd_covid <- get_mean_sd(data, "COVID", col_name)$sd
    mean_control <- get_mean_sd(data, "CONTROL", col_name)$mean
    sd_control <- get_mean_sd(data, "CONTROL", col_name)$sd
    
    
    if (p_value_age_full < 0.05 & p_value_etiv_full < 0.05) {
      # Both age and eTIV are significant, leave the model as is
      df <- append_to_dataframe(df, col_name, mean_covid, sd_covid, mean_control, sd_control, p_value_group_full, f_value_group_full, "Age, eTIV")
      if (p_value_group_full < 0.05) {
        print(paste("Significant effect for", col_name, "in", hemisphere, "hemisphere with ANCOVA on Age and eTIV"))
      }
    } else if (p_value_age_full < 0.05) {
      # Only age is significant, perform ANCOVA correcting only for age
      model_age <- perform_anova(formula_age, data)
      p_value_group_age <- model_age$"Pr(>F)"[2]
      f_value_group_age <- model_age$"F value"[2]
      df <- append_to_dataframe(df, col_name, mean_covid, sd_covid, mean_control, sd_control, p_value_group_age, f_value_group_age, "Age")
      if (p_value_group_age < 0.05) {
        print(paste("Significant effect for", col_name, "in", hemisphere,"hemisphere with ANCOVA on Age"))
        print(model_age)
      }
    } else if (p_value_etiv_full < 0.05) {
      # Only eTIV is significant, perform ANCOVA correcting only for eTIV
      model_etiv <- perform_anova(formula_etiv, data)
      p_value_group_etiv <- model_etiv$"Pr(>F)"[2]
      f_value_group_etiv <- model_etiv$"F value"[2]
      df <- append_to_dataframe(df, col_name, mean_covid, sd_covid, mean_control, sd_control, p_value_group_etiv, f_value_group_etiv, "eTIV")
      if (p_value_group_etiv < 0.05) {
        print(paste("Significant effect for", col_name, "in", hemisphere,"hemisphere with ANCOVA on eTIV"))
        print(model_etiv)
      }
    } else {
      # Neither age nor eTIV is significant, perform ANCOVA correcting only for group
      model_group <- perform_anova(formula_group, data)
      p_value_group_only <- model_group$"Pr(>F)"[2]
      f_value_group_only <- model_group$"F value"[2]
      df <- append_to_dataframe(df, col_name, mean_covid, sd_covid, mean_control, sd_control, p_value_group_only, f_value_group_only, ".")
      if (p_value_group_only < 0.05) {
        print(paste("Significant effect for", col_name, "in", hemisphere,"hemisphere with ANCOVA"))
        print(model_group)
        
      }
    }
  }
  return(list(df = df))
}


cluster_data <- function(data, 
                         free_group = NULL,  # Si se especifica, para los sujetos cuyo Grupo coincida se asigna ese valor en Cluster
                         proteins = c("BDNF", "M6a", "NFL", "TotalGM_CBF"),
                         cognitive_vars = NULL, 
                         centers = 2, 
                         nstart = 25) {
  # Usar todo el dataset sin filtrar por grupo
  data_subset <- data
  
  # Preparar datos para el clustering: quitamos 'Grupo' pero mantenemos 'ID'
  clustering_data <- data_subset %>% dplyr::select(-Grupo)
  
  # Escalar los datos (excluyendo la columna ID)
  clustering_data_scaled <- clustering_data %>%
    dplyr::select(-ID) %>%
    scale() %>%
    as.data.frame()
  
  # Eliminar variables de proteínas y TotalGM_CBF del clustering
  clustering_data_no_protein <- clustering_data_scaled %>%
    dplyr::select(-all_of(proteins)) %>%
    as.data.frame()
  
  # Si se está clusterizando todo el dataset o se clusterizan controles, 
  # se omiten las variables de síntomas cognitivos (no disponibles para controles)
  if (is.null(free_group) || free_group == "CONTROL") {
    clustering_data_no_protein <- clustering_data_no_protein %>%
      dplyr::select(-all_of(cognitive_vars))
  }
  
  # Identificar y eliminar filas con NA en los datos de clustering
  na_positions <- is.na(clustering_data_no_protein)
  clustering_data_clean <- clustering_data_no_protein[!rowSums(na_positions), ]
  
  # Ejecutar k-means en los datos limpios
  set.seed(123)  # Para reproducibilidad
  km <- kmeans(clustering_data_clean, centers = centers, nstart = nstart)
  
  # Identificar los índices de las filas utilizadas en el clustering (sin NA)
  clean_indices <- which(rowSums(is.na(clustering_data_no_protein)) == 0)
  clean_ids <- data_subset$ID[clean_indices]
  
  # Crear un data frame con la asignación de clusters
  cluster_df <- data.frame(ID = clean_ids, Cluster = as.factor(km$cluster))
  
  # Fusionar la información de clusters con la data original
  merged_data <- merge(data_subset, cluster_df, by = "ID", all.x = TRUE)
  
  # Si se especifica free_group, para los sujetos cuyo Grupo sea free_group, se asigna ese valor en Cluster
  if (!is.null(free_group)) {
    merged_data$Cluster <- ifelse(merged_data$Grupo == free_group,
                                  free_group,
                                  merged_data$Cluster)
  }
  
  return(merged_data)
}

cluster_data_by_group <- function(data,
                                  group_filter,    # "CONTROL" o "COVID"
                                  proteins = c("BDNF", "M6a", "NFL", "TotalGM_CBF"),
                                  value_to_exclude = NULL,
                                  centers = 2, 
                                  nstart = 25) {
  # Filtrar el subset del grupo deseado
  data_subset <- data %>% dplyr::filter(Grupo == group_filter)
  
  # Quitar la columna 'Grupo' pero mantener 'ID'
  clustering_data <- data_subset %>% dplyr::select(-Grupo)
  
  # Escalar (excluyendo ID)
  clustering_data_scaled <- clustering_data %>%
    dplyr::select(-ID) %>%
    scale() %>%
    as.data.frame()
  
  # Eliminar variables de proteínas (y TotalGM_CBF) y opcionalmente otra variable
  if (!is.null(value_to_exclude)) {
    clustering_data_no_protein <- clustering_data_scaled %>%
      dplyr::select(-all_of(c(proteins, value_to_exclude))) %>%
      as.data.frame()
  } else {
    clustering_data_no_protein <- clustering_data_scaled %>%
      dplyr::select(-all_of(proteins)) %>%
      as.data.frame()
  }
  
  # Quitar filas con NA
  na_positions <- is.na(clustering_data_no_protein)
  clustering_data_clean <- clustering_data_no_protein[!rowSums(na_positions), ]
  
  # Ejecutar k-means
  set.seed(123)
  km <- kmeans(clustering_data_clean, centers = centers, nstart = nstart)
  
  # Identificar IDs de las filas usadas en el clustering (sin NA)
  clean_indices <- which(rowSums(is.na(clustering_data_no_protein)) == 0)
  clean_ids <- data_subset$ID[clean_indices]
  
  # Data frame con asignación de clusters
  cluster_df <- data.frame(ID = clean_ids, Cluster = as.factor(km$cluster))
  
  # Fusionar la asignación de clusters con el subset original
  merged_subset <- merge(data_subset, cluster_df, by = "ID", all.x = TRUE)
  
  # Devolver una lista con los datos fusionados, el resultado de kmeans y los datos usados para clustering
  return(list(merged_data = merged_subset,
              km_result = km,
              clustering_data_clean = clustering_data_clean))
}

cluster_data_flexible_v2 <- function(data,
                                     group_filter = NULL,    # Si es NULL, usa todo
                                     proteins = c("BDNF", "M6a", "NFL"),
                                     values_to_exclude = c("TotalGM_CBF"),
                                     value_to_exclude = NULL,
                                     centers = 2, 
                                     nstart = 25) {
  
  # Si se pasa un filtro de grupo, aplicarlo
  if (!is.null(group_filter)) {
    data_subset <- data %>% dplyr::filter(Grupo == group_filter)
  } else {
    data_subset <- data
  }
  
  # Quitar la columna 'Grupo' pero mantener 'ID'
  clustering_data <- data_subset %>% dplyr::select(-Grupo)
  
  # Escalar (excluyendo ID)
  clustering_data_scaled <- clustering_data %>%
    dplyr::select(-ID) %>%
    scale() %>%
    as.data.frame()
  
  # Eliminar variables de proteínas (y opcionalmente otra variable)
  if (!is.null(value_to_exclude)) {
    clustering_data_no_protein <- clustering_data_scaled %>%
      dplyr::select(-all_of(c(proteins, value_to_exclude))) %>%
      as.data.frame()
  } else {
    clustering_data_no_protein <- clustering_data_scaled %>%
      dplyr::select(-all_of(proteins)) %>%
      as.data.frame()
  }
  
  # Quitar filas con NA
  na_positions <- is.na(clustering_data_no_protein)
  clustering_data_clean <- clustering_data_no_protein[!rowSums(na_positions), ]
  
  # Ejecutar k-means
  set.seed(123)
  km <- kmeans(clustering_data_clean, centers = centers, nstart = nstart)
  
  # Identificar IDs de las filas usadas en el clustering
  clean_indices <- which(rowSums(is.na(clustering_data_no_protein)) == 0)
  clean_ids <- data_subset$ID[clean_indices]
  
  # Data frame con asignación de clusters
  cluster_df <- data.frame(ID = clean_ids, Cluster = as.factor(km$cluster))
  
  # Fusionar la asignación de clusters con el subset original
  merged_subset <- merge(data_subset, cluster_df, by = "ID", all.x = TRUE)
  
  # Fusionar con TODOS los datos originales
  if (!is.null(group_filter)) {
    # Si filtraste un grupo
    data_full <- data
    data_full$Cluster <- NA
    
    data_full <- data_full %>%
      dplyr::left_join(cluster_df, by = "ID", suffix = c("", ".new")) %>%
      dplyr::mutate(Cluster = dplyr::case_when(
        !is.na(Cluster.new) ~ as.character(Cluster.new),    # Si clusterizó, poner número
        Grupo != group_filter ~ as.character(Grupo),        # Si es de otro grupo, poner nombre del grupo
        TRUE ~ NA_character_                                # Si es del grupo filtrado pero no clusterizó, NA
      )) %>%
      dplyr::select(-Cluster.new)
    
  } else {
    # Si no filtraste, simplemente usás merged_subset
    data_full <- merged_subset
  }
  
  # Devolver una lista completa
  return(list(merged_data = data_full,
              km_result = km,
              clustering_data_clean = clustering_data_clean))
}



crear_boxplot <- function(data, variable, title = NULL, ylim_sup = NULL, mostrar_ids = FALSE) {
  # Armar estética dinámica
  aes_base <- aes_string(x = "Cluster", y = variable, fill = "Cluster")
  
  # Base del gráfico
  p <- ggplot(data, aes_base) +
    geom_boxplot(size = 1.2, outlier.shape = NA) +
    geom_jitter(color = "black", width = 0.2, alpha = 0.5) +
    labs(
      title = ifelse(is.null(title), variable, title),
      x = "Grupo",
      y = variable
    ) +
    theme_bw()
  
  # Agregar etiquetas si se solicita
  if (mostrar_ids) {
    p <- p + geom_text(
      aes(label = ID),
      position = position_jitter(width = 0.2),
      vjust = -0.5, size = 3
    )
  }
  
  # Limitar eje Y si se pasa un valor superior
  if (!is.null(ylim_sup)) {
    p <- p + coord_cartesian(ylim = c(0, ylim_sup))
  }
  
  return(p)
}


estilizar_plot <- function(plot) {
  plot +
    theme_minimal(base_size = 14) +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      strip.text = element_text(size = 16),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 13),
      plot.title = element_text(size = 18, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  # ) +
  # geom_boxplot(size = 1.2, outlier.shape = NA)
}
crear_panel_covid <- function(data, 
                              filename_prefix = "panel", 
                              panel_title = "COVID - Clinical and Protein Markers",
                              mostrar_ids = FALSE) {
  
  # Diccionario de nombres amigables
  titulo_variables <- c(
    "DepresionEscalaUnificada" = "Depresión",
    "AnsiedadEscalaUnificada"  = "Ansiedad",
    "M6a"                      = "M6a",
    "BDNF"                     = "BDNF",
    "NFL"                      = "NFL"
  )
  
  clinicas <- c("DepresionEscalaUnificada", "AnsiedadEscalaUnificada")
  proteinas <- c("M6a", "BDNF", "NFL")
  
  # --- Gráficos clínicos ---
  plots_covid_common <- list()
  for (var in clinicas) {
    p <- crear_boxplot(
      data = data,
      variable = var,
      title = titulo_variables[var],
      ylim_sup = NULL,
      mostrar_ids = mostrar_ids  # 👈 nuevo parámetro
    )
    plots_covid_common[[length(plots_covid_common) + 1]] <- estilizar_plot(p)
  }
  
  # --- Gráficos de proteínas ---
  plots_covid_proteins <- list()
  for (var in proteinas) {
    p <- crear_boxplot(
      data = data,
      variable = var,
      title = titulo_variables[var],
      ylim_sup = NULL,
      mostrar_ids = mostrar_ids  # 👈 nuevo parámetro
    )
    plots_covid_proteins[[length(plots_covid_proteins) + 1]] <- estilizar_plot(p)
  }
  
  # --- Unir en una fila ---
  plots_all <- c(plots_covid_common, plots_covid_proteins)
  
  panel <- wrap_plots(plots_all, ncol = 5) +
    plot_annotation(
      title = panel_title,
      theme = theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
      )
    )
  
  # --- Guardar ---
  ggsave(
    filename = paste0(path_plots_data_analysis, filename_prefix, ".png"),
    plot = panel,
    width = 22, height = 8, dpi = 300
  )
}


crear_panel_cognitivas <- function(data, filename_prefix = "panel_cognitivas", panel_title = "COVID – Variables Cognitivas") {
  # Diccionario de nombres amigables
  titulo_variables <- c(
    "MOCA_perc" = "MoCA (%)",
    "TMT.A_perc" = "TMT-A (%)",
    "FAS" = "FAS"
  )
  
  # Variables cognitivas
  cognitivas <- c("MOCA_perc", "TMT.A_perc", "FAS")
  
  # Lista para los gráficos
  plots_cognitivas <- list()
  for (var in cognitivas) {
    p <- crear_boxplot(
      data = data,
      variable = var,
      title = titulo_variables[var],
      ylim_sup = NULL
    )
    plots_cognitivas[[length(plots_cognitivas) + 1]] <- estilizar_plot(p)
  }
  
  # Crear el panel
  panel <- wrap_plots(plots_cognitivas, ncol = 3) +
    plot_annotation(
      title = panel_title,
      theme = theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
      )
    )
  
  # Guardar imagen
  ggsave(
    filename = paste0(path_plots_data_analysis, filename_prefix, ".png"),
    plot = panel,
    width = 16, height = 6, dpi = 300
  )
}
crear_panel_WMH <- function(data, filename_prefix = "panel_WMH", panel_title = "COVID – WMH Volumes") {
  # Nombres amigables
  titulo_variables <- c(
    "Total_WMH" = "Total WMH",
    "Deep_WMH" = "Deep WMH",
    "PV_WMH" = "Periventricular WMH"
  )
  
  # Variables WMH
  wmh_vars <- c("Total_WMH", "Deep_WMH", "PV_WMH")
  
  # Lista de plots
  plots_wmh <- list()
  for (var in wmh_vars) {
    p <- crear_boxplot(
      data = data,
      variable = var,
      title = titulo_variables[var],
      ylim_sup = NULL
    )
    plots_wmh[[length(plots_wmh) + 1]] <- estilizar_plot(p)
  }
  
  # Crear panel
  panel <- wrap_plots(plots_wmh, ncol = 3) +
    plot_annotation(
      title = panel_title,
      theme = theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
      )
    )
  
  # Guardar imagen
  ggsave(
    filename = paste0(path_plots_data_analysis, filename_prefix, ".png"),
    plot = panel,
    width = 16, height = 6, dpi = 300
  )
}



crear_panel_ASL <- function(data, filename_prefix = "panel_ASL", panel_title = "ASL – Total Gray Matter Metrics") {
  # Nombres amigables
  titulo_variables <- c(
    "TotalGM_sCOV" = "Total GM sCOV",
    "TotalGM_CBF" = "Total GM CBF"
  )
  
  # Variables ASL
  asl_vars <- c("TotalGM_sCOV", "TotalGM_CBF")
  
  # Lista de plots
  plots_asl <- list()
  for (var in asl_vars) {
    p <- crear_boxplot(
      data = data,
      variable = var,
      title = titulo_variables[var],
      ylim_sup = NULL
    )
    plots_asl[[length(plots_asl) + 1]] <- estilizar_plot(p)
  }
  
  # Crear panel
  panel <- wrap_plots(plots_asl, ncol = 2) +
    plot_annotation(
      title = panel_title,
      theme = theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
      )
    )
  
  # Guardar imagen
  ggsave(
    filename = paste0(path_plots_data_analysis, filename_prefix, ".png"),
    plot = panel,
    width = 12, height = 6, dpi = 300
  )
}


