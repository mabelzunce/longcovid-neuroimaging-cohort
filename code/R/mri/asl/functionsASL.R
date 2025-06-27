
# Cargar librerías necesarias
library(mgcv)
library(car)       # Para Anova
library(emmeans)   # Para calcular medias ajustadas
library(dplyr)
library(ppcor)
library(tidyr)

# =========================
# ✅ FUNCIÓN GAM ANALYSIS ✅
# =========================
gam_analysis <- function(df, columns) {
  
  # Covariables fijas (siempre aplicadas)
  covariates <- c("Edad", "Genero")
  
  for (col_name in columns) {
    
    if (!(col_name %in% colnames(df))) {
      next  # Saltar si la columna no existe
    }
    
    # Test de normalidad
    shapiro_test <- shapiro.test(df[[col_name]])
    df$Grupo <- as.factor(df$Grupo)
    
    if (shapiro_test$p.value >= 0.05) {  
      
      cat("\n📊 **Performing ANCOVA for", col_name, "**\n")
      df$Grupo <- as.factor(df$Grupo)
      temp_df <- na.omit(df[, c(col_name, "Grupo", "Edad", "Genero")])
      
  

      if (nrow(temp_df) == 0) {
        cat("\n⚠️ ERROR: No data available after removing NA for", col_name, ". Skipping...\n")
        next
      }

      # Definir modelo ANCOVA con covariables personalizadas
      formula <- as.formula(paste(col_name, "~ Grupo + Edad + Genero"))

      # Ajustar modelo ANCOVA
      ancova_result <- lm(formula, data = temp_df)

      # Aplicar ANCOVA con Type III Sum of Squares
      model <- Anova(ancova_result, type = "III")
      
      p_values <- model$"Pr(>F)"
       
      if (p_values[2] < 0.05) {
        cat("\n🔍 **Number of rows before ANCOVA for", col_name, ":**", nrow(temp_df), "\n")
        cat("\n✅ **Grupo is significant in ANCOVA for", col_name, "** (p =", p_values["Grupo"], ")\n")
        print(model["Grupo", , drop = FALSE])

        em_means <- suppressMessages(emmeans(ancova_result, ~ Grupo))
        em_means_df <- as.data.frame(em_means)

        cat("\n📊 **Adjusted Group Means for", col_name, "**\n")
        print(em_means_df)
      }
      
    } else {  
      cat("\n📊 **Performing GAM for", col_name, "**\n")
      
      # Asegurar que "Edad" siempre se trate como numérica y "Genero" como categórica
      df$Edad <- as.numeric(df$Edad)
      df$Genero <- as.factor(df$Genero)
      
      # Crear la fórmula GAM correctamente
      formula <- as.formula(paste(col_name, "~ Grupo + s(Edad) + Genero"))
      
      # Ajustar modelo GAM
      gam_model <- gam(formula, data = df, family = gaussian())
      gam_summary <- summary(gam_model)
      

      # Extraer p-value de Grupo
      if ("GrupoCOVID" %in% rownames(gam_summary$p.table)) {
        grupo_p_value <- gam_summary$p.table["GrupoCOVID", "Pr(>|t|)"]

      if (grupo_p_value < 0.1) {
          # Calcular medias ajustadas solo si Grupo es significativo
          em_means <- suppressMessages(emmeans(gam_model, ~ Grupo))
          em_means_df <- as.data.frame(em_means)
          
          cat("\n📊 **Adjusted Group Means for", col_name, "**\n")
          print(em_means_df)
        }
      }
      
      
      
      
      
      
      }
  }
}


# ===========================
# ✅ FUNCIÓN MULTIPLE ANCOVA ✅
# ===========================

multiple_ancova <- function(df, dependent_vars) {
  
  for (col_name in dependent_vars) {
    
    if (!(col_name %in% colnames(df))) {
      cat("\n⚠️ ERROR: Column", col_name, "not found in dataframe. Skipping...\n")
      next
    }
    
    cat("\n📊 **Performing ANCOVA for", col_name, "**\n")

    df$Grupo <- as.factor(df$Grupo)
    
    temp_df <- na.omit(df[, c(col_name, "Grupo", "Edad", "Genero")])


    if (nrow(temp_df) == 0) {
      cat("\n⚠️ ERROR: No data available after removing NA for", col_name, ". Skipping...\n")
      next
    }

    # Definir modelo ANCOVA con covariables personalizadas
    formula <- as.formula(paste(col_name, "~ Grupo + Edad + Genero"))

    # Ajustar modelo ANCOVA
    ancova_result <- lm(formula, data = temp_df)

    # Aplicar ANCOVA con Type III Sum of Squares
    model <- Anova(ancova_result, type = "III")

    p_values <- model$"Pr(>F)"
    
    
    if (p_values[2] < 0.05) {

      cat("\n🔍 **Number of rows before ANCOVA for", col_name, ":**", nrow(temp_df), "\n")
      cat("\n✅ **Grupo is significant in ANCOVA for", col_name, "** (p =", p_values["Grupo"], ")\n")
      print(model)

      em_means <- suppressMessages(emmeans(ancova_result, ~ Grupo))
      em_means_df <- as.data.frame(em_means)

      cat("\n📊 **Adjusted Group Means for", col_name, "**\n")
      print(em_means_df)
    }
  }
}


t_test_analysis <- function(df, columns) {
  
  for (col_name in columns) {
    
    if (!(col_name %in% colnames(df))) {
      next  # Saltar si la columna no existe
    }
    
    df$Grupo <- as.factor(df$Grupo)
    
    # Eliminar NAs antes de aplicar el t-test
    temp_df <- na.omit(df[, c(col_name, "Grupo")])
    
    if (nrow(temp_df) == 0) {
      next  # Saltar si no hay datos después de eliminar NA
    }
    
    # Aplicar t-test
    t_test_result <- t.test(temp_df[[col_name]] ~ temp_df$Grupo)
    
    # Mostrar solo si es significativo
    if (t_test_result$p.value < 0.05) {
      cat("\n📊 **Significant t-test results for", col_name, "**\n")
      print(t_test_result)
      
      # Calcular medias y número de datos por grupo
      means <- temp_df %>%
        group_by(Grupo) %>%
        summarise(
          n = n(),  # Número de datos
          Mean = mean(!!sym(col_name), na.rm = TRUE)
        )
      
      cat("\n📊 **Group Means and Sample Size for", col_name, "**\n")
      print(means)
    }
  }
}

kruskal_wallis_analysis <- function(df, columns) {
  
  for (col_name in columns) {
    
    if (!(col_name %in% colnames(df))) {
      next  # Saltar si la columna no existe
    }
    
    df$Grupo <- as.factor(df$Grupo)
    
    # Eliminar NAs antes de aplicar Kruskal-Wallis
    temp_df <- na.omit(df[, c(col_name, "Grupo")])
    
    if (nrow(temp_df) == 0) {
      next  # Saltar si no hay datos después de eliminar NA
    }
    
    # Aplicar Kruskal-Wallis test
    kruskal_test_result <- kruskal.test(temp_df[[col_name]] ~ temp_df$Grupo)
    
    # Mostrar solo si es significativo
    if (kruskal_test_result$p.value < 0.05) {
      cat("\n📊 **Significant Kruskal-Wallis results for", col_name, "**\n")
      print(kruskal_test_result)

      # Calcular medianas y número de datos por grupo
      medians <- temp_df %>%
        group_by(Grupo) %>%
        summarise(
          n = n(),  # Número de datos
          Median = median(!!sym(col_name), na.rm = TRUE)
        )

      cat("\n📊 **Group Medians and Sample Size for", col_name, "**\n")
      print(medians)
    }
    
    # cat("\n📊 **Significant Kruskal-Wallis results for", col_name, "**\n")
    # print(kruskal_test_result)
    # 
    # # Calcular medianas y número de datos por grupo
    # medians <- temp_df %>%
    #   group_by(Grupo) %>%
    #   summarise(
    #     n = n(),  # Número de datos
    #     Median = median(!!sym(col_name), na.rm = TRUE)
    #   )
    # 
    # cat("\n📊 **Group Medians and Sample Size for", col_name, "**\n")
    # print(medians)
  
    }
}

partial_correlation_tests <- function(df, brain_vars, cognitive_vars) {
  
  library(ppcor)
  
  # Check if control variables exist
  control_vars <- c("Edad", "Genero")
  missing_controls <- setdiff(control_vars, colnames(df))
  
  if (length(missing_controls) > 0) {
    cat("\n❌ ERROR: Missing control variables:", paste(missing_controls, collapse = ", "), "\n")
    return()
  }
  
  # Convert Genero to numeric if it is categorical
  df$Genero <- as.numeric(as.factor(df$Genero))
  
  for (brain_col in brain_vars) {
    for (cog_col in cognitive_vars) {
      
      # Check if the columns exist
      if (!(brain_col %in% colnames(df))) {
        cat("\n⚠️ ERROR: Column", brain_col, "not found in dataframe. Skipping...\n")
        next
      }
      if (!(cog_col %in% colnames(df))) {
        cat("\n⚠️ ERROR: Column", cog_col, "not found in dataframe. Skipping...\n")
        next
      }
      
      # Remove NA values
      temp_df <- na.omit(df[, c(brain_col, cog_col, control_vars)])
      
      if (nrow(temp_df) == 0) {
        cat("\n⚠️ ERROR: No data available after removing NA for", brain_col, "-", cog_col, ". Skipping...\n")
        next
      }
      
      # Compute partial correlation
      result <- pcor.test(temp_df[[brain_col]], temp_df[[cog_col]], temp_df[, control_vars])
      
      # Show only significant results (p < 0.05)
      if (result$p.value < 0.05 && abs(result$estimate) > 0.2) {
        cat("\n📊 **Significant Partial Correlation Between", brain_col, "and", cog_col, "**\n")
        cat("🔹 Correlation coefficient:", round(result$estimate, 3), "\n")
        cat("🔹 p-value:", format.pval(result$p.value, digits = 4), "\n")
        print(result)
      }
    }
  }
}

multiple_ancova_old <- function(df, dependent_vars) {
  
  for (col_name in dependent_vars) {
    
    if (!(col_name %in% colnames(df))) {
      cat("\n⚠️ ERROR: Column", col_name, "not found in dataframe. Skipping...\n")
      next
    }
    
    cat("\n📊 **Performing ANCOVA for", col_name, "**\n")
    
    df$Grupo <- as.factor(df$Grupo)
    
    temp_df <- na.omit(df[, c(col_name, "Grupo", "Edad", "Genero")])
    
    df$Edad <- as.numeric(df$Edad)
    if (any(is.na(df$Edad))) {
      cat("Advertencia: Se introdujeron NAs por coerción en 'Edad'.\n")
    }
    
    if (nrow(temp_df) == 0) {
      cat("\n⚠️ ERROR: No data available after removing NA for", col_name, ". Skipping...\n")
      next
    }
    
    # Definir modelo ANCOVA con covariables personalizadas
    formula <- as.formula(paste(col_name, "~ Grupo + Edad + Genero "))
    
    # Ajustar modelo ANCOVA
    ancova_result <- lm(formula, data = temp_df)
    
    # Aplicar ANCOVA con Type III Sum of Squares
    model <- Anova(ancova_result, type = "III")
    
    p_values <- model$"Pr(>F)"
    
    if (p_values[2] < 0.1) {
      
      cat("\n🔍 **Number of rows before ANCOVA for", col_name, ":**", nrow(temp_df), "\n")
       print(model)
      
      em_means <- suppressMessages(emmeans(ancova_result, ~ Grupo))
      em_means_df <- as.data.frame(em_means)
      
      cat("\n📊 **Adjusted Group Means for", col_name, "**\n")
      print(em_means_df)
    }
  }
}


# Function to detect outliers using IQR (only for COVID group)
detect_outliers_covid <- function(df, cols) {
  # df_covid <- df %>% filter(Grupo == "COVID")  # Filter only COVID group
  outlier_rows <- df %>%
    pivot_longer(all_of(cols), names_to = "Region", values_to = "Value") %>%
    group_by(Region) %>%
    mutate(Q1 = quantile(Value, 0.25, na.rm = TRUE),
           Q3 = quantile(Value, 0.75, na.rm = TRUE),
           IQR = Q3 - Q1,
           Lower_Bound = Q1 - 1.5 * IQR,
           Upper_Bound = Q3 + 1.5 * IQR,
           Outlier = (Value < Lower_Bound)) %>%
    filter(Outlier) #%>%

  return(outlier_rows)
}



multiple_ancova_correcting_scov <- function(df, dependent_vars) {
  
  for (col_name in dependent_vars) {
    
    if (!(col_name %in% colnames(df))) {
      cat("\n⚠️ ERROR: Column", col_name, "not found in dataframe. Skipping...\n")
      next
    }
    
    cat("\n📊 **Performing ANCOVA for", col_name, "**\n")
    
    df$Grupo <- as.factor(df$Grupo)
    
    temp_df <- na.omit(df[, c(col_name, "Grupo", "Edad", "Genero", "TotalGM_B_scov")])
    
    df$Edad <- as.numeric(df$Edad)
    if (any(is.na(df$Edad))) {
      cat("Advertencia: Se introdujeron NAs por coerción en 'Edad'.\n")
    }
    
    if (nrow(temp_df) == 0) {
      cat("\n⚠️ ERROR: No data available after removing NA for", col_name, ". Skipping...\n")
      next
    }
    
    # Definir modelo ANCOVA con covariables personalizadas
    formula <- as.formula(paste(col_name, "~ Grupo + Edad + Genero + TotalGM_B_scov"))
    
    # Ajustar modelo ANCOVA
    ancova_result <- lm(formula, data = temp_df)
    
    # Aplicar ANCOVA con Type III Sum of Squares
    model <- Anova(ancova_result, type = "III")
    
    p_values <- model$"Pr(>F)"
    
    if (p_values[2] < 0.05) {
      
      cat("\n🔍 **Number of rows before ANCOVA for", col_name, ":**", nrow(temp_df), "\n")
      print(model)
      
      em_means <- suppressMessages(emmeans(ancova_result, ~ Grupo))
      em_means_df <- as.data.frame(em_means)
      
      cat("\n📊 **Adjusted Group Means for", col_name, "**\n")
      print(em_means_df)
    }
  }
}
