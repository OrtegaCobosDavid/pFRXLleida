#===================================================
# FUNCIONA PERFECTE, NO TOCAR RES
# selecció de millor escala/transformació de dades
# i comparació de resultats de RF
# codi revisat per ChatGpt
# dades: pFR_data.xlsx
# amb recodificacions
# 27/3/2026, 12:00
#===================================================

# ----------------------------
# 0. CONFIGURACIÓ INICIAL
# ----------------------------

# Establir directori de treball
setwd("D:/FRXLleida")
out_dir <- "Resultados"
dir.create(out_dir, showWarnings = FALSE)
set.seed(1234)

options(timeout = 999999)
options(scipen = 999)
setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)

# Carregar llibreries
library(readxl)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Hmisc)
library(corrplot)
library(openxlsx)

# ----------------------------
# 1. IMPORTACIÓ I PREPARACIÓ DE DADES
# ----------------------------

# recodifiem unit i chertType
geologic$unit <- fct_drop(fct_collapse(geologic$unit, "VIL-PIR" = c("VIL", "PIR")))
geologic$chertType<- fct_drop(fct_collapse(geologic$chertType, "A" = c("A", "C", "D")))
archeo$chertType<- fct_drop(fct_collapse(archeo$chertType, "A" = c("A", "C", "D")))
view(geologic)
view(archeo)

#===========================================================
# TRANSFORMACIONES DE VARIABLES (posiciones 9 a 18)
#===========================================================

# Identificar las columnas a transformar (posiciones 9 a 18)
cols_to_transform <- 9:18
vars_to_transform <- names(geologic)[cols_to_transform]

cat("Variables a transformar:\n")
print(vars_to_transform)

# Asegurar que las variables son numéricas
for (col in vars_to_transform) {
  geologic[[col]] <- as.numeric(geologic[[col]])
}

# Verificar que no hay NAs que puedan causar problemas
if (any(is.na(geologic[, vars_to_transform]))) {
  cat("\nADVERTENCIA: Se encontraron NAs en las variables a transformar.\n")
  cat("Filas con NAs:", sum(complete.cases(geologic[, vars_to_transform]) == FALSE), "\n")
}

#========================
# 1) Puntuaciones Z (z-score)
#========================
geologic_zscore <- geologic
geologic_zscore[, vars_to_transform] <- scale(geologic[, vars_to_transform], 
                                              center = TRUE, 
                                              scale = TRUE)

cat("\n✓ Transformación z-score completada\n")

#========================
# 2) RobustScaler (mediana e IQR)
#========================
geologic_robust <- geologic
geologic_robust[, vars_to_transform] <- as.data.frame(
  lapply(geologic[, vars_to_transform], function(x) {
    med <- median(x, na.rm = TRUE)
    iqr <- IQR(x, na.rm = TRUE)
    if (iqr == 0) (x - med) else (x - med) / iqr
  })
)

cat("✓ Transformación RobustScaler completada\n")

#========================
# 3) Logaritmo natural (log)
#========================
geologic_log <- geologic

# Manejar valores <= 0 para log (desplazamiento si es necesario)
shift_applied <- FALSE
for (col in vars_to_transform) {
  min_val <- min(geologic_log[[col]], na.rm = TRUE)
  if (min_val <= 0) {
    shift <- abs(min_val) + 0.001
    geologic_log[[col]] <- geologic_log[[col]] + shift
    shift_applied <- TRUE
    cat("  Desplazamiento aplicado en", col, "(sumando", shift, ")\n")
  }
  geologic_log[[col]] <- log(geologic_log[[col]])
}

if (!shift_applied) {
  geologic_log[, vars_to_transform] <- log(geologic[, vars_to_transform])
}

cat("✓ Transformación logarítmica natural completada\n")

#========================
# 4) Centered log-ratio (clr)
#========================
geologic_clr <- geologic

# Manejar valores <= 0 para clr
for (col in vars_to_transform) {
  min_val <- min(geologic_clr[[col]], na.rm = TRUE)
  if (min_val <= 0) {
    shift <- abs(min_val) + 0.001
    geologic_clr[[col]] <- geologic_clr[[col]] + shift
    cat("  Desplazamiento aplicado en", col, "para clr (sumando", shift, ")\n")
  }
}

# Aplicar clr por fila
clr_matrix <- t(apply(geologic_clr[, vars_to_transform], 1, function(row) {
  gm <- exp(mean(log(row), na.rm = TRUE))
  log(row / gm)
}))

geologic_clr[, vars_to_transform] <- as.data.frame(clr_matrix)
colnames(geologic_clr)[cols_to_transform] <- vars_to_transform

cat("✓ Transformación centered log-ratio (clr) completada\n")


#===========================================================
# RANDOM FOREST - COMPARACIÓN DE TRANSFORMACIONES
# Comparación de F1-score y Accuracy
#===========================================================

#========================
# 1) CONFIGURACIÓN
#========================
suppressPackageStartupMessages({
  library(caret)
  library(randomForest)
  library(pROC)
  library(tidyverse)
})

set.seed(1234)

cfg <- list(
  n_reps = 10,
  p_train = 0.70,
  k_folds = 10,
  ntree = 2000,
  y_name = "unit"
)

# Lista de datasets a evaluar
datasets <- list(
  "Sin_transformar" = geologic,
  "Z-score" = geologic_zscore,
  "RobustScaler" = geologic_robust,
  "Log_natural" = geologic_log,
  "Centered_log_ratio" = geologic_clr
)

# Variables predictoras (posiciones 9 a 18)
vars_predictoras <- names(geologic)[9:18]

cat("Variables predictoras:\n")
print(vars_predictoras)

cat("\nVariable objetivo:", cfg$y_name, "\n")

#========================
# 2) FUNCIONES AUXILIARES
#========================

# Función para calcular métricas por clase
per_class_metrics <- function(conf_mat) {
  cm <- as.matrix(conf_mat)
  classes <- rownames(cm)
  N <- sum(cm)
  
  bind_rows(lapply(classes, function(cl) {
    TP <- cm[cl, cl]
    FN <- sum(cm[cl, ]) - TP
    FP <- sum(cm[, cl]) - TP
    TN <- N - TP - FN - FP
    
    Sensitivity <- if((TP + FN) > 0) TP / (TP + FN) else NA_real_
    Specificity <- if((TN + FP) > 0) TN / (TN + FP) else NA_real_
    Precision <- if((TP + FP) > 0) TP / (TP + FP) else NA_real_
    F1 <- if(!is.na(Precision) && !is.na(Sensitivity) && (Precision + Sensitivity) > 0) {
      2 * Precision * Sensitivity / (Precision + Sensitivity)
    } else NA_real_
    
    tibble(
      Class = cl,
      Support = TP + FN,
      Accuracy = Sensitivity,
      Sensitivity = Sensitivity,
      Specificity = Specificity,
      Precision = Precision,
      F1_Score = F1
    )
  }))
}

# Función para ejecutar una repetición
fit_one_repetition <- function(data, transform_name, rep_id, cfg, vars_predictoras) {
  set.seed(100 + rep_id)
  
  # Eliminar filas con NA en predictoras o objetivo
  data_complete <- data[complete.cases(data[, c(vars_predictoras, cfg$y_name)]), ]
  
  # Asegurar que la variable objetivo es factor
  data_complete[[cfg$y_name]] <- factor(data_complete[[cfg$y_name]])
  
  # Verificar que hay al menos 2 clases
  if (length(unique(data_complete[[cfg$y_name]])) < 2) {
    cat("  Advertencia: Solo una clase en esta repetición\n")
    return(NULL)
  }
  
  # División train/test
  idx_tr <- createDataPartition(data_complete[[cfg$y_name]], p = cfg$p_train, list = FALSE)
  
  # Verificar que train y test tienen todas las clases
  train_classes <- unique(data_complete[[cfg$y_name]][idx_tr])
  test_classes <- unique(data_complete[[cfg$y_name]][-idx_tr])
  
  if (length(train_classes) < 2 || length(test_classes) < 2) {
    cat("  Advertencia: Train o test no tienen suficientes clases\n")
    return(NULL)
  }
  
  train <- data_complete[idx_tr, , drop = FALSE]
  test <- data_complete[-idx_tr, , drop = FALSE]
  
  Xtr <- as.data.frame(train[, vars_predictoras, drop = FALSE])
  Xte <- as.data.frame(test[, vars_predictoras, drop = FALSE])
  ytr <- train[[cfg$y_name]]
  yte <- test[[cfg$y_name]]
  
  # Validación cruzada para mtry
  folds <- tryCatch(createFolds(ytr, k = cfg$k_folds, returnTrain = TRUE),
                    error = function(e) NULL)
  
  if (is.null(folds) || length(unique(ytr)) < 2) {
    # Si falla CV, usar mtry por defecto
    best_mtry <- floor(sqrt(length(vars_predictoras)))
  } else {
    ctrl <- trainControl(method = "cv", number = cfg$k_folds, index = folds)
    grid_mtry <- unique(pmax(1, round(seq(1, length(vars_predictoras), length.out = 7))))
    
    fit_cv <- tryCatch(
      caret::train(
        x = Xtr, y = ytr,
        method = "rf",
        tuneGrid = data.frame(mtry = grid_mtry),
        trControl = ctrl,
        ntree = cfg$ntree
      ),
      error = function(e) NULL
    )
    
    best_mtry <- if(!is.null(fit_cv)) fit_cv$bestTune$mtry else floor(sqrt(length(vars_predictoras)))
  }
  
  # Entrenar modelo final
  set.seed(2000 + rep_id)
  rf <- randomForest(
    x = Xtr, y = ytr,
    ntree = cfg$ntree,
    mtry = best_mtry,
    importance = FALSE
  )
  
  # Predicciones
  pred <- predict(rf, Xte)
  
  # Matriz de confusión
  conf_mat <- table(Real = yte, Pred = pred)
  
  # Métricas globales
  cm <- confusionMatrix(pred, yte)
  accuracy <- cm$overall["Accuracy"]
  
  # Métricas por clase
  class_metrics <- per_class_metrics(conf_mat)
  
  # F1-score promedio (macro)
  f1_macro <- mean(class_metrics$F1_Score, na.rm = TRUE)
  
  # Resultados
  res <- list(
    transform = transform_name,
    repetition = rep_id,
    mtry = best_mtry,
    accuracy = accuracy,
    f1_macro = f1_macro,
    class_metrics = class_metrics,
    confusion = conf_mat
  )
  
  cat(sprintf("  %s - Rep %02d | Acc=%.3f | F1=%.3f\n",
              transform_name, rep_id, accuracy, f1_macro))
  
  res
}

# Función para ejecutar todas las repeticiones de un dataset
run_experiment <- function(data, transform_name, cfg, vars_predictoras) {
  cat("\n", strrep("-", 50), "\n")
  cat("Procesando:", transform_name, "\n")
  cat(strrep("-", 50), "\n")
  
  cat("  Filas totales:", nrow(data), "\n")
  cat("  Clases en unit:", paste(levels(data[[cfg$y_name]]), collapse = ", "), "\n")
  
  results <- list()
  for (i in 1:cfg$n_reps) {
    res <- fit_one_repetition(data, transform_name, i, cfg, vars_predictoras)
    if (!is.null(res)) {
      results[[i]] <- res
    }
  }
  
  # Combinar resultados
  if (length(results) == 0) {
    cat("  ERROR: No se completó ninguna repetición\n")
    return(NULL)
  }
  
  metrics_df <- bind_rows(lapply(results, function(x) {
    tibble(
      Transform = x$transform,
      Repetition = x$repetition,
      Accuracy = x$accuracy,
      F1_Macro = x$f1_macro,
      Mtry = x$mtry
    )
  }))
  
  class_metrics_df <- bind_rows(lapply(results, function(x) {
    x$class_metrics %>% mutate(Repetition = x$repetition, Transform = x$transform)
  }))
  
  list(
    transform_name = transform_name,
    metrics = metrics_df,
    class_metrics = class_metrics_df,
    repetitions_completed = length(results),
    all_results = results
  )
}

#========================
# 3) EJECUTAR EXPERIMENTOS
#========================
cat("\n", strrep("=", 70), "\n")
cat("INICIO DEL EXPERIMENTO - RANDOM FOREST\n")
cat("Comparación de transformaciones de datos\n")
cat(strrep("=", 70), "\n")
cat("Número de repeticiones:", cfg$n_reps, "\n")
cat("Variables predictoras:", length(vars_predictoras), "\n\n")

# Ejecutar para cada dataset
experiment_results <- list()
for (name in names(datasets)) {
  experiment_results[[name]] <- run_experiment(
    datasets[[name]], name, cfg, vars_predictoras
  )
}

#========================
# 4) COMBINAR Y ANALIZAR RESULTADOS
#========================
# Combinar métricas globales de todos los experimentos
all_metrics <- bind_rows(lapply(experiment_results, function(x) {
  if (!is.null(x)) x$metrics
}))

# Combinar métricas por clase
all_class_metrics <- bind_rows(lapply(experiment_results, function(x) {
  if (!is.null(x)) x$class_metrics
}))

# Estadísticas resumen por transformación
summary_stats <- all_metrics %>%
  group_by(Transform) %>%
  summarise(
    Accuracy_mean = mean(Accuracy, na.rm = TRUE),
    Accuracy_sd = sd(Accuracy, na.rm = TRUE),
    F1_Macro_mean = mean(F1_Macro, na.rm = TRUE),
    F1_Macro_sd = sd(F1_Macro, na.rm = TRUE),
    Reps_completadas = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(Accuracy_mean))

# Estadísticas por clase
class_summary <- all_class_metrics %>%
  group_by(Transform, Class) %>%
  summarise(
    F1_mean = mean(F1_Score, na.rm = TRUE),
    F1_sd = sd(F1_Score, na.rm = TRUE),
    Support_mean = mean(Support, na.rm = TRUE),
    .groups = "drop"
  )

#========================
# 5) MOSTRAR RESULTADOS
#========================
cat("\n", strrep("=", 70), "\n")
cat("RESULTADOS DEL EXPERIMENTO\n")
cat(strrep("=", 70), "\n\n")

cat("RESUMEN GLOBAL:\n")
print(summary_stats)

cat("\n\nMEJOR TRANSFORMACIÓN POR ACCURACY:\n")
best_acc <- summary_stats %>% slice_max(Accuracy_mean, n = 1)
print(best_acc)

cat("\nMEJOR TRANSFORMACIÓN POR F1-MACRO:\n")
best_f1 <- summary_stats %>% slice_max(F1_Macro_mean, n = 1)
print(best_f1)

cat("\n\nRESUMEN POR CLASE (F1-score):\n")
print(class_summary %>%
        arrange(Transform, desc(F1_mean)) %>%
        filter(Transform %in% c("Sin_transformar", "Z-score", "RobustScaler", 
                                "Log_natural", "Centered_log_ratio")) %>%
        head(20))

# ANOVA per comparar resultats de les 10 repeticions per cada mètode
resultats_complets <- all_metrics %>%
  filter(Transform %in% c("Log_natural", "RobustScaler", "Sin_transformar", "Z-score", "Centered_log_ratio")) %>%
  select(Transform, Accuracy, F1_Macro, Repetition) %>%
  rename(Metode = Transform)

# ANOVA per Accuracy
anova_accuracy <- aov(Accuracy ~ Metode, data = resultats_complets)
summary_accuracy <- summary(anova_accuracy)

# ANOVA per F1-Macro
anova_f1 <- aov(F1_Macro ~ Metode, data = resultats_complets)
summary_f1 <- summary(anova_f1)

# Crear dataframes amb els resultats dels ANOVA
anova_accuracy_df <- as.data.frame(summary_accuracy[[1]]) %>%
  rownames_to_column("Font") %>%
  mutate(Metrica = "Accuracy")

anova_f1_df <- as.data.frame(summary_f1[[1]]) %>%
  rownames_to_column("Font") %>%
  mutate(Metrica = "F1_Macro")

# Combinar resultats ANOVA
anova_results <- bind_rows(anova_accuracy_df, anova_f1_df) %>%
  select(Metrica, Font, Df, `Sum Sq`, `Mean Sq`, `F value`, `Pr(>F)`)

# Afegir interpretació
anova_results <- anova_results %>%
  mutate(
    Significatiu = ifelse(`Pr(>F)` < 0.05, "Sí (p < 0.05)", "No (p > 0.05)"),
    Interpretacio = case_when(
      `Pr(>F)` < 0.001 ~ "Diferències altament significatives",
      `Pr(>F)` < 0.01 ~ "Diferències significatives",
      `Pr(>F)` < 0.05 ~ "Diferències marginalment significatives",
      TRUE ~ "No hi ha diferències significatives"
    )
  )

# Carregar el workbook existent o crear un de nou
if(file.exists(file.path(out_dir, "resultats_transformacions.xlsx"))) {
  wb_escalados <- loadWorkbook(file.path(out_dir, "resultats_transformacions.xlsx"))
} else {
  wb_escalados <- createWorkbook()
}

# Afegir fulla amb resultats ANOVA
addWorksheet(wb_escalados, "ANOVA_Comparacio")
writeData(wb_escalados, "ANOVA_Comparacio", anova_results)

# Afegir estadístics descriptius per context
descriptius <- all_metrics %>%
  filter(Transform %in% c("Log_natural", "RobustScaler", "Sin_transformar", "Z-score", "Centered_log_ratio")) %>%
  group_by(Transform) %>%
  summarise(
    Accuracy_mitjana = mean(Accuracy),
    Accuracy_sd = sd(Accuracy),
    F1_mitjana = mean(F1_Macro),
    F1_sd = sd(F1_Macro),
    .groups = "drop"
  ) %>%
  arrange(desc(Accuracy_mitjana))

addWorksheet(wb_escalados, "Estadistics_Descriptius")
writeData(wb_escalados, "Estadistics_Descriptius", descriptius)

# Guardar el workbook
saveWorkbook(wb_escalados, file.path(out_dir, "resultats_transformacions.xlsx"), overwrite = TRUE)

cat("\n✓ Resultats ANOVA afegits al fitxer Excel:\n")
cat("  - Full 'ANOVA_Comparacio': Resultats dels tests ANOVA\n")
cat("  - Full 'Estadistics_Descriptius': Mitjanes i desviacions per mètode\n")


#========================
# 6) GUARDAR RESULTADOS
#========================
out_dir <- "Resultados"
dir.create(out_dir, showWarnings = FALSE)

resultados_finales <- list(
  config = cfg,
  variables_predictoras = vars_predictoras,
  summary_stats = summary_stats,
  all_metrics = all_metrics,
  class_summary = class_summary,
  experiment_results = experiment_results,
  date = Sys.time()
)

saveRDS(resultados_finales, file.path(out_dir, "RF_comparacion_transformaciones.rds"))

# Exportar a CSV
write.csv(summary_stats, file.path(out_dir, "resumen_transformaciones.csv"), row.names = FALSE)
write.csv(all_metrics, file.path(out_dir, "metricas_todas_transformaciones.csv"), row.names = FALSE)
write.csv(class_summary, file.path(out_dir, "metricas_por_clase_transformaciones.csv"), row.names = FALSE)

# Crear un nuevo libro de Excel
wb_escalados <- createWorkbook()

# Añadir cada dataframe como una hoja
addWorksheet(wb_escalados, "Resumen Transformaciones")
writeData(wb_escalados, "Resumen Transformaciones", summary_stats)

addWorksheet(wb_escalados, "Métricas Totales")
writeData(wb_escalados, "Métricas Totales", all_metrics)

addWorksheet(wb_escalados, "Métricas por Clase")
writeData(wb_escalados, "Métricas por Clase", class_summary)

# Guardar el libro
saveWorkbook(wb_escalados, file.path(out_dir, "resultados_transformaciones.xlsx"), overwrite = TRUE)

cat("\n", strrep("=", 70), "\n")
cat("ANÁLISIS COMPLETADO\n")
cat(strrep("=", 70), "\n")
cat("Resultados guardados en:", out_dir, "\n")
cat("  - RF_comparacion_transformaciones.rds\n")
cat("  - resumen_transformaciones.csv\n")
cat("  - metricas_todas_transformaciones.csv\n")
cat("  - metricas_por_clase_transformaciones.csv\n")

#========================
# GUARDAR ENTORNO FINAL
#========================
save.image(file = file.path(out_dir, "pFRX.RData"))

cat("\nEntorno guardado correctamente en:\n")
cat(file.path(out_dir, "pFRX.RData"), "\n")
