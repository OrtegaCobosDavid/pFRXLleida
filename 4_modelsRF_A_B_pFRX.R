#===========================================================
# RANDOM FOREST - MODELOS A Y B
# - Modelo A: categóricas + química
# - Modelo B: solo química
# - Guarda resultados para calibración posterior (Platt)
# - Versión corregida y robusta
#===========================================================

suppressPackageStartupMessages({
  library(caret)
  library(randomForest)
  library(pROC)
  library(tidyverse)
  library(openxlsx)
  library(forcats)
})

setwd("D:/FRXLleida")
out_dir <- "Resultados"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(1234)

options(timeout = 999999)
options(scipen = 999)
setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)

#========================
# 1) PARÁMETROS
#========================
n_reps  <- 10
p_train <- 0.70
k_folds <- 10
ntree   <- 2000

class_levels <- c("ALC", "TRC", "CST", "MNT", "TRM", "ULL", "VIL-PIR", "VLL")
y_name <- "unit"

vars_cat <- c("environment", "chertType")
vars_num <- c("Al", "Si", "K", "Ca", "Ti", "Mn", "Fe", "Pb", "U", "Sr")

vars_A <- c(vars_cat, vars_num)
vars_B <- vars_num

#========================
# 2) DATOS Y PREPARACIÓN
#========================
if (!exists("geologic")) {
  stop("No existe el objeto 'geologic' en el entorno.")
}
if (!exists("archeo")) {
  warning("No existe el objeto 'archeo'. El script RF puede ejecutarse igualmente.")
}

required_vars <- unique(c(y_name, vars_A))
missing_vars <- setdiff(required_vars, names(geologic))
if (length(missing_vars) > 0) {
  stop("Faltan variables en 'geologic': ", paste(missing_vars, collapse = ", "))
}

# Recodificación segura
geologic <- as.data.frame(geologic)

geologic[[y_name]] <- as.character(geologic[[y_name]])
geologic[[y_name]][geologic[[y_name]] %in% c("VIL", "PIR")] <- "VIL-PIR"
geologic[[y_name]] <- factor(geologic[[y_name]], levels = class_levels)

geologic$chertType <- as.character(geologic$chertType)
geologic$chertType[geologic$chertType %in% c("A", "C", "D")] <- "A"
geologic$chertType <- factor(geologic$chertType)

if (exists("archeo")) {
  archeo <- as.data.frame(archeo)
  if ("chertType" %in% names(archeo)) {
    archeo$chertType <- as.character(archeo$chertType)
    archeo$chertType[archeo$chertType %in% c("A", "C", "D")] <- "A"
    archeo$chertType <- factor(archeo$chertType)
  }
}

for (v in vars_cat) {
  geologic[[v]] <- factor(geologic[[v]])
}
for (v in vars_num) {
  geologic[[v]] <- as.numeric(geologic[[v]])
}

# Diagnóstico temprano
cat("\n", strrep("=", 60), "\n", sep = "")
cat("INFORMACIÓN DEL DATASET\n")
cat(strrep("=", 60), "\n", sep = "")
cat("Filas:", nrow(geologic), "\n")
cat("Columnas:", ncol(geologic), "\n")
cat("Distribución inicial de la variable objetivo:\n")
print(table(geologic[[y_name]], useNA = "ifany"))

unexpected_classes <- setdiff(unique(as.character(geologic[[y_name]])), c(class_levels, NA))
if (length(unexpected_classes) > 0) {
  cat("\nClases inesperadas detectadas:", paste(unexpected_classes, collapse = ", "), "\n")
}

# Eliminar filas con NA en respuesta
n_na_y <- sum(is.na(geologic[[y_name]]))
if (n_na_y > 0) {
  cat("Se eliminan", n_na_y, "filas con NA en", y_name, "\n")
  geologic <- geologic[!is.na(geologic[[y_name]]), , drop = FALSE]
}

cat("\nDistribución final de la variable objetivo:\n")
print(table(geologic[[y_name]], useNA = "ifany"))

missing_classes <- setdiff(class_levels, unique(as.character(geologic[[y_name]])))
if (length(missing_classes) > 0) {
  cat("\nClases ausentes en los datos:", paste(missing_classes, collapse = ", "), "\n")
}

#========================
# 3) FUNCIONES
#========================
align_factor_levels <- function(train_df, new_df, vars) {
  out <- new_df[, vars, drop = FALSE]
  for (v in vars) {
    if (is.factor(train_df[[v]])) {
      out[[v]] <- factor(as.character(out[[v]]), levels = levels(train_df[[v]]))
    } else {
      out[[v]] <- as.numeric(out[[v]])
    }
  }
  as.data.frame(out)
}

complete_prob_matrix <- function(prob_df, class_levels) {
  prob_df <- as.data.frame(prob_df)
  for (cl in setdiff(class_levels, colnames(prob_df))) prob_df[[cl]] <- 0
  prob_df <- prob_df[, class_levels, drop = FALSE]
  rs <- rowSums(prob_df)
  rs[is.na(rs) | rs == 0] <- 1
  prob_df / rs
}

conf_to_df <- function(conf_mat, class_levels) {
  M <- as.matrix(conf_mat)
  Mfull <- matrix(0, nrow = length(class_levels), ncol = length(class_levels),
                  dimnames = list(class_levels, class_levels))
  rr <- intersect(rownames(M), class_levels)
  cc <- intersect(colnames(M), class_levels)
  if (length(rr) > 0 && length(cc) > 0) Mfull[rr, cc] <- M[rr, cc, drop = FALSE]
  
  df <- as.data.frame(Mfull, stringsAsFactors = FALSE)
  df$Real <- rownames(df)
  df <- df[, c("Real", class_levels), drop = FALSE]
  rownames(df) <- NULL
  df
}

per_class_metrics_from_conf <- function(conf_mat) {
  cm <- as.matrix(conf_mat)   # filas = Real, columnas = Pred
  classes <- rownames(cm)
  N <- sum(cm)
  
  bind_rows(lapply(classes, function(cl) {
    TP <- cm[cl, cl]
    FN <- sum(cm[cl, ]) - TP
    FP <- sum(cm[, cl]) - TP
    TN <- N - TP - FN - FP
    
    Sensitivity <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
    Specificity <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
    Precision   <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
    F1 <- if (!is.na(Sensitivity) && !is.na(Precision) && (Sensitivity + Precision) > 0) {
      2 * Precision * Sensitivity / (Precision + Sensitivity)
    } else NA_real_
    
    tibble(
      Class = cl,
      Support = TP + FN,
      PCC = Sensitivity,
      Sensitivity = Sensitivity,
      Specificity = Specificity,
      Precision = Precision,
      F1_Score = F1
    )
  }))
}

per_class_auc_ovr <- function(truth, probs_df, class_levels) {
  bind_rows(lapply(class_levels, function(cl) {
    y_bin <- as.integer(truth == cl)
    if (length(unique(y_bin)) < 2) {
      auc_val <- NA_real_
    } else {
      roc_obj <- tryCatch(
        pROC::roc(y_bin, probs_df[[cl]], quiet = TRUE),
        error = function(e) NULL
      )
      auc_val <- if (is.null(roc_obj)) NA_real_ else as.numeric(pROC::auc(roc_obj))
    }
    tibble(Class = cl, AUC = auc_val)
  }))
}

#========================
# 4) FUNCIÓN PRINCIPAL
#========================
run_rf_experiment <- function(data, vars, y_name, model_name,
                              n_reps, p_train, k_folds, ntree, class_levels,
                              save_models = FALSE,
                              compute_importance = FALSE,
                              keep_confusion = FALSE,
                              keep_roc = FALSE) {
  
  reps <- vector("list", n_reps)
  
  for (i in seq_len(n_reps)) {
    set.seed(100 + i)
    
    # Mantener solo columnas necesarias para este modelo
    model_df <- data[, c(y_name, vars), drop = FALSE]
    model_df <- model_df[complete.cases(model_df), , drop = FALSE]
    model_df[[y_name]] <- factor(as.character(model_df[[y_name]]), levels = class_levels)
    
    if (nrow(model_df) == 0) {
      stop("No quedan filas completas para el modelo ", model_name)
    }
    if (anyNA(model_df[[y_name]])) {
      stop("Persisten NA en la variable respuesta tras la limpieza para el modelo ", model_name)
    }
    
    idx_tr <- createDataPartition(model_df[[y_name]], p = p_train, list = FALSE)
    idx_te <- setdiff(seq_len(nrow(model_df)), idx_tr)
    
    train <- model_df[idx_tr, , drop = FALSE]
    test  <- model_df[idx_te, , drop = FALSE]
    
    train[[y_name]] <- factor(train[[y_name]], levels = class_levels)
    test[[y_name]]  <- factor(test[[y_name]], levels = class_levels)
    
    Xtr <- train[, vars, drop = FALSE]
    Xte <- align_factor_levels(train, test, vars)
    ytr <- train[[y_name]]
    yte <- test[[y_name]]
    
    folds <- createFolds(ytr, k = k_folds, returnTrain = TRUE)
    ctrl  <- trainControl(method = "cv", number = k_folds, index = folds)
    
    grid_mtry <- unique(pmax(1, round(seq(1, length(vars), length.out = min(7, length(vars))))))
    fit_cv <- tryCatch(
      caret::train(
        x = Xtr, y = ytr,
        method = "rf",
        tuneGrid = data.frame(mtry = grid_mtry),
        trControl = ctrl,
        ntree = ntree
      ),
      error = function(e) NULL
    )
    
    best_mtry <- if (!is.null(fit_cv)) fit_cv$bestTune$mtry else max(1, floor(sqrt(length(vars))))
    
    set.seed(2000 + i)
    rf <- randomForest(
      x = Xtr, y = ytr,
      ntree = ntree,
      mtry = best_mtry,
      importance = compute_importance,
      keep.forest = TRUE
    )
    
    pred <- predict(rf, Xte)
    prob_df <- complete_prob_matrix(predict(rf, Xte, type = "prob"), class_levels)
    
    cm <- confusionMatrix(
      data = factor(pred, levels = class_levels),
      reference = factor(yte, levels = class_levels),
      mode = "everything"
    )
    
    auc <- tryCatch({
      as.numeric(multiclass.roc(yte, as.matrix(prob_df))$auc)
    }, error = function(e) NA_real_)
    
    out_i <- list(
      repetition = i,
      mtry = best_mtry,
      metrics = tibble(
        Repetition = i,
        Model = model_name,
        Accuracy = unname(cm$overall["Accuracy"]),
        Balanced_Accuracy = mean(cm$byClass[, "Balanced Accuracy"], na.rm = TRUE),
        Kappa = unname(cm$overall["Kappa"]),
        AUC = auc
      )
    )
    
    if (keep_roc) {
      out_i$roc_data <- list(truth = yte, probs = prob_df)
    }
    
    if (compute_importance) {
      imp <- as.data.frame(importance(rf))
      imp$Variable <- rownames(imp)
      rownames(imp) <- NULL
      imp$Repetition <- i
      imp$Model <- model_name
      out_i$importance <- imp
    }
    
    if (keep_confusion) {
      out_i$confusion <- table(
        Real = factor(yte, levels = class_levels),
        Pred = factor(pred, levels = class_levels)
      )
    }
    
    if (save_models) {
      out_i$model <- rf
      out_i$y_train <- ytr
      out_i$X_train <- Xtr
      out_i$train_idx <- idx_tr
      out_i$test_idx <- idx_te
      out_i$train_data <- train
      out_i$test_data <- test
    }
    
    reps[[i]] <- out_i
    
    cat(sprintf("%s - Rep %02d/%02d | Acc=%.3f | Bal=%.3f | AUC=%.3f%s\n",
                model_name, i, n_reps,
                out_i$metrics$Accuracy,
                out_i$metrics$Balanced_Accuracy,
                out_i$metrics$AUC,
                if (save_models) " [guardado calibración]" else ""))
  }
  
  out <- list(
    metrics = bind_rows(lapply(reps, `[[`, "metrics")),
    mtry_values = sapply(reps, `[[`, "mtry"),
    model_name = model_name,
    classes = class_levels
  )
  
  if (keep_roc) out$roc_data <- lapply(reps, `[[`, "roc_data")
  if (compute_importance) out$importance <- bind_rows(lapply(reps, `[[`, "importance"))
  if (keep_confusion) out$confusion_list <- lapply(reps, `[[`, "confusion")
  
  if (save_models) {
    out$models <- lapply(reps, `[[`, "model")
    out$y_train <- lapply(reps, `[[`, "y_train")
    out$X_train <- lapply(reps, `[[`, "X_train")
    out$train_idx <- lapply(reps, `[[`, "train_idx")
    out$test_idx <- lapply(reps, `[[`, "test_idx")
    out$train_data <- lapply(reps, `[[`, "train_data")
    out$test_data <- lapply(reps, `[[`, "test_data")
    out$vars_used <- vars
  }
  
  out
}

#========================
# 5) EJECUCIÓN
#========================
cat("\n", strrep("=", 60), "\n", sep = "")
cat("MODELO A: categóricas + química\n")
cat(strrep("=", 60), "\n\n", sep = "")

res_A <- run_rf_experiment(
  data = geologic,
  vars = vars_A,
  y_name = y_name,
  model_name = "A_with_categoricals",
  n_reps = n_reps,
  p_train = p_train,
  k_folds = k_folds,
  ntree = ntree,
  class_levels = class_levels,
  save_models = TRUE,
  compute_importance = TRUE,
  keep_confusion = TRUE,
  keep_roc = TRUE
)

cat("\n", strrep("=", 60), "\n", sep = "")
cat("MODELO B: solo química\n")
cat(strrep("=", 60), "\n\n", sep = "")

res_B <- run_rf_experiment(
  data = geologic,
  vars = vars_B,
  y_name = y_name,
  model_name = "B_chemistry_only",
  n_reps = n_reps,
  p_train = p_train,
  k_folds = k_folds,
  ntree = ntree,
  class_levels = class_levels,
  save_models = FALSE,
  compute_importance = FALSE,
  keep_confusion = FALSE,
  keep_roc = FALSE
)

#========================
# 6) OBJETOS PARA CALIBRACIÓN POSTERIOR
#========================
saveRDS(list(A = res_A, B = res_B), file.path(out_dir, "results_minimal_A_B.rds"))

calibration_data_A <- list(
  models = res_A$models,
  y_train = res_A$y_train,
  X_train = res_A$X_train,
  train_idx = res_A$train_idx,
  test_idx = res_A$test_idx,
  train_data = res_A$train_data,
  test_data = res_A$test_data,
  vars = res_A$vars_used,
  class_levels = class_levels,
  n_reps = n_reps
)
saveRDS(calibration_data_A, file.path(out_dir, "calibration_data_A.rds"))

#========================
# 7) MATRICES DE CONFUSIÓN (A)
#========================
conf_sum_A <- Reduce(`+`, res_A$confusion_list)
conf_avg_A <- conf_sum_A / n_reps
conf_pct_A <- prop.table(conf_sum_A, 1) * 100

for (i in seq_len(n_reps)) {
  write.csv(
    conf_to_df(res_A$confusion_list[[i]], class_levels),
    file.path(out_dir, sprintf("A_confusion_rep_%02d.csv", i)),
    row.names = FALSE
  )
}

write.csv(conf_to_df(conf_sum_A, class_levels),
          file.path(out_dir, "A_confusion_SUM.csv"), row.names = FALSE)
write.csv(conf_to_df(round(conf_avg_A, 2), class_levels),
          file.path(out_dir, "A_confusion_AVG.csv"), row.names = FALSE)
write.csv(conf_to_df(round(conf_pct_A, 1), class_levels),
          file.path(out_dir, "A_confusion_PCT_row.csv"), row.names = FALSE)

#========================
# 8) MÉTRICAS POR CLASE (A)
#========================
metrics_by_class_list <- lapply(seq_len(n_reps), function(i) {
  m_i <- per_class_metrics_from_conf(res_A$confusion_list[[i]]) %>%
    mutate(Repetition = i)
  
  auc_i <- per_class_auc_ovr(
    truth = res_A$roc_data[[i]]$truth,
    probs_df = res_A$roc_data[[i]]$probs,
    class_levels = class_levels
  ) %>%
    mutate(Repetition = i)
  
  left_join(m_i, auc_i, by = c("Class", "Repetition"))
})

res_A$class_metrics <- bind_rows(metrics_by_class_list) %>%
  mutate(Model = res_A$model_name)

res_A$class_metrics_summary <- res_A$class_metrics %>%
  group_by(Model, Class) %>%
  summarise(
    n_mean = mean(Support, na.rm = TRUE),
    PCC_mean = mean(PCC, na.rm = TRUE),
    PCC_sd   = sd(PCC, na.rm = TRUE),
    Sensitivity_mean = mean(Sensitivity, na.rm = TRUE),
    Sensitivity_sd   = sd(Sensitivity, na.rm = TRUE),
    Specificity_mean = mean(Specificity, na.rm = TRUE),
    Specificity_sd   = sd(Specificity, na.rm = TRUE),
    Precision_mean = mean(Precision, na.rm = TRUE),
    Precision_sd   = sd(Precision, na.rm = TRUE),
    F1_mean = mean(F1_Score, na.rm = TRUE),
    F1_sd   = sd(F1_Score, na.rm = TRUE),
    AUC_mean = mean(AUC, na.rm = TRUE),
    AUC_sd   = sd(AUC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Class)

write.csv(res_A$class_metrics,
          file.path(out_dir, "A_class_metrics_detailed.csv"),
          row.names = FALSE)

write.csv(res_A$class_metrics_summary,
          file.path(out_dir, "A_class_metrics_summary.csv"),
          row.names = FALSE)

#========================
# 9) COMPARACIÓN GLOBAL A vs B
#========================
metrics_compare <- bind_rows(res_A$metrics, res_B$metrics)

summary_stats <- metrics_compare %>%
  group_by(Model) %>%
  summarise(
    Accuracy_mean = mean(Accuracy, na.rm = TRUE),
    Accuracy_sd   = sd(Accuracy, na.rm = TRUE),
    Balanced_Accuracy_mean = mean(Balanced_Accuracy, na.rm = TRUE),
    Balanced_Accuracy_sd   = sd(Balanced_Accuracy, na.rm = TRUE),
    Kappa_mean = mean(Kappa, na.rm = TRUE),
    Kappa_sd   = sd(Kappa, na.rm = TRUE),
    AUC_mean   = mean(AUC, na.rm = TRUE),
    AUC_sd     = sd(AUC, na.rm = TRUE),
    .groups = "drop"
  )

wide_metrics <- metrics_compare %>%
  select(Repetition, Model, Accuracy, Balanced_Accuracy, Kappa, AUC) %>%
  pivot_wider(names_from = Model, values_from = c(Accuracy, Balanced_Accuracy, Kappa, AUC))

paired_t <- tibble(metric = c("Accuracy", "Balanced_Accuracy", "Kappa", "AUC")) %>%
  rowwise() %>%
  mutate(
    p_value = {
      a <- wide_metrics[[paste0(metric, "_A_with_categoricals")]]
      b <- wide_metrics[[paste0(metric, "_B_chemistry_only")]]
      tryCatch(t.test(a, b, paired = TRUE)$p.value, error = function(e) NA_real_)
    }
  ) %>%
  ungroup()

#========================
# 10) FIGURAS
#========================
plot_metrics_comparison <- metrics_compare %>%
  pivot_longer(c(Accuracy, Balanced_Accuracy, Kappa, AUC),
               names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.08), alpha = 0.5) +
  facet_wrap(~Metric, scales = "free_y") +
  theme_minimal() +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  labs(title = "Global performance metrics: Model A vs Model B",
       x = "", y = "Value") +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "01_global_metrics_comparison.png"),
       plot_metrics_comparison, width = 12, height = 8, dpi = 300)

if (!is.null(res_A$importance)) {
  impA_plot_data <- res_A$importance %>%
    select(Variable, Repetition, MeanDecreaseAccuracy, MeanDecreaseGini) %>%
    pivot_longer(c(MeanDecreaseAccuracy, MeanDecreaseGini),
                 names_to = "ImportanceMetric", values_to = "Importance") %>%
    group_by(ImportanceMetric, Variable) %>%
    summarise(
      mean_imp = mean(Importance, na.rm = TRUE),
      sd_imp   = sd(Importance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(ImportanceMetric) %>%
    slice_max(mean_imp, n = 10, with_ties = FALSE) %>%
    mutate(Variable = fct_reorder(Variable, mean_imp)) %>%
    ungroup()
  
  plot_impA <- ggplot(impA_plot_data, aes(x = Variable, y = mean_imp)) +
    geom_col(alpha = 0.85) +
    geom_errorbar(aes(ymin = mean_imp - sd_imp, ymax = mean_imp + sd_imp), width = 0.2) +
    coord_flip() +
    facet_wrap(~ImportanceMetric, scales = "free") +
    theme_minimal() +
    labs(title = "Model A - Variable importance (mean ± sd, 10 reps)",
         x = "", y = "Importance")
  
  ggsave(file.path(out_dir, "02_A_importance_MDA_Gini_combined.png"),
         plot_impA, width = 12, height = 8, dpi = 300)
  
  rank_table <- res_A$importance %>%
    group_by(Variable) %>%
    summarise(
      MDA_mean  = mean(MeanDecreaseAccuracy, na.rm = TRUE),
      MDA_sd    = sd(MeanDecreaseAccuracy, na.rm = TRUE),
      Gini_mean = mean(MeanDecreaseGini, na.rm = TRUE),
      Gini_sd   = sd(MeanDecreaseGini, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      rank_MDA = dense_rank(desc(MDA_mean)),
      rank_Gini = dense_rank(desc(Gini_mean)),
      rank_combined = (rank_MDA + rank_Gini) / 2
    ) %>%
    arrange(rank_combined, rank_MDA, rank_Gini)
  
  write.csv(rank_table, file.path(out_dir, "A_importance_ranking.csv"), row.names = FALSE)
}

#========================
# 11) EXPORTACIÓN A EXCEL
#========================
wb <- createWorkbook()

addWorksheet(wb, "Global_Summary")
writeData(wb, "Global_Summary", summary_stats)

addWorksheet(wb, "Paired_ttest")
writeData(wb, "Paired_ttest", paired_t)

addWorksheet(wb, "Metrics_A")
writeData(wb, "Metrics_A", res_A$metrics)

addWorksheet(wb, "Metrics_B")
writeData(wb, "Metrics_B", res_B$metrics)

if (!is.null(res_A$importance)) {
  addWorksheet(wb, "Importance_A")
  writeData(wb, "Importance_A", res_A$importance)
}

addWorksheet(wb, "A_ClassMetrics")
writeData(wb, "A_ClassMetrics", res_A$class_metrics)

addWorksheet(wb, "A_ClassMetrics_Sum")
writeData(wb, "A_ClassMetrics_Sum", res_A$class_metrics_summary)

addWorksheet(wb, "A_Conf_SUM")
writeData(wb, "A_Conf_SUM", conf_to_df(conf_sum_A, class_levels))

addWorksheet(wb, "A_Conf_AVG")
writeData(wb, "A_Conf_AVG", conf_to_df(round(conf_avg_A, 2), class_levels))

addWorksheet(wb, "A_Conf_PCT")
writeData(wb, "A_Conf_PCT", conf_to_df(round(conf_pct_A, 1), class_levels))

for (i in seq_len(n_reps)) {
  sh <- sprintf("A_Conf_rep_%02d", i)
  addWorksheet(wb, sh)
  writeData(wb, sh, conf_to_df(res_A$confusion_list[[i]], class_levels))
}

saveWorkbook(wb, file.path(out_dir, "RF_Comparison_Results.xlsx"), overwrite = TRUE)

#========================
# 12) RESUMEN FINAL
#========================
cat("\n", strrep("=", 60), "\n", sep = "")
cat("ANÁLISIS COMPLETADO\n")
cat(strrep("=", 60), "\n\n", sep = "")

cat("Directorio de salida:", out_dir, "\n\n")
cat("Archivos principales:\n")
cat(" - 01_global_metrics_comparison.png\n")
cat(" - 02_A_importance_MDA_Gini_combined.png\n")
cat(" - A_confusion_rep_01..10.csv\n")
cat(" - A_confusion_SUM.csv | A_confusion_AVG.csv | A_confusion_PCT_row.csv\n")
cat(" - A_class_metrics_detailed.csv | A_class_metrics_summary.csv\n")
cat(" - A_importance_ranking.csv\n")
cat(" - RF_Comparison_Results.xlsx\n")
cat(" - calibration_data_A.rds\n")
cat(" - results_minimal_A_B.rds\n\n")

cat("Resumen global:\n")
print(summary_stats)

cat("\nTest t pareado A vs B:\n")
print(paired_t)

cat("\n", strrep("=", 60), "\n", sep = "")
cat("FIN DEL SCRIPT\n")
cat(strrep("=", 60), "\n", sep = "")

#========================
# GUARDAR ENTORNO FINAL
#========================
save.image(file = file.path(out_dir, "pFRX.RData"))

cat("\nEntorno guardado correctamente en:\n")
cat(file.path(out_dir, "pFRX.RData"), "\n")
