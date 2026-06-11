#=============================================================================
# PLATT SCALING DEL MODELO A
# - Evalúa calibración en geologic
# - Aplica calibración a archeo
# - Requiere: calibration_data_A.rds generado por el script RF
# - Versión corregida y alineada con el último script RF
#=============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
  library(reshape2)
  library(ggplot2)
})

setwd("D:/FRXLleida/")
out_dir <- "Resultados"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(1234)

options(timeout = 999999)
options(scipen = 999)
setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)

#=============================================================================
# 1) CARGA
#=============================================================================

calib_file <- file.path(out_dir, "calibration_data_A.rds")
if (!file.exists(calib_file)) {
  stop("No se encuentra: ", calib_file,
       "\nEjecuta antes el script RF con save_models = TRUE")
}

calib <- readRDS(calib_file)

needed <- c("models", "y_train", "X_train", "train_idx", "test_idx",
            "vars", "class_levels", "n_reps")
miss <- setdiff(needed, names(calib))
if (length(miss) > 0) {
  stop("Faltan elementos en calibration_data_A.rds: ", paste(miss, collapse = ", "))
}

rf_models      <- calib$models
y_train_list   <- calib$y_train
train_idx_list <- calib$train_idx
test_idx_list  <- calib$test_idx
vars           <- calib$vars
class_levels   <- calib$class_levels
n_reps         <- calib$n_reps

# Si existen, usar directamente train/test guardados por el script RF
has_saved_splits <- all(c("train_data", "test_data") %in% names(calib))
if (has_saved_splits) {
  train_data_list <- calib$train_data
  test_data_list  <- calib$test_data
} else {
  train_data_list <- NULL
  test_data_list  <- NULL
}

y_name   <- "unit"
vars_cat <- c("environment", "chertType")
vars_num <- c("Al", "Si", "K", "Ca", "Ti", "Mn", "Fe", "Pb", "U", "Sr")

# Usamos objetos del entorno solo si hacen falta para reconstruir splits
if (!has_saved_splits) {
  if (exists("geologic_recod")) {
    geologic_src <- geologic_recod
    message("Usando geologic_recod")
  } else if (exists("geologic")) {
    geologic_src <- geologic
    message("Usando geologic")
  } else {
    stop("No existe 'geologic' ni 'geologic_recod' y el RDS no contiene train_data/test_data.")
  }
} else {
  geologic_src <- NULL
}

if (exists("archeo_recod")) {
  archeo_src <- archeo_recod
  message("Usando archeo_recod")
} else if (exists("archeo")) {
  archeo_src <- archeo
  message("Usando archeo")
} else {
  stop("No existe 'archeo' ni 'archeo_recod'")
}

#=============================================================================
# 2) FUNCIONES
#=============================================================================

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

fit_platt_models <- function(rf_model, y_train, class_levels) {
  prob_oob <- complete_prob_matrix(rf_model$votes, class_levels)
  
  models <- lapply(class_levels, function(cl) {
    dat <- data.frame(
      prob = pmin(pmax(prob_oob[[cl]], 1e-15), 1 - 1e-15),
      y = as.integer(y_train == cl)
    )
    dat <- dat[complete.cases(dat), , drop = FALSE]
    
    if (length(unique(dat$y)) < 2) {
      glm(y ~ 1, data = dat, family = binomial())
    } else {
      tryCatch(
        glm(y ~ prob, data = dat, family = binomial()),
        error = function(e) glm(y ~ 1, data = dat, family = binomial())
      )
    }
  })
  
  names(models) <- class_levels
  models
}

apply_platt_models <- function(prob_raw, platt_models, class_levels) {
  prob_raw <- complete_prob_matrix(prob_raw, class_levels)
  
  prob_cal <- sapply(class_levels, function(cl) {
    predict(platt_models[[cl]],
            newdata = data.frame(prob = prob_raw[[cl]]),
            type = "response")
  })
  
  complete_prob_matrix(as.data.frame(prob_cal), class_levels)
}

score_probs <- function(prob_df, truth, class_levels) {
  truth <- factor(truth, levels = class_levels)
  
  keep <- !is.na(truth)
  prob_df <- prob_df[keep, , drop = FALSE]
  truth <- truth[keep]
  
  pred <- factor(class_levels[max.col(as.matrix(prob_df), ties.method = "first")],
                 levels = class_levels)
  
  y_onehot <- sapply(class_levels, function(cl) as.integer(truth == cl))
  y_onehot <- as.matrix(y_onehot)
  idx_true <- match(as.character(truth), class_levels)
  
  brier <- mean(rowSums((as.matrix(prob_df) - y_onehot)^2))
  logloss <- -mean(sapply(seq_along(idx_true), function(i) log(max(prob_df[i, idx_true[i]], 1e-15))))
  acc <- mean(pred == truth)
  
  prob_true <- prob_df[cbind(seq_along(idx_true), idx_true)]
  correct <- as.numeric(pred == truth)
  ece <- mean(abs(prob_true - correct))
  
  list(
    metrics = c(Accuracy = acc, Brier = brier, LogLoss = logloss, ECE = ece),
    pred = pred,
    truth = truth
  )
}

sum_confusions <- function(lst, class_levels) {
  out <- matrix(0, nrow = length(class_levels), ncol = length(class_levels),
                dimnames = list(Real = class_levels, Pred = class_levels))
  for (m in lst) out[rownames(m), colnames(m)] <- out[rownames(m), colnames(m)] + m
  out
}

class_metrics <- function(conf_mat) {
  N <- sum(conf_mat)
  
  bind_rows(lapply(rownames(conf_mat), function(cl) {
    TP <- conf_mat[cl, cl]
    FN <- sum(conf_mat[cl, ]) - TP
    FP <- sum(conf_mat[, cl]) - TP
    TN <- N - TP - FN - FP
    
    precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
    recall    <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
    f1        <- if (!is.na(precision) && !is.na(recall) && precision + recall > 0) {
      2 * precision * recall / (precision + recall)
    } else NA_real_
    
    tibble(Class = cl, Precision = precision, Recall = recall, F1 = f1)
  }))
}

predict_calibrated <- function(rf_model, platt_models, newdata, class_levels) {
  prob_raw <- complete_prob_matrix(predict(rf_model, newdata, type = "prob"), class_levels)
  prob_cal <- apply_platt_models(prob_raw, platt_models, class_levels)
  
  list(
    class = factor(class_levels[max.col(as.matrix(prob_cal), ties.method = "first")],
                   levels = class_levels),
    raw_prob = prob_raw,
    cal_prob = prob_cal,
    conf = apply(prob_cal, 1, max)
  )
}

#=============================================================================
# 3) AJUSTE DE MODELOS PLATT
#=============================================================================

platt_list <- lapply(seq_len(n_reps), function(i) {
  fit_platt_models(rf_models[[i]], y_train_list[[i]], class_levels)
})

#=============================================================================
# 4) EVALUACIÓN EN GEOLOGIC
#=============================================================================

test_metrics <- vector("list", n_reps)
conf_raw_list <- vector("list", n_reps)
conf_cal_list <- vector("list", n_reps)

for (i in seq_len(n_reps)) {
  
  if (has_saved_splits) {
    train_data <- train_data_list[[i]]
    test_data  <- test_data_list[[i]]
  } else {
    train_data <- geologic_src[train_idx_list[[i]], , drop = FALSE]
    test_data  <- geologic_src[test_idx_list[[i]], , drop = FALSE]
  }
  
  test_data[[y_name]] <- factor(test_data[[y_name]], levels = class_levels)
  
  X_test <- align_factor_levels(train_data, test_data, vars)
  y_test <- test_data[[y_name]]
  
  prob_raw <- complete_prob_matrix(predict(rf_models[[i]], X_test, type = "prob"), class_levels)
  prob_cal <- apply_platt_models(prob_raw, platt_list[[i]], class_levels)
  
  raw <- score_probs(prob_raw, y_test, class_levels)
  cal <- score_probs(prob_cal, y_test, class_levels)
  
  test_metrics[[i]] <- list(
    repetition = i,
    raw = raw$metrics,
    cal = cal$metrics
  )
  
  conf_raw_list[[i]] <- table(
    Real = factor(raw$truth, levels = class_levels),
    Pred = factor(raw$pred, levels = class_levels)
  )
  
  conf_cal_list[[i]] <- table(
    Real = factor(cal$truth, levels = class_levels),
    Pred = factor(cal$pred, levels = class_levels)
  )
}

metrics_df <- bind_rows(lapply(seq_len(n_reps), function(i) {
  bind_rows(
    data.frame(Repetition = i, Type = "Raw_RF", t(test_metrics[[i]]$raw)),
    data.frame(Repetition = i, Type = "Platt",  t(test_metrics[[i]]$cal))
  )
}))

metrics_summary <- metrics_df %>%
  group_by(Type) %>%
  summarise(
    across(
      Accuracy:ECE,
      list(mean = ~mean(.x, na.rm = TRUE),
           sd   = ~sd(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

paired_tests <- tibble(Metric = c("Accuracy", "Brier", "LogLoss", "ECE")) %>%
  mutate(
    p_value = sapply(Metric, function(m) {
      raw <- metrics_df %>%
        filter(Type == "Raw_RF") %>%
        pull(all_of(m))
      
      cal <- metrics_df %>%
        filter(Type == "Platt") %>%
        pull(all_of(m))
      
      tryCatch(
        t.test(cal, raw, paired = TRUE)$p.value,
        error = function(e) NA_real_
      )
    }),
    
    delta_mean = sapply(Metric, function(m) {
      raw <- metrics_df %>%
        filter(Type == "Raw_RF") %>%
        pull(all_of(m))
      
      cal <- metrics_df %>%
        filter(Type == "Platt") %>%
        pull(all_of(m))
      
      mean(cal - raw, na.rm = TRUE)
    })
  )

conf_raw_aggr <- sum_confusions(conf_raw_list, class_levels)
conf_cal_aggr <- sum_confusions(conf_cal_list, class_levels)

metrics_by_class <- full_join(
  class_metrics(conf_raw_aggr) %>% rename_with(~paste0(.x, "_Raw"), -Class),
  class_metrics(conf_cal_aggr) %>% rename_with(~paste0(.x, "_Cal"), -Class),
  by = "Class"
)

#=============================================================================
# 5) APLICACIÓN A ARCHEO
#=============================================================================

missing_vars <- setdiff(vars, names(archeo_src))
if (length(missing_vars) > 0) {
  stop("Faltan variables en archeo: ", paste(missing_vars, collapse = ", "))
}

archeo_prep <- archeo_src
for (v in intersect(vars_num, names(archeo_prep))) archeo_prep[[v]] <- as.numeric(archeo_prep[[v]])
for (v in intersect(vars_cat, names(archeo_prep))) archeo_prep[[v]] <- as.character(archeo_prep[[v]])

ok <- complete.cases(archeo_prep[, vars, drop = FALSE])
if (!all(ok)) message(sum(!ok), " filas eliminadas de archeo por NA en predictores")
archeo_prep <- archeo_prep[ok, , drop = FALSE]

pred_list <- vector("list", n_reps)
raw_prob_list <- vector("list", n_reps)
cal_prob_list <- vector("list", n_reps)
conf_list <- vector("list", n_reps)

for (i in seq_len(n_reps)) {
  
  if (has_saved_splits) {
    train_data <- train_data_list[[i]]
  } else {
    train_data <- geologic_src[train_idx_list[[i]], , drop = FALSE]
  }
  
  X_arch <- align_factor_levels(train_data, archeo_prep, vars)
  
  pred <- predict_calibrated(rf_models[[i]], platt_list[[i]], X_arch, class_levels)
  
  pred_list[[i]] <- as.character(pred$class)
  raw_prob_list[[i]] <- pred$raw_prob
  cal_prob_list[[i]] <- pred$cal_prob
  conf_list[[i]] <- pred$conf
}

pred_matrix <- do.call(cbind, pred_list)
colnames(pred_matrix) <- paste0("Rep_", seq_len(n_reps))

prob_raw_avg <- complete_prob_matrix(Reduce(`+`, raw_prob_list) / n_reps, class_levels)
prob_cal_avg <- complete_prob_matrix(Reduce(`+`, cal_prob_list) / n_reps, class_levels)
confidence_avg <- rowMeans(do.call(cbind, conf_list))

consensus_class <- apply(pred_matrix, 1, function(x) names(which.max(table(x))))
consensus_strength <- apply(pred_matrix, 1, function(x) max(table(x)) / length(x))
entropy <- apply(prob_cal_avg, 1, function(p) -sum(p * log(p + 1e-15)) / log(length(class_levels)))

top_order <- t(apply(prob_cal_avg, 1, function(x) names(sort(x, decreasing = TRUE))))
prob_sorted <- t(apply(prob_cal_avg, 1, function(x) sort(x, decreasing = TRUE)))

archeo_results <- archeo_prep %>%
  mutate(
    Muestra = seq_len(n()),
    Clase_consenso = factor(consensus_class, levels = class_levels),
    Clase_max_prob = factor(
      class_levels[max.col(as.matrix(prob_cal_avg), ties.method = "first")],
      levels = class_levels
    ),
    Segunda_clase = factor(top_order[,2], levels = class_levels),
    
    Confianza = round(confidence_avg, 4),
    Consenso  = round(consensus_strength, 4),
    Entropia  = round(entropy, 4),
    Diferencia = round(prob_sorted[,1] - prob_sorted[,2], 4),
    
    Prob_Top1 = round(prob_sorted[,1], 4),
    Prob_Top2 = round(prob_sorted[,2], 4),
    Prob_Top3 = round(prob_sorted[,3], 4),
    
    Nivel_certeza = case_when(
      Confianza >= 0.70 & Entropia <= 0.40 & Diferencia >= 0.20 ~ "Alta",
      Confianza < 0.50  | Entropia >= 0.60 | Diferencia < 0.10 ~ "Baja",
      TRUE ~ "Media"
    )
  )

#=============================================================================
# 6) FIGURAS
#=============================================================================

prob_long <- prob_cal_avg
prob_long$Muestra <- seq_len(nrow(prob_long))
prob_melt <- reshape2::melt(prob_long, id.vars = "Muestra",
                            variable.name = "Clase", value.name = "Probabilidad")

p1 <- ggplot(prob_melt, aes(x = Clase, y = Muestra, fill = Probabilidad)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#2c3e50", limits = c(0, 1)) +
  theme_minimal() +
  labs(title = "Probabilidades calibradas - archeo",
       subtitle = paste("Promedio de", n_reps, "repeticiones"))

p2 <- ggplot(archeo_results, aes(x = Clase_consenso, fill = Clase_consenso)) +
  geom_bar() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Distribución de clases por consenso", x = "Clase", y = "Frecuencia")

p3 <- ggplot(archeo_results, aes(x = Consenso, y = Confianza, color = Clase_consenso)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(title = "Consenso vs confianza media",
       x = "Consenso entre repeticiones", y = "Confianza media")

p4 <- ggplot(archeo_results, aes(x = Entropia)) +
  geom_histogram(binwidth = 0.05, fill = "#56B4E9", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  labs(title = "Distribución de incertidumbre",
       x = "Entropía normalizada", y = "Frecuencia")

ggsave(file.path(out_dir, "17_archeo_heatmap.png"), p1, width = 12, height = 8, dpi = 300)
ggsave(file.path(out_dir, "18_archeo_distribucion.png"), p2, width = 10, height = 6, dpi = 300)
ggsave(file.path(out_dir, "19_archeo_confianza_vs_consenso.png"), p3, width = 10, height = 7, dpi = 300)
ggsave(file.path(out_dir, "20_archeo_entropia.png"), p4, width = 10, height = 6, dpi = 300)

#=============================================================================
# 7) EXPORTACIÓN
#=============================================================================

saveRDS(
  list(
    platt_models = platt_list,
    test_metrics = test_metrics,
    metrics_df = metrics_df,
    metrics_summary = metrics_summary,
    paired_tests = paired_tests,
    confusion_raw = conf_raw_aggr,
    confusion_cal = conf_cal_aggr,
    metrics_by_class = metrics_by_class
  ),
  file.path(out_dir, "calibracion_A_resultados.rds")
)

saveRDS(
  list(
    predictions = pred_list,
    probs_raw = raw_prob_list,
    probs_cal = cal_prob_list,
    prob_avg = prob_cal_avg,
    pred_matrix = pred_matrix,
    results = archeo_results
  ),
  file.path(out_dir, "archeo_predictions_calibradas.rds")
)

write.csv(archeo_results,
          file.path(out_dir, "21_archeo_resultados_completos.csv"),
          row.names = FALSE)

write.csv(as.data.frame.matrix(conf_raw_aggr),
          file.path(out_dir, "25_confusion_matrix_raw.csv"))

write.csv(as.data.frame.matrix(conf_cal_aggr),
          file.path(out_dir, "26_confusion_matrix_calibrada.csv"))

write.csv(metrics_by_class,
          file.path(out_dir, "27_metrics_by_class.csv"),
          row.names = FALSE)

if (nrow(high_uncertainty) > 0) {
  write.csv(high_uncertainty,
            file.path(out_dir, "archeo_high_uncertainty.csv"),
            row.names = FALSE)
}

wb <- createWorkbook()

addWorksheet(wb, "Geologic_Metrics")
writeData(wb, "Geologic_Metrics", metrics_df)

addWorksheet(wb, "Geologic_Summary")
writeData(wb, "Geologic_Summary", metrics_summary)

addWorksheet(wb, "Paired_Tests")
writeData(wb, "Paired_Tests", paired_tests)

addWorksheet(wb, "Confusion_Raw")
writeData(wb, "Confusion_Raw", as.data.frame.matrix(conf_raw_aggr), rowNames = TRUE)

addWorksheet(wb, "Confusion_Cal")
writeData(wb, "Confusion_Cal", as.data.frame.matrix(conf_cal_aggr), rowNames = TRUE)

addWorksheet(wb, "Class_Metrics")
writeData(wb, "Class_Metrics", metrics_by_class)

addWorksheet(wb, "Archeo_Results")
writeData(wb, "Archeo_Results", archeo_results)

addWorksheet(wb, "Archeo_Prob_Avg")
writeData(wb, "Archeo_Prob_Avg", cbind(Muestra = seq_len(nrow(prob_cal_avg)), prob_cal_avg))

addWorksheet(wb, "Archeo_Pred_By_Rep")
writeData(wb, "Archeo_Pred_By_Rep", cbind(Muestra = seq_len(nrow(pred_matrix)), as.data.frame(pred_matrix)))

saveWorkbook(wb, file.path(out_dir, "22_archeo_resultados.xlsx"), overwrite = TRUE)

#=============================================================================
# 8) RESUMEN
#=============================================================================

cat("\n", strrep("=", 70), "\n", sep = "")
cat("PROCESO COMPLETADO\n")
cat(strrep("=", 70), "\n\n", sep = "")

cat("Resumen geologic:\n")
print(metrics_summary)

cat("\nTests pareados:\n")
print(paired_tests)

cat("\nResumen archeo:\n")
cat(sprintf("Muestras procesadas: %d\n", nrow(archeo_results)))
cat(sprintf("Clases asignadas: %d\n", length(unique(archeo_results$Clase_consenso))))
cat(sprintf("Confianza media: %.3f\n", mean(archeo_results$Confianza)))
cat(sprintf("Consenso medio: %.3f\n", mean(archeo_results$Consenso)))
cat(sprintf("Entropía media: %.3f\n", mean(archeo_results$Entropia)))
if (nrow(high_uncertainty) > 0) {
  cat(sprintf("Muestras con alta incertidumbre: %d\n", nrow(high_uncertainty)))
}

#========================
# GUARDAR ENTORNO FINAL
#========================
save.image(file = file.path(out_dir, "pFRX.RData"))

cat("\nEntorno guardado correctamente en:\n")
cat(file.path(out_dir, "pFRX.RData"), "\n")

cat("\nFIN DEL SCRIPT\n")

