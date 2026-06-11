#===================================================
# PROVES MULTIVARIANTS DE COMPARACIÓ PER GRUPS
# I PCA EXPLORATORI
# amb puntuacions z i recodificació chert C i D a A,
# i VIL + PIR a VIL-PIR
#===================================================

# ----------------------------
# 0. CONFIGURACIÓ INICIAL
# ----------------------------
setwd("D:/FRXLleida")
out_dir <- "Resultados"
dir.create(out_dir, showWarnings = FALSE)
set.seed(1234)
options(scipen = 999)

options(timeout = 999999)
options(scipen = 999)
setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(psych)
  library(ggplot2)
})

# ----------------------------
# 1. VALIDACIÓ I PREPARACIÓ
# ----------------------------
if (!exists("geologic")) stop("No existeix l'objecte 'geologic' al entorno.")
if (!exists("archeo")) stop("No existeix l'objecte 'archeo' al entorno.")

variables <- c("Al", "Si", "K", "Ca", "Ti", "Mn", "Fe", "Pb", "U", "Sr")
required_geo <- c("unit", "chertType", variables)
required_arch <- c("chertType", variables)

missing_geo <- setdiff(required_geo, names(geologic))
missing_arch <- setdiff(required_arch, names(archeo))

if (length(missing_geo) > 0) {
  stop("Falten variables en geologic: ", paste(missing_geo, collapse = ", "))
}
if (length(missing_arch) > 0) {
  stop("Falten variables en archeo: ", paste(missing_arch, collapse = ", "))
}

# Escalat z només per variables químiques
geologic_z <- geologic %>%
  mutate(across(all_of(variables), ~ as.numeric(scale(.x))))

archeo_z <- archeo %>%
  mutate(across(all_of(variables), ~ as.numeric(scale(.x))))

# Recodificació consistent
recode_chert <- function(x) {
  x <- as.character(x)
  x[x %in% c("A", "C", "D")] <- "A"
  factor(x)
}

recode_unit <- function(x) {
  x <- as.character(x)
  x[x %in% c("VIL", "PIR")] <- "VIL-PIR"
  factor(x)
}

geologic_z$unit <- recode_unit(geologic_z$unit)
geologic_z$chertType <- recode_chert(geologic_z$chertType)
archeo_z$chertType <- recode_chert(archeo_z$chertType)

# Subconjunt net per ANOSIM, MRPP i PCA
geo_mv <- geologic_z %>%
  select(unit, all_of(variables)) %>%
  filter(complete.cases(.))

if (nrow(geo_mv) == 0) {
  stop("No queden files completes a geologic després del filtratge.")
}

cat("Mostres geològiques per a anàlisi multivariant:", nrow(geo_mv), "\n")
cat("Distribució per unitat:\n")
print(table(geo_mv$unit, useNA = "ifany"))

# ----------------------------
# 2. ANOSIM
# ----------------------------
dist_matrix <- dist(geo_mv[, variables], method = "euclidean")

anosim_result <- anosim(
  x = dist_matrix,
  grouping = geo_mv$unit,
  permutations = 9999
)

cat("\n====================\nANOSIM\n====================\n")
print(anosim_result)

# ----------------------------
# 3. MRPP
# ----------------------------
set.seed(123)
mrpp_result <- mrpp(
  dat = dist_matrix,
  grouping = geo_mv$unit,
  permutations = 9999
)

cat("\n====================\nMRPP\n====================\n")
print(mrpp_result)

# ----------------------------
# 4. DIAGNÒSTICS PCA
# ----------------------------
cor_matrix <- cor(geo_mv[, variables], use = "pairwise.complete.obs")

high_cor <- sum(abs(cor_matrix[upper.tri(cor_matrix)]) > 0.3)
total_cor <- length(cor_matrix[upper.tri(cor_matrix)])

cat("\nCorrelaciones > |0.3|:", high_cor, "de", total_cor,
    "(", round(high_cor / total_cor * 100, 1), "%)\n")
cat("¿Suficiente para PCA?:", ifelse(high_cor / total_cor > 0.3, "SÍ", "NO"), "\n\n")

bartlett_test <- cortest.bartlett(cor_matrix, n = nrow(geo_mv))
print(bartlett_test)
cat("p-value:", bartlett_test$p.value, "\n")
cat("¿PCA justificado?:", ifelse(bartlett_test$p.value < 0.05, "SÍ", "NO"), "\n\n")

kmo_result <- KMO(cor_matrix)
cat("KMO:", round(kmo_result$MSA, 3), "\n")
cat("Interpretación:",
    if (kmo_result$MSA >= 0.8) {
      "BUENA"
    } else if (kmo_result$MSA >= 0.7) {
      "ACEPTABLE"
    } else {
      "LIMITADA"
    }, "\n")

# ----------------------------
# 5. PCA
# ----------------------------
pca_result <- prcomp(geo_mv[, variables], scale. = FALSE)

var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)
print(var_explained)

# Scree plot base R
n_vars <- length(var_explained)
eigenvalues <- var_explained * n_vars / 100

plot(1:n_vars, eigenvalues, type = "b", pch = 19,
     main = "Varianza por Componente",
     xlab = "Número de Componente", ylab = "Eigenvalue",
     ylim = c(0, max(eigenvalues) * 1.1))
abline(h = 1, col = "red", lty = 2)
abline(h = mean(eigenvalues), col = "blue", lty = 3)
text(1:n_vars, eigenvalues, round(eigenvalues, 1), pos = 3, cex = 0.8)

# Resultats primers 3 components
var_explained_3 <- summary(pca_result)$importance[2, 1:3]
cum_var_3 <- summary(pca_result)$importance[3, 1:3]
pca_loadings_3 <- pca_result$rotation[, 1:3, drop = FALSE]

cat("\nVarianza explicada 3 primeros componentes:\n")
print(round(var_explained_3 * 100, 1))

cat("\nVarianza acumulada 3 primeros componentes:\n")
print(round(cum_var_3 * 100, 1))

# ----------------------------
# 6. GRÀFIC DE CÀRREGUES PCA (PC1 vs PC2)
# ----------------------------
loadings_df <- data.frame(
  Elemento = rownames(pca_loadings_3),
  PC1 = pca_loadings_3[, 1],
  PC2 = pca_loadings_3[, 2],
  PC3 = pca_loadings_3[, 3],
  row.names = NULL
)

grafico_1 <- ggplot(loadings_df, aes(x = PC1, y = PC2, label = Elemento)) +
  geom_point(size = 3, color = "black") +
  geom_text(size = 3.5, vjust = -0.8, hjust = 0.5, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  labs(
    title = NULL,
    x = paste0("PC1 (", round(var_explained_3[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(var_explained_3[2] * 100, 1), "%)")
  ) +
  coord_fixed(ratio = 1) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "gray70", linewidth = 0.5),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    aspect.ratio = 0.7
  )

print(grafico_1)

ggsave(
  filename = file.path(out_dir, "carga_PCA_1_2_geoquimica.tiff"),
  plot = grafico_1,
  width = 10, height = 8, dpi = 300,
  compression = "lzw", bg = "white", device = "tiff"
)

# ----------------------------
# 7. GRÀFICS DE SCORES PCA
# ----------------------------
pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PC3 = pca_result$x[, 3],
  Unit = geo_mv$unit
)

unit_levels <- sort(unique(as.character(pca_df$Unit)))

palette_vals <- c("blue", "brown", "green", "purple", "orange", "red",
                  "darkgreen", "pink", "cyan")
palette_vals <- rep(palette_vals, length.out = length(unit_levels))
names(palette_vals) <- unit_levels

shape_vals <- seq_along(unit_levels)
names(shape_vals) <- unit_levels

grafico_3 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Unit, shape = Unit)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, alpha = 0.2, linewidth = 0.8, type = "norm") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  scale_color_manual(values = palette_vals) +
  scale_shape_manual(values = shape_vals) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
    color = "Unit",
    shape = "Unit"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "gray70", linewidth = 0.5),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10)
  )

print(grafico_3)

ggsave(
  filename = file.path(out_dir, "pca_1_2_geoquimica.tiff"),
  plot = grafico_3,
  width = 10, height = 8, dpi = 300,
  compression = "lzw", bg = "white", device = "tiff"
)

grafico_4 <- ggplot(pca_df, aes(x = PC1, y = PC3, color = Unit, shape = Unit)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, alpha = 0.2, linewidth = 0.8, type = "norm") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  scale_color_manual(values = palette_vals) +
  scale_shape_manual(values = shape_vals) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC3 (", round(var_explained[3], 1), "%)"),
    shape = "Unit",
    color = "Unit"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "gray70", linewidth = 0.5),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10)
  )

print(grafico_4)

ggsave(
  filename = file.path(out_dir, "pca_1_3_geoquimica.tiff"),
  plot = grafico_4,
  width = 10, height = 8, dpi = 300,
  compression = "lzw", bg = "white", device = "tiff"
)

# ----------------------------
# 8. RESUM FINAL
# ----------------------------
cat("\n", strrep("=", 50), "\n", sep = "")
cat("ANÀLISI MULTIVARIANT COMPLETAT\n")
cat(strrep("=", 50), "\n", sep = "")
cat("Arxius gràfics generats a:", out_dir, "\n")
cat(" - carga_PCA_1_2_geoquimica.tiff\n")
cat(" - pca_1_2_geoquimica.tiff\n")
cat(" - pca_1_3_geoquimica.tiff\n")
cat(strrep("=", 50), "\n", sep = "")

#========================
# GUARDAR ENTORNO FINAL
#========================
save.image(file = file.path(out_dir, "pFRX.RData"))

cat("\nEntorno guardado correctamente en:\n")
cat(file.path(out_dir, "pFRX.RData"), "\n")