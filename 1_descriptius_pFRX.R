# ====================================================
# ANÀLISI DE DADES pFRX DE CHERTS
# ====================================================
# Descripció: Anàlisi descriptiu, boxplots i correlacions
# per a mostres geològiques i arqueològiques
# revisat per ChatGPT: 21/4/26, 9:09
# ====================================================

# ----------------------------
# 0. CONFIGURACIÓ INICIAL
# ----------------------------
setwd("D:/FRXLleida")

out_dir <- "Resultados"
dir.create(out_dir, showWarnings = FALSE)

# TEMP robusto para openxlsx
tmp_dir <- "D:/FRXLleida/temp"
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

Sys.setenv(TMP = tmp_dir)
Sys.setenv(TEMP = tmp_dir)

set.seed(1234)

suppressPackageStartupMessages({
  library(readxl)
  library(tidyverse)
  library(Hmisc)
  library(corrplot)
  library(openxlsx)
  library(grid)
})

# ----------------------------
# 1. IMPORTACIÓ I PREPARACIÓ
# ----------------------------
elements <- c("Al", "Si", "K", "Ca", "Ti", "Mn", "Fe", "Pb", "U", "Sr")
order_unit <- c("ALC", "TRC", "CST", "MNT", "TRM", "ULL", "PIR", "VIL", "VLL")
order_chertType <- c("A", "B", "C", "D", "E", "F")
order_element <- c("Si", "Ca", "Sr", "Al", "K", "Ti", "Fe", "Mn", "Pb", "U")

round_2 <- function(x) if (is.numeric(x)) round(x, 2) else x

pFRX_data <- read_excel(
  "D:/FRXLleida/pFRX_data2.xls",
  col_types = c("text", "text", "text", "text", "text", "text", "text", "text",
                rep("numeric", 10))
) %>%
  mutate(
    across(c(collection, site, unit, environment, chertType), as.factor)
  )

archeo <- pFRX_data %>% filter(collection == "arch")
geologic <- pFRX_data %>% filter(collection == "geo")

if ("idLab" %in% names(archeo)) {
  geologic <- archeo %>%
    filter(site == "Montvell", idLab != 100) %>%
    bind_rows(geologic)
} else {
  warning("La columna 'idLab' no existe. No es poden filtrar casos de Montvell.")
}

# ----------------------------
# 2. FUNCIONS AUXILIARS
# ----------------------------
make_global_stats <- function(df, vars) {
  df %>%
    select(all_of(vars)) %>%
    summarise(across(
      everything(),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd = ~sd(.x, na.rm = TRUE),
        median = ~median(.x, na.rm = TRUE),
        min = ~min(.x, na.rm = TRUE),
        max = ~max(.x, na.rm = TRUE),
        q25 = ~quantile(.x, 0.25, na.rm = TRUE),
        q75 = ~quantile(.x, 0.75, na.rm = TRUE),
        n = ~sum(!is.na(.x))
      ),
      .names = "{.col}_{.fn}"
    )) %>%
    pivot_longer(everything(), names_to = "Estadístic", values_to = "Valor") %>%
    mutate(Valor = round_2(Valor))
}

make_group_stats <- function(df, group_var, vars) {
  df %>%
    group_by(.data[[group_var]]) %>%
    summarise(across(
      all_of(vars),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd = ~sd(.x, na.rm = TRUE),
        median = ~median(.x, na.rm = TRUE),
        min = ~min(.x, na.rm = TRUE),
        max = ~max(.x, na.rm = TRUE),
        n = ~sum(!is.na(.x))
      ),
      .names = "{.col}_{.fn}"
    ), .groups = "drop") %>%
    arrange(.data[[group_var]]) %>%
    mutate(across(-all_of(group_var), round_2))
}

make_publication_table <- function(df, group_var, vars) {
  tmp <- df %>%
    group_by(.data[[group_var]]) %>%
    summarise(across(
      all_of(vars),
      list(mitjana = ~mean(.x, na.rm = TRUE),
           sd = ~sd(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ), .groups = "drop") %>%
    arrange(.data[[group_var]])
  
  out <- tmp %>% select(all_of(group_var))
  
  for (el in vars) {
    out[[el]] <- paste0(
      format(round(tmp[[paste0(el, "_mitjana")]], 2), nsmall = 2),
      " ± ",
      format(round(tmp[[paste0(el, "_sd")]], 2), nsmall = 2)
    )
  }
  out
}

make_corr_outputs <- function(df, vars, order_vars) {
  rc <- rcorr(as.matrix(df %>% select(all_of(vars))))
  cor_mat <- rc$r
  p_mat <- rc$P
  
  present <- intersect(order_vars, rownames(cor_mat))
  cor_mat <- cor_mat[present, present, drop = FALSE]
  p_mat <- p_mat[present, present, drop = FALSE]
  p_mat[is.na(p_mat)] <- 1
  
  cor_df <- as.data.frame(cor_mat) %>%
    rownames_to_column("Element") %>%
    mutate(across(where(is.numeric), round_2))
  
  p_df <- as.data.frame(p_mat) %>%
    rownames_to_column("Element") %>%
    mutate(across(where(is.numeric), round_2))
  
  list(cor = cor_mat, p = p_mat, cor_df = cor_df, p_df = p_df)
}

make_boxplot <- function(df, group_var, group_levels, vars, var_levels) {
  plot_df <- df %>%
    select(all_of(c(group_var, vars))) %>%
    pivot_longer(cols = all_of(vars), names_to = "variable", values_to = "value")
  
  plot_df[[group_var]] <- factor(plot_df[[group_var]], levels = group_levels)
  plot_df$variable <- factor(plot_df$variable, levels = var_levels)
  
  ggplot(plot_df, aes(x = .data[[group_var]], y = value, fill = .data[[group_var]])) +
    geom_boxplot(alpha = 0.8, outlier.size = 1.5) +
    facet_wrap(~ variable, scales = "free_y", ncol = 2) +
    scale_fill_brewer(palette = "Set3") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "none",
      panel.spacing = unit(1, "lines"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    labs(title = NULL)
}

save_corrplot_tiff <- function(file, cor_mat, p_mat = NULL, significant = FALSE) {
  tiff(file, width = 8, height = 8, units = "in", res = 300, compression = "lzw")
  if (significant) {
    corrplot(
      cor_mat,
      type = "upper",
      order = "original",
      tl.col = "black",
      tl.srt = 45,
      p.mat = p_mat,
      insig = "label_sig",
      sig.level = c(0.001, 0.01, 0.05),
      pch.cex = 1.2,
      pch.col = "white",
      diag = TRUE,
      title = NULL,
      mar = c(0, 0, 2, 0)
    )
  } else {
    corrplot(
      cor_mat,
      type = "upper",
      order = "original",
      tl.col = "black",
      tl.srt = 45,
      diag = TRUE,
      title = NULL,
      mar = c(0, 0, 2, 0)
    )
  }
  dev.off()
}

# ----------------------------
# 3. TAULES
# ----------------------------

tabla_unit_site <- geologic %>%
  count(unit, site) %>%
  group_by(unit) %>%
  mutate(total_unit = sum(n)) %>%
  arrange(unit, desc(n)) %>%
  ungroup() %>%
  mutate(across(where(is.numeric), round_2))

tabla_typeChert_site <- archeo %>%
  count(chertType, site) %>%
  group_by(site) %>%
  mutate(total_site = sum(n)) %>%
  arrange(site, desc(n)) %>%
  ungroup() %>%
  mutate(across(where(is.numeric), round_2))

geologic_global <- make_global_stats(geologic, elements)
geologic_per_unit <- make_group_stats(geologic, "unit", elements)
publication_table <- make_publication_table(geologic, "unit", elements)

archeo_global <- make_global_stats(archeo, elements)
archeo_per_type <- make_group_stats(archeo, "chertType", elements)
publication_table_archeo <- make_publication_table(archeo, "chertType", elements)

corr_geo <- make_corr_outputs(geologic, elements, order_element)
corr_archeo <- make_corr_outputs(archeo, elements, order_element)

# ----------------------------
# 4. EXPORTACIÓ A EXCEL
# ----------------------------

llista_taules <- list(
  "Resum_geologic_unitats" = tabla_unit_site,
  "Resum_archeo_tipus" = tabla_typeChert_site,
  "Estad_globals_geologic" = geologic_global,
  "Estad_per_unitat_geologic" = geologic_per_unit,
  "Publicacio_geologic" = publication_table,
  "Estad_globals_archeo" = archeo_global,
  "Estad_per_tipus_archeo" = archeo_per_type,
  "Publicacio_archeo" = publication_table_archeo,
  "Correlacions_geo" = corr_geo$cor_df,
  "Pvalors_geo" = corr_geo$p_df,
  "Correlacions_archeo" = corr_archeo$cor_df,
  "Pvalors_archeo" = corr_archeo$p_df
)

excel_file <- file.path(out_dir, "resultats_complets.xlsx")
write.xlsx(llista_taules, file = excel_file, overwrite = TRUE)

cat("\n✅ Taules exportades a:", excel_file, "\n")
cat("Nombre de fulls creats:", length(llista_taules), "\n")

# ----------------------------
# 5. GRÀFICS
# ----------------------------

plot_element <- make_boxplot(
  df = geologic,
  group_var = "unit",
  group_levels = order_unit,
  vars = elements,
  var_levels = order_element
)

ggsave(
  file.path(out_dir, "boxplot_geologic.tiff"),
  plot = plot_element,
  width = 12, height = 10, dpi = 300,
  compression = "lzw", bg = "white"
)

plot_tipus_chert <- make_boxplot(
  df = archeo,
  group_var = "chertType",
  group_levels = order_chertType,
  vars = elements,
  var_levels = order_element
)

ggsave(
  file.path(out_dir, "boxplot_archeo.tiff"),
  plot = plot_tipus_chert,
  width = 12, height = 10, dpi = 300,
  compression = "lzw", bg = "white"
)

save_corrplot_tiff(
  file.path(out_dir, "correlacions_geologic_ordenat.tiff"),
  corr_geo$cor,
  significant = FALSE
)

save_corrplot_tiff(
  file.path(out_dir, "correlacions_geologic_significatives_ordenat.tiff"),
  corr_geo$cor,
  p_mat = corr_geo$p,
  significant = TRUE
)

save_corrplot_tiff(
  file.path(out_dir, "correlacions_archeo_ordenat.tiff"),
  corr_archeo$cor,
  significant = FALSE
)

save_corrplot_tiff(
  file.path(out_dir, "correlacions_archeo_significatives_ordenat.tiff"),
  corr_archeo$cor,
  p_mat = corr_archeo$p,
  significant = TRUE
)

# ----------------------------
# 6. RESUM FINAL
# ----------------------------

cat("\n", strrep("=", 50), "\n", sep = "")
cat("EXECUCIÓ COMPLETADA AMB ÈXIT\n")
cat(strrep("=", 50), "\n", sep = "")
cat("Arxius generats a:", out_dir, "\n")
cat("\n-- FITXER EXCEL (amb 2 decimals) --\n")
cat("  - resultats_complets.xlsx\n")
cat("\n-- GRÀFICS TIFF (300 ppp) --\n")
cat("  - boxplot_geologic.tiff\n")
cat("  - boxplot_archeo.tiff\n")
cat("  - correlacions_geologic_ordenat.tiff\n")
cat("  - correlacions_geologic_significatives_ordenat.tiff\n")
cat("  - correlacions_archeo_ordenat.tiff\n")
cat("  - correlacions_archeo_significatives_ordenat.tiff\n")
cat(strrep("=", 50), "\n", sep = "")

#========================
# GUARDAR ENTORNO FINAL
#========================
save.image(file = file.path(out_dir, "pFRX.RData"))

cat("\nEntorno guardado correctamente en:\n")
cat(file.path(out_dir, "pFRX.RData"), "\n")
