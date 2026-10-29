#!/usr/bin/env Rscript
# ============================================================
# Generate Phosphorus Interaction Plots for Alpha Diversity GLMMs
# Both Alternatives requested by supervisor:
#   1. Linetype (Solid for P0 vs Dashed for P3)
#   2. Duplicated figure:
#      - Grid 2x2 (facet_grid: Rows = Phosphorus, Cols = Pop)
#      - Separate standalone figures (P0 and P3)
# ============================================================

suppressPackageStartupMessages({
  library(glmmTMB)
  library(emmeans)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  # Default paths if run directly in the repository
  shannon_file  <- "glmm_inputs/shannon_metadata.tsv"
  observed_file <- "glmm_inputs/observed_features_metadata.tsv"
  faithpd_file  <- "glmm_inputs/faith_pd_metadata.tsv"
  metadata_file <- "00-Ampullaceana_balthica_raw_data/Ampullaceana_balthica_metadata.tsv"
  out_dir       <- "glmm_outputs"
} else {
  shannon_file  <- args[1]
  observed_file <- args[2]
  faithpd_file  <- args[3]
  metadata_file <- args[4]
  out_dir       <- args[5]
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading metadata and diversity data...\n")
metadata <- read.table(metadata_file, sep="\t", header=TRUE, check.names=FALSE)
if (metadata[1,1] == "#q2:types") metadata <- metadata[-1,]

shannon_data  <- read.table(shannon_file, sep="\t", header=TRUE, check.names=FALSE)
observed_data <- read.table(observed_file, sep="\t", header=TRUE, check.names=FALSE)
faithpd_data  <- read.table(faithpd_file,  sep="\t", header=TRUE, check.names=FALSE)

colnames(shannon_data)[1]  <- "SampleID"
colnames(observed_data)[1] <- "SampleID"
colnames(faithpd_data)[1]  <- "SampleID"
colnames(metadata)[1]      <- "SampleID"

# Merge
df <- metadata %>%
  inner_join(shannon_data %>% select(SampleID, shannon_entropy), by="SampleID") %>%
  inner_join(observed_data %>% select(SampleID, observed_features), by="SampleID") %>%
  inner_join(faithpd_data %>% select(SampleID, faith_pd), by="SampleID")

df <- df %>%
  mutate(
    Temp = as.factor(Temp),
    Diet = as.factor(Diet),
    Phosphorus = as.factor(Phosphorus),
    Pop = as.factor(Pop),
    Box = as.factor(Box)
  )

# Calculate Initial_Size and Final_Size vectors via PCA
ini_vars <- c("Ini_weight", "Shell_IL", "Shell_IAA")
fin_vars <- c("Fin_weight", "Shell_FL", "Shell_FAA")
other_vars <- c("N_clutches_per day", "Av_stick_height", "Growth_rate")
all_morpho_vars <- intersect(unique(c(ini_vars, fin_vars, other_vars)), colnames(df))

if (length(all_morpho_vars) >= 2) {
  for (var in all_morpho_vars) {
    df[[var]] <- as.numeric(df[[var]])
    df[[var]][is.na(df[[var]])] <- mean(df[[var]], na.rm = TRUE)
  }
  
  if (all(ini_vars %in% colnames(df))) {
    pca_ini <- prcomp(df[, ini_vars], center = TRUE, scale. = TRUE)
    df$Initial_Size <- pca_ini$x[, 1]
    if (cor(df$Initial_Size, df$Ini_weight) < 0) df$Initial_Size <- -df$Initial_Size
  }
  
  if (all(fin_vars %in% colnames(df))) {
    pca_fin <- prcomp(df[, fin_vars], center = TRUE, scale. = TRUE)
    df$Final_Size <- pca_fin$x[, 1]
    if (cor(df$Final_Size, df$Fin_weight) < 0) df$Final_Size <- -df$Final_Size
  }
  
  pca_vars <- intersect(c("Initial_Size", "Final_Size", other_vars), colnames(df))
  pca_res <- prcomp(df[, pca_vars], center = TRUE, scale. = TRUE)
  df$Size_PC1 <- df$Initial_Size # Covariável dos GLMMs: Tamanho Inicial (pré-tratamento)
  formula_rhs <- "Pop * Temp * Diet * Phosphorus + Size_PC1 + (1|Box)"
} else {
  formula_rhs <- "Pop * Temp * Diet * Phosphorus + (1|Box)"
}

diet_colors <- c("A" = "#740000", "M" = "#DAA520", "P" = "#4F734E")

clean_theme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "italic", size = 12),
    legend.position = "right"
  )

generate_phosphorus_plots <- function(metric_col, ylab_name) {
  cat("Fitting GLMM for", metric_col, "...\n")
  sub_df <- df[!is.na(df[[metric_col]]), ]
  sub_df[[metric_col]] <- as.numeric(sub_df[[metric_col]])
  
  fml <- as.formula(paste0(metric_col, " ~ ", formula_rhs))
  model <- glmmTMB(fml, data = sub_df)
  
  cat("Calculating Estimated Marginal Means (emmeans)...\n")
  emm <- as.data.frame(emmeans(model, ~ Temp * Diet * Pop * Phosphorus, type = "response"))
  
  # Standardize confidence limit column names (asymp.LCL vs lower.CL)
  if ("asymp.LCL" %in% colnames(emm)) {
    emm$ymin <- emm$asymp.LCL
    emm$ymax <- emm$asymp.UCL
  } else if ("lower.CL" %in% colnames(emm)) {
    emm$ymin <- emm$lower.CL
    emm$ymax <- emm$upper.CL
  } else {
    emm$ymin <- emm$emmean - 1.96 * emm$SE
    emm$ymax <- emm$emmean + 1.96 * emm$SE
  }
  
  emm$Temp <- factor(emm$Temp)
  emm$Diet <- factor(emm$Diet, levels = c("A", "M", "P"))
  emm$Pop  <- factor(emm$Pop)
  emm$Phosphorus <- factor(emm$Phosphorus, levels = c("0", "3"))
  
  # -------------------------------------------------------------------------
  # ALTERNATIVA 1: Linetype para Fósforo (0 = solid, 3 = dashed)
  # -------------------------------------------------------------------------
  cat("  -> Saving Alternativa 1 (Linetype)...\n")
  dodge_line <- position_dodge(width = 0.6)
  p_line <- ggplot(emm, aes(x = Temp, y = emmean, color = Diet, linetype = Phosphorus, shape = Phosphorus,
                           group = interaction(Diet, Phosphorus))) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.3, position = dodge_line, linewidth = 0.8) +
    geom_line(position = dodge_line, linewidth = 1) +
    geom_point(position = dodge_line, size = 3) +
    scale_color_manual(name = "Diet", values = diet_colors) +
    scale_linetype_manual(name = "Phosphorus", values = c("0" = "solid", "3" = "dashed"),
                          labels = c("0" = "0 (P0)", "3" = "3 (P3)")) +
    scale_shape_manual(name = "Phosphorus", values = c("0" = 16, "3" = 17),
                       labels = c("0" = "0 (P0)", "3" = "3 (P3)")) +
    ylab(ylab_name) +
    xlab("Temp") +
    facet_wrap(~ Pop, labeller = ggplot2::label_both) +
    clean_theme
  
  ggsave(file.path(out_dir, paste0(metric_col, "_interaction_Temp_vs_Diet_by_Pop_Phosphorus_linetype.png")),
         plot = p_line, width = 10, height = 7, dpi = 300)
  
  # -------------------------------------------------------------------------
  # ALTERNATIVA 2A: Grelha Duplicada 2x2 (Linhas = Phosphorus, Colunas = Pop)
  # -------------------------------------------------------------------------
  cat("  -> Saving Alternativa 2A (Grid 2x2)...\n")
  dodge_grid <- position_dodge(width = 0.5)
  p_grid <- ggplot(emm, aes(x = Temp, y = emmean, color = Diet, group = Diet)) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.3, position = dodge_grid, linewidth = 0.8) +
    geom_line(position = dodge_grid, linewidth = 1) +
    geom_point(position = dodge_grid, size = 3) +
    scale_color_manual(name = "Diet", values = diet_colors) +
    ylab(ylab_name) +
    xlab("Temp") +
    facet_grid(Phosphorus ~ Pop, labeller = ggplot2::label_both) +
    clean_theme
  
  ggsave(file.path(out_dir, paste0(metric_col, "_interaction_Temp_vs_Diet_by_Pop_Phosphorus_grid.png")),
         plot = p_grid, width = 10, height = 9, dpi = 300)
  
  # -------------------------------------------------------------------------
  # ALTERNATIVA 2B: Duas figuras separadas (uma para P0 e outra para P3)
  # -------------------------------------------------------------------------
  cat("  -> Saving Alternativa 2B (Separate P0 & P3)...\n")
  for (p_lvl in c("0", "3")) {
    sub_emm <- emm %>% filter(Phosphorus == p_lvl)
    p_single <- ggplot(sub_emm, aes(x = Temp, y = emmean, color = Diet, group = Diet)) +
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.3, position = dodge_grid, linewidth = 0.8) +
      geom_line(position = dodge_grid, linewidth = 1) +
      geom_point(position = dodge_grid, size = 3) +
      scale_color_manual(name = "Diet", values = diet_colors) +
      ylab(ylab_name) +
      xlab("Temp") +
      ggtitle(paste0("Phosphorus = ", p_lvl)) +
      facet_wrap(~ Pop, labeller = ggplot2::label_both) +
      clean_theme
    
    ggsave(file.path(out_dir, paste0(metric_col, "_interaction_Temp_vs_Diet_by_Pop_P", p_lvl, ".png")),
           plot = p_single, width = 10, height = 7, dpi = 300)
  }
}

generate_phosphorus_plots("observed_features", "Observed Features (ASVs)")
generate_phosphorus_plots("shannon_entropy", "Shannon Diversity")
generate_phosphorus_plots("faith_pd", "Faith's Phylogenetic Diversity")

cat("\nAll Phosphorus interaction figures generated successfully in:", out_dir, "\n")
