#!/usr/bin/env Rscript
# =============================================================================
# generate_temp_phosphorus_plots.R
# Generate Interaction Plots for Temperatura : Fósforo across Alpha Diversity metrics:
#   - Shannon Entropy (Chisq = 7.21, p = 0.0072)
#   - Observed Features (Chisq = 6.58, p = 0.0103)
#   - Faith's Phylogenetic Diversity (Chisq = 6.87, p = 0.0088)
#
# Generates:
#   1. Standalone plots for each metric (PNG 300dpi & vector PDF)
#   2. Population-stratified plots (Portugal vs Suécia)
#   3. A combined 3-panel publication-ready figure (A, B, C)
# =============================================================================

suppressPackageStartupMessages({
  library(glmmTMB)
  library(emmeans)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 5) {
  shannon_file  <- args[1]
  observed_file <- args[2]
  faithpd_file  <- args[3]
  metadata_file <- args[4]
  out_dir       <- args[5]
} else {
  shannon_file  <- "glmm_inputs/shannon_metadata.tsv"
  observed_file <- "glmm_inputs/observed_features_metadata.tsv"
  faithpd_file  <- "glmm_inputs/faith_pd_metadata.tsv"
  metadata_file <- "00-Ampullaceana_balthica_raw_data/Ampullaceana_balthica_metadata.tsv"
  out_dir       <- "Outputs"
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading metadata and alpha diversity inputs...\n")
metadata <- read.table(metadata_file, sep="\t", header=TRUE, check.names=FALSE)
if (metadata[1,1] == "#q2:types") metadata <- metadata[-1,]

shannon_data  <- read.table(shannon_file, sep="\t", header=TRUE, check.names=FALSE)
observed_data <- read.table(observed_file, sep="\t", header=TRUE, check.names=FALSE)
faithpd_data  <- read.table(faithpd_file,  sep="\t", header=TRUE, check.names=FALSE)

colnames(shannon_data)[1]  <- "SampleID"
colnames(observed_data)[1] <- "SampleID"
colnames(faithpd_data)[1]  <- "SampleID"
colnames(metadata)[1]      <- "SampleID"

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

# Covariável Tamanho Inicial (PC1)
ini_vars <- c("Ini_weight", "Shell_IL", "Shell_IAA")
fin_vars <- c("Fin_weight", "Shell_FL", "Shell_FAA")
other_vars <- c("N_clutches_per day", "Cons_rate", "Growth_rate")
all_morpho_vars <- intersect(unique(c(ini_vars, fin_vars, other_vars)), colnames(df))

for (var in all_morpho_vars) {
  df[[var]] <- as.numeric(df[[var]])
  df[[var]][is.na(df[[var]])] <- mean(df[[var]], na.rm = TRUE)
}

pca_ini <- prcomp(df[, ini_vars], center = TRUE, scale. = TRUE)
df$Initial_Size <- pca_ini$x[, 1]
if (cor(df$Initial_Size, df$Ini_weight) < 0) df$Initial_Size <- -df$Initial_Size

df$Size_PC1 <- df$Initial_Size
formula_rhs <- "Pop * Temp * Diet * Phosphorus + Size_PC1 + (1|Box)"

clean_theme <- theme_bw() +
  theme(
    panel.grid.major = element_line(color = "grey92", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_blank(),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    legend.background = element_rect(color = "grey80", fill = "white", linewidth = 0.4),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 6)),
    plot.subtitle = element_text(size = 11, face = "italic", hjust = 0.5, margin = margin(b = 10))
  )

phos_colors <- c("0" = "#4A5568", "3" = "#9400D3")
phos_labels <- c("0" = "P0 (Ausente)", "3" = "P3 (+3% Fósforo)")
phos_shapes <- c("0" = 16, "3" = 17)
phos_lines  <- c("0" = "solid", "3" = "dashed")
dodge_val   <- position_dodge(width = 0.25)

generate_metric_plots <- function(metric_col, ylab_pt, title_pt, panel_letter, p_stat_str) {
  cat("Fitting GLMM for", metric_col, "...\n")
  sub_df <- df[!is.na(df[[metric_col]]), ]
  sub_df[[metric_col]] <- as.numeric(sub_df[[metric_col]])
  fml <- as.formula(paste0(metric_col, " ~ ", formula_rhs))
  model <- glmmTMB(fml, data = sub_df)
  
  emm <- as.data.frame(emmeans(model, ~ Temp * Phosphorus))
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
  
  emm$Temp_Label <- factor(ifelse(emm$Temp == "14", "14°C", "20°C"), levels = c("14°C", "20°C"))
  emm$Phosphorus <- factor(emm$Phosphorus, levels = c("0", "3"))
  
  # Standalone single plot
  p_standalone <- ggplot(emm, aes(x = Temp_Label, y = emmean, color = Phosphorus, group = Phosphorus,
                                  linetype = Phosphorus, shape = Phosphorus)) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.15, linewidth = 0.8, position = dodge_val) +
    geom_line(linewidth = 1.1, position = dodge_val) +
    geom_point(size = 3.8, position = dodge_val) +
    scale_color_manual(name = "Fósforo", values = phos_colors, labels = phos_labels) +
    scale_linetype_manual(name = "Fósforo", values = phos_lines, labels = phos_labels) +
    scale_shape_manual(name = "Fósforo", values = phos_shapes, labels = phos_labels) +
    labs(
      title = title_pt,
      subtitle = paste0("Interação Temperatura × Fósforo (GLMM ANOVA: ", p_stat_str, ")"),
      x = "Temperatura",
      y = ylab_pt
    ) +
    clean_theme
  
  # Plot adapted for 3-panel composite figure
  p_panel <- ggplot(emm, aes(x = Temp_Label, y = emmean, color = Phosphorus, group = Phosphorus,
                             linetype = Phosphorus, shape = Phosphorus)) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.15, linewidth = 0.8, position = dodge_val) +
    geom_line(linewidth = 1.1, position = dodge_val) +
    geom_point(size = 3.6, position = dodge_val) +
    scale_color_manual(name = "Fósforo", values = phos_colors, labels = phos_labels) +
    scale_linetype_manual(name = "Fósforo", values = phos_lines, labels = phos_labels) +
    scale_shape_manual(name = "Fósforo", values = phos_shapes, labels = phos_labels) +
    labs(
      title = paste0(panel_letter, "  ", title_pt),
      subtitle = p_stat_str,
      x = "Temperatura",
      y = ylab_pt
    ) +
    clean_theme +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 10, face = "italic", hjust = 0, color = "grey30")
    )
  
  # Population breakdown
  emm_pop <- as.data.frame(emmeans(model, ~ Temp * Phosphorus * Pop))
  if ("asymp.LCL" %in% colnames(emm_pop)) {
    emm_pop$ymin <- emm_pop$asymp.LCL
    emm_pop$ymax <- emm_pop$asymp.UCL
  } else if ("lower.CL" %in% colnames(emm_pop)) {
    emm_pop$ymin <- emm_pop$lower.CL
    emm_pop$ymax <- emm_pop$upper.CL
  } else {
    emm_pop$ymin <- emm_pop$emmean - 1.96 * emm_pop$SE
    emm_pop$ymax <- emm_pop$emmean + 1.96 * emm_pop$SE
  }
  emm_pop$Temp_Label <- factor(ifelse(emm_pop$Temp == "14", "14°C", "20°C"), levels = c("14°C", "20°C"))
  emm_pop$Phosphorus <- factor(emm_pop$Phosphorus, levels = c("0", "3"))
  emm_pop$Pop_Label <- factor(ifelse(emm_pop$Pop == "PT", "Portugal (PT)", "Suécia (SE)"),
                              levels = c("Portugal (PT)", "Suécia (SE)"))
  
  p_pop <- ggplot(emm_pop, aes(x = Temp_Label, y = emmean, color = Phosphorus, group = Phosphorus,
                               linetype = Phosphorus, shape = Phosphorus)) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.15, linewidth = 0.8, position = dodge_val) +
    geom_line(linewidth = 1.1, position = dodge_val) +
    geom_point(size = 3.8, position = dodge_val) +
    scale_color_manual(name = "Fósforo", values = phos_colors, labels = phos_labels) +
    scale_linetype_manual(name = "Fósforo", values = phos_lines, labels = phos_labels) +
    scale_shape_manual(name = "Fósforo", values = phos_shapes, labels = phos_labels) +
    labs(
      title = paste0(title_pt, " por População"),
      subtitle = "Médias marginais estimadas (emmeans ± 95% IC)",
      x = "Temperatura",
      y = ylab_pt
    ) +
    facet_wrap(~ Pop_Label) +
    clean_theme +
    theme(strip.background = element_rect(fill = "grey95", color = "black"),
          strip.text = element_text(face = "bold", size = 12))
  
  # Save individual
  out_png <- file.path(out_dir, paste0(metric_col, "_interaction_Temp_vs_Phosphorus.png"))
  out_pdf <- file.path(out_dir, paste0(metric_col, "_interaction_Temp_vs_Phosphorus.pdf"))
  ggsave(out_png, plot = p_standalone, width = 8, height = 6, dpi = 300)
  ggsave(out_pdf, plot = p_standalone, width = 8, height = 6)
  
  out_pop_png <- file.path(out_dir, paste0(metric_col, "_interaction_Temp_vs_Phosphorus_by_Pop.png"))
  out_pop_pdf <- file.path(out_dir, paste0(metric_col, "_interaction_Temp_vs_Phosphorus_by_Pop.pdf"))
  ggsave(out_pop_png, plot = p_pop, width = 10, height = 6, dpi = 300)
  ggsave(out_pop_pdf, plot = p_pop, width = 10, height = 6)
  
  return(list(standalone = p_standalone, panel = p_panel, pop = p_pop))
}

cat("Generating plots for all 3 metrics...\n")
res_shannon  <- generate_metric_plots("shannon_entropy",  "Índice de Shannon", "Entropia de Shannon", "(A)", "Chisq = 7.21, p = 0.0072")
res_observed <- generate_metric_plots("observed_features", "Riqueza de ASVs", "Features Observadas (ASVs)", "(B)", "Chisq = 6.58, p = 0.0103")
res_faith    <- generate_metric_plots("faith_pd",          "Diversidade de Faith (PD)", "Diversidade Filogenética (Faith)", "(C)", "Chisq = 6.87, p = 0.0088")

# Combine 3 panels with shared legend at bottom
combined_fig <- (res_shannon$panel + res_observed$panel + res_faith$panel) +
  plot_layout(ncol = 3, guides = "collect") +
  plot_annotation(
    title = "Efeito da Interação Temperatura × Fósforo na Diversidade Alfa",
    subtitle = "Médias marginais estimadas (emmeans) ± 95% IC a partir dos GLMMs (Pop * Temp * Diet * Phosphorus + Size_PC1 + (1|Box))"
  )

comb_png <- file.path(out_dir, "AlphaDiversity_interaction_Temp_vs_Phosphorus_combined.png")
comb_pdf <- file.path(out_dir, "AlphaDiversity_interaction_Temp_vs_Phosphorus_combined.pdf")

ggsave(comb_png, plot = combined_fig, width = 16, height = 6.5, dpi = 300)
ggsave(comb_pdf, plot = combined_fig, width = 16, height = 6.5)

cat("\n=======================================================\n")
cat("SUCCESS! All Temperatura:Fósforo plots generated:\n")
cat("  1. Individual plots:\n")
cat("     -", file.path(out_dir, "shannon_entropy_interaction_Temp_vs_Phosphorus.png"), "\n")
cat("     -", file.path(out_dir, "observed_features_interaction_Temp_vs_Phosphorus.png"), "\n")
cat("     -", file.path(out_dir, "faith_pd_interaction_Temp_vs_Phosphorus.png"), "\n")
cat("  2. Stratified by Population (Portugal vs Suécia):\n")
cat("     -", file.path(out_dir, "shannon_entropy_interaction_Temp_vs_Phosphorus_by_Pop.png"), "\n")
cat("     -", file.path(out_dir, "observed_features_interaction_Temp_vs_Phosphorus_by_Pop.png"), "\n")
cat("     -", file.path(out_dir, "faith_pd_interaction_Temp_vs_Phosphorus_by_Pop.png"), "\n")
cat("  3. 3-Panel publication composite:\n")
cat("     -", comb_png, "\n")
cat("     -", comb_pdf, "\n")
cat("=======================================================\n")
