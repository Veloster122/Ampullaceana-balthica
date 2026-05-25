#!/usr/bin/env Rscript
# ============================================================
# Advanced Microbiome GLMMs Analysis
# Adapted from Sara & Julia (ce3c) standard script
# ============================================================
# ---- Auto-Install and Load Libraries ----
packages <- c("glmmTMB", "lmodel2", "car", "effects", "emmeans", 
              "DHARMa", "ggeffects", "ggplot2", "patchwork", "interactions", 
              "grid", "readr", "plyr", "jtools", "sjPlot", "dplyr", "rlang")
for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "http://cran.us.r-project.org")
  }
  library(p, character.only = TRUE)
}
# ---- Arguments ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript GLMMs_automatized.R <shannon_data> <observed_data> <faithpd_data> <metadata> <out_dir>")
}
shannon_file  <- args[1]
observed_file <- args[2]
faithpd_file  <- args[3]
metadata_file <- args[4]
out_dir       <- args[5]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# ---- Load Data ----
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
# ---- Factor and Numeric Conversions ----
df <- df %>%
  mutate(
    Temp = as.factor(Temp),
    Diet = as.factor(Diet),
    Phosphorus = as.factor(Phosphorus),
    Pop = as.factor(Pop),
    Box = as.factor(Box)
  )
# ---- PCA for Size Variables (Bonus Quest) ----
# Supervisor's request: Combine Ini_weight, Shell_IL, and Shell_IAA into a single Size_PC1
size_vars <- c("Ini_weight", "Shell_IL", "Shell_IAA")
for (var in size_vars) {
  df[[var]] <- as.numeric(df[[var]])
  df[[var]][is.na(df[[var]])] <- mean(df[[var]], na.rm = TRUE) # Impute NAs to prevent PCA failure
}
pca_res <- prcomp(df[, size_vars], center = TRUE, scale. = TRUE)
df$Size_PC1 <- pca_res$x[, 1]
# Save PCA summary to file
pca_log <- file.path(out_dir, "PCA_Size_summary.txt")
sink(pca_log)
cat("=== PCA Summary for Size Variables ===\n")
cat("Variables included: Ini_weight, Shell_IL, Shell_IAA\n\n")
print(summary(pca_res))
cat("\n=== PCA Loadings (Eigenvectors) ===\n")
print(pca_res$rotation)
sink()
# Save PCA biplot (Advanced ggplot2 version using Binned Metadata)
binned_meta_path <- file.path(dirname(metadata_file), "Ampullaceana_balthica_metadata_binned.tsv")
if (file.exists(binned_meta_path)) {
  binned_df <- read.table(binned_meta_path, sep="\t", header=TRUE, check.names=FALSE)
  if (binned_df[1,1] == "#q2:types") binned_df <- binned_df[-1,]
  
  pca_data <- data.frame(SampleID = df$SampleID, PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2])
  pca_data <- merge(pca_data, binned_df, by.x="SampleID", by.y=names(binned_df)[1], all.x=TRUE)
  
  loadings <- as.data.frame(pca_res$rotation)
  loadings$var <- rownames(loadings)
  # Scale arrows to fit plot
  mult <- min(max(pca_data$PC1, na.rm=T)/max(abs(loadings$PC1)), max(pca_data$PC2, na.rm=T)/max(abs(loadings$PC2))) * 0.8
  
  pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = Pop), alpha = 0.8, size = 3) +
    geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1 * mult, yend = PC2 * mult), 
                 arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 1) +
    geom_text(data = loadings, aes(x = PC1 * mult * 1.15, y = PC2 * mult * 1.15, label = var), 
              color = "darkred", size = 5, fontface="bold") +
    scale_color_manual(values = c("PT" = "#F8766D", "SE" = "#00BFC4")) +
    theme_bw() +
    labs(title = "PCA Biplot of Initial Size Variables", x = "PC1", y = "PC2", color = "Population") +
    theme(panel.grid.minor = element_blank())
  
  ggsave(file.path(out_dir, "PCA_Size_biplot.png"), plot = pca_plot, width = 8, height = 6, dpi = 300)
}
# ---- Helper Functions ----
main_effects1 <- function(model, a, metric_name, log_file) {
  fml <- reformulate(a)
  effects.a <- emmeans(model, fml)
  
  cat("\n--- Main Effect:", a, "---\n", file=log_file, append=TRUE)
  sink(log_file, append=TRUE)
  print(effects.a)
  cat("------------------------------------\n")
  print(test(pairs(effects.a), joint = TRUE))
  cat("------------------------------------\n")
  print(pairs(effects.a, adjust='bonferroni'))
  sink()
}
main_effects2 <- function(model, a, b, metric_name, log_file) {
  fml <- as.formula(paste("~", a, "|", b))
  effects.a.b <- emmeans(model, fml)
  
  cat("\n--- Effect of", a, "at each level of", b, "---\n", file=log_file, append=TRUE)
  sink(log_file, append=TRUE)
  print(effects.a.b)
  cat("------------------------------------\n")
  print(test(pairs(effects.a.b), joint =TRUE))
  cat("------------------------------------\n")
  print(pairs(effects.a.b, adjust='bonferroni'))
  sink()
}
create_cat_plot <- function(model, pred_var, modx_var, mod2_var, ylab, metric_name, data, png_path) {
  p <- interactions::cat_plot(
    model = model,
    data = data,
    pred = !!sym(pred_var),
    modx = !!sym(modx_var),
    mod2 = !!sym(mod2_var),
    geom = "line",
    size = 1,
    error.width = 1,
    dodge.width = 0.5,
    panel = TRUE
  ) +
    ylab(ylab) +
    xlab(pred_var) +
    scale_color_manual(name = modx_var, values = c("#740000", "#DAA520", "#4F734E", "#112233", "#445566")) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "italic", size = 12)
    ) +
    facet_wrap(as.formula(paste("~", mod2_var)), labeller = ggplot2::label_both)
  
  ggsave(filename = png_path, plot = p, width = 10, height = 7, dpi=300)
}
# ---- Core Runner ----
run_analysis <- function(metric_col, ylab_name) {
  print(paste("Running GLMM for", metric_col))
  
  # Remove NAs
  sub_df <- df[!is.na(df[[metric_col]]), ]
  
  # Ensure target is numeric
  sub_df[[metric_col]] <- as.numeric(sub_df[[metric_col]])
  
  # Log file
  log_file <- file.path(out_dir, paste0(metric_col, "_summary_results.txt"))
  cat("GLMM Analysis for:", metric_col, "\n", file=log_file)
  
  # 1. Build Model: Pop * Temp * Diet * Phosphorus + Size_PC1 + (1|Box)
  formula_str <- paste0(metric_col, " ~ Pop * Temp * Diet * Phosphorus + Size_PC1 + (1|Box)")
  best_model <- glmmTMB(as.formula(formula_str), data=sub_df)
  
  # 2. ANOVA & Summary
  sink(log_file, append=TRUE)
  cat("\n____________________________________________________\nANOVA Results (Type II):\n")
  print(Anova(best_model, type = 2))
  cat("\n____________________________________________________\nGLMM Coefficients Summary:\n")
  print(summary(best_model))
  sink()
  
  # 3. DHARMa Residuals
  png(filename = file.path(out_dir, paste0(metric_col, "_DHARMa_residuals.png")), width = 1000, height = 800)
  par(mfrow = c(2, 2))
  residuals_gr <- simulateResiduals(fittedModel = best_model, quantreg = TRUE)
  # Save residual tests to log
  sink(log_file, append=TRUE)
  cat("\n____________________________________________________\nResiduals Test:\n")
  print(testResiduals(residuals_gr, plot=FALSE))
  sink()
  # Plot residuals vs factors
  plotResiduals(residuals_gr, form = sub_df$Temp, main = "Residuals vs Temp") 
  plotResiduals(residuals_gr, form = sub_df$Diet, main = "Residuals vs Diet") 
  plotResiduals(residuals_gr, form = sub_df$Pop, main = "Residuals vs Pop") 
  plotResiduals(residuals_gr, form = sub_df$Phosphorus, main = "Residuals vs Phosphorus")
  plotResiduals(residuals_gr, form = sub_df$Size_PC1, main = "Residuals vs Size_PC1")
  dev.off()
  
  # 4. EMMEANS Post-Hoc Comparisons
  # Main effects
  main_effects1(best_model, "Pop", metric_col, log_file)
  main_effects1(best_model, "Temp", metric_col, log_file)
  main_effects1(best_model, "Diet", metric_col, log_file)
  
  # Interactions
  main_effects2(best_model, "Diet", "Temp", metric_col, log_file)
  main_effects2(best_model, "Diet", "Pop", metric_col, log_file)
  main_effects2(best_model, "Temp", "Pop", metric_col, log_file)
  
  # 5. Interaction Plots (cat_plot)
  # Pop vs Temp separated by Diet
  png1_path <- file.path(out_dir, paste0(metric_col, "_interaction_Pop_vs_Temp_by_Diet.png"))
  create_cat_plot(best_model, "Pop", "Temp", "Diet", ylab_name, metric_col, sub_df, png1_path)
  
  # Temp vs Pop separated by Diet
  png2_path <- file.path(out_dir, paste0(metric_col, "_interaction_Temp_vs_Diet_by_Pop.png"))
  create_cat_plot(best_model, "Temp", "Diet", "Pop", ylab_name, metric_col, sub_df, png2_path)
  
  print(paste("Done for", metric_col))
}
# ---- Execute Pipeline ----
run_analysis("shannon_entropy", "Shannon Diversity")
run_analysis("observed_features", "Observed Features (ASVs)")
run_analysis("faith_pd", "Faith's Phylogenetic Diversity")
# 5. Plotting
plot_metric <- function(metric_name, title, output_file) {
  p <- ggplot(df, aes_string(x="Temp", y=metric_name, fill="Diet")) +
    geom_boxplot(outlier.shape = NA, alpha=0.7, width=0.7, color="black") +
    geom_jitter(aes(color=Phosphorus, group=Diet), position=position_jitterdodge(jitter.width=0.15, dodge.width=0.7), 
                alpha=0.8, size=1.5) +
    facet_wrap(~Pop) +
    scale_fill_manual(values = c("A" = "#F8766D", "M" = "#619CFF", "P" = "#00BA38")) +
    scale_color_manual(values = c("0" = "grey40", "3" = "#9400D3"), labels = c("Ausente (0)", "Presente (3)")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill="grey95"),
          legend.position = "right",
          plot.title = element_text(hjust = 0, size = 14)) +
    labs(title=title, x="Temperature (°C)", y=metric_name, fill="Diet")
  
  ggsave(output_file, p, width=10, height=7, dpi=300)
}
plot_metric("shannon_entropy",   "Shannon Diversity by Temp, Diet, and Population",  file.path(out_dir, "Shannon_plot.png"))
plot_metric("observed_features", "Observed Features by Temp, Diet, and Population", file.path(out_dir, "Observed_plot.png"))
plot_metric("faith_pd",          "Faith's Phylogenetic Diversity by Temp, Diet, and Population", file.path(out_dir, "FaithPD_plot.png"))
print("All advanced GLMM analyses completed successfully.")
