#!/usr/bin/env Rscript
# scripts/GLMMs_automatized.R
# Generalized Linear Mixed Models for Ampullaceana balthica 16S diversity

# Load libraries
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# Get arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  stop("Usage: Rscript GLMMs_automatized.R <shannon_meta> <observed_meta> <faithpd_meta> <metadata> <shannon_plot> <observed_plot> <pvalues_tsv>")
}

shannon_file  <- args[1]
observed_file <- args[2]
faithpd_file  <- args[3]
metadata_file <- args[4]
shannon_plot  <- args[5]
observed_plot <- args[6]
pvalues_file  <- args[7]

# 1. Load data
metadata <- read.table(metadata_file, sep="\t", header=TRUE, check.names=FALSE)
# Remove QIIME2 header row if present (#q2:types)
if (metadata[1,1] == "#q2:types") {
  metadata <- metadata[-1,]
}

shannon_data  <- read.table(shannon_file, sep="\t", header=TRUE, check.names=FALSE)
observed_data <- read.table(observed_file, sep="\t", header=TRUE, check.names=FALSE)
faithpd_data  <- read.table(faithpd_file,  sep="\t", header=TRUE, check.names=FALSE)

# 2. Merge data
# Use the first column (Sample ID) for merging
colnames(shannon_data)[1]  <- "SampleID"
colnames(observed_data)[1] <- "SampleID"
colnames(faithpd_data)[1]  <- "SampleID"
colnames(metadata)[1]      <- "SampleID"

# Clean up data
shannon_data  <- shannon_data  %>% select(SampleID, shannon_entropy)
observed_data <- observed_data %>% select(SampleID, observed_features)
faithpd_data  <- faithpd_data  %>% select(SampleID, faith_pd)

df <- metadata %>%
  inner_join(shannon_data, by="SampleID") %>%
  inner_join(observed_data, by="SampleID") %>%
  inner_join(faithpd_data, by="SampleID")

# Convert factors
df$Temp <- as.factor(df$Temp)
df$Diet <- as.factor(df$Diet)
df$Pop  <- as.factor(df$Pop)
df$Box  <- as.factor(df$Box)
df$Phosphorus <- as.numeric(as.character(df$Phosphorus))

# 3. Model Function
run_glmm <- function(metric_name, data) {
  # Formula based on Ampullaceana balthica study design
  # metric ~ Temp * Diet + Phosphorus + Pop + (1|Box)
  formula_str <- paste(metric_name, "~ Temp * Diet + Phosphorus + Pop + (1|Box)")
  model <- lmer(as.formula(formula_str), data = data)
  
  # Get p-values
  summary_model <- summary(model)
  coefs <- as.data.frame(summary_model$coefficients)
  coefs$Variable <- rownames(coefs)
  coefs$Metric <- metric_name
  
  return(list(model=model, coefs=coefs))
}

# Run models
shannon_res  <- run_glmm("shannon_entropy", df)
observed_res <- run_glmm("observed_features", df)
faithpd_res  <- run_glmm("faith_pd", df)

# 4. FDR correction (Benjamini-Hochberg), applied per metric separately
adjust_fdr <- function(coefs) {
  p_col <- grep("Pr", colnames(coefs), value = TRUE)[1]
  coefs$p_fdr <- p.adjust(coefs[[p_col]], method = "fdr")
  return(coefs)
}

shannon_res$coefs  <- adjust_fdr(shannon_res$coefs)
observed_res$coefs <- adjust_fdr(observed_res$coefs)
faithpd_res$coefs  <- adjust_fdr(faithpd_res$coefs)

# 5. Save P-values
pvalues <- rbind(shannon_res$coefs, observed_res$coefs, faithpd_res$coefs)
write.table(pvalues, pvalues_file, sep="\t", row.names=FALSE, quote=FALSE)

# 5. Plotting
plot_metric <- function(metric_name, title, output_file) {
  p <- ggplot(df, aes_string(x="Temp", y=metric_name, fill="Diet")) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position=position_jitterdodge(), alpha=0.3, size=1) +
    facet_wrap(~Pop) +
    theme_bw() +
    labs(title=title, x="Temperature (°C)", y=metric_name) +
    scale_fill_brewer(palette="Set1")
  
  ggsave(output_file, p, width=8, height=6)
}

plot_metric("shannon_entropy",   "Shannon Diversity by Temp, Diet, and Population",  shannon_plot)
plot_metric("observed_features", "Observed Features by Temp, Diet, and Population", observed_plot)

message("GLMM analysis complete. Results saved to glmm_outputs/")
