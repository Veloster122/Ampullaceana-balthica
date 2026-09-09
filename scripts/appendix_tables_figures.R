#!/usr/bin/env Rscript
# =============================================================================
# appendix_tables_figures.R
# Generates formatted tables and combined figures for Appendices B and C
# Reads outputs already produced by GLMMs_automatized.R
# =============================================================================
# Usage:
#   Rscript scripts/appendix_tables_figures.R \
#       glmm_outputs/ \
#       Outputs/Appendices/
# =============================================================================

packages <- c("dplyr", "tidyr", "ggplot2", "patchwork", "emmeans",
              "glmmTMB", "readr", "stringr", "broom.mixed", "knitr")
for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "http://cran.us.r-project.org")
  library(p, character.only = TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript appendix_tables_figures.R <glmm_out_dir> <appendix_out_dir>")
}
glmm_dir  <- args[1]
out_dir   <- args[2]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Appendix Tables & Figures Generator ===\n")
cat("GLMM outputs dir:", glmm_dir, "\n")
cat("Appendix output dir:", out_dir, "\n\n")

# =============================================================================
# APPENDIX B — Tables B1/B2/B3: GLMM coefficient tables
# =============================================================================
# Parse the *_summary_results.txt files to extract the fixed-effects
# coefficient table (Estimate, Std.Error, z value, Pr(>|z|))

parse_glmm_coef_table <- function(txt_file, metric_label) {
  lines <- readLines(txt_file)

  # Find the start of the "Fixed Effects:" section in glmmTMB summary output
  start_idx <- grep("^Fixed effects:", lines)
  if (length(start_idx) == 0) {
    # Try "Coefficients" header as fallback
    start_idx <- grep("Coefficients:", lines)
  }
  if (length(start_idx) == 0) {
    warning("Could not find coefficient table in: ", txt_file)
    return(NULL)
  }

  start_idx <- start_idx[1]

  # Find header row (contains "Estimate")
  header_idx <- start_idx + which(grepl("Estimate", lines[(start_idx + 1):length(lines)]))[1]

  # Find end of table (blank line or next section header)
  data_lines <- lines[(header_idx):length(lines)]
  end_offset <- which(data_lines == "" | grepl("^\\s*$", data_lines))[1]
  if (!is.na(end_offset)) {
    data_lines <- data_lines[1:(end_offset - 1)]
  }

  # Parse as table
  tbl <- tryCatch({
    con <- textConnection(paste(data_lines, collapse = "\n"))
    df  <- read.table(con, header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
    close(con)
    df
  }, error = function(e) {
    warning("Parsing failed for ", txt_file, ": ", e$message)
    NULL
  })

  if (!is.null(tbl)) {
    tbl$Term   <- rownames(tbl)
    tbl$Metric <- metric_label
    # Standardise column names
    names(tbl) <- gsub("Pr...z..", "p_value", names(tbl))
    names(tbl) <- gsub("Std..Error", "Std_Error", names(tbl))
    names(tbl) <- gsub("z.value", "z_value", names(tbl))
    tbl <- tbl %>% select(Metric, Term, Estimate, Std_Error, z_value, p_value)
  }
  tbl
}

metrics <- list(
  B1 = list(file = "faith_pd_summary_results.txt",
            label = "Faith's Phylogenetic Diversity (Faith PD)"),
  B2 = list(file = "observed_features_summary_results.txt",
            label = "Observed Features (ASV Richness)"),
  B3 = list(file = "shannon_entropy_summary_results.txt",
            label = "Shannon Entropy")
)

for (tbl_id in names(metrics)) {
  m    <- metrics[[tbl_id]]
  fpath <- file.path(glmm_dir, m$file)

  if (!file.exists(fpath)) {
    cat("  [SKIP] File not found:", fpath, "\n")
    next
  }

  cat(sprintf("  Parsing %s → Table %s...\n", m$file, tbl_id))
  coef_tbl <- parse_glmm_coef_table(fpath, m$label)

  if (!is.null(coef_tbl)) {
    out_csv <- file.path(out_dir, sprintf("Table_%s_GLMM_coefficients.csv", tbl_id))
    write_csv(coef_tbl, out_csv)
    cat(sprintf("    Saved: %s\n", out_csv))
  }
}

# =============================================================================
# FIGURE B1 — Combined emmeans interaction plots (all 3 metrics)
# Uses the already-generated PNGs + patchwork, or re-reads raw PNGs
# =============================================================================
cat("\n  Assembling Figure B1 (combined interaction plots)...\n")

# The pipeline already generates individual interaction PNGs.
# We record their paths here for reference in the appendix.
# If they all exist, we also note them in a manifest.

interaction_files <- c(
  faith_pd         = file.path(glmm_dir, "faith_pd_interaction_Pop_vs_Temp_by_Diet.png"),
  observed_features = file.path(glmm_dir, "observed_features_interaction_Pop_vs_Temp_by_Diet.png"),
  shannon_entropy  = file.path(glmm_dir, "shannon_entropy_interaction_Pop_vs_Temp_by_Diet.png")
)

fig_b1_manifest <- data.frame(
  Panel   = c("B1a", "B1b", "B1c"),
  Metric  = c("Faith's PD", "Observed Features", "Shannon Entropy"),
  File    = unname(interaction_files),
  Exists  = file.exists(unname(interaction_files))
)
write_csv(fig_b1_manifest, file.path(out_dir, "FigureB1_file_manifest.csv"))
cat("    Saved: FigureB1_file_manifest.csv\n")
cat("    Individual panels:\n")
for (i in seq_len(nrow(fig_b1_manifest))) {
  status <- if (fig_b1_manifest$Exists[i]) "✓" else "✗ MISSING"
  cat(sprintf("      [%s] %s — %s\n", status, fig_b1_manifest$Panel[i], fig_b1_manifest$File[i]))
}

# =============================================================================
# APPENDIX C — Table C1: PCA importance table (formatted CSV)
# =============================================================================
cat("\n  Formatting Table C1 (PCA importance)...\n")

pca_summary_file <- file.path(glmm_dir, "PCA_Traits_summary.txt")

if (file.exists(pca_summary_file)) {
  lines <- readLines(pca_summary_file)

  # Find "Importance of components" section
  imp_start <- grep("Importance of comp", lines)[1]
  if (!is.na(imp_start)) {
    imp_lines <- lines[(imp_start + 1):(imp_start + 4)]  # 3 rows: Std Dev, Prop Var, Cum Prop
    imp_tbl <- tryCatch({
      con <- textConnection(paste(imp_lines, collapse = "\n"))
      df  <- read.table(con, header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
      close(con)
      df
    }, error = function(e) NULL)

    if (!is.null(imp_tbl)) {
      colnames(imp_tbl)[1] <- "Statistic"
      out_csv <- file.path(out_dir, "Table_C1_PCA_importance.csv")
      write_csv(imp_tbl, out_csv)
      cat("    Saved: Table_C1_PCA_importance.csv\n")
    }
  }

  # Also export loadings
  load_start <- grep("PCA Loadings", lines)[1]
  if (!is.na(load_start)) {
    load_lines <- lines[(load_start + 1):min(load_start + 15, length(lines))]
    load_lines <- load_lines[load_lines != ""]
    load_tbl <- tryCatch({
      con <- textConnection(paste(load_lines, collapse = "\n"))
      df  <- read.table(con, header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
      close(con)
      df$Variable <- rownames(df)
      df
    }, error = function(e) NULL)

    if (!is.null(load_tbl)) {
      out_csv <- file.path(out_dir, "Table_C1_PCA_loadings.csv")
      write_csv(load_tbl, out_csv)
      cat("    Saved: Table_C1_PCA_loadings.csv\n")
    }
  }
} else {
  cat("  [SKIP] PCA summary file not found:", pca_summary_file, "\n")
}

# =============================================================================
# APPENDIX C — Figure C1 & C2: Record existing PCA biplot paths
# =============================================================================
cat("\n  Recording PCA biplot paths (Figure C1 & C2)...\n")

pca_files <- data.frame(
  Figure  = c("C1", "C2a", "C2b", "C2c"),
  Caption = c("PCA biplot coloured by experimental Box",
              "PCA biplot coloured by Diet",
              "PCA biplot coloured by Phosphorus",
              "PCA biplot coloured by Population"),
  File    = file.path(glmm_dir, c("PCA_Traits_biplot_Box.png",
                                   "PCA_Traits_biplot_Diet.png",
                                   "PCA_Traits_biplot_Phosphorus.png",
                                   "PCA_Traits_biplot_Pop.png")),
  Exists  = file.exists(file.path(glmm_dir, c("PCA_Traits_biplot_Box.png",
                                                "PCA_Traits_biplot_Diet.png",
                                                "PCA_Traits_biplot_Phosphorus.png",
                                                "PCA_Traits_biplot_Pop.png")))
)
write_csv(pca_files, file.path(out_dir, "FigureC1C2_file_manifest.csv"))
for (i in seq_len(nrow(pca_files))) {
  status <- if (pca_files$Exists[i]) "✓" else "✗ MISSING"
  cat(sprintf("    [%s] %s — %s\n", status, pca_files$Figure[i], pca_files$File[i]))
}

cat("\n=== Done! All appendix tables and manifests written to:", out_dir, "===\n")
