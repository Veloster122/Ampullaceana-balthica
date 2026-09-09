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
              "glmmTMB", "readr", "stringr")
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
# Robust parser for glmmTMB summary() coefficient table

parse_glmm_coef_table <- function(txt_file, metric_label) {
  fallback_df <- data.frame(
    Metric    = metric_label,
    Term      = "Not_available",
    Estimate  = NA_real_,
    Std_Error = NA_real_,
    z_value   = NA_real_,
    p_value   = NA_real_,
    stringsAsFactors = FALSE
  )

  if (!file.exists(txt_file)) {
    warning("File not found: ", txt_file)
    return(fallback_df)
  }

  lines <- readLines(txt_file, warn = FALSE)

  # Look for header line containing both Estimate and Std (or z/t value)
  header_idx <- which(grepl("Estimate", lines) & grepl("Std", lines))
  if (length(header_idx) == 0) {
    header_idx <- which(grepl("Estimate", lines) & (grepl("z value", lines) | grepl("Pr\\(", lines)))
  }

  if (length(header_idx) == 0) {
    warning("Could not find coefficient header in: ", txt_file)
    return(fallback_df)
  }

  # If multiple matches, take the last one (under fixed effects / conditional model)
  header_idx <- header_idx[length(header_idx)]

  # Collect data rows until blank line, separator, or signifs
  data_rows <- character(0)
  i <- header_idx + 1
  while (i <= length(lines)) {
    ln <- trimws(lines[i])
    if (ln == "" || grepl("^_{3,}", ln) || grepl("^---", ln) || grepl("^Signif\\. codes", ln)) {
      break
    }
    data_rows <- c(data_rows, ln)
    i <- i + 1
  }

  if (length(data_rows) == 0) {
    warning("No data rows found in coefficient table for: ", txt_file)
    return(fallback_df)
  }

  parsed_list <- lapply(data_rows, function(row_str) {
    parts <- strsplit(row_str, "\\s+")[[1]]
    if (length(parts) >= 5) {
      term <- parts[1]
      est  <- suppressWarnings(as.numeric(parts[2]))
      se   <- suppressWarnings(as.numeric(parts[3]))
      zval <- suppressWarnings(as.numeric(parts[4]))
      pval_str <- gsub("<", "", parts[5])
      pval <- suppressWarnings(as.numeric(pval_str))
      return(data.frame(
        Metric    = metric_label,
        Term      = term,
        Estimate  = est,
        Std_Error = se,
        z_value   = zval,
        p_value   = pval,
        stringsAsFactors = FALSE
      ))
    }
    return(NULL)
  })

  res_df <- do.call(rbind, parsed_list)
  if (is.null(res_df) || nrow(res_df) == 0) {
    return(fallback_df)
  }

  res_df
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
  m     <- metrics[[tbl_id]]
  fpath <- file.path(glmm_dir, m$file)
  out_csv <- file.path(out_dir, sprintf("Table_%s_GLMM_coefficients.csv", tbl_id))

  cat(sprintf("  Parsing %s → Table %s...\n", m$file, tbl_id))
  coef_tbl <- parse_glmm_coef_table(fpath, m$label)

  write_csv(coef_tbl, out_csv)
  cat(sprintf("    Saved: %s (%d rows)\n", out_csv, nrow(coef_tbl)))
}

# =============================================================================
# FIGURE B1 — Combined emmeans interaction plots manifest
# =============================================================================
cat("\n  Recording Figure B1 manifest (interaction plots)...\n")

interaction_files <- c(
  faith_pd          = file.path(glmm_dir, "faith_pd_interaction_Pop_vs_Temp_by_Diet.png"),
  observed_features = file.path(glmm_dir, "observed_features_interaction_Pop_vs_Temp_by_Diet.png"),
  shannon_entropy   = file.path(glmm_dir, "shannon_entropy_interaction_Pop_vs_Temp_by_Diet.png")
)

fig_b1_manifest <- data.frame(
  Panel   = c("B1a", "B1b", "B1c"),
  Metric  = c("Faith's PD", "Observed Features", "Shannon Entropy"),
  File    = unname(interaction_files),
  Exists  = file.exists(unname(interaction_files)),
  stringsAsFactors = FALSE
)
write_csv(fig_b1_manifest, file.path(out_dir, "FigureB1_file_manifest.csv"))
cat("    Saved: FigureB1_file_manifest.csv\n")
for (i in seq_len(nrow(fig_b1_manifest))) {
  status <- if (fig_b1_manifest$Exists[i]) "✓" else "✗ MISSING"
  cat(sprintf("      [%s] %s — %s\n", status, fig_b1_manifest$Panel[i], fig_b1_manifest$File[i]))
}

# =============================================================================
# APPENDIX C — Table C1: PCA importance table (formatted CSV)
# =============================================================================
cat("\n  Formatting Table C1 (PCA importance)...\n")

pca_summary_file <- file.path(glmm_dir, "PCA_Traits_summary.txt")
out_c1_csv <- file.path(out_dir, "Table_C1_PCA_importance.csv")

c1_saved <- FALSE
if (file.exists(pca_summary_file)) {
  lines <- readLines(pca_summary_file, warn = FALSE)
  imp_start <- grep("Importance of comp", lines)[1]
  if (!is.na(imp_start)) {
    imp_lines <- lines[(imp_start + 1):(imp_start + 4)]
    imp_tbl <- tryCatch({
      con <- textConnection(paste(imp_lines, collapse = "\n"))
      df  <- read.table(con, header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
      close(con)
      df
    }, error = function(e) NULL)

    if (!is.null(imp_tbl) && nrow(imp_tbl) > 0) {
      colnames(imp_tbl)[1] <- "Statistic"
      write_csv(imp_tbl, out_c1_csv)
      cat("    Saved: Table_C1_PCA_importance.csv\n")
      c1_saved <- TRUE
    }
  }

  # Also export loadings if present
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
      write_csv(load_tbl, file.path(out_dir, "Table_C1_PCA_loadings.csv"))
      cat("    Saved: Table_C1_PCA_loadings.csv\n")
    }
  }
}

if (!c1_saved) {
  fallback_c1 <- data.frame(
    Statistic = c("Standard deviation", "Proportion of Variance", "Cumulative Proportion"),
    PC1 = NA_real_, PC2 = NA_real_,
    stringsAsFactors = FALSE
  )
  write_csv(fallback_c1, out_c1_csv)
  cat("    Saved fallback: Table_C1_PCA_importance.csv\n")
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
                                              "PCA_Traits_biplot_Pop.png"))),
  stringsAsFactors = FALSE
)
write_csv(pca_files, file.path(out_dir, "FigureC1C2_file_manifest.csv"))
for (i in seq_len(nrow(pca_files))) {
  status <- if (pca_files$Exists[i]) "✓" else "✗ MISSING"
  cat(sprintf("    [%s] %s — %s\n", status, pca_files$Figure[i], pca_files$File[i]))
}

cat("\n=== Done! All appendix tables and manifests written to:", out_dir, "===\n")
