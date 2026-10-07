# scripts/HeatTrees_metacoder.R
# ============================================================
# Heat Trees (metacoder) with differential abundance tested by
# per-taxon GLMMs (default) or Wilcoxon-Mann-Whitney (legacy).
#
# GLMM approach (method = "glmm"):
#   For EVERY node of the taxonomic tree (domain -> ASV level), one GLMM:
#       abundance ~ Temp + Pop + Diet + Phosphorus + (1 | Box)
#   Several candidate distributions are fitted and the one with the best
#   residuals (DHARMa) is retained:
#     Counts (offset = log(total bacterial reads of the sample)):
#       - Poisson, Negative binomial (nbinom1, nbinom2), Zero-inflated NB2
#     Relative abundance (continuous):
#       - Tweedie (log link; handles exact zeros)
#       - Gamma   (log link; only fitted if the taxon has no zeros)
#       - Log-normal (Gaussian on log(prop + pseudocount))
#   Selection rule:
#     1. keep models that converged (positive-definite Hessian)
#     2. prefer models whose DHARMa tests (uniformity, dispersion,
#        zero-inflation) are all non-significant (p >= 0.05)
#     3. among those, lowest AIC (count models compared with each other;
#        continuous models compared on the proportion scale, with the
#        Jacobian correction for the log-normal). Count models are preferred
#        over continuous ones when both pass, because they use the raw data.
#     4. if no model passes DHARMa, the one with the best residuals is kept
#        and flagged (residuals_ok = FALSE) in the output tables.
#   Log2 fold change = model coefficient (log link) / log(2).
#   P-values are corrected with Benjamini-Hochberg (FDR) per comparison.
#
# Usage:
#   Rscript HeatTrees_metacoder.R <feature_table.qza> <taxonomy.qza> <metadata.tsv> <out_dir> [glmm|wilcoxon] [cores]
# ============================================================
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

if (!require("metacoder", quietly = TRUE)) {
  cat("\n=======================================================\n")
  cat("ERROR: The 'metacoder' package is missing.\n")
  cat("Please install it by running this in your terminal:\n")
  cat("conda install -c conda-forge -c bioconda r-metacoder r-taxa\n")
  cat("=======================================================\n")
  quit(status=1)
}

if (!require("qiime2R", quietly = TRUE)) {
  stop("The 'qiime2R' package is missing.")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript HeatTrees_metacoder.R <feature_table.qza> <taxonomy.qza> <metadata.tsv> <out_dir> [glmm|wilcoxon] [cores]")
}

feature_file <- args[1]
tax_file <- args[2]
meta_file <- args[3]
out_dir <- args[4]
method  <- if (length(args) >= 5) tolower(args[5]) else "glmm"
n_cores <- if (length(args) >= 6) as.integer(args[6]) else 1L
if (!method %in% c("glmm", "wilcoxon")) stop("method must be 'glmm' or 'wilcoxon'")

if (method == "glmm") {
  for (p in c("glmmTMB", "DHARMa")) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "http://cran.us.r-project.org")
  }
  suppressPackageStartupMessages({
    library(glmmTMB)
    library(DHARMa)
  })
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Reading QIIME2 artifacts...\n")
features <- read_qza(feature_file)$data
taxonomy <- read_qza(tax_file)$data
metadata <- read_tsv(meta_file, show_col_types = FALSE)

if (metadata[1, 1] == "#q2:types") {
  metadata <- metadata[-1, ]
}

# Rename the first column to 'SampleID' (in case it is named 'ID' or '#SampleID')
colnames(metadata)[1] <- "SampleID"

cat("Converting to Taxmap object...\n")
tax_abund <- as.data.frame(features)
tax_abund$Feature.ID <- rownames(tax_abund)
tax_abund <- merge(taxonomy, tax_abund, by = "Feature.ID")

# Remove unclassified or empty Taxon entries
tax_abund <- tax_abund[tax_abund$Taxon != "" & !is.na(tax_abund$Taxon), ]
tax_abund <- tax_abund[!grepl("^Unassigned", tax_abund$Taxon, ignore.case = TRUE), ]

# Filter ONLY Bacteria (remove Archaea, Eukaryota, and Host DNA)
tax_abund <- tax_abund[grepl("d__Bacteria", tax_abund$Taxon, ignore.case = TRUE), ]

# Parse taxonomy for metacoder
obj <- parse_tax_data(tax_abund,
                      class_cols = "Taxon",
                      class_sep = ";",
                      class_regex = "^([a-z0-9_]+)__?(.*)$",
                      class_key = c(tax_rank = "info", tax_name = "taxon_name"))

# Filter to include only samples in metadata
sample_cols <- intersect(metadata$SampleID, colnames(obj$data$tax_data))
obj$data$tax_data <- obj$data$tax_data[, c("taxon_id", "Feature.ID", "Taxon", "Confidence", sample_cols)]
metadata <- metadata[metadata$SampleID %in% sample_cols, ]

# Subsample or aggregate abundance
obj$data$tax_abund <- calc_taxon_abund(obj, "tax_data", cols = sample_cols)
obj$data$tax_abund_prop <- calc_obs_props(obj, "tax_abund")

# ============================================================
# Differential abundance: per-taxon GLMMs
# ============================================================
# Comparisons to extract from each model: column, reference level, treatment level
comparisons <- list(
  Temp       = c(ref = "14", treat = "20"),
  Pop        = c(ref = "SE", treat = "PT"),
  Phosphorus = c(ref = "0",  treat = "3")
)
# Fixed-effects structure (additive: a full 4-way interaction per taxon is over-parameterised)
glmm_rhs <- "Temp + Pop + Diet + Phosphorus"
# Minimum prevalence for a taxon to be modelled (others are left grey in the tree)
min_prev <- max(4, ceiling(0.10 * length(sample_cols)))

run_taxon_glmms <- function() {
  meta <- as.data.frame(metadata[match(sample_cols, metadata$SampleID), ])
  for (v in names(comparisons)) {
    meta[[v]] <- factor(as.character(meta[[v]]),
                        levels = c(comparisons[[v]][["ref"]],
                                   setdiff(unique(as.character(meta[[v]])), comparisons[[v]][["ref"]])))
  }
  meta$Diet <- factor(meta$Diet)
  meta$Box  <- factor(meta$Box)
  meta$lib  <- colSums(obj$data$tax_data[, sample_cols])   # total bacterial reads per sample

  counts <- as.data.frame(obj$data$tax_abund)
  # Pseudocount for log-normal: half of the smallest non-zero proportion in the dataset
  all_props <- as.matrix(counts[, sample_cols]) / matrix(meta$lib, nrow(counts), length(sample_cols), byrow = TRUE)
  pseudo <- min(all_props[all_props > 0]) / 2

  terms <- setNames(paste0(names(comparisons), sapply(comparisons, `[[`, "treat")), names(comparisons))
  fe    <- paste("~", glmm_rhs, "+ (1 | Box)")

  candidates <- list(
    poisson   = list(type = "count", resp = "y",      fam = poisson(),                 zi = ~0, offset = TRUE),
    nbinom1   = list(type = "count", resp = "y",      fam = nbinom1(),                 zi = ~0, offset = TRUE),
    nbinom2   = list(type = "count", resp = "y",      fam = nbinom2(),                 zi = ~0, offset = TRUE),
    zinbinom2 = list(type = "count", resp = "y",      fam = nbinom2(),                 zi = ~1, offset = TRUE),
    tweedie   = list(type = "cont",  resp = "prop",   fam = tweedie(link = "log"),     zi = ~0, offset = FALSE),
    gamma     = list(type = "cont",  resp = "prop",   fam = Gamma(link = "log"),       zi = ~0, offset = FALSE),
    lognormal = list(type = "cont",  resp = "logp",   fam = gaussian(),                zi = ~0, offset = FALSE)
  )

  fit_taxon <- function(i) {
    d <- meta
    d$y    <- as.numeric(counts[i, sample_cols])
    d$prop <- d$y / d$lib
    d$logp <- log(d$prop + pseudo)
    empty  <- data.frame(taxon_id = counts$taxon_id[i], family = NA, residuals_ok = NA,
                         dharma_unif_p = NA, dharma_disp_p = NA, dharma_zi_p = NA)
    for (v in names(terms)) { empty[[paste0(v, "_est")]] <- NA; empty[[paste0(v, "_p")]] <- NA }
    if (sum(d$y > 0) < min_prev) return(empty)

    # Step 1: fit every candidate distribution (cheap) and keep converged ones
    fits <- list()
    for (nm in names(candidates)) {
      cand <- candidates[[nm]]
      if (nm == "gamma" && any(d$y == 0)) next
      f <- paste(cand$resp, fe)
      if (cand$offset) f <- paste(f, "+ offset(log(lib))")
      m <- tryCatch(suppressWarnings(glmmTMB(as.formula(f), data = d, family = cand$fam, ziformula = cand$zi)),
                    error = function(e) NULL)
      if (is.null(m) || !isTRUE(m$sdr$pdHess)) next
      co <- tryCatch(summary(m)$coefficients$cond, error = function(e) NULL)
      if (is.null(co) || any(!terms %in% rownames(co)) || any(is.na(co[terms, "Pr(>|z|)"]))) next
      aic <- AIC(m)
      if (nm == "lognormal") aic <- aic + 2 * sum(d$logp)   # Jacobian: put on the proportion scale
      fits[[nm]] <- list(model = m, type = cand$type, aic = aic, co = co)
    }
    if (length(fits) == 0) return(empty)

    # Step 2: DHARMa residual checks (expensive) in order of preference:
    # count models by AIC, then continuous models by AIC. Stop at the first that passes.
    ord <- names(fits)[order(sapply(fits, `[[`, "type") != "count", sapply(fits, `[[`, "aic"))]
    best_name <- NULL; best_dh <- NULL; fallback <- NULL
    for (nm in ord) {
      res <- tryCatch(simulateResiduals(fits[[nm]]$model, n = 250, plot = FALSE), error = function(e) NULL)
      if (is.null(res)) next
      dh <- c(unif = tryCatch(testUniformity(res, plot = FALSE)$p.value, error = function(e) NA),
              disp = tryCatch(testDispersion(res, plot = FALSE)$p.value, error = function(e) NA),
              zi   = if (fits[[nm]]$type == "count")
                       tryCatch(testZeroInflation(res, plot = FALSE)$p.value, error = function(e) NA) else NA)
      worst <- suppressWarnings(min(dh, na.rm = TRUE))
      if (all(dh[!is.na(dh)] >= 0.05)) { best_name <- nm; best_dh <- dh; break }
      if (is.null(fallback) || worst > fallback$worst) fallback <- list(name = nm, dh = dh, worst = worst)
    }
    ok <- !is.null(best_name)
    if (!ok) {
      if (is.null(fallback)) return(empty)
      best_name <- fallback$name; best_dh <- fallback$dh     # best residuals, flagged
    }
    best <- fits[[best_name]]

    out <- data.frame(taxon_id = counts$taxon_id[i], family = best_name, residuals_ok = ok,
                      dharma_unif_p = best_dh[["unif"]], dharma_disp_p = best_dh[["disp"]],
                      dharma_zi_p = best_dh[["zi"]])
    for (v in names(terms)) {
      out[[paste0(v, "_est")]] <- best$co[terms[[v]], "Estimate"]
      out[[paste0(v, "_p")]]   <- best$co[terms[[v]], "Pr(>|z|)"]
    }
    out
  }

  # Robust CSV output written progressively by workers: immune to fork memory crashes
  results_tsv <- file.path(out_dir, "GLMM_raw_taxa_results.tsv")
  res_header <- c("taxon_id", "family", "residuals_ok", "dharma_unif_p", "dharma_disp_p", "dharma_zi_p",
                  paste0(rep(names(terms), each=2), c("_est", "_p")))

  # Check if models were already computed in a previous run (smart resume)
  already_computed <- file.exists(results_tsv) && length(readLines(results_tsv, warn = FALSE)) >= nrow(counts)

  if (already_computed) {
    cat("Found existing complete GLMM results at:", results_tsv, "\n")
    cat("Skipping recalculation and proceeding directly to aggregation and Heat Trees!\n")
  } else {
    writeLines(paste(res_header, collapse="\t"), results_tsv)

    # Progress log
    progress_file <- file.path(out_dir, "GLMM_taxa_progress.log")
    cat("", file = progress_file)

    fit_and_write_taxon <- function(i) {
      out <- tryCatch(fit_taxon(i), error = function(e) NULL)
      
      if (is.null(out) || is.na(out$family)) {
        cat(sprintf("%s\t%d\tskipped\n", format(Sys.time(), "%H:%M:%S"), i),
            file = progress_file, append = TRUE)
        empty_row <- c(counts$taxon_id[i], rep("NA", length(res_header) - 1))
        cat(paste(empty_row, collapse="\t"), "\n", file = results_tsv, append = TRUE)
      } else {
        cat(sprintf("%s\t%d\t%s\n", format(Sys.time(), "%H:%M:%S"), i, out$family),
            file = progress_file, append = TRUE)
        row_vals <- sapply(res_header, function(col) if (is.null(out[[col]])) "NA" else as.character(out[[col]]))
        cat(paste(row_vals, collapse="\t"), "\n", file = results_tsv, append = TRUE)
      }
      NULL
    }

    n_modelled <- sum(rowSums(as.matrix(counts[, sample_cols]) > 0) >= min_prev)
    cat("Fitting per-taxon GLMMs for", nrow(counts), "taxa;", n_modelled,
        "pass the prevalence filter (>=", min_prev, "samples) ...\n")
    cat("Follow progress with:  wc -l", progress_file, "\n")
    idx <- seq_len(nrow(counts))

    if (n_cores > 1 && .Platform$OS.type == "unix") {
      parallel::mclapply(idx, fit_and_write_taxon, mc.cores = n_cores, mc.preschedule = FALSE)
    } else {
      lapply(idx, function(i) { if (i %% 50 == 0) cat("  ", i, "/", length(idx), "\n"); fit_and_write_taxon(i) })
    }
  }

  cat("\nAggregating results from disk...\n")
  res <- read_tsv(results_tsv, show_col_types = FALSE)

  # Deduplicate in case of parallel writes
  res <- res %>% distinct(taxon_id, .keep_all = TRUE)

  # Taxon names / ranks for the output tables
  cat("Extracting taxonomy annotations...\n")
  t_names <- taxon_names(obj)
  t_ranks <- obj$taxon_ranks()
  res$taxon_name <- t_names[as.character(res$taxon_id)]
  res$taxon_rank <- t_ranks[as.character(res$taxon_id)]

  # Convert columns to numeric
  num_cols <- c("dharma_unif_p", "dharma_disp_p", "dharma_zi_p",
                paste0(rep(names(terms), each=2), c("_est", "_p")))
  for (col in num_cols) {
    res[[col]] <- as.numeric(res[[col]])
  }

  # Log2 FC and BH-FDR per comparison
  for (v in names(comparisons)) {
    res[[paste0(v, "_log2FC")]] <- res[[paste0(v, "_est")]] / log(2)
    res[[paste0(v, "_p_adj")]]  <- p.adjust(res[[paste0(v, "_p")]], method = "BH")
  }

  cat("Writing GLMM differential abundance tables...\n")
  write_csv(res, file.path(out_dir, "GLMM_taxa_differential_abundance.csv"))

  # Summary of distribution choice and residual quality
  summary_file <- file.path(out_dir, "GLMM_taxa_summary.txt")
  sink(summary_file)
  cat("=== Per-taxon GLMM differential abundance ===\n")
  cat("Model: abundance ~", glmm_rhs, "+ (1 | Box)\n")
  cat("Taxa in tree:", nrow(res), " | modelled:", sum(!is.na(res$family)),
      " | skipped (prevalence <", min_prev, "or no convergence):", sum(is.na(res$family)), "\n\n")
  cat("Selected distribution per taxon:\n"); print(table(res$family, useNA = "ifany"))
  cat("\nTaxa whose best model still failed DHARMa tests:", sum(res$residuals_ok == FALSE, na.rm = TRUE), "\n\n")
  for (v in names(comparisons)) {
    cat(sprintf("%s (%s vs %s): %d taxa with FDR < 0.05 (%d raw p < 0.05)\n", v,
                comparisons[[v]][["treat"]], comparisons[[v]][["ref"]],
                sum(res[[paste0(v, "_p_adj")]] < 0.05, na.rm = TRUE),
                sum(res[[paste0(v, "_p")]] < 0.05, na.rm = TRUE)))
  }
  sink()
  cat("Summary saved to:", summary_file, "\n")

  res
}

glmm_res <- if (method == "glmm") run_taxon_glmms() else NULL

# Custom Heat Tree function
# treat_label / ref_label: optional display strings for titles/axis (e.g. "20°C")
# treat_level / ref_level: actual values in the metadata column (e.g. "20")
make_heat_tree <- function(group_col, ref_level, treat_level, file_name,
                           ref_label = NULL, treat_label = NULL) {

  # Fall back to metadata values if no display labels provided
  if (is.null(ref_label))   ref_label   <- ref_level
  if (is.null(treat_label)) treat_label <- treat_level

  cat("Generating Heat Tree for", group_col, ":", treat_label, "vs", ref_label, "...\n")

  if (method == "glmm") {
    # Significant (FDR < 0.05) GLMM log2 fold changes; everything else grey (0)
    diff_sub <- data.frame(taxon_id = glmm_res$taxon_id,
                           log2_fold_change = glmm_res[[paste0(group_col, "_log2FC")]],
                           p_adj = glmm_res[[paste0(group_col, "_p_adj")]])
    diff_sub$log2_fold_change[is.na(diff_sub$p_adj) | diff_sub$p_adj >= 0.05] <- 0
    test_label <- "GLMM, FDR < 0.05"
  } else {
    # Remove NAs in the group column
    valid_samples <- metadata$SampleID[!is.na(metadata[[group_col]])]

    diff_table <- compare_groups(obj, data = "tax_abund_prop",
                                 cols = valid_samples,
                                 groups = metadata[[group_col]][metadata$SampleID %in% valid_samples])

    if ("treatment_1" %in% colnames(diff_table)) {
      diff_table <- diff_table %>% rename(treat1 = treatment_1, treat2 = treatment_2)
    }

    # Filter using actual metadata values (not display labels)
    diff_sub <- diff_table %>%
      filter((treat1 == treat_level & treat2 == ref_level) | (treat1 == ref_level & treat2 == treat_level))

    if (nrow(diff_sub) == 0) {
      cat("Warning: No differences found for", file_name, "\n")
      return()
    }

    # In newer metacoder versions, the fold change column is called log2_median_ratio
    if (!"log2_fold_change" %in% colnames(diff_sub) && "log2_median_ratio" %in% colnames(diff_sub)) {
      diff_sub <- diff_sub %>% rename(log2_fold_change = log2_median_ratio)
    }

    # Adjust sign so treat is always compared to ref
    diff_sub <- diff_sub %>%
      mutate(log2_fold_change = ifelse(treat1 == ref_level, -log2_fold_change, log2_fold_change))

    # Add significant filter (Wilcox p < 0.05) to mask non-significant changes
    if ("wilcox_p_value" %in% colnames(diff_sub)) {
      diff_sub <- diff_sub %>%
        mutate(log2_fold_change = ifelse(is.na(wilcox_p_value) | wilcox_p_value >= 0.05, 0, log2_fold_change))
    }
    test_label <- "Wilcoxon, p < 0.05"
  }

  # Merge back
  obj$data$diff <- diff_sub

  set.seed(50) # For reproducible layout — update if a different seed was used previously
  p <- heat_tree(obj,
                 node_label = taxon_names,
                 node_size = n_obs,
                 node_color = log2_fold_change,
                 node_label_size_range = c(0.0035, 0.035),
                 node_color_range = c("darkolivegreen4", "gray80", "firebrick"),
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_range = c("darkolivegreen4", "gray80", "firebrick"),
                 node_size_axis_label = "ASVs count",
                 node_color_axis_label = paste0("Log2 FC (", treat_label, " vs ", ref_label, ")"),
                 layout = "fruchterman-reingold") +
    ggtitle(paste0("Heat Tree: ", group_col, " (", treat_label, " vs ", ref_label, ")")) +
    labs(subtitle = paste0("Vermelho: Mais abundante em ", treat_label, "  |  Verde: Mais abundante em ", ref_label,
                           "  (", test_label, ")")) +
    theme(plot.title    = element_text(size = 34, face = "bold",   hjust = 0.5, margin = margin(b = 10)),
          plot.subtitle = element_text(size = 24, face = "italic", hjust = 0.5, margin = margin(b = 20)))

  # Save ultra-high resolution PNG
  ggsave(file.path(out_dir, file_name), plot = p, width = 22, height = 22, dpi = 600)

  # Save vector PDF (infinite resolution)
  pdf_name <- sub("\\.png$", ".pdf", file_name)
  ggsave(file.path(out_dir, pdf_name), plot = p, width = 22, height = 22, dpi = 600)
}

# Run trees
# treat_label / ref_label control display text; treat_level / ref_level must match metadata exactly
make_heat_tree("Temp", ref_level = "14", treat_level = "20",
               file_name = "HeatTree_Temp_20_vs_14.png",
               ref_label = "14°C", treat_label = "20°C")
make_heat_tree("Pop", "SE", "PT", "HeatTree_Pop_PT_vs_SE.png")
make_heat_tree("Phosphorus", "0", "3", "HeatTree_Phosphorus_3_vs_0.png")

cat("All Heat Trees generated successfully!\n")
