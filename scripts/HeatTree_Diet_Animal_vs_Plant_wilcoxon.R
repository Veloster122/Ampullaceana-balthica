#!/usr/bin/env Rscript
# =============================================================================
# HeatTree_Diet_Animal_vs_Plant_wilcoxon.R
# Differential abundance between Animal (A) vs Vegetal (P) diets
# using the original Wilcoxon-Mann-Whitney test in metacoder.
#
# Usage:
#   Rscript scripts/HeatTree_Diet_Animal_vs_Plant_wilcoxon.R \
#       [feature_table.qza] [taxonomy.qza] [metadata.tsv] [out_dir]
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

if (!require("metacoder", quietly = TRUE)) {
  stop("The 'metacoder' package is missing. Install with: conda install -c conda-forge -c bioconda r-metacoder r-taxa")
}

if (!require("qiime2R", quietly = TRUE)) {
  stop("The 'qiime2R' package is missing.")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 4) {
  feature_file <- args[1]
  tax_file     <- args[2]
  meta_file    <- args[3]
  out_dir      <- args[4]
} else {
  # Default paths within the repository
  feature_file <- if (file.exists("artifact_Ampullaceana/freq_taxo_tbl.qza")) {
    "artifact_Ampullaceana/freq_taxo_tbl.qza"
  } else {
    "artifact_Ampullaceana/Ampullaceana_freq_table.qza"
  }
  tax_file  <- "artifact_Ampullaceana/Ampullaceana_taxonomy.qza"
  meta_file <- "00-Ampullaceana_balthica_raw_data/Ampullaceana_balthica_metadata.tsv"
  out_dir   <- "Outputs"
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("===================================================================\n")
cat("Heat Tree: Diet Contrast — Animal (A) vs Vegetal/Plant (P)\n")
cat("Method: Wilcoxon-Mann-Whitney rank-sum test (p < 0.05)\n")
cat("Feature table:", feature_file, "\n")
cat("Taxonomy:     ", tax_file, "\n")
cat("Metadata:     ", meta_file, "\n")
cat("Output dir:   ", out_dir, "\n")
cat("===================================================================\n\n")

cat("Reading QIIME2 artifacts...\n")
features <- read_qza(feature_file)$data
taxonomy <- read_qza(tax_file)$data
metadata <- read_tsv(meta_file, show_col_types = FALSE)

if (metadata[1, 1] == "#q2:types") {
  metadata <- metadata[-1, ]
}
colnames(metadata)[1] <- "SampleID"

# Filter metadata to strictly Animal (A) and Vegetal/Plant (P) diets
metadata_diet <- metadata %>%
  filter(Diet %in% c("A", "P"))

cat(sprintf("Samples selected: %d (Animal A: %d | Vegetal P: %d)\n",
            nrow(metadata_diet),
            sum(metadata_diet$Diet == "A"),
            sum(metadata_diet$Diet == "P")))

cat("Converting to Taxmap object...\n")
tax_abund <- as.data.frame(features)
tax_abund$Feature.ID <- rownames(tax_abund)
tax_abund <- merge(taxonomy, tax_abund, by = "Feature.ID")

# Clean unclassified and filter ONLY Bacteria
tax_abund <- tax_abund[tax_abund$Taxon != "" & !is.na(tax_abund$Taxon), ]
tax_abund <- tax_abund[!grepl("^Unassigned", tax_abund$Taxon, ignore.case = TRUE), ]
tax_abund <- tax_abund[grepl("d__Bacteria", tax_abund$Taxon, ignore.case = TRUE), ]

# Parse taxonomy for metacoder
obj <- parse_tax_data(tax_abund,
                      class_cols = "Taxon",
                      class_sep = ";",
                      class_regex = "^([a-z0-9_]+)__?(.*)$",
                      class_key = c(tax_rank = "info", tax_name = "taxon_name"))

# Intersect samples
diet_samples <- intersect(metadata_diet$SampleID, colnames(obj$data$tax_data))
obj$data$tax_data <- obj$data$tax_data[, c("taxon_id", "Feature.ID", "Taxon", "Confidence", diet_samples)]
metadata_diet <- metadata_diet[metadata_diet$SampleID %in% diet_samples, ]

# Calculate relative abundance
obj$data$tax_abund <- calc_taxon_abund(obj, "tax_data", cols = diet_samples)
obj$data$tax_abund_prop <- calc_obs_props(obj, "tax_abund")

cat("Running Wilcoxon-Mann-Whitney test between Animal (A) and Vegetal (P)...\n")
diff_table <- compare_groups(obj,
                             data = "tax_abund_prop",
                             cols = diet_samples,
                             groups = metadata_diet$Diet[match(diet_samples, metadata_diet$SampleID)])

if ("treatment_1" %in% colnames(diff_table)) {
  diff_table <- diff_table %>% rename(treat1 = treatment_1, treat2 = treatment_2)
}

# Standardize column name for fold change
if (!"log2_fold_change" %in% colnames(diff_table) && "log2_median_ratio" %in% colnames(diff_table)) {
  diff_table <- diff_table %>% rename(log2_fold_change = log2_median_ratio)
}

# Filter for A vs P
diff_sub <- diff_table %>%
  filter((treat1 == "A" & treat2 == "P") | (treat1 == "P" & treat2 == "A"))

# Adjust sign so Animal (A) is treatment (numerator) and Vegetal (P) is reference (denominator)
diff_sub <- diff_sub %>%
  mutate(log2_fold_change = ifelse(treat1 == "P", -log2_fold_change, log2_fold_change))

# Mask non-significant taxa (p >= 0.05) or NAs to 0 (gray in Heat Tree)
pval_col <- if ("wilcox_p_value" %in% colnames(diff_sub)) "wilcox_p_value" else "p_value"
diff_sub <- diff_sub %>%
  mutate(is_sig = !is.na(.data[[pval_col]]) & .data[[pval_col]] < 0.05,
         log2_fold_change = ifelse(is_sig, log2_fold_change, 0))

diff_sub$log2_fold_change[is.na(diff_sub$log2_fold_change)] <- 0

# Align with all taxa in the tree to prevent NA node colors
all_taxa <- taxon_ids(obj)
match_idx <- match(all_taxa, as.character(diff_sub$taxon_id))
aligned_fc <- diff_sub$log2_fold_change[match_idx]
aligned_fc[is.na(aligned_fc)] <- 0
aligned_p  <- diff_sub[[pval_col]][match_idx]

obj$data$diff <- data.frame(
  taxon_id = all_taxa,
  log2_fold_change = aligned_fc,
  p_value = ifelse(is.na(aligned_p), 1, aligned_p),
  stringsAsFactors = FALSE
)

# Export detailed summary table
t_names <- taxon_names(obj)
t_ranks <- taxon_ranks(obj)
export_table <- diff_sub %>%
  mutate(taxon_name = t_names[as.character(taxon_id)],
         taxon_rank = t_ranks[as.character(taxon_id)],
         Enriched_in = case_when(
           is_sig & log2_fold_change > 0 ~ "Animal (A)",
           is_sig & log2_fold_change < 0 ~ "Vegetal (P)",
           TRUE ~ "Not Significant"
         )) %>%
  select(taxon_id, taxon_name, taxon_rank, log2_fold_change, all_of(pval_col), Enriched_in)

out_tsv <- file.path(out_dir, "HeatTree_Diet_Animal_vs_Plant_taxa.tsv")
write_tsv(export_table, out_tsv)
cat("Taxa statistical table saved to:", out_tsv, "\n")

n_animal <- sum(export_table$Enriched_in == "Animal (A)")
n_plant  <- sum(export_table$Enriched_in == "Vegetal (P)")
cat(sprintf("Significant taxa (Wilcoxon p < 0.05): %d enriched in Animal | %d enriched in Vegetal\n",
            n_animal, n_plant))

# Plot Heat Tree
cat("Rendering Heat Tree...\n")
set.seed(50)
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
               node_color_axis_label = "Log2 FC (Animal vs Vegetal)",
               layout = "fruchterman-reingold") +
  ggtitle("Heat Tree: Dieta (Animal vs Vegetal)") +
  labs(subtitle = "Vermelho: Mais abundante em Dieta Animal (A)  |  Verde: Mais abundante em Dieta Vegetal (P)  (Wilcoxon, p < 0.05)") +
  theme(plot.title    = element_text(size = 34, face = "bold",   hjust = 0.5, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 24, face = "italic", hjust = 0.5, margin = margin(b = 20)))

png_file <- file.path(out_dir, "HeatTree_Diet_Animal_vs_Plant.png")
pdf_file <- file.path(out_dir, "HeatTree_Diet_Animal_vs_Plant.pdf")

cat("Saving ultra-high resolution PNG (600 dpi)...\n")
ggsave(png_file, plot = p, width = 22, height = 22, dpi = 600)

cat("Saving vector PDF (infinite resolution)...\n")
ggsave(pdf_file, plot = p, width = 22, height = 22, dpi = 600)

cat("Successfully generated:\n")
cat("  -", png_file, "\n")
cat("  -", pdf_file, "\n")
cat("===================================================================\n")
