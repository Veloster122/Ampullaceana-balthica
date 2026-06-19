# scripts/HeatTrees_metacoder.R
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
  stop("Usage: Rscript HeatTrees_metacoder.R <feature_table.qza> <taxonomy.qza> <metadata.tsv> <out_dir>")
}
feature_file <- args[1]
tax_file <- args[2]
meta_file <- args[3]
out_dir <- args[4]
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
# Custom Heat Tree function
make_heat_tree <- function(group_col, ref_level, treat_level, file_name) {
  cat("Generating Heat Tree for", group_col, ":", treat_level, "vs", ref_level, "...\n")
  
  # Remove NAs in the group column
  valid_samples <- metadata$SampleID[!is.na(metadata[[group_col]])]
  
  diff_table <- compare_groups(obj, data = "tax_abund_prop",
                               cols = valid_samples,
                               groups = metadata[[group_col]][metadata$SampleID %in% valid_samples])
  
  if ("treatment_1" %in% colnames(diff_table)) {
    diff_table <- diff_table %>% rename(treat1 = treatment_1, treat2 = treatment_2)
  }
  
  # Filter diff_table for the specific comparison
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
  
  # Remove taxa with 0 fold change to clean up tree (optional, keeping for structure)
  
  # Merge back
  obj$data$diff <- diff_sub 
  
  set.seed(50) # For reproducible layout
  p <- heat_tree(obj,
                 node_label = taxon_names,
                 node_size = n_obs,
                 node_color = log2_fold_change,
                 node_label_size_range = c(0.0035, 0.035), 
                 node_color_range = c("darkolivegreen4", "gray80", "firebrick"), # Cyan to Red
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_range = c("darkolivegreen4", "gray80", "firebrick"),
                 node_size_axis_label = "ASVs count",
                 node_color_axis_label = paste("Log2 FC (", treat_level, " vs ", ref_level, ")", sep=""),
                 layout = "fruchterman-reingold") +
    ggtitle(paste("Heat Tree: ", group_col, " (", treat_level, " vs ", ref_level, ")", sep="")) +
    labs(subtitle = paste("Vermelho: Mais abundante em", treat_level, "  |  Verde: Mais abundante em", ref_level)) +
    theme(plot.title = element_text(size = 34, face = "bold", hjust = 0.5, margin = margin(b = 10)),
          plot.subtitle = element_text(size = 24, face = "italic", hjust = 0.5, margin = margin(b = 20)))
  
  # Save ultra-high resolution PNG
  ggsave(file.path(out_dir, file_name), plot = p, width = 22, height = 22, dpi = 600)
  
  # Save vector PDF (infinite resolution)
  pdf_name <- sub("\\.png$", ".pdf", file_name)
  ggsave(file.path(out_dir, pdf_name), plot = p, width = 22, height = 22, dpi = 600)
}
# Run trees
make_heat_tree("Temp", "14", "20", "HeatTree_Temp_20_vs_14.png")
make_heat_tree("Pop", "SE", "PT", "HeatTree_Pop_PT_vs_SE.png")
make_heat_tree("Phosphorus", "0", "3", "HeatTree_Phosphorus_3_vs_0.png")
cat("All Heat Trees generated successfully!\n")
