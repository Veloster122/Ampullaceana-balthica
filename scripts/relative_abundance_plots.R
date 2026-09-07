# scripts/relative_abundance_plots.R
# Generates stacked bar plots of the top 5 bacterial phyla and orders
# for each experimental factor: Diet, Temp, Phosphorus, Pop
#
# Usage (called by Snakemake rule):
#   Rscript scripts/relative_abundance_plots.R \
#       <feature_table.qza> <taxonomy.qza> <metadata.tsv> <out_dir>

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# ── Install qiime2R silently if missing ──────────────────────────────────────
if (!requireNamespace("qiime2R", quietly = TRUE)) {
  message("Installing qiime2R from GitHub...")
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  remotes::install_github("jbisanz/qiime2R", quiet = TRUE)
}
library(qiime2R)

# ── Parse command-line arguments ─────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript relative_abundance_plots.R <feature_table.qza> <taxonomy.qza> <metadata.tsv> <out_dir>")
}

feat_file <- args[1]
tax_file  <- args[2]
meta_file <- args[3]
out_dir   <- args[4]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load data ─────────────────────────────────────────────────────────────────
cat("Reading QIIME 2 artifacts...\n")
feat_mat <- as.data.frame(read_qza(feat_file)$data) # rows = ASVs, cols = samples
taxonomy  <- read_qza(tax_file)$data                 # Feature.ID + Taxon
metadata  <- read_tsv(meta_file, show_col_types = FALSE)

# Remove QIIME 2 types row if present
if (!is.na(metadata[1, 1]) && metadata[1, 1] == "#q2:types") {
  metadata <- metadata[-1, ]
}
colnames(metadata)[1] <- "SampleID"

# ── Parse taxonomy string into ranks ─────────────────────────────────────────
cat("Parsing taxonomy...\n")
tax_parsed <- taxonomy %>%
  separate(Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";\\s*", fill = "right", extra = "drop") %>%
  mutate(across(Domain:Species, ~gsub("^[a-z]__", "", .x)),
         across(Domain:Species, ~ifelse(.x == "" | is.na(.x), "Unclassified", .x)),
         # Keep only Bacteria; exclude host contamination
         Phylum = ifelse(Domain != "Bacteria", "Other", Phylum),
         Order  = ifelse(Domain != "Bacteria", "Other", Order)) %>%
  select(Feature.ID, Phylum, Order)

# ── Merge feature table with taxonomy ────────────────────────────────────────
feat_mat$Feature.ID <- rownames(feat_mat)
merged <- feat_mat %>%
  inner_join(tax_parsed, by = "Feature.ID") %>%
  filter(!Phylum %in% c("Unclassified", "", "Other", "Chloroplast", "Mitochondria"))

# ── Compute relative abundance per sample ────────────────────────────────────
sample_cols <- intersect(metadata$SampleID, colnames(merged))

# Normalise each sample to relative abundance (0–100 %)
abund_mat <- merged[, sample_cols] %>%
  mutate(across(everything(), as.numeric))

col_totals <- colSums(abund_mat, na.rm = TRUE)
col_totals[col_totals == 0] <- 1           # avoid division by zero
rel_mat <- sweep(abund_mat, 2, col_totals, FUN = "/") * 100

rel_mat$Phylum <- merged$Phylum
rel_mat$Order  <- merged$Order

# ── Helper: long table aggregated by taxon level ─────────────────────────────
# Add combined interaction columns to metadata (mirrors ancombc2_pop logic)
metadata <- metadata %>%
  mutate(
    Pop_x_Temp = paste(Pop, Temp, sep = "_"),
    Pop_x_Diet = paste(Pop, Diet, sep = "_")
  )

make_long <- function(df, rank_col) {
  df %>%
    group_by(.data[[rank_col]]) %>%
    summarise(across(all_of(sample_cols), sum, na.rm = TRUE), .groups = "drop") %>%
    pivot_longer(-all_of(rank_col), names_to = "SampleID", values_to = "RelAbund") %>%
    rename(Taxon = all_of(rank_col)) %>%
    left_join(metadata %>% select(SampleID, Diet, Temp, Phosphorus, Pop,
                                  Pop_x_Temp, Pop_x_Diet), by = "SampleID")
}

long_phylum <- make_long(rel_mat, "Phylum")
long_order  <- make_long(rel_mat, "Order")

# ── Helper: collapse to top 10 + "Other" ──────────────────────────────────────
top10_and_other <- function(long_df) {
  top10 <- long_df %>%
    group_by(Taxon) %>%
    summarise(TotalAbund = sum(RelAbund, na.rm = TRUE), .groups = "drop") %>%
    slice_max(TotalAbund, n = 10) %>%
    pull(Taxon)

  long_df %>%
    mutate(Taxon = ifelse(Taxon %in% top10, Taxon, "Other"))
}

phylum_top10 <- top10_and_other(long_phylum)
order_top10  <- top10_and_other(long_order)

# ── Colour palettes ───────────────────────────────────────────────────────────
build_palette <- function(taxa_vec, other_col = "grey80") {
  taxa <- setdiff(unique(taxa_vec), "Other")
  n    <- length(taxa)
  cols <- c(
    "#4E79A7", "#F28E2B", "#8B1A1A", "#76B7B2", "#B07AA1",
    "#59A14F", "#EDC948", "#FF9DA7", "#9C755F", "#BAB0AC"
  )[seq_len(n)]
  setNames(c(cols, other_col), c(taxa, "Other"))
}

phy_pal <- build_palette(phylum_top10$Taxon)
ord_pal <- build_palette(order_top10$Taxon)

# ── Helper: one stacked bar plot ─────────────────────────────────────────────
stacked_bar <- function(long_df, x_var, x_lab, fill_lab, palette, n_top = 10) {

  # Re-collapse top 10 within this specific factor
  top10_local <- long_df %>%
    group_by(Taxon) %>%
    summarise(tot = sum(RelAbund, na.rm = TRUE), .groups = "drop") %>%
    slice_max(tot, n = n_top) %>%
    pull(Taxon)

  plot_df <- long_df %>%
    mutate(Taxon = ifelse(Taxon %in% top10_local, Taxon, "Other")) %>%
    # Step 1: sum all "Other" taxa within each sample first
    group_by(.data[[x_var]], SampleID, Taxon) %>%
    summarise(RelAbund = sum(RelAbund, na.rm = TRUE), .groups = "drop") %>%
    # Step 2: then average across samples within each group
    group_by(.data[[x_var]], Taxon) %>%
    summarise(RelAbund = mean(RelAbund, na.rm = TRUE), .groups = "drop")

  # Order taxa by total mean abundance (descending), Other always last
  taxon_order <- plot_df %>%
    group_by(Taxon) %>%
    summarise(tot = sum(RelAbund), .groups = "drop") %>%
    filter(Taxon != "Other") %>%
    arrange(desc(tot)) %>%
    pull(Taxon)
  taxa_ordered <- c(taxon_order, "Other")

  # Factor levels: bottom-to-top for stacking (least abundant at bottom of bar = reversed)
  # Legend: most abundant at top → use guides(fill = guide_legend(reverse = FALSE))
  # We want bars stacked with most abundant at BOTTOM (visually prominent),
  # and legend showing most abundant at TOP — so factor levels ascending for stacking,
  # then reverse = TRUE in guide to flip legend only.
  plot_df$Taxon <- factor(plot_df$Taxon, levels = rev(taxa_ordered))

  local_pal <- build_palette(plot_df$Taxon)

  ggplot(plot_df, aes(x = as.factor(.data[[x_var]]),
                      y = RelAbund,
                      fill = Taxon)) +
    geom_bar(stat = "identity", colour = "white", linewidth = 0.25) +
    scale_fill_manual(values = local_pal, name = fill_lab) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
    labs(x = x_lab, y = "Relative Abundance (%)") +
    guides(fill = guide_legend(reverse = TRUE)) +   # most abundant at top of legend
    theme_classic(base_size = 12) +
    theme(
      legend.key.size  = unit(0.45, "cm"),
      legend.text      = element_text(size = 9),
      legend.title     = element_text(face = "bold", size = 10),
      axis.title       = element_text(face = "bold"),
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 11),
      panel.grid.major = element_blank()
    )
}

# ── Main loop: one figure per experimental factor ─────────────────────────────
factors <- list(
  Diet       = list(col = "Diet",       label = "Diet"),
  Temp       = list(col = "Temp",       label = "Temperature"),
  Phosphorus = list(col = "Phosphorus", label = "Phosphorus"),
  Pop        = list(col = "Pop",        label = "Population"),
  Pop_x_Temp = list(col = "Pop_x_Temp", label = "Population × Temperature"),
  Pop_x_Diet = list(col = "Pop_x_Diet", label = "Population × Diet")
)

for (factor_name in names(factors)) {
  fac   <- factors[[factor_name]]
  x_col <- fac$col
  x_lab <- fac$label

  cat(sprintf("  Generating plot for factor: %s\n", factor_name))

  # Panel A – Phyla
  pA <- stacked_bar(phylum_top10, x_col, x_lab, "Bacterial Phylum", phy_pal) +
    ggtitle(sprintf("Top 10 Bacterial Phyla Relative Abundance by %s", x_lab)) +
    labs(tag = "A")

  # Panel B – Orders
  pB <- stacked_bar(order_top10, x_col, x_lab, "Bacterial Order", ord_pal) +
    ggtitle(sprintf("Top 10 Bacterial Order Relative Abundance by %s", x_lab)) +
    labs(tag = "B")

  # Combine panels A + B with patchwork
  combined <- pA / pB +
    plot_layout(heights = c(1, 1))

  out_png <- file.path(out_dir, sprintf("RelativeAbundance_%s.png", factor_name))
  out_pdf <- file.path(out_dir, sprintf("RelativeAbundance_%s.pdf", factor_name))

  ggsave(out_png, combined, width = 8, height = 10, dpi = 300)
  ggsave(out_pdf, combined, width = 8, height = 10)

  cat(sprintf("  Saved: %s\n", out_png))
}

cat("All relative abundance plots generated successfully!\n")
