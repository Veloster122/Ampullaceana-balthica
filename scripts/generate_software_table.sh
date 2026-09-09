#!/usr/bin/env bash
# =============================================================================
# generate_software_table.sh
# Generates Table D1: exact software versions used in the pipeline
# Run on the university server inside the conda environment:
#   conda activate qiime2-amplicon-2026.1
#   bash scripts/generate_software_table.sh > Outputs/Appendices/Table_D1_software_versions.tsv
# =============================================================================

set -euo pipefail

TSV_FILE="${1:-Outputs/Appendices/Table_D1_software_versions.tsv}"
mkdir -p "$(dirname "$TSV_FILE")"

# Header
echo -e "Software\tVersion\tSource\tRole" | tee "$TSV_FILE"

# ── Operating System ──────────────────────────────────────────────────────────
OS_VER=$(lsb_release -d 2>/dev/null | awk -F'\t' '{print $2}' || echo "Ubuntu (unknown)")
echo -e "Ubuntu\t${OS_VER}\tOS\tHost operating system" | tee -a "$TSV_FILE"

# ── Python ────────────────────────────────────────────────────────────────────
PY_VER=$(python --version 2>&1 | awk '{print $2}')
echo -e "Python\t${PY_VER}\tconda\tScripting language" | tee -a "$TSV_FILE"

# ── Snakemake ─────────────────────────────────────────────────────────────────
SM_VER=$(snakemake --version 2>/dev/null || echo "N/A")
echo -e "Snakemake\t${SM_VER}\tconda\tWorkflow manager" | tee -a "$TSV_FILE"

# ── QIIME2 ────────────────────────────────────────────────────────────────────
Q2_VER=$(python -c "import qiime2; print(qiime2.__version__)" 2>/dev/null || echo "N/A")
echo -e "QIIME2\t${Q2_VER}\tconda\tAmplicon microbiome analysis platform" | tee -a "$TSV_FILE"

# ── DADA2 (via QIIME2) ────────────────────────────────────────────────────────
DADA2_VER=$(python -c "
import subprocess, re
r = subprocess.run(['qiime', 'info'], capture_output=True, text=True)
m = re.search(r'dada2\s+([\d.]+)', r.stdout, re.IGNORECASE)
print(m.group(1) if m else 'N/A')
" 2>/dev/null || echo "N/A")
echo -e "DADA2\t${DADA2_VER}\tQIIME2 plugin\tASV denoising" | tee -a "$TSV_FILE"

# ── fastp ─────────────────────────────────────────────────────────────────────
FASTP_VER=$(fastp --version 2>&1 | head -1 | awk '{print $2}' || echo "N/A")
echo -e "fastp\t${FASTP_VER}\tconda\tRead quality control and adapter trimming" | tee -a "$TSV_FILE"

# ── cutadapt ─────────────────────────────────────────────────────────────────
CUTADAPT_VER=$(cutadapt --version 2>/dev/null || echo "N/A")
echo -e "cutadapt\t${CUTADAPT_VER}\tconda\tPrimer removal" | tee -a "$TSV_FILE"

# ── MAFFT ─────────────────────────────────────────────────────────────────────
MAFFT_VER=$(mafft --version 2>&1 | head -1 | awk '{print $2}' || echo "N/A")
echo -e "MAFFT\t${MAFFT_VER}\tQIIME2\tMultiple sequence alignment" | tee -a "$TSV_FILE"

# ── FastTree ──────────────────────────────────────────────────────────────────
FT_VER=$(FastTree 2>&1 | head -1 | grep -oP '(?<=version )\S+' || echo "N/A")
echo -e "FastTree\t${FT_VER}\tQIIME2\tPhylogenetic tree inference" | tee -a "$TSV_FILE"

# ── SILVA ─────────────────────────────────────────────────────────────────────
echo -e "SILVA\t138.2\tmanual download\t16S rRNA taxonomy reference database" | tee -a "$TSV_FILE"

# ── R packages (via Rscript) ──────────────────────────────────────────────────
Rscript - <<'REOF' | grep -v "^>" | tee -a "$TSV_FILE"
pkgs <- c(
  "glmmTMB"   , "Mixed models for alpha diversity",
  "emmeans"   , "Estimated marginal means and post-hoc contrasts",
  "DHARMa"    , "Residual diagnostics for GLMMs",
  "metacoder" , "Heat tree visualisation of taxonomic differences",
  "ANCOMBC"   , "Differential abundance (ANCOM-BC2)",
  "vegan"     , "Multivariate community ecology",
  "ggplot2"   , "Publication-quality graphics",
  "patchwork" , "Multi-panel figure composition",
  "qiime2R"   , "Import QIIME2 artifacts into R",
  "dplyr"     , "Data wrangling"
)
mat <- matrix(pkgs, ncol = 2, byrow = TRUE)
for (i in seq_len(nrow(mat))) {
  pkg  <- mat[i, 1]
  role <- mat[i, 2]
  ver  <- tryCatch(as.character(packageVersion(pkg)), error = function(e) "not installed")
  cat(paste(pkg, ver, "R CRAN/Bioconductor", role, sep = "\t"), "\n")
}
REOF

echo ""
echo "Table D1 written to: $TSV_FILE"
