<<<<<<< HEAD
# Ampullaceana balthica — 16S Gut Microbiome Snakemake Pipeline

> **Species:** *Ampullaceana balthica*  
> **Study:** Effect of temperature (14°C vs 20°C) on gut microbiome composition  
> **Populations:** Portugal (PT) and Sweden (SE)  
> **Primers:** 341F (`CCTACGGGNGGCWGCAG`) / 806R (`GACTACHVGGGTATCTAATCC`)

---

## Project Structure

```
Ampullaceana balthica/
├── config.yaml              ← Edit this file (and only this file) between runs
├── Snakefile_phase1         ← Phase 1: QC → Feature Tables → Taxonomy
├── Snakefile_phase2         ← Phase 2: Diversity → GLMMs
├── rules/
│   ├── fastp_import.smk     ← fastp QC + QIIME2 import
│   ├── qc_reports.smk       ← FastQC + MultiQC reports
│   ├── feature_table.smk    ← DADA2 denoising + feature tables
│   ├── taxonomy.smk         ← SILVA classifier + taxonomy + bar plots
│   ├── phylo_diversity.smk  ← Phylogeny + alpha/beta diversity + ANCOMBC2
│   └── glmm.smk             ← Shannon extraction + GLMM in R
├── 00-Ampullaceana_balthica_raw_data/
│   ├── fastqs/              ← Raw paired-end fastq.gz files
│   └── Ampullaceana_balthica_metadata.tsv
└── manifest.tsv             ← Already generated
```

---

## Requirements

```bash
conda activate qiime2-amplicon-2024.10
conda install -c bioconda snakemake fastp fastqc multiqc
```

---

## How to Run

### Phase 1 — From raw reads to taxonomy

**Step 1:** Open `config.yaml` and verify all paths are correct.

**Step 2:** Run Phase 1:
```bash
snakemake -s Snakefile_phase1 --configfile config.yaml --cores 24
```

**⚠️ CHECKPOINT 1:** Inspect the demux visualization:
```
visualizations_Ampullaceana/ampullaceana_demux.qzv
```
Open it at [view.qiime2.org](https://view.qiime2.org) → "Interactive Quality Plot" tab.  
Find where quality drops and update `config.yaml`:
```yaml
denoise:
  trunc_f: ??  # ← set your value here
  trunc_r: ??  # ← set your value here
```
Then re-run Phase 1 — Snakemake will skip already-done steps and only re-run DADA2 onwards.

---

### Phase 2 — Diversity and statistics

**Step 3:** After Phase 1 completes, inspect the taxonomy frequency table:
```
visualizations_Ampullaceana/freq_taxo_tbl.qzv
```
Note the **Median Frequency** and update `config.yaml`:
```yaml
diversity:
  max_depth: ??    # ← Median Frequency from taxa freq table
  raref_depth: ??  # ← will be set after viewing rarefaction curve
```

**Step 4:** Start Phase 2 (only `phylo_tree` will run until you set `raref_depth`):
```bash
snakemake -s Snakefile_phase2 --configfile config.yaml --cores 24 --until phylo_tree
```

**⚠️ CHECKPOINT 2:** Inspect the rarefaction curve:
```
diversity_results_viz_Ampullaceana/alpha_rarefaction_curves.qzv
```
Find where Shannon diversity stabilizes and update `config.yaml`:
```yaml
diversity:
  raref_depth: ??  # ← where the curve stabilizes
```

**Step 5:** Run the full Phase 2:
```bash
snakemake -s Snakefile_phase2 --configfile config.yaml --cores 24
```

---

## Dry-run (validate without executing)

```bash
# Phase 1
snakemake -s Snakefile_phase1 --configfile config.yaml -n

# Phase 2
snakemake -s Snakefile_phase2 --configfile config.yaml -n

# Visualize the DAG (requires graphviz)
snakemake -s Snakefile_phase1 --configfile config.yaml --dag | dot -Tpdf > dag_phase1.pdf
```

---

## Key outputs

| File | What to look for |
|---|---|
| `qc_reports/raw/multiqc_report.html` | Read quality before filtering |
| `qc_reports/filtered/multiqc_report.html` | Read quality after fastp |
| `visualizations_Ampullaceana/ampullaceana_demux.qzv` | **Set trunc_f/trunc_r** |
| `visualizations_Ampullaceana/taxa_bar_plots.qzv` | Taxonomic composition |
| `diversity_results_viz_Ampullaceana/alpha_rarefaction_curves.qzv` | **Set raref_depth** |
| `diversity_results_viz_Ampullaceana/faiths_pd_anova.qzv` | Faith PD ~ Temp * Diet |
| `diversity_results_viz_Ampullaceana/da_barplot_temp.qzv` | Diff. abundance by Temperature |
| `glmm_outputs/shannon_glmm_results.pdf` | GLMM Shannon diversity results |

---

## Citation

Pipeline adapted from:  
> DanielaDeodato. (2025). HonoMystes/R_temporaria_Metagenomics: v1.1.0. Zenodo. https://doi.org/10.5281/zenodo.18257526
=======
# Ampullaceana-balthica
>>>>>>> 4a337291ad035b54fde31b4aee6db200234e83ad
