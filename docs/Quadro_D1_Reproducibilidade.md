# Quadro D1 — Infraestrutura Computacional e Reprodutibilidade

## Repositório GitHub

| Campo | Valor |
|---|---|
| **URL** | https://github.com/Veloster122/Ampullaceana-balthica |
| **Branch principal** | `main` |
| **Licença** | — |
| **Última versão usada** | ver commit no servidor com `git log --oneline -5` |

---

## Estrutura do Repositório

```
Ampullaceana-balthica/
│
├── Snakefile_phase1              # Fase 1: QC → DADA2 → Taxonomia
├── Snakefile_phase2              # Fase 2: Filogenia → Diversidade → Análises R
├── config.yaml                   # Parâmetros configuráveis (paths, profundidades, etc.)
│
├── rules/                        # Regras Snakemake modulares
│   ├── phylo_diversity.smk       #   Filogenia, core metrics, ANCOMBC2
│   ├── glmm.smk                  #   GLMMs e PCA morfológico
│   ├── heattrees.smk             #   Heat Trees (metacoder)
│   └── relative_abundance.smk   #   Gráficos de abundância relativa
│
├── scripts/                      # Scripts de análise
│   ├── GLMMs_automatized.R       #   GLMMs + PCA + boxplots
│   ├── HeatTrees_metacoder.R     #   Heat Trees
│   ├── relative_abundance_plots.R #  Barras de abundância relativa (Top 10)
│   ├── appendix_tables_figures.R #  Tabelas e figuras dos apêndices B e C
│   └── generate_software_table.sh #  Tabela D1 de versões de software
│
├── 00-Ampullaceana_balthica_raw_data/
│   └── Ampullaceana_balthica_metadata.tsv   # Metadados das amostras
│
├── artifact_Ampullaceana/        # Artefactos QIIME2 (gerado pelo pipeline)
├── diversity_core_metrics/       # Métricas de core diversity (rarefied)
├── diversity_results_viz_Ampullaceana/  # Visualizações .qzv QIIME2
├── glmm_outputs/                 # Outputs dos GLMMs e PCA
├── Outputs/                      # Figuras finais (.png / .pdf)
│   └── Appendices/               #   Tabelas dos apêndices
└── .snakemake/                   # Metadados internos do Snakemake (não rastrear)
```

---

## Ficheiro `config.yaml` — Parâmetros Principais

| Parâmetro | Descrição |
|---|---|
| `directory_name.artifact` | Pasta dos artefactos QIIME2 (`artifact_Ampullaceana`) |
| `directory_name.phylogeny` | Pasta da árvore filogenética (`phylogeny_Ampullaceana`) |
| `directory_name.viz_dir_div` | Pasta das visualizações de diversidade |
| `tables.taxa_freq` | Tabela de frequências filtrada (`freq_taxo_tbl.qza`) |
| `taxonomy.taxo_data` | Artefacto de taxonomia SILVA (`Ampullaceana_taxonomy.qza`) |
| `diversity.raref_depth` | Profundidade de rarefação (definir após inspecção da curva) |
| `diversity.max_depth` | Profundidade máxima para a curva de rarefação |
| `metadata_columns.fixed_effects` | Efeitos fixos dos GLMMs (`Temp + Diet + Phosphorus + Pop`) |
| `metadata_columns.random_effect` | Efeito aleatório (`Box`) |

---

## Comandos de Execução

### Pré-requisitos

```bash
# 1. Clonar o repositório
git clone https://github.com/Veloster122/Ampullaceana-balthica.git
cd Ampullaceana-balthica

# 2. Activar o ambiente conda
conda activate qiime2-amplicon-2026.1
```

### Fase 1 — Controlo de Qualidade, DADA2 e Taxonomia

```bash
snakemake -s Snakefile_phase1 \
    --configfile config.yaml \
    --cores 20
```

> Após a Fase 1, inspecionar `diversity_results_viz_Ampullaceana/alpha_rarefaction_curves.qzv`  
> em [view.qiime2.org](https://view.qiime2.org) e definir `diversity.raref_depth` no `config.yaml`.

### Fase 2 — Diversidade, GLMMs, Heat Trees, Abundância Relativa

```bash
snakemake -s Snakefile_phase2 \
    --configfile config.yaml \
    --cores 20
```

### Gerar Tabelas dos Apêndices B e C

```bash
Rscript scripts/appendix_tables_figures.R \
    glmm_outputs/ \
    Outputs/Appendices/
```

### Gerar Tabela D1 — Versões de Software

```bash
conda activate qiime2-amplicon-2026.1
bash scripts/generate_software_table.sh \
    Outputs/Appendices/Table_D1_software_versions.tsv
```

---

## Nota sobre Reprodutibilidade

- Todos os parâmetros de análise estão centralizados em `config.yaml` — nenhum valor está *hardcoded* nos scripts.
- O Snakemake garante que apenas as regras cujos outputs estão em falta (ou cujo código mudou) são executadas novamente.
- Os artefactos QIIME2 (`.qza`/`.qzv`) são autocontidos e incluem os parâmetros exactos de cada passo.
- O ambiente conda `qiime2-amplicon-2026.1` pode ser exportado e reinstalado com:
  ```bash
  conda env export > environment.yml
  conda env create -f environment.yml
  ```
