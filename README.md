# Ampullaceana balthica — 16S Gut Microbiome Snakemake Pipeline

[![Snakemake](https://img.shields.io/badge/Snakemake-7.32.4-brightgreen.svg)](https://snakemake.github.io)
[![QIIME 2](https://img.shields.io/badge/QIIME_2-amplicon--2026.1-blue.svg)](https://qiime2.org)
[![R](https://img.shields.io/badge/R-4.3.3-blue.svg)](https://www.r-project.org)
[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org)

Pipeline bioinformático automatizado, rigoroso e reproduzível em **Snakemake** para análise de dados de sequenciação do gene **16S rRNA** (região V3–V4) em gastrópodes de água doce (*Ampullaceana balthica*).

---

## 🔬 Contexto Experimental

* **Espécie Hospedeira:** *Ampullaceana balthica* (Lymnaeidae)
* **Fatores Ambientais Testados:**
  * **Temperatura:** 14°C vs. 20°C (Stress térmico)
  * **Dieta:** Macroalgas / perifíton (Qualidade nutricional e estequiometria)
  * **Fósforo:** Limitação e enriquecimento estequiométrico
  * **Populações Geográficas:** Portugal (PT) e Suécia (SE)
* **Gene Alvo:** 16S rRNA (região hipervariável V3–V4)
* **Primers:**
  * Forward: `341F` (`CCTACGGGNGGCWGCAG`)
  * Reverse: `806R` (`GACTACHVGGGTATCTAATCC`)

---

## 📁 Estrutura do Repositório

```text
Ampullaceana-balthica/
│
├── config.yaml                    # Ficheiro central de parâmetros do pipeline
├── Snakefile_phase1               # Fase 1: Leituras Brutas ➔ DADA2 ➔ Classificação SILVA
├── Snakefile_phase2               # Fase 2: Filogenia ➔ Diversidade ➔ GLMMs ➔ Heat Trees
│
├── rules/                         # Regras modulares do Snakemake (.smk)
│   ├── fastp_import.smk           # QC e filtragem inicial com fastp
│   ├── cutadapt.smk               # Remoção precisa dos primers 341F e 806R
│   ├── qc_reports.smk             # Relatórios consolidados com MultiQC
│   ├── feature_table.smk          # Denoising DADA2, geração de ASVs e filtragem
│   ├── taxonomy.smk               # Classificador Naive Bayes SILVA 138.2 e filtro
│   ├── phylo_diversity.smk        # Filogenia (FastTree), rarefação e diversidade
│   ├── heattrees.smk              # Abundância diferencial com metacoder
│   └── glmm.smk                   # Integração dos scripts em R para GLMMs
│
├── scripts/                       # Scripts analíticos complementares
│   ├── GLMMs_automatized.R        # Modelos Lineares Mistos (glmmTMB) e PCA morfológico
│   └── HeatTrees_metacoder.R      # Geração dos cladogramas (metacoder + ANCOMBC2)
│
├── 00-Ampullaceana_balthica_raw_data/  # Metadados experimentais e pasta de FASTQs
│   ├── Ampullaceana_balthica_metadata.tsv
│   └── fastqs/
│
├── Outputs/                       # Figuras e resultados finais (PDFs, PNGs, tabelas)
│
│   # --- Componentes Opcionais (GUI Extra) ---
├── app.py                         # Interface gráfica opcional (Streamlit)
├── pipeline_runner.py             # Módulo de suporte à GUI e sincronização
├── iniciar_pipeline.bat           # Lançador de 1 clique para Windows
├── iniciar_pipeline.sh            # Lançador para Linux / Ubuntu
└── README.md
```

---

## 🛠️ Requisitos e Instalação

O pipeline foi desenhado para correr dentro do ambiente Conda do QIIME 2:

```bash
conda activate qiime2-amplicon-2026.1
conda install -c bioconda -c conda-forge snakemake fastp fastqc multiqc
```

---

## 🚀 Como Executar o Pipeline (Linha de Comandos)

O workflow está estruturado em **duas fases sequenciais** para permitir a inspeção de pontos de controlo (*checkpoints*) de qualidade:

### Fase 1: Do Pré-processamento à Taxonomia

1. **Configuração:** Verifique os caminhos e parâmetros no ficheiro `config.yaml`.
2. **Execução:**
   ```bash
   snakemake -s Snakefile_phase1 --configfile config.yaml --cores 20
   ```

> [!TIP]
> **⚠️ Ponto de Controlo 1 (Inspeção de Qualidade):**  
> Abra o ficheiro `visualizations_Ampullaceana/ampullaceana_demux.qzv` em [view.qiime2.org](https://view.qiime2.org).  
> Verifique onde a qualidade das leituras desce para atualizar `trunc_f` e `trunc_r` no `config.yaml`, se necessário.

---

### Fase 2: Filogenia, Diversidade, GLMMs e Heat Trees

1. Inspecione a tabela taxonómica gerada na Fase 1 (`visualizations_Ampullaceana/freq_taxo_tbl.qzv`) e confirme a **Mediana** de leituras para definir `max_depth` e `raref_depth` no `config.yaml`.
2. **Executar a Fase 2 completa:**
   ```bash
   snakemake -s Snakefile_phase2 --configfile config.yaml --cores 20
   ```

---

### Simulação Prévia (Dry-run)

Para testar a lógica do grafo de dependências (DAG) sem executar comandos pesados:
```bash
# Validar Fase 1
snakemake -s Snakefile_phase1 --configfile config.yaml -n

# Validar Fase 2
snakemake -s Snakefile_phase2 --configfile config.yaml -n
```

---

## 📊 Principais Resultados Gerados

| Ficheiro / Pasta | Descrição |
| :--- | :--- |
| `qc_reports/filtered/multiqc_report.html` | Relatório consolidado da qualidade das leituras após o Fastp. |
| `visualizations_Ampullaceana/dada2_stats.qzv` | Tabela detalhada da retenção de sequências ao longo do DADA2. |
| `Outputs/Shannon_plot.png`, `FaithPD_plot.png` | Boxplots de Diversidade Alfa comparando os tratamentos. |
| `Outputs/HeatTree_*.pdf` | Cladogramas de abundância diferencial gerados com o `metacoder` e `ANCOM-BC2`. |
| `Outputs/PCA_Traits_biplot_*.png` | Biplots do PCA aos traços morfológicos do hospedeiro (`Size_PC1`). |
| `glmm_pvalues.tsv`, `*_summary_results.txt` | Resumos estatísticos e tabelas ANOVA Tipo II dos modelos GLMM (`glmmTMB`). |

---

## ✨ Extra: Interface Gráfica Interativa (Streamlit Dashboard)

Como funcionalidade complementar e opcional para utilizadores que prefiram uma abordagem visual para explorar os resultados ou orquestrar o pipeline:

* **No Windows:** Dê um duplo-clique no ficheiro `iniciar_pipeline.bat`.
* **No Linux / Ubuntu (Servidor):** Execute:
  ```bash
  chmod +x iniciar_pipeline.sh
  ./iniciar_pipeline.sh
  ```
  *(Se executar no servidor, pode aceder a partir do navegador de qualquer computador da rede em `http://<IP_DO_SERVIDOR>:8501`).*

A aplicação abre no navegador (`http://localhost:8501`) e permite:
- Inspecionar a tabela de metadados e os ficheiros FASTQ;
- Ajustar os parâmetros do `config.yaml` com controlos visuais;
- Executar e monitorizar os comandos do Snakemake em tempo real;
- Visualizar interativamente todas as *Heat Trees*, boxplots de diversidade e relatórios estatísticos, com suporte a download dos dados compactados num ficheiro `.zip`.

*(Nota: Caso utilize um servidor remoto, as credenciais privadas podem ser configuradas no ficheiro `.env` sem nunca serem enviadas para o GitHub).*

---

## 📚 Citações e Referências

* **Snakemake:** Mölder, F., et al. (2021). *Sustainable data analysis with Snakemake*. F1000Research, 10:33. [doi:10.12688/f1000research.29032.2](https://doi.org/10.12688/f1000research.29032.2).
* **QIIME 2:** Bolyen, E., et al. (2019). *Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2*. Nature Biotechnology, 37(8), 852–860. [doi:10.1038/s41587-019-0209-9](https://doi.org/10.1038/s41587-019-0209-9).
* **SILVA Database:** Chuvochina, M., et al. (2026). *SILVA in 2026: a global core biodata resource for rRNA within the DSMZ digital diversity*. Nucleic Acids Research, 54(D1), D1–D8. [doi:10.1093/nar/gkaf1247](https://doi.org/10.1093/nar/gkaf1247).
* **DADA2:** Callahan, B. J., et al. (2016). *DADA2: High-resolution sample inference from Illumina amplicon data*. Nature Methods, 13(7), 581–583. [doi:10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869).
* **Metacoder:** Foster, Z. S. L., Sharpton, T. J., & Grünwald, N. J. (2017). *Metacoder: An R package for visualization and manipulation of community taxonomic diversity data*. PLOS Computational Biology, 13(2), e1005404. [doi:10.1371/journal.pcbi.1005404](https://doi.org/10.1371/journal.pcbi.1005404).
* **glmmTMB:** Brooks, M. E., et al. (2017). *glmmTMB balances speed and flexibility among packages for zero-inflated generalized linear mixed modeling*. The R Journal, 9(2), 378–400.
* **Fastp:** Chen, S., et al. (2018). *fastp: an ultra-fast all-in-one FASTQ preprocessor*. Bioinformatics, 34(17), i884–i890.
* **Pipeline Base:** DanielaDeodato. (2025). *HonoMystes/R_temporaria_Metagenomics: v1.1.0*. Zenodo. [doi:10.5281/zenodo.18257526](https://doi.org/10.5281/zenodo.18257526).

---

## 👥 Autoria

Desenvolvido no âmbito do estágio e dissertação de mestrado no **Centro de Ecologia, Evolução e Alterações Ambientais (cE3c)**.
