# Apêndice C – Infraestrutura Computacional, Ferramentas e Reprodutibilidade

Este apêndice detalha a infraestrutura computacional, ambiente de execução, versões rigorosas de software e o controlo de qualidade das leituras brutas e processadas do estudo do microbioma de *Ampullaceana balthica*.

---

## Tabela C1: Versões exatas do sistema operativo, plataformas e pacotes de software utilizados

| Categoria | Software / Pacote | Versão | Fonte / Repositório | Função Principal no Pipeline | Citação / Referência |
|---|---|---|---|---|---|
| **Sistema Operativo** | Ubuntu Linux | 22.04.4 LTS (x86_64) | Canonical | Ambiente de computação de alto desempenho (HPC) | — |
| **Gestor de Workflow** | Snakemake | 7.32.4 | Conda / Bioconda | Orquestração reproduzível e automatizada dos DAGs das Fases 1 e 2 | Köster & Rahmann (2012) |
| **Linguagem Base** | Python | 3.10.14 | Conda-forge | Ambiente de execução de plugins e scripts auxiliares | Van Rossum & Drake (2009) |
| **Linguagem Estatística**| R | 4.3.3 / 4.4.0 | CRAN | Modelação GLMM, testes post-hoc, PCA e visualizações gráficas | R Core Team (2024) |
| **Plataforma Central** | QIIME 2 | 2026.1 (amplicon) | QIIME 2 Conda | Plataforma bioinformática central de microbioma | Bolyen et al. (2019) |
| **Controlo de Qualidade** | fastp | 0.23.4 | Bioconda | Filtragem de qualidade por janela móvel (Phred Q > 20) e remoção de adaptadores | Chen et al. (2018) |
| **Controlo de Qualidade** | FastQC | 0.12.1 | Babraham Inst. | Avaliação da qualidade por base das leituras brutas e filtradas | Andrews (2010) |
| **Relatório Integrado** | MultiQC | 1.19 / 1.21 | PyPI / Bioconda | Agregação e visualização interativa dos relatórios FastQC | Ewels et al. (2016) |
| **Remoção de Primers** | Cutadapt | 4.4 | Bioconda | Remoção precisa dos primers 341F e 806R nas orientações forward e reverse | Martin (2011) |
| **Inferência de ASVs** | DADA2 (`q2-dada2`) | 1.28.0 | Bioconductor / QIIME 2 | Modelação de erros de amplificação, desreplicação e resolução de ASVs | Callahan et al. (2016) |
| **Alinhamento Múltiplo** | MAFFT | 7.505 | Bioconda / QIIME 2 | Alinhamento múltiplo das sequências representativas (ASVs) | Katoh & Standley (2013) |
| **Árvore Filogenética** | FastTree | 2.1.11 | Bioconda / QIIME 2 | Inferência filogenética de máxima verossimilhança aproximada | Price et al. (2010) |
| **Base de Dados** | SILVA SSU NR99 | 138.2 | SILVA / DSMZ | Base taxonómica de referência para o gene 16S rRNA | Chuvochina et al. (2026) |
| **Classificador** | Naïve Bayes (`q2-feature-classifier`) | 2026.1 | QIIME 2 | Treino e atribuição taxonómica supervisionada de ASVs | Bokulich et al. (2018) |
| **Modelos Lineares Mistos**| `glmmTMB` | 1.1.14 | CRAN (R) | Modelação GLMM com famílias gaussianas e efeitos aleatórios (Box) | Brooks et al. (2017) |
| **Diagnóstico GLMM** | `DHARMa` | 0.4.6 | CRAN (R) | Diagnóstico e testes de resíduos quantílicos simulados | Hartig (2022) |
| **Médias Marginais** | `emmeans` | 1.10.0 | CRAN (R) | Estimativa de médias marginais (EMMs) e contrastes post-hoc com Bonferroni | Lenth (2024) |
| **Ecologia Comunitária** | `vegan` | 2.6-4 | CRAN (R) | Análises multivariadas e PERMANOVA | Oksanen et al. (2022) |
| **Abundância Diferencial**| `ANCOMBC` (ANCOM-BC2) | 2.4.0 | Bioconductor | Modelo log-ratio de abundância diferencial com controlo de viés | Lin & Peddada (2020, 2023) |
| **Árvores de Calor** | `metacoder` | 0.3.9 | CRAN (R) | Análise taxonómica diferencial com árvores de calor (Heat Trees) | Foster et al. (2017) |
| **Integração QIIME2-R** | `qiime2R` | 0.99.6 | GitHub | Importação direta de artefactos `.qza` e `.qzv` para estruturas do R | Bisanz (2018) |
| **Visualização Gráfica** | `ggplot2` | 3.5.0 | CRAN (R) | Construção de gráficos de barras, boxplots e PCA | Wickham (2016) |
| **Composição de Painéis** | `patchwork` | 1.2.0 | CRAN (R) | Arranjo multi-painel de figuras (Painéis A, B, C) | Pedersen (2024) |

---

## Figura C1: Perfil de qualidade das leituras antes e depois do controlo e filtragem

**Legenda para a Tese / Relatório:**
> **Figura C1: Perfil de qualidade médio por base (Phred Quality Score) das leituras das 180 amostras de *Ampullaceana balthica*.**  
> **(A)** Leituras brutas (*Raw reads*, R1 e R2) antes de qualquer intervenção, evidenciando a presença de adaptadores de sequenciação e o declínio típico da qualidade nas posições terminais (especialmente em R2).  
> **(B)** Leituras após controlo de qualidade e filtragem com `fastp` (janela móvel, remoção de adaptadores Illumina e descartes de leituras com $Q < 20$) e remoção dos primers 341F/806R com `Cutadapt`. Observa-se a retenção uniforme de qualidade acima de $Q30$ ao longo de quase toda a extensão dos amplicões antes da entrada no algoritmo de *denoising* DADA2.

### Como aceder e exportar os gráficos originais no servidor:

Os relatórios completos com os gráficos interativos de FastQC/MultiQC antes e depois da filtragem já foram gerados pela regra `qc_reports` da Fase 1 e estão disponíveis nos seguintes caminhos no servidor:

1. **Relatório MultiQC das Leituras Brutas (Antes da Filtragem):**
   ```bash
   qc_reports/raw/multiqc_report.html
   ```
   *Gráfico a exportar:* Secção **FastQC: Mean Quality Scores** (ou `qc_reports/raw/multiqc_data/mqc_fastqc_per_base_sequence_quality_plot_1.png` caso tenhas a pasta de dados exportada).

2. **Relatório MultiQC das Leituras Filtradas (Depois da Filtragem):**
   ```bash
   qc_reports/filtered/multiqc_report.html
   ```
   *Gráfico a exportar:* Secção **FastQC: Mean Quality Scores** (mostra a qualidade recuperada em verde $> Q30$).

3. **Visualizador QIIME2 do Demultiplexing:**
   ```bash
   visualizations_Ampullaceana/ampullaceana_demux.qzv
   ```
   *Abra em [view.qiime2.org](https://view.qiime2.org)*: O separador **Interactive Quality Plot** contém exatamente os boxplots interativos de $Q$-score por posição para Forward (1–300 bp) e Reverse (1–300 bp), permitindo exportar diretamente o gráfico em formato PDF ou PNG de alta resolução.
