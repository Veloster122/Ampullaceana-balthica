# rules/phylo_diversity.smk
# Phase 2 — Phylogeny, diversity metrics, and differential abundance
# Refactored into independent rules so a failure in one step does NOT
# cause the others to be deleted and rerun from scratch.
#
# Dependency chain:
#   phylo_tree
#       └── core_diversity
#               ├── alpha_significance
#               ├── beta_significance
#               └── ancombc2_prep
#                       ├── ancombc2_main
#                       └── ancombc2_pop
# !! Phase 2 rules: run ONLY after setting raref_depth in config.yaml !!

# ── Shortcuts ────────────────────────────────────────────────────────────────
_ART  = config["directory_name"]["artifact"]
_VIZ  = config["directory_name"]["viz_dir_div"]
_PHY  = config["directory_name"]["phylogeny"]


# ── Rule 1: Phylogenetic tree + alpha rarefaction curve ──────────────────────
rule phylo_tree:
    """Build phylogenetic tree and compute alpha rarefaction curves."""
    input:
        taxo_seqs = _ART + "/" + config["tables"]["taxa_seqs"],
        taxo_tbl  = _ART + "/" + config["tables"]["taxa_freq"],
        metadata  = config["raw"]["metadata"]
    output:
        tree     = _PHY + "/tree.qza",
        rarefact = _VIZ + "/alpha_rarefaction_curves.qzv"
    params:
        max_depth = config["diversity"]["max_depth"]
    shell:
        """
        mkdir -p {_PHY} {_VIZ}

        qiime phylogeny align-to-tree-mafft-fasttree \
            --i-sequences {input.taxo_seqs} \
            --o-alignment {_PHY}/aligned-rep-seqs.qza \
            --o-masked-alignment {_PHY}/masked-aligned-rep-seqs.qza \
            --o-tree {_PHY}/unrooted-tree.qza \
            --o-rooted-tree {output.tree}

        qiime diversity alpha-rarefaction \
            --i-table {input.taxo_tbl} \
            --m-metadata-file {input.metadata} \
            --p-min-depth 1 \
            --p-max-depth {params.max_depth} \
            --o-visualization {output.rarefact}
        """


# ── Rule 2: Core phylogenetic diversity metrics ───────────────────────────────
rule core_diversity:
    """Rarefy and compute core phylogenetic diversity metrics (shannon, faith_pd, unifrac…).
    !! Requires raref_depth to be set in config.yaml !!"""
    input:
        taxo_tbl = _ART + "/" + config["tables"]["taxa_freq"],
        tree     = _PHY + "/tree.qza",
        metadata = config["raw"]["metadata"]
    output:
        core_flag = "diversity_core_metrics/.done",
        shannon   = "diversity_core_metrics/shannon_vector.qza",
        observed  = "diversity_core_metrics/observed_features_vector.qza"
    params:
        raref_depth = config["diversity"]["raref_depth"]
    shell:
        """
        [ -d diversity_core_metrics ] && rm -rf diversity_core_metrics

        qiime diversity core-metrics-phylogenetic \
            --i-table {input.taxo_tbl} \
            --i-phylogeny {input.tree} \
            --m-metadata-file {input.metadata} \
            --p-sampling-depth {params.raref_depth} \
            --output-dir diversity_core_metrics

        touch {output.core_flag}
        """


# ── Rule 3: Alpha diversity significance ─────────────────────────────────────
rule alpha_significance:
    """Kruskal-Wallis tests on faith_pd and evenness; ANOVA faith_pd ~ Temp * Diet."""
    input:
        core_flag = "diversity_core_metrics/.done",
        metadata  = config["raw"]["metadata"]
    output:
        faithpd_stat  = _VIZ + "/faiths_pd_statistics.qzv",
        evenness_stat = _VIZ + "/evenness_statistics.qzv",
        faithpd_anova = _VIZ + "/faiths_pd_anova.qzv"
    shell:
        """
        mkdir -p {_VIZ}

        qiime diversity alpha-group-significance \
            --i-alpha-diversity diversity_core_metrics/faith_pd_vector.qza \
            --m-metadata-file {input.metadata} \
            --o-visualization {output.faithpd_stat}

        qiime diversity alpha-group-significance \
            --i-alpha-diversity diversity_core_metrics/evenness_vector.qza \
            --m-metadata-file {input.metadata} \
            --o-visualization {output.evenness_stat}

        qiime longitudinal anova \
            --m-metadata-file diversity_core_metrics/faith_pd_vector.qza \
            --m-metadata-file {input.metadata} \
            --p-formula 'faith_pd ~ Temp * Diet' \
            --o-visualization {output.faithpd_anova}
        """


# ── Rule 4: Beta diversity significance ──────────────────────────────────────
rule beta_significance:
    """PERMANOVA tests on UniFrac distance matrices for Temp and Diet."""
    input:
        core_flag = "diversity_core_metrics/.done",
        metadata  = config["raw"]["metadata"]
    output:
        uw_temp = _VIZ + "/unweighted_unifrac_temp_significance.qzv",
        w_temp  = _VIZ + "/weighted_unifrac_temp_significance.qzv",
        uw_diet = _VIZ + "/unweighted_unifrac_diet_significance.qzv",
        w_diet  = _VIZ + "/weighted_unifrac_diet_significance.qzv"
    shell:
        """
        mkdir -p {_VIZ}

        qiime diversity beta-group-significance \
            --i-distance-matrix diversity_core_metrics/unweighted_unifrac_distance_matrix.qza \
            --m-metadata-file {input.metadata} \
            --m-metadata-column Temp \
            --o-visualization {output.uw_temp}

        qiime diversity beta-group-significance \
            --i-distance-matrix diversity_core_metrics/weighted_unifrac_distance_matrix.qza \
            --m-metadata-file {input.metadata} \
            --m-metadata-column Temp \
            --o-visualization {output.w_temp}

        qiime diversity beta-group-significance \
            --i-distance-matrix diversity_core_metrics/unweighted_unifrac_distance_matrix.qza \
            --m-metadata-file {input.metadata} \
            --m-metadata-column Diet \
            --o-visualization {output.uw_diet}

        qiime diversity beta-group-significance \
            --i-distance-matrix diversity_core_metrics/weighted_unifrac_distance_matrix.qza \
            --m-metadata-file {input.metadata} \
            --m-metadata-column Diet \
            --o-visualization {output.w_diet}
        """


# ── Rule 5: ANCOMBC2 — prepare collapsed table ───────────────────────────────
rule ancombc2_prep:
    """Filter rare features and collapse taxonomy to order level (level 4)
    as pre-processing step shared by all ANCOMBC2 models."""
    input:
        taxo_tbl = _ART + "/" + config["tables"]["taxa_freq"],
        taxonomy = _ART + "/" + config["taxonomy"]["taxo_data"]
    output:
        abund_tbl     = _ART + "/table_abund.qza",
        collapsed_tbl = _ART + "/table_abund_collapsed.qza"
    shell:
        """
        qiime feature-table filter-features \
            --i-table {input.taxo_tbl} \
            --p-min-frequency 50 \
            --p-min-samples 4 \
            --o-filtered-table {output.abund_tbl}

        qiime taxa collapse \
            --i-table {output.abund_tbl} \
            --i-taxonomy {input.taxonomy} \
            --p-level 4 \
            --o-collapsed-table {output.collapsed_tbl}
        """


# ── Rule 6a: ANCOMBC2 — main environmental effects ──────────────────────────
rule ancombc2_main:
    """Differential abundance for Temp, Diet, Phosphorus (individual and combined)."""
    input:
        collapsed_tbl = _ART + "/table_abund_collapsed.qza",
        metadata      = config["raw"]["metadata"]
    output:
        da_temp     = _VIZ + "/da_barplot_temp.qzv",
        da_temp_ref = _VIZ + "/da_barplot_temp_reference_20.qzv",
        da_diet     = _VIZ + "/da_barplot_diet.qzv",
        da_p        = _VIZ + "/da_barplot_phosphorus.qzv",
        da_combined = _VIZ + "/da_barplot_temp_diet_phosphorus.qzv"
    shell:
        """
        mkdir -p {_VIZ}

        # Temperature effect (reference = 14°C, default)
        qiime composition ancombc2 \
            --i-table {input.collapsed_tbl} \
            --m-metadata-file {input.metadata} \
            --p-fixed-effects-formula 'Temp' \
            --o-ancombc2-output {_ART}/ancombc2_temp.qza
        qiime composition ancombc2-visualizer \
            --i-data {_ART}/ancombc2_temp.qza \
            --o-visualization {output.da_temp}

        # Temperature effect (reference = 20°C)
        qiime composition ancombc2 \
            --i-table {input.collapsed_tbl} \
            --m-metadata-file {input.metadata} \
            --p-fixed-effects-formula 'Temp' \
            --p-reference-levels Temp::20 \
            --o-ancombc2-output {_ART}/ancombc2_temp_ref20.qza
        qiime composition ancombc2-visualizer \
            --i-data {_ART}/ancombc2_temp_ref20.qza \
            --o-visualization {output.da_temp_ref}

        # Diet effect
        qiime composition ancombc2 \
            --i-table {input.collapsed_tbl} \
            --m-metadata-file {input.metadata} \
            --p-fixed-effects-formula 'Diet' \
            --o-ancombc2-output {_ART}/ancombc2_diet.qza
        qiime composition ancombc2-visualizer \
            --i-data {_ART}/ancombc2_diet.qza \
            --o-visualization {output.da_diet}

        # Phosphorus effect
        qiime composition ancombc2 \
            --i-table {input.collapsed_tbl} \
            --m-metadata-file {input.metadata} \
            --p-fixed-effects-formula 'Phosphorus' \
            --o-ancombc2-output {_ART}/ancombc2_p.qza
        qiime composition ancombc2-visualizer \
            --i-data {_ART}/ancombc2_p.qza \
            --o-visualization {output.da_p}

        # Combined: Temp + Diet + Phosphorus
        qiime composition ancombc2 \
            --i-table {input.collapsed_tbl} \
            --m-metadata-file {input.metadata} \
            --p-fixed-effects-formula 'Temp + Diet + Phosphorus' \
            --o-ancombc2-output {_ART}/ancombc2_combined.qza
        qiime composition ancombc2-visualizer \
            --i-data {_ART}/ancombc2_combined.qza \
            --o-visualization {output.da_combined}
        """


# ── Rule 6b: ANCOMBC2 — population interaction effects ──────────────────────
rule ancombc2_pop:
    """Differential abundance for Pop (additive) and Pop x Temp / Pop x Diet interactions.
    QIIME2 does not support R-style interaction syntax (Pop:Temp) in formula strings.
    Workaround: create combined categorical columns (Pop_x_Temp, Pop_x_Diet) in a
    temporary metadata file, then use each as a single variable in ANCOMBC2.
    This is statistically equivalent and captures all group contrasts."""
    input:
        collapsed_tbl = _ART + "/table_abund_collapsed.qza",
        metadata      = config["raw"]["metadata"]
    output:
        meta_inter    = _ART + "/metadata_interactions.tsv",
        da_pop_temp   = _VIZ + "/da_barplot_pop_temp.qzv",
        da_pop_x_temp = _VIZ + "/da_barplot_pop_x_temp.qzv",
        da_pop_x_diet = _VIZ + "/da_barplot_pop_x_diet.qzv"
    shell:
        """
        mkdir -p {_VIZ}

        # Pop + Temp (additive — main effects only)
        qiime composition ancombc2 \
            --i-table {input.collapsed_tbl} \
            --m-metadata-file {input.metadata} \
            --p-fixed-effects-formula 'Pop + Temp' \
            --o-ancombc2-output {_ART}/ancombc2_pop_temp.qza
        qiime composition ancombc2-visualizer \
            --i-data {_ART}/ancombc2_pop_temp.qza \
            --o-visualization {output.da_pop_temp}

        # Create metadata with combined interaction columns
        # Pop_x_Temp: PT_14, PT_20, SE_14, SE_20
        # Pop_x_Diet: PT_A,  PT_M,  PT_P,  SE_A,  SE_M,  SE_P
        python3 - <<'PYEOF'
import pandas as pd

meta = pd.read_csv("{input.metadata}", sep="\\t")

# Detect and preserve QIIME2 #q2:types directive row
has_types = str(meta.iloc[0, 0]).strip() == "#q2:types"
if has_types:
    types_row = meta.iloc[[0]].copy()
    meta      = meta.iloc[1:].copy()

meta["Pop_x_Temp"] = meta["Pop"].astype(str) + "_" + meta["Temp"].astype(str)
meta["Pop_x_Diet"] = meta["Pop"].astype(str) + "_" + meta["Diet"].astype(str)

if has_types:
    types_row["Pop_x_Temp"] = "categorical"
    types_row["Pop_x_Diet"] = "categorical"
    meta = pd.concat([types_row, meta], ignore_index=True)

meta.to_csv("{output.meta_inter}", sep="\\t", index=False)
PYEOF

        # Pop x Temp: each coefficient contrasts one Pop_x_Temp group vs reference (PT_14)
        qiime composition ancombc2 \
            --i-table {input.collapsed_tbl} \
            --m-metadata-file {output.meta_inter} \
            --p-fixed-effects-formula 'Pop_x_Temp' \
            --p-reference-levels 'Pop_x_Temp::PT_14' \
            --o-ancombc2-output {_ART}/ancombc2_pop_x_temp.qza
        qiime composition ancombc2-visualizer \
            --i-data {_ART}/ancombc2_pop_x_temp.qza \
            --o-visualization {output.da_pop_x_temp}

        # Pop x Diet: each coefficient contrasts one Pop_x_Diet group vs reference (PT_A)
        qiime composition ancombc2 \
            --i-table {input.collapsed_tbl} \
            --m-metadata-file {output.meta_inter} \
            --p-fixed-effects-formula 'Pop_x_Diet' \
            --p-reference-levels 'Pop_x_Diet::PT_A' \
            --o-ancombc2-output {_ART}/ancombc2_pop_x_diet.qza
        qiime composition ancombc2-visualizer \
            --i-data {_ART}/ancombc2_pop_x_diet.qza \
            --o-visualization {output.da_pop_x_diet}
        """
