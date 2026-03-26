# rules/phylo_diversity.smk
# Rules 5 & 6: Phylogeny, alpha rarefaction, and full diversity analysis
# !! Phase 2 rules: run ONLY after setting raref_depth in config.yaml !!

rule phylo_tree:
    """Build phylogenetic tree and compute alpha rarefaction curves."""
    input:
        taxo_seqs = config["directory_name"]["artifact"] + "/" + config["tables"]["taxa_seqs"],
        taxo_tbl  = config["directory_name"]["artifact"] + "/" + config["tables"]["taxa_freq"],
        metadata  = config["raw"]["metadata"]
    output:
        tree     = config["directory_name"]["phylogeny"] + "/tree.qza",
        rarefact = config["directory_name"]["viz_dir_div"] + "/alpha_rarefaction_curves.qzv"
    params:
        max_depth = config["diversity"]["max_depth"]
    shell:
        """
        mkdir -p {config[directory_name][phylogeny]} {config[directory_name][viz_dir_div]}

        # Build alignment and tree
        qiime phylogeny align-to-tree-mafft-fasttree \
            --i-sequences {input.taxo_seqs} \
            --o-alignment {config[directory_name][phylogeny]}/aligned-rep-seqs.qza \
            --o-masked-alignment {config[directory_name][phylogeny]}/masked-aligned-rep-seqs.qza \
            --o-tree {config[directory_name][phylogeny]}/unrooted-tree.qza \
            --o-rooted-tree {output.tree}

        # Alpha rarefaction (inspect this to set raref_depth!)
        qiime diversity alpha-rarefaction \
            --i-table {input.taxo_tbl} \
            --m-metadata-file {input.metadata} \
            --p-min-depth 1 \
            --p-max-depth {params.max_depth} \
            --o-visualization {output.rarefact}
        """


rule diversity_analysis:
    """Run core phylogenetic diversity metrics, alpha/beta significance tests,
    and ANCOMBC2 differential abundance for all key metadata variables.
    !! Requires raref_depth to be set in config.yaml !!
    """
    input:
        taxo_tbl = config["directory_name"]["artifact"] + "/" + config["tables"]["taxa_freq"],
        taxonomy = config["directory_name"]["artifact"] + "/" + config["taxonomy"]["taxo_data"],
        tree     = config["directory_name"]["phylogeny"] + "/tree.qza",
        metadata = config["raw"]["metadata"]
    output:
        # Core metrics directory sentinel
        core_flag    = "diversity_core_metrics/.done",
        shannon      = "diversity_core_metrics/shannon_vector.qza",
        observed     = "diversity_core_metrics/observed_features_vector.qza",
        # Alpha significance
        faithpd_stat = config["directory_name"]["viz_dir_div"] + "/faiths_pd_statistics.qzv",
        evenness_stat= config["directory_name"]["viz_dir_div"] + "/evenness_statistics.qzv",
        faithpd_anova= config["directory_name"]["viz_dir_div"] + "/faiths_pd_anova.qzv",
        # Beta significance (temp + diet)
        uw_temp      = config["directory_name"]["viz_dir_div"] + "/unweighted_unifrac_temp_significance.qzv",
        w_temp       = config["directory_name"]["viz_dir_div"] + "/weighted_unifrac_temp_significance.qzv",
        uw_diet      = config["directory_name"]["viz_dir_div"] + "/unweighted_unifrac_diet_significance.qzv",
        w_diet       = config["directory_name"]["viz_dir_div"] + "/weighted_unifrac_diet_significance.qzv",
        # Differential abundance (ANCOMBC2)
        da_temp      = config["directory_name"]["viz_dir_div"] + "/da_barplot_temp.qzv",
        da_temp_ref  = config["directory_name"]["viz_dir_div"] + "/da_barplot_temp_reference_20.qzv",
        da_diet      = config["directory_name"]["viz_dir_div"] + "/da_barplot_diet.qzv",
        da_p         = config["directory_name"]["viz_dir_div"] + "/da_barplot_phosphorus.qzv",
        da_combined  = config["directory_name"]["viz_dir_div"] + "/da_barplot_temp_diet_phosphorus.qzv"
    params:
        raref_depth = config["diversity"]["raref_depth"],
        artifact    = config["directory_name"]["artifact"],
        taxo_tbl    = config["tables"]["taxa_freq"],
        taxonomy    = config["taxonomy"]["taxo_data"]
    shell:
        """
        mkdir -p {config[directory_name][viz_dir_div]}

        # Core phylogenetic diversity metrics
        [ -d diversity_core_metrics ] && rm -rf diversity_core_metrics
        qiime diversity core-metrics-phylogenetic \
            --i-table {input.taxo_tbl} \
            --i-phylogeny {input.tree} \
            --m-metadata-file {input.metadata} \
            --p-sampling-depth {params.raref_depth} \
            --output-dir diversity_core_metrics
        touch {output.core_flag}

        # Alpha diversity significance
        qiime diversity alpha-group-significance \
            --i-alpha-diversity diversity_core_metrics/faith_pd_vector.qza \
            --m-metadata-file {input.metadata} \
            --o-visualization {output.faithpd_stat}

        qiime diversity alpha-group-significance \
            --i-alpha-diversity diversity_core_metrics/evenness_vector.qza \
            --m-metadata-file {input.metadata} \
            --o-visualization {output.evenness_stat}

        # ANOVA for faith_pd ~ Temp * Diet
        qiime longitudinal anova \
            --m-metadata-file diversity_core_metrics/faith_pd_vector.qza \
            --m-metadata-file {input.metadata} \
            --p-formula 'faith_pd ~ Temp * Diet' \
            --o-visualization {output.faithpd_anova}

        # Beta diversity significance (Temp)
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

        # Beta diversity significance (Diet)
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

        # --- Differential Abundance (ANCOMBC2) ---
        # Filter to features present in >= 4 samples with >= 50 reads
        qiime feature-table filter-features \
            --i-table {input.taxo_tbl} \
            --p-min-frequency 50 \
            --p-min-samples 4 \
            --o-filtered-table {params.artifact}/table_abund.qza

        # Collapse to order level (level 4)
        qiime taxa collapse \
            --i-table {params.artifact}/table_abund.qza \
            --i-taxonomy {input.taxonomy} \
            --p-level 4 \
            --o-collapsed-table {params.artifact}/table_abund_collapsed.qza

        # Temperature effect
        qiime composition ancombc2 \
            --i-table {params.artifact}/table_abund_collapsed.qza \
            --m-metadata-file {input.metadata} \
            --p-fixed-effects-formula 'Temp' \
            --o-ancombc2-output {params.artifact}/ancombc2_temp.qza
        qiime composition ancombc2-visualizer \
            --i-data {params.artifact}/ancombc2_temp.qza \
            --o-visualization {output.da_temp}

        # Temperature effect (reference = 20°C)
        qiime composition ancombc2 \
            --i-table {params.artifact}/table_abund_collapsed.qza \
            --m-metadata-file {input.metadata} \
            --p-fixed-effects-formula 'Temp' \
            --p-reference-levels Temp::20 \
            --o-ancombc2-output {params.artifact}/ancombc2_temp_ref20.qza
        qiime composition ancombc2-visualizer \
            --i-data {params.artifact}/ancombc2_temp_ref20.qza \
            --o-visualization {output.da_temp_ref}

        # Diet effect
        qiime composition ancombc2 \
            --i-table {params.artifact}/table_abund_collapsed.qza \
            --m-metadata-file {input.metadata} \
            --p-fixed-effects-formula 'Diet' \
            --o-ancombc2-output {params.artifact}/ancombc2_diet.qza
        qiime composition ancombc2-visualizer \
            --i-data {params.artifact}/ancombc2_diet.qza \
            --o-visualization {output.da_diet}

        # Phosphorus effect
        qiime composition ancombc2 \
            --i-table {params.artifact}/table_abund_collapsed.qza \
            --m-metadata-file {input.metadata} \
            --p-fixed-effects-formula 'Phosphorus' \
            --o-ancombc2-output {params.artifact}/ancombc2_p.qza
        qiime composition ancombc2-visualizer \
            --i-data {params.artifact}/ancombc2_p.qza \
            --o-visualization {output.da_p}

        # Combined: Temp + Diet + Phosphorus
        qiime composition ancombc2 \
            --i-table {params.artifact}/table_abund_collapsed.qza \
            --m-metadata-file {input.metadata} \
            --p-fixed-effects-formula 'Temp + Diet + Phosphorus' \
            --o-ancombc2-output {params.artifact}/ancombc2_combined.qza
        qiime composition ancombc2-visualizer \
            --i-data {params.artifact}/ancombc2_combined.qza \
            --o-visualization {output.da_combined}
        """
