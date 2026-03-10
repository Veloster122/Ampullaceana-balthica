# rules/feature_table.smk
# Rule 3: DADA2 denoising and feature table creation
# !! Run ONLY after setting trunc_f / trunc_r in config.yaml !!

rule filter_samples:
    """Filter out samples below the minimum sequence count threshold."""
    input:
        demux  = config["directory_name"]["artifact"] + "/ampullaceana_demux.qza",
        counts = config["directory_name"]["demux_export"] + "/per-sample-fastq-counts.tsv"
    output:
        filtered = config["directory_name"]["artifact"] + "/filtered_demux.qza"
    params:
        min_seq = config["denoise"]["num_min_seq"]
    shell:
        """
        qiime demux filter-samples \
            --i-demux {input.demux} \
            --m-metadata-file {input.counts} \
            --p-where "CAST([forward sequence count] AS INT) > {params.min_seq}" \
            --o-filtered-demux {output.filtered}
        """


rule dada2_denoise:
    """Denoise, merge, and chimera-filter with DADA2."""
    input:
        filtered = config["directory_name"]["artifact"] + "/filtered_demux.qza"
    output:
        freq_tbl         = config["directory_name"]["artifact"] + "/" + config["tables"]["freq_tbl"],
        seqs_rep         = config["directory_name"]["artifact"] + "/" + config["tables"]["seqs_rep"],
        stats_qza        = config["directory_name"]["artifact"] + "/dada2_stats.qza",
        base_trans_stats = config["directory_name"]["artifact"] + "/base_transition_stats.qza",
        stats_qzv        = config["directory_name"]["visualizations"] + "/dada2_stats.qzv"
    resources:
        mem_intensive = 1  # prevents parallel execution with train_classifier
    params:
        trim_f    = config["denoise"]["trim_f"],
        trim_r    = config["denoise"]["trim_r"],
        trunc_f   = config["denoise"]["trunc_f"],
        trunc_r   = config["denoise"]["trunc_r"],
        chimera   = config["denoise"]["chimera_method"],
        tucker    = config["denoise"]["chimeric_parent_over_abundance"],
        overlap   = config["denoise"]["min_overlap"],
        threads   = config["raw"]["threads"]
    shell:
        """
        qiime dada2 denoise-paired \
            --i-demultiplexed-seqs {input.filtered} \
            --p-trim-left-f {params.trim_f} \
            --p-trim-left-r {params.trim_r} \
            --p-trunc-len-f {params.trunc_f} \
            --p-trunc-len-r {params.trunc_r} \
            --p-chimera-method {params.chimera} \
            --p-min-fold-parent-over-abundance {params.tucker} \
            --p-min-overlap {params.overlap} \
            --p-n-threads {params.threads} \
            --o-table {output.freq_tbl} \
            --o-representative-sequences {output.seqs_rep} \
            --o-denoising-stats {output.stats_qza} \
            --o-base-transition-stats {output.base_trans_stats}

        qiime metadata tabulate \
            --m-input-file {output.stats_qza} \
            --o-visualization {output.stats_qzv}
        """


rule summarize_feature_table:
    """Summarize frequency table and representative sequences for visualization."""
    input:
        freq_tbl  = config["directory_name"]["artifact"] + "/" + config["tables"]["freq_tbl"],
        seqs_rep  = config["directory_name"]["artifact"] + "/" + config["tables"]["seqs_rep"],
        metadata  = config["raw"]["metadata"]
    output:
        freq_viz  = config["directory_name"]["visualizations"] + "/" + config["tables"]["freq_tbl_viz"],
        seqs_viz  = config["directory_name"]["visualizations"] + "/" + config["tables"]["seqs_rep_viz"]
    shell:
        """
        qiime feature-table summarize \
            --i-table {input.freq_tbl} \
            --m-sample-metadata-file {input.metadata} \
            --o-visualization {output.freq_viz}

        qiime feature-table tabulate-seqs \
            --i-data {input.seqs_rep} \
            --o-visualization {output.seqs_viz}
        """
