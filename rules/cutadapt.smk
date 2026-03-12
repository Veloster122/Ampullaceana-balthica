# rules/cutadapt.smk
# Rule: Remove primers with cutadapt before DADA2
# Uses the actual primer sequences (with IUPAC ambiguity codes) for accurate trimming.
# This is more reliable than trimming a fixed number of bases.

rule cutadapt_trim:
    """Trim primers from paired-end reads using cutadapt via QIIME2.
    Uses the actual 341F / 806R primer sequences so trimming is
    sequence-aware and handles degenerate IUPAC bases correctly.
    """
    input:
        demux = config["directory_name"]["artifact"] + "/ampullaceana_demux.qza"
    output:
        trimmed = config["directory_name"]["artifact"] + "/ampullaceana_trimmed.qza",
        trimmed_viz = config["directory_name"]["visualizations"] + "/ampullaceana_trimmed.qzv"
    params:
        primer_f = config["illumina"]["primer_f"],
        primer_r = config["illumina"]["primer_r"],
        threads   = config["raw"]["threads"]
    shell:
        """
        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences {input.demux} \
            --p-front-f {params.primer_f} \
            --p-front-r {params.primer_r} \
            --p-cores {params.threads} \
            --o-trimmed-sequences {output.trimmed} \
            --verbose

        qiime demux summarize \
            --i-data {output.trimmed} \
            --o-visualization {output.trimmed_viz}
        """


rule export_trimmed_demux:
    """Export post-cutadapt demux visualization to get per-sample read counts.
    These counts reflect reads AFTER primer removal, so samples emptied by
    cutadapt are correctly detected and removed by filter_samples.
    """
    input:
        qzv = config["directory_name"]["visualizations"] + "/ampullaceana_trimmed.qzv"
    output:
        counts = config["directory_name"]["trimmed_export"] + "/per-sample-fastq-counts.tsv"
    shell:
        """
        rm -rf {config[directory_name][trimmed_export]}
        qiime tools export \
            --input-path {input.qzv} \
            --output-path {config[directory_name][trimmed_export]}
        """
