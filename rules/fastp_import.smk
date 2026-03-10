# rules/fastp_import.smk
# Rule 1: Quality filter with fastp and import into QIIME2

rule fastp_qc:
    """Run fastp quality filtering on all raw paired-end reads."""
    input:
        r1 = config["raw"]["data_directory"] + "/{sample}" + config["illumina"]["paired_end_tag_1"] + ".fastq.gz",
        r2 = config["raw"]["data_directory"] + "/{sample}" + config["illumina"]["paired_end_tag_2"] + ".fastq.gz"
    output:
        r1 = config["directory_name"]["quality"] + "/out_{sample}" + config["illumina"]["paired_end_tag_1"] + ".fastq.gz",
        r2 = config["directory_name"]["quality"] + "/out_{sample}" + config["illumina"]["paired_end_tag_2"] + ".fastq.gz",
        html = config["directory_name"]["fastp"] + "/fastp_{sample}.html",
        json = config["directory_name"]["fastp"] + "/fastp_{sample}.json"
    threads: min(config["raw"]["threads"], 16)  # fastp supports max 16 threads
    shell:
        """
        mkdir -p {config[directory_name][quality]} {config[directory_name][fastp]}
        fastp \
            -i {input.r1} -I {input.r2} \
            -w {threads} \
            -h {output.html} -j {output.json} \
            -o {output.r1} -O {output.r2}
        """


rule make_manifest:
    """Build the manifest TSV from fastp-filtered reads."""
    input:
        expand(
            config["directory_name"]["quality"] + "/out_{sample}" + config["illumina"]["paired_end_tag_1"] + ".fastq.gz",
            sample=SAMPLES
        )
    output:
        manifest = config["raw"]["manifest"]
    run:
        import os
        quality_dir = os.path.abspath(config["directory_name"]["quality"])
        tag_r1 = config["illumina"]["paired_end_tag_1"]
        tag_r2 = config["illumina"]["paired_end_tag_2"]
        with open(output.manifest, "w") as f:
            f.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
            for s in SAMPLES:
                r1 = os.path.join(quality_dir, f"out_{s}{tag_r1}.fastq.gz")
                r2 = os.path.join(quality_dir, f"out_{s}{tag_r2}.fastq.gz")
                f.write(f"{s}\t{r1}\t{r2}\n")


rule import_demux:
    """Import fastp-filtered reads into a QIIME2 artifact."""
    input:
        manifest = config["raw"]["manifest"]
    output:
        qza = config["directory_name"]["artifact"] + "/ampullaceana_demux.qza",
        qzv = config["directory_name"]["visualizations"] + "/ampullaceana_demux.qzv"
    shell:
        """
        mkdir -p {config[directory_name][artifact]} {config[directory_name][visualizations]}
        qiime tools import \
            --type 'SampleData[PairedEndSequencesWithQuality]' \
            --input-path {input.manifest} \
            --output-path {output.qza} \
            --input-format PairedEndFastqManifestPhred33V2

        qiime demux summarize \
            --i-data {output.qza} \
            --o-visualization {output.qzv}
        """


rule export_demux:
    """Export the demux visualization to get per-sample-fastq-counts.tsv.
    This file is required by filter_samples to remove low-coverage samples.
    NOTE: If ampullaceana_demux.qzv already exists in the project root,
    copy it to the visualizations dir first so Snakemake can track it.
    """
    input:
        qzv = config["directory_name"]["visualizations"] + "/ampullaceana_demux.qzv"
    output:
        counts = config["directory_name"]["demux_export"] + "/per-sample-fastq-counts.tsv"
    shell:
        """
        rm -rf {config[directory_name][demux_export]}
        qiime tools export \
            --input-path {input.qzv} \
            --output-path {config[directory_name][demux_export]}
        """
