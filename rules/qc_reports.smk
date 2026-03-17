# rules/qc_reports.smk
# Rule 2: FastQC + MultiQC reports before and after fastp

rule qc_reports:
    """Run FastQC + MultiQC on raw reads and fastp-filtered reads for comparison."""
    input:
        raw_r1 = expand(
            config["raw"]["data_directory"] + "/{sample}" + config["illumina"]["paired_end_tag_1"] + ".fastq.gz",
            sample=SAMPLES
        ),
        filtered_r1 = expand(
            config["directory_name"]["quality"] + "/out_{sample}" + config["illumina"]["paired_end_tag_1"] + ".fastq.gz",
            sample=SAMPLES
        )
    output:
        multiqc_raw = "qc_reports/raw/multiqc_report.html",
        multiqc_filtered = "qc_reports/filtered/multiqc_report.html"
    threads: config["raw"]["threads"]
    shell:
        """
        mkdir -p qc_reports/raw qc_reports/filtered

        # Before QC
        fastqc {config[raw][data_directory]}/*.fastq.gz \
            -t {threads} \
            -o qc_reports/raw/
        multiqc --force qc_reports/raw/ -o qc_reports/raw/

        # After QC
        fastqc {config[directory_name][quality]}/*.fastq.gz \
            -t {threads} \
            -o qc_reports/filtered/
        multiqc --force qc_reports/filtered/ -o qc_reports/filtered/
        """
