# rules/feature_table.smk
# Rule 3: DADA2 denoising and feature table creation
# !! Run ONLY after setting trunc_f / trunc_r in config.yaml !!

rule filter_samples:
    """Filter out samples below the minimum sequence count threshold.
    Uses the cutadapt-trimmed artifact (primers already removed).
    """
    input:
        demux  = config["directory_name"]["artifact"] + "/ampullaceana_trimmed.qza",
        counts = config["directory_name"]["trimmed_export"] + "/per-sample-fastq-counts.tsv"
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


rule ensure_nonempty_demux:
    """Inspect filtered_demux.qza for empty FASTQ files and remove them.
    Definitive workaround for DADA2 1.30.0 R bug: mcmapply(filterAndTrim)
    crashes with a names-length mismatch when ANY sample FASTQ is empty.
    Directly opens the QZA ZIP to detect empties, then exports/filters/reimports.
    """
    input:
        demux = config["directory_name"]["artifact"] + "/filtered_demux.qza"
    output:
        demux_clean = config["directory_name"]["artifact"] + "/filtered_demux_clean.qza",
        report      = config["directory_name"]["artifact"] + "/empty_samples_removed.txt"
    run:
        import zipfile, os, subprocess, shutil, re, tempfile
        EMPTY_MAX = 100  # empty gzip stream is ~20-50 bytes on disk
        # Step 1: quick peek — are there any suspiciously small FASTQ.gz files?
        has_empty = False
        with zipfile.ZipFile(str(input.demux), 'r') as zf:
            for info in zf.infolist():
                if '/data/' in info.filename and info.filename.endswith('.fastq.gz'):
                    if info.compress_size < EMPTY_MAX:
                        has_empty = True
                        break
        if not has_empty:
            print("No empty FASTQ files detected — fast-copying artifact.")
            shutil.copy(str(input.demux), str(output.demux_clean))
            with open(str(output.report), 'w') as f:
                f.write("No empty samples detected.\n")
        else:
            # Step 2: full export → filter → reimport
            export_dir = tempfile.mkdtemp(prefix='demux_nonempty_')
            try:
                subprocess.run(["qiime", "tools", "export",
                                "--input-path",  str(input.demux),
                                "--output-path", export_dir], check=True)
                fwd_re = re.compile(r'^(.+?)(_S\d+_L\d+)?_R1(_\d+)?\.fastq\.gz$')
                lines  = ["sample-id\tforward-absolute-filepath\treverse-absolute-filepath"]
                removed = []
                for fname in sorted(os.listdir(export_dir)):
                    m = fwd_re.match(fname)
                    if not m:
                        continue
                    fwd = os.path.abspath(os.path.join(export_dir, fname))
                    rev = os.path.abspath(os.path.join(export_dir, fname.replace('_R1','_R2',1)))
                    if not os.path.exists(rev):
                        continue
                    sid = m.group(1)
                    if os.path.getsize(fwd) < EMPTY_MAX or os.path.getsize(rev) < EMPTY_MAX:
                        removed.append(sid)
                    else:
                        lines.append(f"{sid}\t{fwd}\t{rev}")
                manifest = os.path.join(export_dir, "manifest.tsv")
                with open(manifest, 'w') as f:
                    f.write('\n'.join(lines) + '\n')
                with open(str(output.report), 'w') as f:
                    f.write(f"Removed {len(removed)} empty samples:\n" + '\n'.join(removed) + '\n')
                print(f"Removed {len(removed)} empty samples. Keeping {len(lines)-1}.")
                subprocess.run(["qiime", "tools", "import",
                                "--type",         "SampleData[PairedEndSequencesWithQuality]",
                                "--input-path",   manifest,
                                "--output-path",  str(output.demux_clean),
                                "--input-format", "PairedEndFastqManifestPhred33V2"],
                               check=True)
            finally:
                shutil.rmtree(export_dir)


rule dada2_denoise:
    """Denoise, merge, and chimera-filter with DADA2."""
    input:
        filtered = config["directory_name"]["artifact"] + "/filtered_demux_clean.qza"
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
        max_ee_f  = config["denoise"]["max_ee_f"],
        max_ee_r  = config["denoise"]["max_ee_r"],
        chimera   = config["denoise"]["chimera_method"],
        tucker    = config["denoise"]["chimeric_parent_over_abundance"],
        overlap   = config["denoise"]["min_overlap"],
        threads   = config["denoise"]["dada2_threads"]
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
            --p-max-ee-f {params.max_ee_f} \
            --p-max-ee-r {params.max_ee_r} \
            --p-n-threads {params.threads} \
            --p-no-retain-all-samples \
            --verbose \
            --o-table {output.freq_tbl} \
            --o-representative-sequences {output.seqs_rep} \
            --o-denoising-stats {output.stats_qza} \
            --o-base-transition-stats {output.base_trans_stats}

        qiime metadata tabulate \
            --m-input-file {output.stats_qza} \
            --o-visualization {output.stats_qzv}
        """


rule summarize_feature_table:
    """Summarize frequency table and representative sequences for visualization.
    NOTE: In QIIME 2 2026.1, 'feature-table summarize' became a pipeline.
    The old --o-visualization is now --o-summary (visualisation output).
    We also unset R_HOME to prevent rpy2 from loading the system R library.
    """
    input:
        freq_tbl  = config["directory_name"]["artifact"] + "/" + config["tables"]["freq_tbl"],
        seqs_rep  = config["directory_name"]["artifact"] + "/" + config["tables"]["seqs_rep"],
        metadata  = config["raw"]["metadata"]
    output:
        freq_viz  = config["directory_name"]["visualizations"] + "/" + config["tables"]["freq_tbl_viz"],
        seqs_viz  = config["directory_name"]["visualizations"] + "/" + config["tables"]["seqs_rep_viz"],
        feat_freq = config["directory_name"]["artifact"] + "/feature_frequencies.qza",
        samp_freq = config["directory_name"]["artifact"] + "/sample_frequencies.qza"
    shell:
        """
        unset R_HOME R_LIBS_USER R_LIBS_SITE LD_LIBRARY_PATH

        qiime feature-table summarize \
            --i-table {input.freq_tbl} \
            --m-metadata {input.metadata} \
            --o-summary {output.freq_viz} \
            --o-feature-frequencies {output.feat_freq} \
            --o-sample-frequencies {output.samp_freq}

        qiime feature-table tabulate-seqs \
            --i-data {input.seqs_rep} \
            --o-visualization {output.seqs_viz}
        """
