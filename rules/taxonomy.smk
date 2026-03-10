# rules/taxonomy.smk
# Rule 4: Train SILVA classifier and classify taxonomy

rule get_silva:
    """Download and prepare SILVA reference database for 16S classification."""
    output:
        rna_seqs  = "silva_refs/rna_origin_silva_ref_seqs.qza",
        taxa      = "silva_refs/origin_silva_ref_taxa.qza"
    params:
        version = config["taxonomy"]["version"],
        target  = config["taxonomy"]["target"]
    shell:
        """
        mkdir -p silva_refs
        qiime rescript get-silva-data \
            --p-version "{params.version}" \
            --p-target "{params.target}" \
            --p-include-species-labels \
            --p-no-rank-propagation \
            --o-silva-sequences {output.rna_seqs} \
            --o-silva-taxonomy {output.taxa}
        """


rule prepare_silva:
    """Reverse transcribe, cull, filter, and dereplicate SILVA reference."""
    input:
        rna_seqs = "silva_refs/rna_origin_silva_ref_seqs.qza",
        taxa     = "silva_refs/origin_silva_ref_taxa.qza"
    output:
        ref_seqs = "silva_refs/silva_ref_seqs.qza",
        ref_taxa = "silva_refs/silva_ref_taxa.qza"
    shell:
        """
        # Reverse transcribe RNA -> DNA
        qiime rescript reverse-transcribe \
            --i-rna-sequences {input.rna_seqs} \
            --o-dna-sequences silva_refs/origin_silva_ref_seqs.qza

        # Remove low-quality sequences
        qiime rescript cull-seqs \
            --i-sequences silva_refs/origin_silva_ref_seqs.qza \
            --o-clean-sequences silva_refs/cleaned_silva_ref_seqs.qza

        # Filter by length by taxon group
        qiime rescript filter-seqs-length-by-taxon \
            --i-sequences silva_refs/cleaned_silva_ref_seqs.qza \
            --i-taxonomy {input.taxa} \
            --p-labels Archaea Bacteria Eukaryota \
            --p-min-lens 900 1200 1400 \
            --o-filtered-seqs silva_refs/filt_silva_ref_seqs.qza \
            --o-discarded-seqs silva_refs/discard_silva_ref_seqs.qza

        # Dereplicate full-length seqs
        qiime rescript dereplicate \
            --i-sequences silva_refs/filt_silva_ref_seqs.qza \
            --i-taxa {input.taxa} \
            --p-mode uniq \
            --o-dereplicated-sequences silva_refs/derep_uniq_silva_ref_seqs.qza \
            --o-dereplicated-taxa silva_refs/derep_uniq_silva_ref_taxa.qza

        # Extract amplicon region (341F / 806R)
        qiime feature-classifier extract-reads \
            --i-sequences silva_refs/derep_uniq_silva_ref_seqs.qza \
            --p-f-primer {config[illumina][primer_f]} \
            --p-r-primer {config[illumina][primer_r]} \
            --p-n-jobs {config[raw][threads]} \
            --p-read-orientation forward \
            --o-reads silva_refs/derep_uniq_silva_ref_amplicon.qza

        # Dereplicate amplicon seqs
        qiime rescript dereplicate \
            --i-sequences silva_refs/derep_uniq_silva_ref_amplicon.qza \
            --i-taxa silva_refs/derep_uniq_silva_ref_taxa.qza \
            --p-mode uniq \
            --o-dereplicated-sequences {output.ref_seqs} \
            --o-dereplicated-taxa {output.ref_taxa}
        """


rule train_classifier:
    """Train a Naive Bayes taxonomy classifier on the SILVA amplicon region."""
    input:
        ref_seqs = "silva_refs/silva_ref_seqs.qza",
        ref_taxa = "silva_refs/silva_ref_taxa.qza"
    output:
        classifier = config["directory_name"]["artifact"] + "/" + config["taxonomy"]["classifier"]
    resources:
        mem_intensive = 1  # prevents parallel execution with dada2_denoise
    shell:
        """
        qiime feature-classifier fit-classifier-naive-bayes \
            --i-reference-reads {input.ref_seqs} \
            --i-reference-taxonomy {input.ref_taxa} \
            --o-classifier {output.classifier}
        """


rule classify_taxonomy:
    """Classify ASVs, filter tables by taxonomy, and produce bar plots."""
    input:
        classifier = config["directory_name"]["artifact"] + "/" + config["taxonomy"]["classifier"],
        seqs_rep   = config["directory_name"]["artifact"] + "/" + config["tables"]["seqs_rep"],
        freq_tbl   = config["directory_name"]["artifact"] + "/" + config["tables"]["freq_tbl"],
        metadata   = config["raw"]["metadata"]
    output:
        taxonomy   = config["directory_name"]["artifact"] + "/" + config["taxonomy"]["taxo_data"],
        taxo_tbl   = config["directory_name"]["artifact"] + "/" + config["tables"]["taxa_freq"],
        taxo_seqs  = config["directory_name"]["artifact"] + "/" + config["tables"]["taxa_seqs"],
        barplot    = config["directory_name"]["visualizations"] + "/taxa_bar_plots.qzv",
        taxo_viz   = config["directory_name"]["visualizations"] + "/taxonomy.qzv",
        taxo_tbl_viz = config["directory_name"]["visualizations"] + "/" + config["tables"]["taxa_freq_viz"],
        taxo_seqs_viz = config["directory_name"]["visualizations"] + "/" + config["tables"]["taxa_seqs_viz"]
    shell:
        """
        # Classify ASVs
        qiime feature-classifier classify-sklearn \
            --i-classifier {input.classifier} \
            --i-reads {input.seqs_rep} \
            --o-classification {output.taxonomy}

        qiime metadata tabulate \
            --m-input-file {output.taxonomy} \
            --o-visualization {output.taxo_viz}

        # Filter frequency table: keep bacteria/archaea, exclude chloroplasts/mitochondria
        qiime taxa filter-table \
            --i-table {input.freq_tbl} \
            --i-taxonomy {output.taxonomy} \
            --p-mode contains \
            --p-include p__ \
            --p-exclude 'p__;,Chloroplast,Mitochondria' \
            --o-filtered-table {output.taxo_tbl}

        # Filter representative sequences to match filtered table
        qiime feature-table filter-seqs \
            --i-data {input.seqs_rep} \
            --i-table {output.taxo_tbl} \
            --o-filtered-data {output.taxo_seqs}

        # Summarize filtered tables
        qiime feature-table summarize \
            --i-table {output.taxo_tbl} \
            --m-sample-metadata-file {input.metadata} \
            --o-visualization {output.taxo_tbl_viz}

        qiime feature-table tabulate-seqs \
            --i-data {output.taxo_seqs} \
            --o-visualization {output.taxo_seqs_viz}

        # Taxa bar plot
        qiime taxa barplot \
            --i-table {output.taxo_tbl} \
            --i-taxonomy {output.taxonomy} \
            --m-metadata-file {input.metadata} \
            --o-visualization {output.barplot}
        """
