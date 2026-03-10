# rules/glmm.smk
# Rule 7: Shannon extraction + GLMM statistical analysis in R

rule shannon_extract:
    """Extract Shannon and Observed Features metrics from alpha rarefaction output.
    Requires the alpha rarefaction .qza to be unzipped first."""
    input:
        rarefact_dir = "diversity_core_metrics/shannon_vector.qza"
    output:
        shannon_meta = "glmm_inputs/shannon_metadata.tsv",
        observed_meta = "glmm_inputs/observed_features_metadata.tsv"
    shell:
        """
        mkdir -p glmm_inputs
        # Unzip the shannon vector artifact
        mkdir -p tmp_shannon
        unzip -o {input.rarefact_dir} -d tmp_shannon/

        # Extract the data CSV from the QIIME2 artifact
        find tmp_shannon/ -name "*.tsv" | head -1 | xargs -I{{}} cp {{}} {output.shannon_meta}

        # Same for observed features
        mkdir -p tmp_observed
        unzip -o diversity_core_metrics/observed_features_vector.qza -d tmp_observed/
        find tmp_observed/ -name "*.tsv" | head -1 | xargs -I{{}} cp {{}} {output.observed_meta}

        rm -rf tmp_shannon tmp_observed
        """


rule glmm_analysis:
    """Run Generalized Linear Mixed Models on alpha diversity metrics.
    Uses Temp, Diet, Phosphorus, Pop as fixed effects and Box as random effect."""
    input:
        shannon_meta  = "glmm_inputs/shannon_metadata.tsv",
        observed_meta = "glmm_inputs/observed_features_metadata.tsv",
        metadata      = config["raw"]["metadata"],
        glmm_script   = "scripts/GLMMs_automatized.R"
    output:
        shannon_plot  = "glmm_outputs/shannon_glmm_results.pdf",
        observed_plot = "glmm_outputs/observed_features_glmm_results.pdf",
        pvalues_tsv   = "glmm_outputs/glmm_pvalues.tsv"
    shell:
        """
        mkdir -p glmm_outputs
        Rscript {input.glmm_script} \
            {input.shannon_meta} \
            {input.observed_meta} \
            {input.metadata} \
            {output.shannon_plot} \
            {output.observed_plot} \
            {output.pvalues_tsv}
        """
