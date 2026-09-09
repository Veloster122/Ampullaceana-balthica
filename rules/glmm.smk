# rules/glmm.smk
# Rule 7: Shannon extraction + GLMM statistical analysis in R

rule shannon_extract:
    """Extract Shannon, Observed Features, and Faith's PD metrics from alpha diversity output."""
    input:
        rarefact_dir = "diversity_core_metrics/shannon_vector.qza"
    output:
        shannon_meta  = "glmm_inputs/shannon_metadata.tsv",
        observed_meta = "glmm_inputs/observed_features_metadata.tsv",
        faithpd_meta  = "glmm_inputs/faith_pd_metadata.tsv"
    shell:
        """
        mkdir -p glmm_inputs

        # Unzip shannon vector artifact
        mkdir -p tmp_shannon
        unzip -o {input.rarefact_dir} -d tmp_shannon/
        find tmp_shannon/ -name "*.tsv" | head -1 | xargs -I{{}} cp {{}} {output.shannon_meta}
        rm -rf tmp_shannon

        # Unzip observed features vector artifact
        mkdir -p tmp_observed
        unzip -o diversity_core_metrics/observed_features_vector.qza -d tmp_observed/
        find tmp_observed/ -name "*.tsv" | head -1 | xargs -I{{}} cp {{}} {output.observed_meta}
        rm -rf tmp_observed

        # Unzip Faith's PD vector artifact
        mkdir -p tmp_faithpd
        unzip -o diversity_core_metrics/faith_pd_vector.qza -d tmp_faithpd/
        find tmp_faithpd/ -name "*.tsv" | head -1 | xargs -I{{}} cp {{}} {output.faithpd_meta}
        rm -rf tmp_faithpd
        """


rule glmm_analysis:
    """Run Generalized Linear Mixed Models on alpha diversity metrics.
    Uses Temp, Diet, Phosphorus, Pop as fixed effects and Box as random effect."""
    input:
        shannon_meta  = "glmm_inputs/shannon_metadata.tsv",
        observed_meta = "glmm_inputs/observed_features_metadata.tsv",
        faithpd_meta  = "glmm_inputs/faith_pd_metadata.tsv",
        metadata      = config["raw"]["metadata"],
        glmm_script   = "scripts/GLMMs_automatized.R"
    output:
        shannon_summary  = "glmm_outputs/shannon_entropy_summary_results.txt",
        observed_summary = "glmm_outputs/observed_features_summary_results.txt",
        faithpd_summary  = "glmm_outputs/faith_pd_summary_results.txt",
        shannon_dharma   = "glmm_outputs/shannon_entropy_DHARMa_residuals.png",
        shannon_plot     = "glmm_outputs/shannon_entropy_interaction_Pop_vs_Temp_by_Diet.png",
        pca_temp         = "glmm_outputs/PCA_Traits_biplot_Temp.png",
        pca_diet         = "glmm_outputs/PCA_Traits_biplot_Diet.png"
    shell:
        """
        mkdir -p glmm_outputs
        Rscript {input.glmm_script} \\
            {input.shannon_meta} \\
            {input.observed_meta} \\
            {input.faithpd_meta} \\
            {input.metadata} \\
            glmm_outputs
        """


rule appendix_tables:
    """Format GLMM outputs into clean appendix tables (B1/B2/B3) and
    record PCA figure paths (C1/C2). Runs after glmm_analysis."""
    input:
        faithpd_summary  = "glmm_outputs/faith_pd_summary_results.txt",
        observed_summary = "glmm_outputs/observed_features_summary_results.txt",
        shannon_summary  = "glmm_outputs/shannon_entropy_summary_results.txt",
        pca_summary      = "glmm_outputs/PCA_Traits_summary.txt",
        script           = "scripts/appendix_tables_figures.R"
    output:
        tbl_b1 = "Outputs/Appendices/Table_B1_GLMM_coefficients.csv",
        tbl_b2 = "Outputs/Appendices/Table_B2_GLMM_coefficients.csv",
        tbl_b3 = "Outputs/Appendices/Table_B3_GLMM_coefficients.csv",
        tbl_c1 = "Outputs/Appendices/Table_C1_PCA_importance.csv"
    shell:
        """
        mkdir -p Outputs/Appendices
        Rscript {input.script} glmm_outputs/ Outputs/Appendices/
        """
