# rules/heattrees.smk
# Rule: Differential Abundance via metacoder Heat Trees in R

rule heattrees_metacoder:
    """Generate taxonomic Heat Trees using metacoder to visualize differential abundance."""
    input:
        feature_tbl = config["directory_name"]["artifact"] + "/Ampullaceana_freq_table.qza",
        taxonomy    = config["directory_name"]["artifact"] + "/Ampullaceana_taxonomy.qza",
        metadata    = config["raw"]["metadata"],
        script      = "scripts/HeatTrees_metacoder.R"
    output:
        out_temp_png = "Outputs/HeatTree_Temp_20_vs_14.png",
        out_temp_pdf = "Outputs/HeatTree_Temp_20_vs_14.pdf",
        out_pop_png  = "Outputs/HeatTree_Pop_PT_vs_SE.png",
        out_pop_pdf  = "Outputs/HeatTree_Pop_PT_vs_SE.pdf",
        out_phos_png = "Outputs/HeatTree_Phosphorus_3_vs_0.png",
        out_phos_pdf = "Outputs/HeatTree_Phosphorus_3_vs_0.pdf"
    shell:
        """
        Rscript {input.script} {input.feature_tbl} {input.taxonomy} {input.metadata} Outputs/
        """
