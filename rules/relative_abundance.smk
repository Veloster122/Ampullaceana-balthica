# rules/relative_abundance.smk
# Rule: Stacked bar plots of top 5 bacterial phyla and orders
# for each experimental factor (Diet, Temp, Phosphorus, Pop)

rule relative_abundance_plots:
    """Generate top-5 phylum and order relative abundance stacked bar charts
    for Diet, Temperature, Phosphorus and Population comparisons."""
    input:
        feature_tbl = config["directory_name"]["artifact"] + "/" + config["tables"]["taxa_freq"],
        taxonomy    = config["directory_name"]["artifact"] + "/" + config["taxonomy"]["taxo_data"],
        metadata    = config["raw"]["metadata"],
        script      = "scripts/relative_abundance_plots.R"
    output:
        diet_png        = "Outputs/RelativeAbundance_Diet.png",
        diet_pdf        = "Outputs/RelativeAbundance_Diet.pdf",
        temp_png        = "Outputs/RelativeAbundance_Temp.png",
        temp_pdf        = "Outputs/RelativeAbundance_Temp.pdf",
        phos_png        = "Outputs/RelativeAbundance_Phosphorus.png",
        phos_pdf        = "Outputs/RelativeAbundance_Phosphorus.pdf",
        pop_png         = "Outputs/RelativeAbundance_Pop.png",
        pop_pdf         = "Outputs/RelativeAbundance_Pop.pdf",
        pop_x_temp_png  = "Outputs/RelativeAbundance_Pop_x_Temp.png",
        pop_x_temp_pdf  = "Outputs/RelativeAbundance_Pop_x_Temp.pdf",
        pop_x_diet_png  = "Outputs/RelativeAbundance_Pop_x_Diet.png",
        pop_x_diet_pdf  = "Outputs/RelativeAbundance_Pop_x_Diet.pdf"
    shell:
        """
        Rscript {input.script} \
            {input.feature_tbl} \
            {input.taxonomy} \
            {input.metadata} \
            Outputs/
        """
