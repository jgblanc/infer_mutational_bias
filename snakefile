
rule parse_alphas_funct:
    input:
        "data/47UKBB/{trait}.funct.alphas.parquet",
	"data/47UKBB/bolt_337K_unrelStringentBrit_MAF0.001_v3.{trait}.bgen.stats"
    output:
        "output/47UKBB/parsed_alphas/{trait}_alphas_parsed.txt"
    script:
        "code/parse_alphas_47UKBB.R"

rule add_evo_info:
    input:
        "output/47UKBB/parsed_alphas/{trait}_alphas_parsed.txt"
    output:
        "output/47UKBB/evo_info_added/{trait}_evo_added.txt"
    script:
        "code/get_evolutionary_information_from_1kg_rsID.py"

rule add_RAF_info:
    input:
        "output/47UKBB/evo_info_added/{trait}_evo_added.txt"
    output:
        "output/47UKBB/RAF_info_added/{trait}_RAF.txt"
    script:
        "code/get_RAF.R"
     	