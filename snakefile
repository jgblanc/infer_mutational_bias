
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

rule parse_gwas_atlas:
    input:
        "data/GWAS_ATLAS/{trait}.txt"
    output:
        "output/GWAS_ATLAS/parsed_gwas/{trait}_parsed.txt"
    script:
        "code/parse_gwas_atlas/{wildcards.trait}.R"

rule parse_gwas_atlas_RSID:
    input:
        "data/GWAS_ATLAS/{trait}_RS.txt"
    output:
        "output/GWAS_ATLAS/parsed_gwas/{trait}.RS_parsed.txt"
    script:
        "code/parse_gwas_atlas/{wildcards.trait}_RS.R"

rule add_evo_atlas:
    input:
        "output/GWAS_ATLAS/parsed_gwas/{trait}_parsed.txt"
    output:
        "output/GWAS_ATLAS/evo_added/{trait}_evo.txt"
    script:
        "code/get_evolutionary_information_from_1kg_GWAS_ATLAS.py"

rule add_evo_atlas_RSID:
    input:
        "output/GWAS_ATLAS/parsed_gwas/{trait}.RS_parsed.txt"
    output:
        "output/GWAS_ATLAS/evo_added/{trait}_evo.txt"
    script:
        "code/get_evolutionary_information_from_1kg_GWAS_ATLAS_RSID.py"

rule parse_bbj:
    input:
        "data/BBJ/{trait}.txt"
    output:
        "output/BBJ/parsed_gwas/{trait}_parsed.txt"
    script:
        "code/parse_BBJ/{wildcards.trait}.R"

rule parse_neal:
    input:
        "data/UKBB/{trait}.txt"
    output:
        "output/UKBB/parsed_gwas/{trait}_parsed.txt"
    script:
        "code/Neal_UKBB.R"

rule LD_clumping_GWAS_atlas:
    input:
        "output/GWAS_ATLAS/parsed_gwas/{trait}_parsed.txt"
    output:
        "output/GWAS_ATLAS/clumped/{trait}.clumped"
    shell:
        "code/plink \
        --noweb \
        --bfile data/1000G_20101123_v3_GIANT_chr1_23_minimacnamesifnotRS_CEU_MAF0.01/1000G_20101123_v3_GIANT_chr1_23_minimacnamesifnotRS_CEU_MAF0.01_VARID \
        --clump {input} \
        --clump-field P \
        --clump-p1 1 \
        --clump-p2 1 \
        --clump-r2 0.5 \
        --clump-kb 250 \
        --out output/GWAS_ATLAS/clumped/{wildcards.trait}"

rule Extract_clumped_SNPs_GWAS_ATLAS:
    input:
        "output/GWAS_ATLAS/clumped/{trait}.clumped"
    output:
        "output/GWAS_ATLAS/clumped/{trait}_SNPs.txt"
    shell:
        "awk '{{ print $3}}' {input} > {output}"