# pull out a list of snps from summary stats - extract snps from 1kg - prune those snps - now pull out pruned snps from summary stats 

cut -f 1 -d" " ../output/GWAS_ATLAS/parsed_gwas/T2D_4085_0.95_parsed.txt > all_ss.snps | ./plink --bfile ../data/1000G_20101123_v3_GIANT_chr1_23_minimacnamesifnotRS_CEU_MAF0.01/1000G_20101123_v3_GIANT_chr1_23_minimacnamesifnotRS_CEU_MAF0.01_VARID --extract all_ss.snps --make-bed --out all_ss_plink --noweb

#./plink --bfile all_ss_plink --indep-pairwise 50 5 0.5 --noweb --out pruned_snps