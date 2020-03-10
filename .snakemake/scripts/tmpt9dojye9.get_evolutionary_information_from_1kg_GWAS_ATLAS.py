
######## Snakemake header ########
import sys; sys.path.extend(["/Users/jenniferblanc/miniconda3/envs/GWAS_catalog/lib/python3.6/site-packages", "/Users/jenniferblanc/infer_mutational_bias/code"]); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05X1\x00\x00\x00output/GWAS_ATLAS/parsed_gwas/BIP_4039_parsed.txtq\x06a}q\x07X\x06\x00\x00\x00_namesq\x08}q\tsbX\x06\x00\x00\x00outputq\ncsnakemake.io\nOutputFiles\nq\x0b)\x81q\x0cX,\x00\x00\x00output/GWAS_ATLAS/evo_added/BIP_4039_evo.txtq\ra}q\x0eh\x08}q\x0fsbX\x06\x00\x00\x00paramsq\x10csnakemake.io\nParams\nq\x11)\x81q\x12}q\x13h\x08}q\x14sbX\t\x00\x00\x00wildcardsq\x15csnakemake.io\nWildcards\nq\x16)\x81q\x17X\x08\x00\x00\x00BIP_4039q\x18a}q\x19(h\x08}q\x1aX\x05\x00\x00\x00traitq\x1bK\x00N\x86q\x1csh\x1bh\x18ubX\x07\x00\x00\x00threadsq\x1dK\x01X\t\x00\x00\x00resourcesq\x1ecsnakemake.io\nResources\nq\x1f)\x81q (K\x01K\x01e}q!(h\x08}q"(X\x06\x00\x00\x00_coresq#K\x00N\x86q$X\x06\x00\x00\x00_nodesq%K\x01N\x86q&uh#K\x01h%K\x01ubX\x03\x00\x00\x00logq\'csnakemake.io\nLog\nq()\x81q)}q*h\x08}q+sbX\x06\x00\x00\x00configq,}q-X\x04\x00\x00\x00ruleq.X\r\x00\x00\x00add_evo_atlasq/X\x0f\x00\x00\x00bench_iterationq0NX\t\x00\x00\x00scriptdirq1X/\x00\x00\x00/Users/jenniferblanc/infer_mutational_bias/codeq2ub.'); from snakemake.logging import logger; logger.printshellcmds = False; __real_file__ = __file__; __file__ = '/Users/jenniferblanc/infer_mutational_bias/code/get_evolutionary_information_from_1kg_GWAS_ATLAS.py';
######## Original script #########
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 11:24:13 2019

@author: jenniferblanc
"""

print("Im running")
# Open file with D/A info and create a dictionary that has key=chrom:pos and key_value=derived_allele
# Open file with D/A info and create a dictionary that has key=rsID and key_value=derived_allelee
filepath = "data/1kg_phase3_snps.tsv"
da_dict = {}
daf_dict = {}
with open(filepath) as fp:
    line = fp.readline()
    while line:
        line = line.strip()
        line_list = line.split()
        if line_list[11] == "derived_ancestral":
            key = line_list[0] + ":" + line_list[1] 
            key_value_allele = line_list[9]
            key_value_daf = line_list[19]
            da_dict[key] = key_value_allele
            daf_dict[key] = key_value_daf
        line = fp.readline()
        
print("Made Dictionary") 

# Open GWAS file and pull out derived allele for each SNP based on CHRO:POS matching in directory
filepath_GWAS = snakemake.input[0]
new_file_lines = []
with open(filepath_GWAS) as fp: 
    line = fp.readline()
    while line:
        line = line.strip()
        line_list = line.split()
        SNP = line_list[0]
        print(SNP)
        EA = line_list[1]
        RA = line_list[11]
        EAF = line_list[8]
        MAF = line_list[9]
        RAF = line_list[10]
        try: 
            DA = da_dict[SNP]
            DAF = float(daf_dict[SNP])
            if MAF != "NA" and EAF == "NA":
                print("GOOD")
                if DAF < 0.5:
                    if DA == EA:
                        line_list[8] = MAF
                    else:
                        line_list[8] = str(1 - float(MAF)) 
                if DAF >= 0.5:
                    if DA == EA:
                        line_list[8]= str(1 - float(MAF)) 
                    else:
                        line_list[8] = MAF
                if RA == EA:
                    line_list[10] = str(line_list[8])
                else:
                    line_list[10] = str(1 - float(line_list[8]))
                EAF = line_list[8]
                RAF = line_list[10]
            if DA == EA: 
                new_entry = "T"
                if EAF == "NA":
                    line_list[8] = str(DAF) 
                daf = str(line_list[8])
            else: 
                new_entry = "F" 
                if EAF == "NA":
                    line_list[8] = str(1 - DAF)
                daf = str(1 - float(line_list[8]))
            if DA == RA:
                r_derived = "T"
                if RAF == "NA":
                    line_list[10] = str(DAF)
            else:
                r_derived = "F"
                if RAF == "NA":
                    line_list[10] = str(1 - DAF)
            if MAF == "NA":
                if DAF < 0.5:
                    line_list[9] = str(DAF)
                else:
                    line_list[9] = str(1-DAF)
        except KeyError:
            new_entry = "NA" 
            daf = str("NA")
            r_derived = "NA"
        line_list.append(new_entry)
        line_list.append(daf)
        line_list.append(r_derived)
        new_line = ",".join(line_list)
        print(new_line)
        new_file_lines.append(new_line)
        line = fp.readline()     
print("Matched GWAS")       
    
        
# Re-write GWAS summary statistics table with column that says if effect allele is derieved
outpath = snakemake.output[0]
with open(outpath, "w") as doc:
    doc.writelines("%s\n" % place for place in new_file_lines)
      
print("Wrote new file")      
      
      
      
      
      
      
      
      
      
      
      
      
