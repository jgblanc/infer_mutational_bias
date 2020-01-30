
######## Snakemake header ########
import sys; sys.path.extend(["/Users/jenniferblanc/miniconda3/envs/GWAS_catalog/lib/python3.6/site-packages", "/Users/jenniferblanc/infer_mutational_bias/code"]); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05X;\x00\x00\x00output/GWAS_ATLAS/parsed_gwas/schizophrenia_3982_parsed.txtq\x06a}q\x07X\x06\x00\x00\x00_namesq\x08}q\tsbX\x06\x00\x00\x00outputq\ncsnakemake.io\nOutputFiles\nq\x0b)\x81q\x0cX6\x00\x00\x00output/GWAS_ATLAS/evo_added/schizophrenia_3982_evo.txtq\ra}q\x0eh\x08}q\x0fsbX\x06\x00\x00\x00paramsq\x10csnakemake.io\nParams\nq\x11)\x81q\x12}q\x13h\x08}q\x14sbX\t\x00\x00\x00wildcardsq\x15csnakemake.io\nWildcards\nq\x16)\x81q\x17X\x12\x00\x00\x00schizophrenia_3982q\x18a}q\x19(h\x08}q\x1aX\x05\x00\x00\x00traitq\x1bK\x00N\x86q\x1csh\x1bh\x18ubX\x07\x00\x00\x00threadsq\x1dK\x01X\t\x00\x00\x00resourcesq\x1ecsnakemake.io\nResources\nq\x1f)\x81q (K\x01K\x01e}q!(h\x08}q"(X\x06\x00\x00\x00_coresq#K\x00N\x86q$X\x06\x00\x00\x00_nodesq%K\x01N\x86q&uh#K\x01h%K\x01ubX\x03\x00\x00\x00logq\'csnakemake.io\nLog\nq()\x81q)}q*h\x08}q+sbX\x06\x00\x00\x00configq,}q-X\x04\x00\x00\x00ruleq.X\r\x00\x00\x00add_evo_atlasq/X\x0f\x00\x00\x00bench_iterationq0NX\t\x00\x00\x00scriptdirq1X/\x00\x00\x00/Users/jenniferblanc/infer_mutational_bias/codeq2ub.'); from snakemake.logging import logger; logger.printshellcmds = False; __real_file__ = __file__; __file__ = '/Users/jenniferblanc/infer_mutational_bias/code/get_evolutionary_information_from_1kg_GWAS_ATLAS.py';
######## Original script #########
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 11:24:13 2019

@author: jenniferblanc
"""
import json

print("Im running")
# Open file with D/A info and create a dictionary that has key=chrom:pos and key_value=derived_allele
da_dict = {}
with open('data/DA_Dict.json') as f:
    da_dict = json.load(f)
print("Loaded Dictionary")

# Open GWAS file and pull out derived allele for each SNP based on CHRO:POS matching in directory
filepath_GWAS = snakemake.input[0]
new_file_lines = []
with open(filepath_GWAS) as fp: 
   line = fp.readline()
   while line:
      line = line.strip()
      line_list = line.split()
      SNP = line_list[0]
      EA = line_list[1]
      try: 
          DA = da_dict[SNP]
          if DA == EA: 
              new_entry = "T"
          else: 
              new_entry = "F"
      except KeyError:
          new_entry = "NA" 
      line_list.append(new_entry)
      new_line = ",".join(line_list)
      new_file_lines.append(new_line)
      line = fp.readline()
      
print("Matched GWAS")    
    
        
# Re-write GWAS summary statistics table with column that says if effect allele is derieved
#outpath = "../output/evo_info_added/t2D_evo.jbb.txt"
outpath = snakemake.output[0]
with open(outpath, "w") as doc:
    doc.writelines("%s\n" % place for place in new_file_lines)
      
print("Wrote new file")      
      
      
      
      
      
      
      
      
      
      
      
      
