#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 15:27:31 2019

@author: jenniferblanc
"""

print("Im running")
# Open file with D/A info and create a dictionary that has key=rsID and key_value=derived_allelee
filepath = "data/1kg_phase3_snps.tsv"
da_dict = {}
with open(filepath) as fp:
   line = fp.readline()
   while line:
      line = line.strip()
      line_list = line.split()
      key = line_list[2] 
      key_value = line_list[9]
      da_dict[key] = key_value
      line = fp.readline()

print("Made Dictionary")


# Open GWAS file and pull out derived allele for each SNP based on rsID matching in directory
filepath_GWAS = snakemake.input[0]
new_file_lines = []
with open(filepath_GWAS) as fp: 
   line = fp.readline()
   while line:
      line = line.strip()
      line_list = line.split()
      SNP = line_list[0]
      EA = line_list[6]
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



