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
snp_dict = {}
with open(filepath) as fp:
    line = fp.readline()
#    counter = 1
    while line:
        line = line.strip()
        line_list = line.split()
        if line_list[11] == "derived_ancestral":
            key = line_list[2]
            key_value_snp = line_list[0] + ":" + line_list[1] 
            key_value_allele = line_list[9]
            key_value_daf = line_list[19]
            da_dict[key] = key_value_allele
            daf_dict[key] = key_value_daf
            snp_dict[key] = key_value_snp
            #counter = counter + 1
        line = fp.readline()
        #if counter % 100000 == 0: 
        #    print(counter)
        
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
            line_list[0] = snp_dict[SNP]
            chr_bp = line_list[0].split(":")
            line_list[3] = chr_bp[0]
            line_list[4] = chr_bp[1]
        except KeyError:
            new_entry = "NA" 
            daf = str("NA")
            r_derived = "NA"
            #line_lits[0] = snp_dict[SNP]
            #chr_bp = line_list[0].split(":")
            #line_list[3] = chr_bp[0]
            #line_list[4] = chr_bp[1]
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
      
      
      
      
      
      
      
      
      
      
      
      
