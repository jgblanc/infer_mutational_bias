print("Hello")

print("Im running")
# Open file with D/A info and create a dictionary that has key=chrom:pos and key_value=derived_allelee                                                                                 
filepath = "../data/1kg_phase3_snps.tsv"
da_dict = {"1":"A", "2":"T"}
with open(filepath) as fp:
   line = fp.readline()
   while line:
      line = line.strip()
      line_list = line.split()
      key = line_list[0] + ":" + line_list[1]
      key_value = line_list[9]
      da_dict[key] = key_value
      line = fp.readline()


print("Made Dictionary")

f = open("DA_Dict.txt","w")
f.write( str(da_dict) )
f.close()
print("Wrote File")


