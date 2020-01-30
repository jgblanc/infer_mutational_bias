# Parse schizophrenia_EAS.txt
# It was very difficult to tell which allele was what in this GWAS
# The downloads page does not exist, I was able to access the scz gwas's at this link (https://www.med.unc.edu/pgc/results-and-downloads__trashed/scz/)
# I downloaded the EAS samples from EAS GWAS (211a46eb0d5b8642a47df4d131ec43f6)
# The paper is https://www.nature.com/articles/s41588-019-0512-x?draft=marketing#MOESM1 - I can't find a way to tell polarization from the paper
# From the supplement table 2 of that paper https://static-content.springer.com/esm/art%3A10.1038%2Fnature13595/MediaObjects/41586_2014_BFnature13595_MOESM75_ESM.pdf,
# they say that A1 is the effect allele for the frequency and the effect columns - I check a couple snps to make sure that they have the same OR as in my data



# Load libraries
library(data.table)
library(dplyr)

# Check that script is running
print("working")
print(snakemake@input[[1]])

# Load data
#gwas_data <- fread(snakemake@input[[1]], fill = T)
gwas_data <- fread("~/infer_mutational_bias/data/GWAS_ATLAS/schizophrenia_EAS.txt", fill= T)
col_nasmes <- colnames(gwas_data)

# Filter to P-values smaller than 5e-8
data <- gwas_data %>% filter(P <= 5e-8)

# Delete uneeded columns
data <- data[, -c(8,12,13,14,15,16,17,18)]

# Calculate the EAF - use weighted average of case/control freqs
data$EAF <- ((22778/58140)*data$FRQ_A_22778) + ((35362/58140)*data$FRQ_U_35362)

# Remove frequency columns
data <- data[, -c(6,7)]

# Change SNP ID's to be CHR:BP
for (i in 1:nrow(data)) {
  data[i,2] <- paste0(data[i,1],":", data[i,3])
}

# Re-name columns
colnames(data) <- c("CHR", "SNP", "BP", "EA", "NEA", "OR", "SE", "P", "EAF")

# Re-arrange columns
data <- data[,c("SNP","EA" ,"NEA", "CHR" ,"BP","OR" ,"SE" ,"P", "EAF")]

# Calculate MAF
minor <- subset(data, data$EAF <= 0.5)
minor$MAF <- minor$EAF
major <- subset(data, data$EAF > 0.5)
major$MAF <- 1 - major$EAF
data <- rbind(major, minor)

# Calculate RAF
risk <- subset(data,data$OR >= 1)
risk$RAF <- risk$EAF
risk$RISK_ALLELE <- risk$EA
risk$RISK <- T
protective <- subset(data, data$OR < 1)
protective$RAF <- 1 - protective$EAF
protective$RISK_ALLELE <- protective$NEA
protective$RISK <- F
data <- rbind(risk, protective)

## Write Results to Table
write.table(data, snakemake@output[[1]], quote = F, row.names = F)
print("Wrote table")
