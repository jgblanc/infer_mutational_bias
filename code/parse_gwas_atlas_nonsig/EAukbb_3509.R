# Parse EA_3409.txt
# This is UKBB data
# README for data (https://ctg.cncr.nl/documents/p1651/README_atlas_ukb2_sumstats)


# Load libraries
library(data.table)
library(dplyr)

# Check that script is running
print("working")
print(snakemake@input[[1]])

# Load data
gwas_data <- fread(snakemake@input[[1]])
#gwas_data <- fread("~/infer_mutational_bias/data/GWAS_ATLAS/myopia_3539.txt")
col_names <- colnames(gwas_data)

# Filter to P-values smaller than 5e-8
data <- gwas_data %>% filter(P <= 5e-8)

# Delete uneeded columns
data <- data[, -c(5,6,9,10)]
data <- data[, -c(7,11,12,13,14,15,16)]

# Change SNP ID's to be CHR:BP
for (i in 1:nrow(data)) {
  data[i,1] <- paste0(data[i,2],":", data[i,3])
}

# Re-name columns
colnames(data) <- c("SNP", "CHR", "BP", "EA", "OR", "SE", "P", "NEA", "MAF")

# Make dummy EAF column
data$EAF <- NA

# Document risk allele - only MAF info, can't get RAF
risk <- subset(data,data$OR >= 1)
risk$RAF <- NA
risk$RISK_ALLELE <- risk$EA
risk$RISK <- T
protective <- subset(data, data$OR < 1)
protective$RAF <- NA
protective$RISK_ALLELE <- protective$NEA
protective$RISK <- F
data <- rbind(risk, protective)

# Re-arrange columns
data <- data[,c("SNP","EA" ,"NEA", "CHR" ,"BP","OR" ,"SE" ,"P", "EAF", "MAF", "RAF", "RISK_ALLELE", "RISK")]

## Write Results to Table
write.table(data, snakemake@output[[1]], quote = F, row.names = F)
print("Wrote table")
