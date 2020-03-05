# Parse MS_24076602
# The Allels were very nicely labeled EA and NEA
# This was not from the GWAS ATLAS but I got the full summary stats from the GWAS catalog (https://www.ebi.ac.uk/gwas/studies/GCST005531)
# Results from this paper https://www.ncbi.nlm.nih.gov/pubmed/?term=24076602

# Load libraries
library(data.table)
library(dplyr)

# Check that script is running
print("working")
print(snakemake@input[[1]])

# Load data
gwas_data <- fread(snakemake@input[[1]])
gwas_data <- fread("~/infer_mutational_bias/data/GWAS_ATLAS/MS_24076602.txt")

# Filter to P-values smaller than 5e-8
data <- gwas_data %>% filter(p <= 5e-8)

# Delete uneeded columns
data <- data[, -c(9,10,11)]

# Change SNP ID's to be CHR:BP
for (i in 1:nrow(data)) {
  data[i,3] <- paste0(data[i,1],":", data[i,2])
}

# Re-name columns
colnames(data) <- c("CHR", "BP", "SNP", "NEA","EA", "P", "BETA", "SE")


# Document risk allele - only MAF info, can't get RAF
risk <- subset(data,data$BETA >= 0)
risk$RAF <- NA
risk$RISK_ALLELE <- risk$EA
risk$RISK <- T
protective <- subset(data, data$BETA < 0)
protective$RAF <-NA
protective$RISK_ALLELE <- protective$NEA
protective$RISK <- F
data <- rbind(risk, protective)

# Dummy Cols
data$MAF <- NA
data$EAF <- NA

# Re-arrange columns
data <- data[,c("SNP","EA" ,"NEA", "CHR" ,"BP","BETA" ,"SE" ,"P", "EAF", "MAF", "RAF", "RISK_ALLELE", "RISK")]

## Write Results to Table
write.table(data, snakemake@output[[1]], quote = F, row.names = F)
print("Wrote table")
