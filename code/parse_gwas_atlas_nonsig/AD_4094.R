# Parse AD_4094.txt
# The EA and EAF were given column headings

# Load libraries
library(data.table)
library(dplyr)

# Check that script is running
print("working")
print(snakemake@input[[1]])

# Load data
gwas_data <- fread(snakemake@input[[1]])
#gwas_data <- fread("~/infer_mutational_bias/data/GWAS_ATLAS/AD_4094.txt")

# Filter to P-values smaller than 5e-8
data <- gwas_data %>% filter(P <= 5e-8)

# Delete uneeded columns
data <- data[, -c(6,7,9,10,11)]

# Change SNP ID's to be CHR:BP
for (i in 1:nrow(data)) {
  data[i,1] <- paste0(data[i,2],":", data[i,3])
}

# Re-name columns
colnames(data) <- c("SNP", "CHR", "BP", "EA","NEA", "P", "EAF", "BETA", "SE")


# Document risk allele - only MAF info, can't get RAF
risk <- subset(data,data$BETA >= 0)
risk$RAF <- risk$EAF
risk$RISK_ALLELE <- risk$EA
risk$RISK <- T
protective <- subset(data, data$BETA < 0)
protective$RAF <- 1 - protective$EAF
protective$RISK_ALLELE <- protective$NEA
protective$RISK <- F
data <- rbind(risk, protective)

# Calculate MAF
minor <- subset(data, data$EAF <= 0.5)
minor$MAF <- minor$EAF
major <- subset(data, data$EAF > 0.5)
major$MAF <- 1 - major$EAF
data <- rbind(major, minor)

# Re-arrange columns
data <- data[,c("SNP","EA" ,"NEA", "CHR" ,"BP","BETA" ,"SE" ,"P", "EAF", "MAF", "RAF", "RISK_ALLELE", "RISK")]

## Write Results to Table
write.table(data, snakemake@output[[1]], quote = F, row.names = F)
print("Wrote table")
