# Parse depression_4293_RS.txt
# README for data (https://datashare.is.ed.ac.uk/handle/10283/3203)
# A1 is the effect allele


# Load libraries
library(data.table)
library(dplyr)

# Check that script is running
print("working")
print(snakemake@input[[1]])

# Load data
gwas_data <- fread(snakemake@input[[1]])
#gwas_data <- fread("~/infer_mutational_bias/data/GWAS_ATLAS/depression_4293_RS.txt")
col_names <- colnames(gwas_data)

# Filter to P-values smaller than 5e-8
data <- gwas_data %>% filter(P <= 5e-8)

# Add dummy columns to fill in later
data$CHR <- NA
data$BP <- NA

# Re-name columns
colnames(data) <- c("SNP", "EA", "NEA", "EAF", "LogOR", "SE", "P", "CHR", "BP")


# Document risk allele - only MAF info, can't get RAF
risk <- subset(data,data$LogOR >= 0)
risk$RAF <- risk$EAF
risk$RISK_ALLELE <- risk$EA
risk$RISK <- T
protective <- subset(data, data$LogOR < 0)
protective$RAF <- 1 - protective$EAF
protective$RISK_ALLELE <- protective$NEA
protective$RISK <- F
data <- rbind(risk, protective)

# Get MAF
minor <- subset(data, data$EAF <= 0.5)
minor$MAF <- minor$EAF
major <- subset(data, data$EAF > 0.5)
major$MAF <- 1 - major$EAF
data <- rbind(major, minor)

# Re-arrange columns
data <- data[,c("SNP","EA" ,"NEA", "CHR" ,"BP","LogOR" ,"SE" ,"P", "EAF", "MAF", "RAF", "RISK_ALLELE", "RISK")]

## Write Results to Table
write.table(data, snakemake@output[[1]], quote = F, row.names = F)
print("Wrote table")
