# Parse parkinson_4167.txt
# From the table in at the bottom of the preprint (https://www.biorxiv.org/content/10.1101/388165v3.full.pdf)
# I figured out that A1 is the EA and freq refers to the EAF

# Load libraries
library(data.table)
library(dplyr)

# Check that script is running
print("working")
print(snakemake@input[[1]])

# Load data
gwas_data <- fread(snakemake@input[[1]])
#gwas_data <- fread("~/infer_mutational_bias/data/GWAS_ATLAS/parkinson_4167.tab")

# Filter to P-values smaller than 5e-8
data <- gwas_data %>% filter(p <= 5e-8)

# Delete uneeded columns
data <- data[, -c(8,9)]

# Change SNP ID's to be CHR:BP
data$CHR <- 0
data$BP <- 0
for (i in 1:nrow(data)) {
  print(i)
  l <- strsplit(data[i,1], ":")
  data[i,9] <- as.numeric(l[[1]][2])
  k <- strsplit(l[[1]][1], "r")
  data[i,8] <- as.numeric(k[[1]][2])
  data[i,1] <- paste0(as.numeric(k[[1]][2]),":", as.numeric(l[[1]][2]))
}

# Re-name columns
colnames(data) <- c("SNP", "EA","NEA", "EAF", "BETA", "SE", "P", "CHR", "BP")


# Document risk allele
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
