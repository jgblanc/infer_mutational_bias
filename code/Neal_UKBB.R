# Parse all Neale Lab UKBB stats

# Load libraries
library(data.table)
library(dplyr)

# Check that script is running
print("working")
print(snakemake@input[[1]])

# Load data
gwas_data <- fread(snakemake@input[[1]])
#gwas_data <- fread("~/infer_mutational_bias/data/UKBB/Height_50.txt")
col_names <- colnames(gwas_data)

# Filter to P-values smaller than 5e-8
data <- gwas_data %>% filter(pval <= 5e-8)

# Delete uneeded columns
data <- data[, -c(4,5,6,7,10)]

# Change SNP ID's to be CHR:BP
data$CHR <- NA
data$BP <- NA
data$EA <- NA
data$NEA <- NA
for (i in 1:nrow(data)) {
  splits <- strsplit(data[i,1], ":")
  #print(i)
  data[i,7] <- splits[[1]][1]
  data[i,8] <- splits[[1]][2]
  data[i,10] <- splits[[1]][3]
  data[i,9] <- splits[[1]][4]
  data[i,1] <- paste0(data[i,7],":", data[i,8])
}

# Get EAF
data$EAF <- NA
for (i in 1:nrow(data)) {
  if (data[i,9] == data[i,2]) {
    data[i,11] <- data[i,3]
  } else {
    data[i,11] <- 1- data[i,3]
  }
}

data <- data[,-2]

# Re-name columns
colnames(data) <- c("SNP", "MAF", "BETA", "SE", "P", "CHR", "BP", "EA", "NEA", "EAF")


# Calculate RAF - no freq info
risk <- subset(data,data$BETA >= 0)
risk$RAF <- risk$EAF
risk$RISK_ALLELE <- risk$EA
risk$RISK <- T
protective <- subset(data, data$BETA < 0)
protective$RAF <- 1 - protective$EAF
protective$RISK_ALLELE <- protective$NEA
protective$RISK <- F
data <- rbind(risk, protective)

# Re-arrange columns
data <- data[,c("SNP","EA" ,"NEA", "CHR" ,"BP","BETA" ,"SE" ,"P", "EAF", "MAF", "RAF", "RISK_ALLELE", "RISK")]
colnames(data) <- c("SNP","EA" ,"NEA", "CHR" ,"BP","BETA" ,"SE" ,"P", "EAF", "MAF", "IAF", "INCREASING_ALLELE", "INCREASE")

## Write Results to Table
write.table(data, snakemake@output[[1]], quote = F, row.names = F)
print("Wrote table")
