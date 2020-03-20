# Parse PBC_3822.txt
# It was difficult to actually find the summary statitics, I was able to get the full summary statistics from
# (https://www.ebi.ac.uk/gwas/studies/GCST003129)
# The column headings tell us which is the effect allele

# Load libraries
library(data.table)
library(dplyr)

# Check that script is running
print("working")
print(snakemake@input[[1]])

# Load data
gwas_data <- fread(snakemake@input[[1]])
#gwas_data <- fread("~/infer_mutational_bias/data/GWAS_ATLAS/PBC_3822.txt")
col_names <- colnames(gwas_data)

# Filter to P-values smaller than 5e-8
data <- gwas_data %>% filter(p <= 5e-8)

# Delete column with total sample
data <- data[, -c(9,10,11)]

# Change SNP ID's to be CHR:BP
for (i in 1:nrow(data)) {
  data[i,3] <- paste0(data[i,1],":", data[i,2])
}

# Make dummy EAF column
data$EAF <- NA

# Re-name columns
colnames(data) <- c("CHR", "BP", "SNP", "NEA", "EA", "P", "BETA", "SE", "EAF")

# Re-arrange columns
data <- data[,c("SNP","EA" ,"NEA", "CHR" ,"BP","BETA" ,"SE" ,"P", "EAF")]

# Calculate MAF
data$MAF <- NA


# Calculate RAF - no freq info
risk <- subset(data,data$BETA >= 0)
risk$RAF <- NA
risk$RISK_ALLELE <- risk$EA
risk$RISK <- T
protective <- subset(data, data$BETA < 0)
protective$RAF <- NA
protective$RISK_ALLELE <- protective$NEA
protective$RISK <- F
data <- rbind(risk, protective)

## Write Results to Table
write.table(data, snakemake@output[[1]], quote = F, row.names = F)
print("Wrote table")
