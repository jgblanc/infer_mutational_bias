# Parse CAD_3925.txt

# Load libraries
library(data.table)
library(dplyr)

# Check that script is running
print("working")
print(snakemake@input[[1]])
print(snakemake@input[[2]])

# Get P-value threshold
#thresh <- fread("data/GWAS_ATLAS/pval_thresholds/0.9.txt")
thresh <- fread(snakemake@input[[2]])

# Load data
gwas_data <- fread(snakemake@input[[1]])
gwas_data <- fread("~/infer_mutational_bias/data/GWAS_ATLAS/CAD_3925.txt")
col_names <- colnames(gwas_data)

# Filter to P-values smaller than 5e-8
data <- gwas_data %>% filter(`p-value_gc` >= thresh$V1)

# Get Rid of SNP ID's
data <- data[, -2]

# Change SNP ID's to be CHR:BP
#for (i in 1:nrow(data)) {
#  print(i)
#  data[i,1] <- paste0(data[i,2],":", data[i,3])
#}
data$Markername <- paste0(data$chr,":", data$bp_hg19)

# Get rid of n_sample, exome, info_ukbb columnes
data <- data[, -c(10,11,12)]

# Re-name columns
colnames(data) <- c("SNP", "CHR", "BP", "EA", "NEA", "EAF", "Log(OR)", "SE", "P")

# Re-arrange columns
data <- data[,c("SNP","EA" ,"NEA", "CHR" ,"BP","Log(OR)" ,"SE" ,"P", "EAF")]

# Calculate MAF
minor <- subset(data, data$EAF <= 0.5)
minor$MAF <- minor$EAF
major <- subset(data, data$EAF > 0.5)
major$MAF <- 1 - major$EAF
data <- rbind(major, minor)

# Calculate RAF - it seems like all log(OR) are greater than 0, meaning all risk SNPs?
risk <- subset(data,data$`Log(OR)` >= 0)
risk$RAF <- risk$EAF
risk$RISK_ALLELE <- risk$EA
risk$RISK <- T
protective <- subset(data, data$`Log(OR)` < 0)
protective$RAF <- 1 - protective$EAF
protective$RISK_ALLELE <- protective$NEA
data <- rbind(risk, protective)

## Write Results to Table
write.table(data, snakemake@output[[1]], quote = F, row.names = F)
print("Wrote table")
