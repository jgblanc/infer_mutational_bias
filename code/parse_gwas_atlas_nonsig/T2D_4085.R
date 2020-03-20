# Parse T2D_4085.txt
# From the website and the column headings I know which allele is the effect allele and th EAF
# More info http://diagram-consortium.org/downloads/Mahajan.et.al.2018b.European.GWAS.readme.pdf

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
#gwas_data <- fread("~/infer_mutational_bias/data/GWAS_ATLAS/T2D_4085.txt")
col_names <- colnames(gwas_data)

# Filter to P-values smaller than 5e-8
gwas_data$Pvalue <- as.numeric(gwas_data$Pvalue)
#data <- gwas_data %>% filter(`Pvalue` >= thresh$V1)
data <- gwas_data %>% filter(abs(Beta) <= 0.0001)

# Delete column with total sample
data <- data[, -10]

# Re-name columns
colnames(data) <- c("SNP", "CHR", "BP", "EA", "NEA", "EAF", "BETA", "SE", "P")

# Re-arrange columns
data <- data[,c("SNP","EA" ,"NEA", "CHR" ,"BP","BETA" ,"SE" ,"P", "EAF")]

# Calculate MAF
minor <- subset(data, data$EAF <= 0.5)
minor$MAF <- minor$EAF
major <- subset(data, data$EAF > 0.5)
major$MAF <- 1 - major$EAF
data <- rbind(major, minor)


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

## Write Results to Table
write.table(data, snakemake@output[[1]], quote = F, row.names = F)
print("Wrote table")
