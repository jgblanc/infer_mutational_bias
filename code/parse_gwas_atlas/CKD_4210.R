# Parse CKD_4210.txt
# From the website (http://ckdgen.imbi.uni-freiburg.de) I know that Allele1 is the effect allele

# Load libraries
library(data.table)
library(dplyr)

# Check that script is running
print("working")
print(snakemake@input[[1]])

# Load data
gwas_data <- fread(snakemake@input[[1]])
#gwas_data <- fread("~/infer_mutational_bias/data/GWAS_ATLAS/CKD_4210.txt")
col_names <- colnames(gwas_data)

# Filter to P-values smaller than 5e-8
data <- gwas_data %>% filter(`P-value` <= 5e-8)


# Change SNP ID's to be CHR:BP
for (i in 1:nrow(data)) {
  data[i,3] <- paste0(data[i,1],":", data[i,2])
}

# Delete column with total sample
data <- data[, -10]

# Convert alleles to upper case
data <- data.frame(lapply(data, function(v) {
  if (is.character(v)) return(toupper(v))
  else return(v)
}))

# Re-name columns
colnames(data) <- c("CHR", "BP", "SNP", "EA", "NEA", "EAF", "BETA", "SE", "P")

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
