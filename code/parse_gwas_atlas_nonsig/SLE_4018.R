# Parse SLE_4018.txt
# From table 1 in paper (https://arthritis-research.biomedcentral.com/articles/10.1186/s13075-018-1604-1/tables/1)
# I figured out that a2 was the effect allele

# Load libraries
library(data.table)
library(dplyr)

# Check that script is running
print("working")
print(snakemake@input[[1]])

# Load data
gwas_data <- fread(snakemake@input[[1]])
#gwas_data <- fread("~/infer_mutational_bias/data/GWAS_ATLAS/SLE_4018.txt")
col_names <- colnames(gwas_data)

# Filter to P-values smaller than 5e-8
data <- gwas_data %>% filter(`P-value` <= 5e-8)


# Change SNP ID's to be CHR:BP
for (i in 1:nrow(data)) {
  data[i,2] <- paste0(data[i,1],":", data[i,3])
}

# Convert alleles to upper case
data <- data.frame(lapply(data, function(v) {
  if (is.character(v)) return(toupper(v))
  else return(v)
}))

# Create dummy SE column so it matches other GWAS pages
data$SE <- NA
data$EAF <- NA

# Re-name columns
colnames(data) <- c("CHR", "SNP", "BP", "NEA", "EA", "Z_score", "P", "SE", "EAF")

# Re-arrange columns
data <- data[,c("SNP","EA" ,"NEA", "CHR" ,"BP","Z_score" ,"SE" ,"P", "EAF")]

# No MAF info
data$MAF <- NA

# Calculate RAF - no freq info
risk <- subset(data,data$Z_score >= 0)
risk$RAF <- NA
risk$RISK_ALLELE <- risk$EA
risk$RISK <- T
protective <- subset(data, data$Z_score < 0)
protective$RAF <- NA
protective$RISK_ALLELE <- protective$NEA
protective$RISK <- F
data <- rbind(risk, protective)

## Write Results to Table
write.table(data, snakemake@output[[1]], quote = F, row.names = F)
print("Wrote table")
