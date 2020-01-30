# Parse vitiliago_3832.txt
# The summary statistics are here (ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/JinY_27723757_GCST004785)
# They were separated by chromosome and I concatened them together into vitiliago_3832.txt
# I figured out that the effect allele was allele1 from this table (https://www.nature.com/articles/ng.3680/tables/1)

# Load libraries
library(data.table)
library(dplyr)

# Check that script is running
print("working")
print(snakemake@input[[1]])

# Load data
gwas_data <- read.table(snakemake@input[[1]], header = F, fill=T)
print("loaded data")
col_names <- c("CHR","SNP","BP","A1","A2","MAF"	,"CHISQ","P","ORX","SE",	"L95",	"U95")
colnames(gwas_data) <- col_names

# Filter to P-values smaller than 5e-8
data <- gwas_data %>% filter(P <= 5e-8)

# Delete column with total sample
data <- data[, -c(7,11,12)]

# Change SNP ID's to be CHR:BP
data$CHR <- as.integer(data$CHR)
data$BP <- as.integer(data$BP)
data$SNP <- as.character(data$SNP)
for (i in 1:nrow(data)) {
  data[i,2] <- paste0(data[i,1],":", data[i,3])
}

# Re-name columns
colnames(data) <- c("CHR", "SNP", "BP", "EA", "NEA", "MAF", "P", "OR", "SE")

# Add EAF column
data$EAF <- NA

# Calculate RAF - no freq info
risk <- subset(data,data$OR >= 1)
risk$RAF <- NA
risk$RISK_ALLELE <- risk$EA
risk$RISK <- T
protective <- subset(data, data$OR < 1)
protective$RAF <- NA
protective$RISK_ALLELE <- protective$NEA
protective$RISK <- F
data <- rbind(risk, protective)

# Re-arrange columns
data <- data[,c("SNP","EA" ,"NEA", "CHR" ,"BP","OR" ,"SE" ,"P", "EAF", "MAF", "RAF", "RISK_ALLELE", "RISK")]

## Write Results to Table
write.table(data, snakemake@output[[1]], quote = F, row.names = F)
print("Wrote table")
