# Parse UC_2030.txt
# I figured out the EA is allele 2 from comparing SNPs in table 1 from the paper (https://www.nature.com/articles/ng.3760/tables/1)

# Load libraries
library(data.table)
library(dplyr)

# Check that script is running
print("working")
print(snakemake@input[[1]])

# Load data
gwas_data <- fread(snakemake@input[[1]])
#gwas_data <- fread("~/infer_mutational_bias/data/GWAS_ATLAS/UC_2030.txt")
col_names <- colnames(gwas_data)

# Filter to P-values smaller than 5e-8
data <- gwas_data %>% filter(P.value <= 5e-8)

# Convert alleles to upper case
data <- data.frame(lapply(data, function(v) {
  if (is.character(v)) return(toupper(v))
  else return(v)
}))

# Delete uneeded columns
data <- data[, -c(7,8,9,10,11,12,13,14,15)]

# Change SNP ID's to be CHR:BP and add CHR and BP cols
data$CHR <- as.integer(rep(0,nrow(data)))
data$BP <- as.integer(rep(0,nrow(data)))
data$MarkerName <- as.character(data$MarkerName)
for (i in 1:nrow(data)) {
  print(i)
  splits <- base::strsplit(data[i,1],  "_")
  loc <- base::strsplit(splits[[1]][1],  ":")
  data[i,1] <- paste0(loc[[1]][1],":", loc[[1]][2])
  data[i, 7] <- loc[[1]][1]
  data[i, 8] <- loc[[1]][2]
}

# Re-name columns
colnames(data) <- c("SNP", "NEA", "EA", "BETA", "SE", "P", "CHR", "BP")

# Make dummy EAF/MAF column
data$EAF <- NA
data$MAF <- NA

# Document risk allele - can't get RAF, no freq info
risk <- subset(data,data$BETA >= 0)
risk$RAF <- NA
risk$RISK <- T
risk$RISK_ALLELE <- risk$EA
protective <- subset(data, data$BETA < 0)
protective$RAF <- NA
protective$RISK_ALLELE <- protective$NEA
protective$RISK <- F
data <- rbind(risk, protective)

# Re-arrange columns
data <- data[,c("SNP","EA" ,"NEA", "CHR" ,"BP","BETA" ,"SE" ,"P", "EAF", "MAF", "RAF", "RISK_ALLELE", "RISK")]

## Write Results to Table
write.table(data, snakemake@output[[1]], quote = F, row.names = F)
print("Wrote table")
