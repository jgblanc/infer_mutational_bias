# Parse schizophrenia_3982.txt

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
#gwas_data <- fread("~/infer_mutational_bias/data/GWAS_ATLAS/schizophrenia_3982.txt")
high_quality_snps <- fread("data/GWAS_ATLAS/clozuk_pgc2.meta.sumstats.info9.snplist.txt")
col_names <- colnames(gwas_data)

# Filter to P-values smaller than 5e-8
data <- gwas_data %>% filter(P >= thresh$V1)

# Get only high quality SNPs
data_hq <- inner_join(high_quality_snps, data)

# Re-write colume names
colnames(data_hq) <- c("SNP","EA" ,"NEA", "EAF", "CHR" ,"BP","OR" ,"SE" ,"P")

# Change SNP ID's to be CHR:BP
for (i in 1:nrow(data_hq)) {
  #print(i)
  data_hq[i,1] <- paste0(data_hq[i,5],":", data_hq[i,6])
}

# Re-arrange columns
data <- data_hq[,c("SNP","EA" ,"NEA", "CHR" ,"BP","OR" ,"SE" ,"P", "EAF")]

# Calculate MAF
minor <- subset(data, data$EAF <= 0.5)
minor$MAF <- minor$EAF
major <- subset(data, data$EAF > 0.5)
major$MAF <- 1 - major$EAF
data <- rbind(major, minor)

# Calculate RAF
risk <- subset(data,data$OR >= 1)
risk$RAF <- risk$EAF
risk$RISK_ALLELE <- risk$EA
risk$RISK <- T
protective <- subset(data, data$OR < 1)
protective$RAF <- 1 - protective$EAF
protective$RISK_ALLELE <- protective$NEA
protective$RISK <- F
data <- rbind(risk, protective)


## Write Results to Table
write.table(data, snakemake@output[[1]], quote = F, row.names = F)
print("Wrote table")


