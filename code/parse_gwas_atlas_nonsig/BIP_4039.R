# Parse BIP_4039.txt
# I based the polarization on the previous PGC gwas I parsed, I double checked using the SNPs in this table http://pgcdata.med.unc.edu/bipolar_disorder/BIP%202018%20readme.pdf
# and the allele freqs in the GWAS



# Load libraries
library(data.table)
library(dplyr)

# Check that script is running
print("working")
print(snakemake@input[[1]])

# Load data
gwas_data <- fread(snakemake@input[[1]], fill = T)
gwas_data <- fread("~/infer_mutational_bias/data/GWAS_ATLAS/BIP_4039.txt", fill= T)
col_nasmes <- colnames(gwas_data)

# Filter to P-values smaller than 5e-8
data <- gwas_data %>% filter(P <= 5e-8)

# Delete uneeded columns
data <- data[, -c(8,12,13,14,15,16,17,18,19)]

# Calculate the EAF - use weighted average of case/control freqs
data$EAF <- ((20352/51710)*data$FRQ_A_20352) + ((31358/51710)*data$FRQ_U_31358)

# Remove frequency columns
data <- data[, -c(6,7)]

# Change SNP ID's to be CHR:BP
for (i in 1:nrow(data)) {
  data[i,2] <- paste0(data[i,1],":", data[i,3])
}

# Re-name columns
colnames(data) <- c("CHR", "SNP", "BP", "EA", "NEA", "OR", "SE", "P", "EAF")

# Re-arrange columns
data <- data[,c("SNP","EA" ,"NEA", "CHR" ,"BP","OR" ,"SE" ,"P", "EAF")]

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

