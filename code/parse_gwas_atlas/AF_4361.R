# Parse AF_4361.txt

# Load libraries
library(data.table)
library(dplyr)


# Check that script is running
print("working")
print(snakemake@input[[1]])

# Load each split data
gwas_data <- fread(file = "~/infer_mutational_bias/data/GWAS_ATLAS/AF_split_aa")
gwas_data$Pvalue <- as.numeric(gwas_data$Pvalue)
data1 <- gwas_data %>% filter(`Pvalue` <= 5e-8)
col_names <- colnames(data1)
print(head(data1))

gwas_data <- fread(file = "~/infer_mutational_bias/data/GWAS_ATLAS/AF_split_ab")
gwas_data$V10 <- as.numeric(gwas_data$V10)
data2 <- gwas_data %>% filter(V10 <= 5e-8)
colnames(data2) <- col_names
print(head(data2))

gwas_data <- fread(file = "~/infer_mutational_bias/data/GWAS_ATLAS/AF_split_ac")
gwas_data$V10 <- as.numeric(gwas_data$V10)
data3 <- gwas_data %>% filter(V10 <= 5e-8)
colnames(data3) <- col_names
print(head(data3))

gwas_data <- fread(file = "~/infer_mutational_bias/data/GWAS_ATLAS/AF_split_ad")
gwas_data$V10 <- as.numeric(gwas_data$V10)
data4 <- gwas_data %>% filter(V10 <= 5e-8)
colnames(data4) <- col_names
print(head(data4))

data <- rbind(data1, data2, data3, data4)


# Get Rid of SNP ID's
data <- data[, -2]

# Change SNP ID's to be CHR:BP
for (i in 1:nrow(data)) {
  data[i,1] <- paste0(data[i,2],":", data[i,3])
}

# Re-name columns
colnames(data) <- c("SNP", "CHR", "BP", "NEA", "EA", "EAF", "EFFECT", "SE", "P")

# Re-arrange columns
data <- data[,c("SNP","EA" ,"NEA", "CHR" ,"BP","EFFECT" ,"SE" ,"P", "EAF")]

# Calculate MAF
minor <- subset(data, data$EAF <= 0.5)
minor$MAF <- minor$EAF
major <- subset(data, data$EAF > 0.5)
major$MAF <- 1 - major$EAF
data <- rbind(major, minor)

# Calculate RAF
risk <- subset(data,data$EFFECT >= 0)
risk$RAF <- risk$EAF
risk$RISK_ALLELE <- risk$EA
risk$RISK <- T
protective <- subset(data, data$EFFECT < 0)
protective$RAF <- 1 - protective$EAF
protective$RISK_ALLELE <- protective$NEA
protective$RISK <- F
data <- rbind(risk, protective)

## Write Results to Table
write.table(data, snakemake@output[[1]], quote = F, row.names = F)
print("Wrote table")
