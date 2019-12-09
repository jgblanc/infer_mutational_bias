library(arrow)
library(data.table)
library(dplyr)

# Check that script is running
print("working")
print(snakemake@input[[1]])
print(snakemake@input[[2]])

## Read in .parquet file
alpha_data <- read_parquet(snakemake@input[[1]])
#alpha_data <- read_parquet("~/infer_mutational_bias/data/47UKBB/body_HEIGHTz.funct.alphas.parquet")
print("Read in alpha file")

## Read in GWAS/Finemapping results
new_gwas_data <- fread(snakemake@input[[2]])
#gwas_data  <- fread("~/infer_mutational_bias/data/47UKBB/biochemistry_HDLcholesterol.txt")
#new_gwas_data <- fread("~/infer_mutational_bias/data/47UKBB/bolt_337K_unrelStringentBrit_MAF0.001_v3.biochemistry_HDLcholesterol.bgen.stats")
colnames(new_gwas_data) <- c("SNP", "CHR", "BP", "GENPOS" , "A1",  "A2", "A1FREQ" ,
                             "INFO", "CHISQ_LINREG" ,"P_LINREG",  "BETA" ,"SE","CHISQ_BOLT_LMM_INF" ,
                             "P_BOLT_LMM_INF","CHISQ_BOLT_LMM" , "P_BOLT_LMM")
print("Read in GWAS file")
print(head(new_gwas_data))

## Set up df to collect SNPs in credible sets
chr <- seq(1,22,1)
col_titles <- c("ALPHA1","ALPHA2", "ALPHA3", "ALPHA4", "ALPHA5", "ALPHA6","ALPHA7", "ALPHA8","ALPHA9","ALPHA10")
collect_data <- as.tbl(alpha_data[1,c(1,2,17,6,3,4,5)]) # Use existing table to preserve structure
collect_data$UNIQ_ID <- "test:test:test"
colnames(collect_data) <- c("SNP","CHR","REGION_START", "ALPHA", "BP", "A1", "A2","UNIQ_ID")
collect_data <- collect_data[,c(1,2,3,4,8,5,6,7)]

## Loop through all regions on all chromosomes and label all 95% CS for each regions
for (i in 1:22) { # Loop through chromosomes
  print(paste0("Working on CHR", i))
  chr_dat <- subset(alpha_data, alpha_data$CHR == chr[i])
  regions_start <- unique(chr_dat$REGION_START)
  for (j in 1:length(regions_start)) { # Loop through all regions
    region_dat <- subset(chr_dat, chr_dat$REGION_START == regions_start[j])
    for (k in 1:10) { # Loop through all ten columns of alpha values
      alpha_sorted <- arrange(region_dat, desc(region_dat[[col_titles[k]]])) # Sort by alpha values
      alpha_sorted_cumsum <- mutate(alpha_sorted, cumsum = cumsum(alpha_sorted[[col_titles[k]]])) # Calculate cumulative sum
      if (max(alpha_sorted_cumsum$cumsum) >= 0.95) { # If the cumulative sum reaches 0.95 there is a CS
        cutoff <- which(alpha_sorted_cumsum$cumsum >= 0.95)[1]
        target_snps <- alpha_sorted_cumsum[1:cutoff,c(1,2,17,k+5,3,4,5)] # Pull out rows of SNP in CS
        target_snps$UNIQ_ID <- paste0(i, ":", j, ":", k)
        target_snps <- target_snps[, c(1,2,3,4,8,5,6,7)]
        colnames(target_snps) <- c("SNP","CHR","REGION_START", "ALPHA" ,"UNIQ_ID", "BP", "A1", "A2")
        collect_data <- rbind(collect_data, target_snps) # Add new SNPs
      }
    }
  }
}
collect_data <- collect_data[-1,] # Remove initialization row

## For all SNPs in credible sets extract from FINEMAPP file
collect_data$CHR <- as.integer(collect_data$CHR)
collect_data$BP <- as.integer(collect_data$BP)
new <- left_join(collect_data, new_gwas_data, by = c("SNP", "CHR", "BP", "A1", "A2"))
print("Joined table")

## Write Results to Table
write.table(new, snakemake@output[[1]], quote = F, row.names = F)
print("Wrote table")
