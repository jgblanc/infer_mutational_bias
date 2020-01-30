
######## Snakemake header ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character",
        bench_iteration = "numeric",
        scriptdir = "character",
        source = "function"
    )
)
snakemake <- Snakemake(
    input = list('data/GWAS_ATLAS/AF_4361.txt'),
    output = list('output/GWAS_ATLAS/parsed_gwas/AF_4361_parsed.txt'),
    params = list(),
    wildcards = list('AF_4361', "trait" = 'AF_4361'),
    threads = 1,
    log = list(),
    resources = list(),
    config = list(),
    rule = 'parse_gwas_atlas',
    bench_iteration = as.numeric(NA),
    scriptdir = '/Users/jenniferblanc/infer_mutational_bias/code/parse_gwas_atlas',
    source = function(...){
        wd <- getwd()
        setwd(snakemake@scriptdir)
        source(...)
        setwd(wd)
    }
)

######## Original script #########
# Parse CAD_3925.txt

# Load libraries
library(data.table)
library(dplyr)

# Check that script is running
print("working")
print(snakemake@input[[1]])

# Load data
gwas_data <- fread(snakemake@input[[1]])
#gwas_data <- fread("~/infer_mutational_bias/data/GWAS_ATLAS/small.txt")
col_names <- colnames(gwas_data)

# Filter to P-values smaller than 5e-8
data <- gwas_data %>% filter(`Pvalue` <= 5e-8)

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
protective <- subset(data, data$EFFECT < 0)
protective$RAF <- 1 - protective$EAF
protective$RISK_ALLELE <- protective$NEA
data <- rbind(risk, protective)

## Write Results to Table
write.table(data, snakemake@output[[1]], quote = F, row.names = F)
print("Wrote table")
