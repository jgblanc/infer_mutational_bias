library(dplyr)
library(data.table)

# Check that script is running
print("working")
print(snakemake@input[[1]])


# Read in data and adjust column names
data <- fread(snakemake@input[[1]])
#data <- fread("~/infer_mutational_bias/output/47UKBB/evo_info_added/body_HEIGHTz_evo_added.txt")
col_names <- colnames(data)
col_names[10] <- "EAF"
col_names[20] <- "EA_DERIVED"
colnames(data) <- col_names


# Create new columns
data$RISK <- NA
data$RAF <- NA
data$DAF <- NA
data$RISK_DERIVED <- NA
data$MAF <- NA
data$ES <- NA


# Tabulate Risk vs Protective SNPs
data_risk <- data %>% filter(BETA > 0) %>% mutate(RISK = T) %>% mutate(RAF = EAF)
data_protective <- data %>% filter(BETA < 0) %>% mutate(RISK = F) %>% mutate(RAF = 1 - EAF)
data_no_gwas <- data %>% filter(is.na(BETA))
rp_data <- rbind(data_risk, data_protective, data_no_gwas)

# Calculate DAF
D_derived <- rp_data %>% filter(EA_DERIVED == "T") %>% mutate(DAF = EAF)
D_ancestral <- rp_data %>% filter(EA_DERIVED == "F") %>% mutate(DAF = 1 - EAF)
D_no_info <- rp_data %>% filter(is.na(EA_DERIVED))
data <- rbind(D_derived, D_ancestral, D_no_info)

# Get RISK derived
derived_risk <- data %>% filter((EA_DERIVED == "T" & RISK == T) | (EA_DERIVED == "F" & RISK == F)) %>% mutate(RISK_DERIVED = T)
derived_proc <- data %>% filter((EA_DERIVED == "T" & RISK == F) | (EA_DERIVED == "F" & RISK == T)) %>% mutate(RISK_DERIVED = F)
derived_no_info <- data %>% filter(is.na(EA_DERIVED) | is.na(RISK))
data <- rbind(derived_risk, derived_proc, derived_no_info)

# Get Minor Allele Frequency
minor_effect <- data %>% filter(EAF > 0.5) %>% mutate(MAF = 1 - EAF)
major_effect <- data %>% filter(EAF <= 0.5) %>% mutate(MAF = EAF)
data <- rbind(minor_effect, major_effect)

# Transform beta cofficients into meaning effect sizes B=a*sqrt(2p(1-p))
data$ES <- data$BETA / sqrt(2 * data$MAF * (1 - data$MAF))

# Split Uniq-ID column
splits <- strsplit(data$UNIQ_ID, ":")
data$REGION_ID <- sapply(splits,function(x) x[2])
data$CS_ID <- sapply(splits,function(x) x[3])

# Write to table
write.table(data, snakemake@output[[1]], quote = F, row.names = F)
print("Wrote table")
