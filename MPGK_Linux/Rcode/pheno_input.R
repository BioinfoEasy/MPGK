library(data.table)
library(dplyr)

args <- commandArgs(TRUE)
wd <- args[1]
pheno <- args[2]

setwd(wd)

fam <- read.table("dis.fam", sep = " ")
phe <- read.table(pheno, sep = " ")
colnames(fam) <- c("FID", "IID", "PID", "MID", "Sex", "Phe")
colnames(phe) <- c("FID", "IID", "GEN", "PHE")
temp <- merge(fam, phe, by = c("FID", "IID"))
temp_fam <- subset(temp, select = c("FID", "IID", "PID", "MID", "GEN", "PHE"))
temp_fam <- arrange(temp_fam, temp_fam$FID)
write.table(temp_fam, "dis.fam", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
rm(list = ls())