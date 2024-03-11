library(data.table)
library(openxlsx)
library(dplyr)
library(ggplot2)

args <- commandArgs(TRUE)
wd <- args[1]
phen <- args[2]

setwd(wd)

a <- read.table(phen, sep = " ", header = FALSE)
b <- data.frame(FID = a$V1, IID = a$V2, PHE = a$V4)
c <- read.table("prs_phe.txt", sep = " ", header = TRUE)
d <- inner_join(b, c, by = c("FID", "IID"))
d$PHE <- sub(1, "case", d$PHE)
d$PHE <- sub(2, "control", d$PHE)
df <- d[c("PHE", "PRS")]
df$PHE <- factor(df$PHE)
plt <- ggplot(df, aes(x = PRS, fill = PHE, alpha = 0.5)) + geom_density(position = "identity")
ggsave(plot = plt, "PRS_distribution.png", dpi = 300, path = wd)

rm(list = ls())