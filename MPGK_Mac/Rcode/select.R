library(data.table)
library(tidyverse)

args <- commandArgs(TRUE)
wd <- args[1]

setwd(wd)

a <- fread("phe.txt", na.strings = "-9")
b <- a %>% drop_na(V3)

row_n <- nrow(b)
row_target <- ceiling(row_n * 0.3)

set.seed(150)

re1 <- b %>% sample_n(row_target) %>% select(V1, V2)
re2 <- b %>% filter(!V2 %in% re1$V2) %>% select(V1, V2)

fwrite(re1, "name_target.txt", col.names = F, sep = " ", quote = F)
fwrite(re2, "name_base.txt", col.names = F, sep = " ", quote = F)

rm(list = ls())