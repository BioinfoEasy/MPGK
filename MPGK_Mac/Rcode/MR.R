library(TwoSampleMR)
library(ggplot2)

args <- commandArgs(TRUE)
workdir <- args[1]
exposure_id <- args[2]
outcome_id <- args[3]
p_cutoff <- args[4]
p_cutoff <- as.numeric(p_cutoff)

setwd(workdir)

df_exp <- read.table(exposure_id, header = TRUE, sep = " ")
df_out <- read.table(outcome_id, header = TRUE, sep = " ")

#relevance assumption
exp_filter <- subset(df_exp, pval < p_cutoff) 
write.csv(exp_filter, file = "exp_filter.csv") 
exp_clumped <- read_exposure_data(filename = "exp_filter.csv", sep = ",", snp_col = "SNP", beta_col = "beta", eaf_col = "eaf", effect_allele_col = "ref", other_allele_col = "alt", pval_col = "pval")

#exclusivity assumption
out_merge <- merge(exp_clumped, df_out, by.x = "SNP", by.y = "SNP")
write.csv(out_merge, file = "out_merged.csv")
out_merged <- read_outcome_data(snps = exp_clumped$SNP, filename = "out_merged.csv", sep = ",", snp_col = "SNP", beta_col = "beta", eaf_col = "eaf", effect_allele_col = "ref", other_allele_col = "alt", pval_col = "pval")

#Mendelian randomization analysis
har_data <- harmonise_data(exposure_dat = exp_clumped, outcome_dat = out_merged)
res_mr <- mr(har_data)
res_mr <- generate_odds_ratios(mr_res = mr(har_data)) # calculate OR for dichotomous variables

#sensitivity analysis
res_sen <- mr_heterogeneity(har_data)
#multivalency analysis
res_ple <- mr_pleiotropy_test(har_data) 

write.csv(res_mr, "./mr_results.csv")
write.csv(res_sen, "./mr_sensitivity.csv")
write.csv(res_ple, "./mr_pleiotropy.csv")

#visualization
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(har_data))
ggsave("./mr_leaveoneout.png", plot = last_plot(), dpi = 300)
mr_scatter_plot(res_mr, har_data)
ggsave("./mr_scatter.png", plot = last_plot(), dpi = 300)
mr_forest_plot(singlesnp_results = mr_singlesnp(har_data))
ggsave("./mr_forest.png", plot = last_plot(), dpi = 300)
mr_funnel_plot(singlesnp_results = mr_singlesnp(har_data))
ggsave("./mr_funnel.png", plot = last_plot(), dpi = 300)

rm(list = ls())
