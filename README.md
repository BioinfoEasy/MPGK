MPGK is a tool that enables MR, PRS, GO and KEGG analysis. MPGK is built on R and some command line languages. You can run the MPGK programme from the command line.

# Installation of the MPGK programme
MPGK is currently available for Linux and Mac_OS, and you can download it from https://github.com/BioinfoEasy/MPGK . Take MPGK_Mac as an example, the folder contains Rcode, lib, PRSice_mac, plink, PRSice.R, and MPGK.sh. Among them, MPGK.sh is the main runtime program.

# 

The MPGK program includes the following parameters: M for MR, P for PRS, G for GO, and K for KEGG. Additionally, it requires the following input files: an exposure file for MR (specified with the letter e), an outcome file for MR (specified with the letter o), a VCF file of the gene sequencing (specified with the letter v), a phenotype file (specified with the letter f), and GWAS summary data of the disease for which GO and KEGG analyses were performed (specified with the letter d). Finally, a p-value threshold must be specified with the letter p. For assistance, type 'bash MPGK.sh -h' at the shell prompt.

These include TwoSampleMR, ieugwasr, clusterProfiler, org.Hs.eg.db, tidyverse, stringr, dplyr, openxlsx, Hmisc, data.table, ggplot2, enrichplot and RcolorBrewer. Before using the MPGK program, ensure that the required R packages have been installed. It is important to note that the MPGK program relies on these packages and will not function properly without them.

# demo data
demo1.txt: GWAS summary data of the psoriasis.
demo2.txt: GWAS summary data of the diabetes.
demo3.vcf: Gene sequencing results of the diabetes.
demo3_phe.txt: Phenotype data of diabetes.
demo4.txt: GWAS summary data of the psoriasis. The demo4 is accessible at this site: https://pan.baidu.com/s/1RQXH3mO8GqV9bcrEkel4lg?pwd=rtvx extraction code: rtvx.
The GWAS summary file should include the following columns: chr, pos, ref, alt, beta, pval, eaf, and SNP. The phenotype file should consist of four columns: FID, IID, gender, and phenotype.

# software
Computer workstations using MacOS, Unix or Linux operating systems
R (http://cran.r-project.org/) for data analysis and graphing
PLINK software (http://pngu.mgh.harvard.edu/~purcell/plink/download.shtml)
PRSice-2 software (http://prsice.info)
