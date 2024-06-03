![image](https://github.com/BioinfoEasy/MPGK/assets/162677537/540fd00b-62dd-4587-8fc9-b05374cb5ddf)![image](https://github.com/BioinfoEasy/MPGK/assets/162677537/6e16701f-027d-4ae5-8e31-24cd247b8683)MPGK is a tool that enables MR, PRS, GO and KEGG analysis. MPGK is built on R and some command line languages. You can run the MPGK program from the command line.

# Installation of the MPGK programme
MPGK is currently available for Linux and Mac_OS, and you can download it from https://github.com/BioinfoEasy/MPGK . Take MPGK_Mac as an example, the folder contains Rcode, lib, PRSice_mac, plink, PRSice.R, and MPGK.sh. Among them, MPGK.sh is the main runtime program. Before using MPGK, please make sure that R language is installed in your computer and test whether you can run PRSice, plink programme normally. 

MPGK calls the following R package: TwoSampleMR, ieugwasr, clusterProfiler, org.Hs.eg.db, tidyverse, stringr, dplyr, openxlsx, Hmisc, data.table, ggplot2, enrichplot and RcolorBrewer. Before using the MPGK program, ensure that the required R packages have been installed. It is important to note that the MPGK program relies on these packages and will not function properly without them.

# Data format transformation
MPGK requires the following input files: a GWAS summary data file of the exposure in txt format and  a GWAS summary data file of the outcome in txt format for MR; a gene sequencing file in vcf format and a phenotype file in txt format for PRS; a GWAS summary data file of the disease in txt format for GO and KEGG analyses.

The GWAS summary data file should include the following columns: chr, pos, ref, alt, beta, pval, eaf, and SNP.
chr pos ref alt beta se pval eaf samplesize N SNP
1 49298 T C 9.20724e-05 0.000392468 0.809999964789546 0.623763 462933 5314 rs10399793
1 54676 C T -7.52885e-05 0.000388812 0.849999949672056 0.400401 462933 5314 rs2462492
1 86028 T C -0.000465849 0.000621635 0.450000503808422 0.103556 462933 5314 rs114608975
1 91536 G T -9.69521e-05 0.00038284 0.800000023961726 0.456851 462933 5314 rs6702460

The phenotype file should consist of four columns: FID, IID, gender, and phenotype.
1 N0025 1 2
2 N0037 2 0
3 N0026 1 0
4 N0038 1 0
5 N0027 2 2

# Running MPGK
The MPGK program includes the following parameters: M for MR, P for PRS, G for GO, and K for KEGG. Additionally, it requires the following input files: an exposure file for MR (specified with the letter e), an outcome file for MR (specified with the letter o), a VCF file of the gene sequencing (specified with the letter v), a phenotype file (specified with the letter f), and GWAS summary data of the disease for which GO and KEGG analyses were performed (specified with the letter d). Finally, a p-value threshold must be specified with the letter p. For assistance, type 'bash MPGK.sh -h' at the shell prompt.

## Help
% bash MPGK.sh -h
usage: MPGK.sh -M<MR> -P<PRS> -G<GO> -K<KEGG>
-e	GWAS summary txt file of exposure.
-o	GWAS summary txt file of outcome.
-p	threshold of p-value.
-d	GWAS summary txt file of disease.
-f	txt file of phenotype.
-v	vcf file of disease.

## MR
% bash MPGK.sh -M -e exposure.txt -o outcome.txt -p 5e-8
MPGK automatically generates a folder called MR, which includes the funnel plot, the forest plot, the scatter plot, and the leave-one-out plot, as well as the results of MR tests, sensitivity, and pleiotropy tests.

## PRS
% bash MPGK.sh -P -v diease.vcf -f phenotype.txt -p 5e-8
MPGK automatically generates a folder called PRS, which includes PRSice.best, PRSice.prsice, PRSice.summary, bar graphs, high-resolution plots, and probability distributions of PRS scores for the case and control groups.

## GO and KEGG
% bash MPGK.sh -G -d disease.txt -p 5e-8
% bash MPGK.sh -K -d disease.txt -p 5e-8
MPGK generates folders named GO and KEGG respectively. the GO folder contains bar plot, dot plot, DAG plot and map plot for BP, CC and MF. the KEGG folder contains bar plot and dot plot.

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
