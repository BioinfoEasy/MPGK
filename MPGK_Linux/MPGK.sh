#!/bin/bash

workspace="$(pwd)"

exposure="Unknown"
outcome="Unknown"
p_value="Unknown"
disease="Unknown"
phenotype="Unknown"
vcf="Unknown"
M=0
P=0
G=0
K=0

while getopts "hMPGKe:o:p:d:f:v:" arg
do
case $arg in
	h)
	printf "usage: ${0##*/} -M<MR> -P<PRS> -G<GO> -K<KEGG>\n"
	printf "%s\t%s\n" "-e" "GWAS summary txt file of exposure."
	printf "%s\t%s\n" "-o" "GWAS summary txt file of outcome."
	printf "%s\t%s\n" "-p" "threshold of p-value."
	printf "%s\t%s\n" "-d" "GWAS summary txt file of disease."
	printf "%s\t%s\n" "-f" "txt file of phenotype."
	printf "%s\t%s\n" "-v" "vcf file of disease."
	exit 0;;
	M)
	M=1;;
	P)
	P=1;;
	G)
	G=1;;
	K)
	K=1;;
	e)
	exposure=$OPTARG;;
	o)
	outcome=$OPTARG;;
	p)
	p_value=$OPTARG;;
	d)
	disease=$OPTARG;;
	f)
	phenotype=$OPTARG;;
	v)
	vcf=$OPTARG;;
	?)
	printf "error: unknown argument\nuse -h option to see help information.\n"
	exit 1;;
esac
done

if [ $M -eq 1 ]
then
	Rscript ./Rcode/MR.R $workspace $exposure $outcome $p_value
	mkdir ./MR
	rm -f exp_filter.csv out_merged.csv Rplots.pdf
	mv mr_* ./MR
fi

if [ $P -eq 1 ]
then
	./plink --vcf $vcf --allow-extra-chr --make-bed -out dis
	Rscript ./Rcode/pheno_input.R $workspace $phenotype 
	./plink --bfile dis --recode --out d
	./plink --file d --mind 0.1 --geno 0.1 --hwe ${p_value} --maf 0.01 --recode --out dd
	awk '{print $1,$2,$6}' dd.ped > phe.txt
	Rscript ./Rcode/select.R $workspace
	./plink --file dd --keep name_base.txt --recode --out base_data
	./plink --file dd --keep name_target.txt --recode --out target_data
	awk '{print $1,$2,$5}' base_data.ped > cov1.txt
	awk '{print $1,$2,$6}' base_data.ped > phe1.txt
	./plink --file base_data --pca 3
	awk '{print $3,$4,$5}' plink.eigenvec > pca1.txt
	paste cov1.txt pca1.txt | sed 's/\s\+/ /g' > covar1.txt
	./plink --file base_data --pheno phe1.txt --logistic --covar covar1.txt --out re1 --hide-covar 
	awk '{print $1,$2,$5}' target_data.ped > cov2.txt
	awk '{print $1,$2,$6}' target_data.ped > phe2.txt
	./plink --file target_data --pca 3
	awk '{print $3,$4,$5}' plink.eigenvec > pca2.txt
	paste cov2.txt pca2.txt | sed 's/\s\+/ /g' > covar2.txt
	awk '{$6=-9;print $0}' target_data.ped > tmp.ped 
	mv tmp.ped target_data.ped
	./plink --file target_data --make-bed --out target_bi
	Rscript PRSice.R --dir . --prsice ./PRSice_linux --base re1.assoc.logistic --target target_bi --thread 1 --stat OR --binary-target T --pheno phe2.txt
	Rscript PRSice.R --dir . --prsice ./PRSice_linux --base re1.assoc.logistic --target target_bi --thread 1 --stat OR --binary-target T --pheno phe2.txt --extract PRSice.valid	
	awk '{print $1,$2,$4}' PRSice.best > prs_phe.txt
	Rscript ./Rcode/prs_distribute.R $workspace $phenotype 
	mkdir ./PRS
	mv PRSice.best PRSice.prsice PRSice.summary PRS*.png ./PRS
	rm -f dis.* d.* dd.* re1.* base_data.* target_data.* target_bi.* plink.*
	rm -f cov1.txt cov2.txt covar1.txt covar2.txt pca1.txt pca2.txt phe.txt phe1.txt phe2.txt
	rm -f name_base.txt name_target.txt prs_phe.txt
	rm -f PRSice.log PRSice.valid PRSice.mismatch
fi

if [ $G -eq 1 ]
then
	mkdir ./GO
	Rscript ./Rcode/GO_KEGG.R $workspace $disease $p_value
	mv CC_*.png BP_*.png MF_*.png ./GO
	if [ $K -eq 1 ]
	then
		mkdir ./KEGG
		mv KEGG_*.png ./KEGG
	fi
fi

if [ $K -eq 1 ] && [ $G -eq 0 ]
then
	mkdir ./KEGG
	Rscript ./Rcode/GO_KEGG.R $workspace $disease $p_value
	mv KEGG_*.png ./KEGG
fi
