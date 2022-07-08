# NEED INPUT FILES:
# plink, plink2 in path
# betas_all being the PGS scoring file from Chris

# TODO
# betas_all=...
# ga_cleaned.csv=...
GT=/mnt/archive/MOBAGENETICS/genotypes-base/imputed/all/plink/

# create maternal list from the GA phenotype file
awk -F';' -v OFS='\t' 'NR>1{print $11, $11}' ga_cleaned.csv > maternal_ids.fam

# NOTE: old plink (1.90b6.26) was used for this basic merging here
for chr in {1..22}
do
	plink --bfile ${GT}/${chr} \
		--keep maternal_ids.fam \
		--maf 0.01 \
		--make-bed \
		--out chr${chr}-maf1-M
	echo chr${chr}-maf1-M >> bfilelist.txt
done
plink --bfile ${GT}/X \
       --keep maternal_ids.fam \
       --maf 0.01 \
       --make-bed \
       --out chrX-maf1-M
echo chrX-maf1-M >> bfilelist.txt

plink --merge-list bfilelist.txt --extract betas_all.txt --make-bed --out all-pgs-snps


# plink2 below is 2.00 as indicated in the paper
awk '$3 < -0.001{print $1, $2, 1}' betas_all.txt > burden_below0.001.txt
plink2 --bfile all-pgs-snps \
	--max-maf 0.10 --score burden_below0.001.txt \
	--out res_burden_below0.001_rare
plink2 --bfile all-pgs-snps \
	--autosome --score burden_below0.001.txt recessive \
	--out res_burden_below0.001_rec
plink2 --bfile all-pgs-snps \
	--indep-pairwise 100 10 0.2 \
	--out all-pgs-pruning
plink2 --bfile all-pgs-snps \
	--extract all-pgs-pruning.prune.in \
	--pca allele-wts approx --threads 35 \
	--out pca


