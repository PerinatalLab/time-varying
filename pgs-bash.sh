# Script for creating various summary statistics for PGS analysis.
# NEED: plink, plink2 in path
set -e
GITDIR=$(dirname $(realpath "$0"))
SNPLIST=${GITDIR}/snplists/snplist-alltop.txt  # provided in git

WORKDIR=/mnt/work2/jjuod/tmp/plinktests/
cd ${WORKDIR}

GT=/mnt/archive/MOBAGENETICS/genotypes-base/imputed/all/plink/
MFRFILE=../ga_cleaned.csv
BETAS=../betas_all.txt   # this is the PGS scoring file from Pol's meta

# create maternal list from the GA phenotype file
awk -F';' -v OFS='\t' 'NR>1{print $11, $11}' ${MFRFILE} > maternal_ids.fam

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

plink --merge-list bfilelist.txt --extract ${BETAS} --make-bed --out all-pgs-snps


# plink2 below is 2.00 as indicated in the paper

# rare deleterious allele and genotype counts
awk '$3 < -0.001{print $1, $2, 1}' ${BETAS} > burden_below0.001.txt
plink2 --bfile all-pgs-snps \
	--max-maf 0.10 --score burden_below0.001.txt \
	--out res_burden_below0.001_rare
plink2 --bfile all-pgs-snps \
	--autosome --score burden_below0.001.txt recessive \
	--out res_burden_below0.001_rec

# pca
plink2 --bfile all-pgs-snps \
	--indep-pairwise 100 10 0.2 \
	--out all-pgs-pruning
plink2 --bfile all-pgs-snps \
	--extract all-pgs-pruning.prune.in \
	--pca allele-wts approx --threads 35 \
	--out pca

# do the extra PGS variance calculation
rm -f snpranges.txt
while read t_snp t_chr t_bp
do
	awk -v t_bp=$t_bp -v t_chr=$t_chr -v ORS="," 'start=="" && $1==t_chr && $4>t_bp-50000{start=$2} start!="" && $4>t_bp+50000{print start "-" $2; start=""; exit}' all-pgs-snps.bim >> snpranges.txt
done < ${SNPLIST}
plink2 --bfile all-pgs-snps --snps "$(<snpranges.txt)" --score ${BETAS} sum --out res_top_regions

