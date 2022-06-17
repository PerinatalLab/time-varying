#!/bin/bash
set -e

## extracts requested SNPs from MOBA30k.
## All samples will be used - not filtered by QC or by mothers/children.

## USAGE: extract-snps.sh INFILE.txt analysis_name
## INFILE.txt: no header, space separated, RSID - CHR - POS
## remember that MOBA uses GRCh37 coordinates!
## remember that X and 1-22 must NOT be mixed in one snplist!

if [[ "$#" -ne 2 ]]; then
	echo "wrong number of parameters, exiting"
	echo "USAGE: extract-snps.sh INFILE.txt analysis_name"
	exit
fi

export PATH=$PATH:/media/local-disk2/jjuod/erc-genotypes/bin/
GEND=/mnt/archive/MOBAGENETICS/genotypes-base/imputed/all/vcf
WD=/mnt/work2/jjuod/tmp
NAME=$2

# Warning: deletes any previous files with that name!
[ -f ${WD}/$NAME-moba30k-extracted.vcf ] && rm ${WD}/$NAME-moba30k-extracted.vcf

# create one long vcf with the extracted snps
while read -ra snp
do
	echo "looking for snp ${snp}"
	chr=${snp[1]}
	pos=${snp[2]}
	bcftools filter ${GEND}/${chr}.vcf.gz -r ${chr}:${pos} -Ov >> ${WD}/$NAME-moba30k-extracted.vcf
done < $1

# drop extra comment lines
awk '$0!~/^#/{h=1; print; next} h!=1{print}' ${WD}/$NAME-moba30k-extracted.vcf > ${WD}/$NAME-moba30k-extracted.vcf2

mv ${WD}/$NAME-moba30k-extracted.vcf2 ${WD}/$NAME-moba30k-extracted.vcf

# to produce a dosage table:
# first, make the headers
bcftools query -f 'ID CHROM POS REF ALT INFO [%SAMPLE ]\n' ${WD}/$NAME-moba30k-extracted.vcf | head -n1 > ${WD}/$NAME-moba30k-dosage.csv

# now extract the DS field (dosage)
bcftools query -f '%ID %CHROM %POS %REF %ALT %INFO [%DS ]\n' ${WD}/$NAME-moba30k-extracted.vcf >> ${WD}/$NAME-moba30k-dosage.csv

gzip ${WD}/$NAME-moba30k-dosage.csv
