# In R, use snakemake@input[[1]] etc to get the variables
# on harvest:
# create snplist-alltop.txt, snplist-xtop.txt from meta results
# /mnt/work2/jjuod/scripts/extract-snps-moba30k.sh snplist-alltop.txt top1
# /mnt/work2/jjuod/scripts/extract-snps-moba30k.sh snplist-xtop.txt top1x
# gzip both

rule all:
	"Create all target files."
	input:
		expand("/home/julius/Documents/results/tv/report_tvmodels{i}.pdf", i=range(1,23)),
		expand("/home/julius/Documents/results/tv/report_tvmodelsX{i}.pdf", i=range(1,3)),
		"/home/julius/Documents/results/tv/null-pvals.RData"

rule prep_pheno:
	"Clean MFR, create a censoring indicator, export selected columns."
	input:
		"/mnt/HARVEST/PDB1724_MFR_541_v12.csv",
		"/mnt/HARVEST/parental_ID_to_PREG_ID.csv",
		"/mnt/HARVEST/linkage_Mother_PDB1724.csv",
		"/mnt/HARVEST/mobagen-flaglist-n99259.txt"
	output:
		"/mnt/HARVEST/ga_cleaned.csv"
	script:
		"clean-pheno.R"
		
rule analyze_tv:
	"Run the actual analyses of top SNPs in TV models (autosomal SNPs only)."
	input:
		mfr="/mnt/HARVEST/ga_cleaned.csv",
		gt="/mnt/HARVEST/top1-moba30k-dosage.csv.gz",
		moba="/mnt/HARVEST/topsnps_moba_summaries.txt"
	params:
		outstem="/home/julius/Documents/results/tv/report_tvmodels"
	output:
		expand("/home/julius/Documents/results/tv/report_tvmodels{i}.pdf", i=range(1,23))
	shell:  # This is daft but very important to enforce serial Rmd knitting!!!
		""" Rscript -e "for(i in 1:22){{ rmarkdown::render('run-tvmodels.Rmd', params=list(mfrfile='{input.mfr}', gtfile='{input.gt}', mobaresfile='{input.moba}', i=i), output_format='pdf_document', output_file=paste0('{params.outstem}', i, '.pdf')) }}"
		"""

rule analyze_tvX:
	"Run the actual analyses of top SNPs in TV models (X chr. SNPs only)."
	input:
		mfr="/mnt/HARVEST/ga_cleaned.csv",
		gt="/mnt/HARVEST/top1x-moba30k-dosage.csv.gz",  # note: this uses different sample size
		moba="/mnt/HARVEST/topsnps_moba_summaries.txt"
	params:
		outstem="/home/julius/Documents/results/tv/report_tvmodelsX"
	output:
		expand("/home/julius/Documents/results/tv/report_tvmodelsX{i}.pdf", i=range(1,3))
	shell:  # This is daft but very important to enforce serial Rmd knitting!!!
		""" Rscript -e "for(i in 1:2){{ rmarkdown::render('run-tvmodels.Rmd', params=list(mfrfile='{input.mfr}', gtfile='{input.gt}', mobaresfile='{input.moba}', i=i), output_format='pdf_document', output_file=paste0('{params.outstem}', i, '.pdf')) }}"
		"""

rule simulate_null:
	"Check if the pam tests maintain alpha by bootstrapping."
	input:
		mfr="/mnt/HARVEST/ga_cleaned.csv"
	output:
		"/home/julius/Documents/results/tv/null-pvals.RData"
	script:
		"simulate-null.R"
