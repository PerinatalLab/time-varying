# In R, use snakemake@input[[1]] etc to get the variables
# on harvest:
# create snplist-alltop.txt, snplist-xtop.txt from meta results
# /mnt/work2/jjuod/scripts/extract-snps-moba30k.sh snplist-alltop.txt top1
# /mnt/work2/jjuod/scripts/extract-snps-moba30k.sh snplist-xtop.txt top1x
# gzip both

rule all:
	"Create all target files."
	input:
		"/home/julius/Documents/results/tv/table_main.tsv",
		"/home/julius/Documents/results/tv/plot_allmain.png"

rule preliminary:
	"Create files for initial analyses. No real need to run this."
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
		"/mnt/HARVEST/mobagen-flaglist-n99259.txt",
		"/mnt/HARVEST/PDB1724_Q1_v12.csv"
	output:
		"/mnt/HARVEST/ga_cleaned.csv"
	script:
		"clean-pheno.R"

rule analyze_all_main:
	output:
		maintable="/home/julius/Documents/results/tv/table_main.tsv",
		mainplot="/home/julius/Documents/results/tv/plot_allmain.png"
	input:
		mfr="/mnt/HARVEST/ga_cleaned.csv",
		gt="/mnt/HARVEST/top1-moba30k-dosage.csv.gz",
		gtX="/mnt/HARVEST/top1x-moba30k-dosage.csv.gz",
		mobares="/mnt/HARVEST/topsnps_moba_summaries.txt"
	script:
		"run-tvmodels-all.R"
		
rule prel_analyze_tv:
	"Run the preliminary reports for top SNPs in TV models (autosomal SNPs only)."
	input:
		mfr="/mnt/HARVEST/ga_cleaned.csv",
		gt="/mnt/HARVEST/top1-moba30k-dosage.csv.gz",
		moba="/mnt/HARVEST/topsnps_moba_summaries.txt"
	params:
		outstem="/home/julius/Documents/results/tv/report_tvmodels"
	output:
		expand("/home/julius/Documents/results/tv/report_tvmodels{i}.pdf", i=range(1,24))
	shell:  # This is daft but very important to enforce serial Rmd knitting!!!
		""" Rscript -e "for(i in 1:23){{ rmarkdown::render('run-tvmodels.Rmd', params=list(mfrfile='{input.mfr}', gtfile='{input.gt}', mobaresfile='{input.moba}', i=i), output_format='pdf_document', output_file=paste0('{params.outstem}', i, '.pdf')) }}"
		"""

rule prel_analyze_tvX:
	"Run the preliminary reports for top SNPs in TV models (X chr. SNPs only)."
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

# THIS IS RUN ON HARVEST ONLY
# rule extract_gt:
# 	"Extract the selected SNPs from complete genotyping data."
# 	input:
# 		inauto="../snplist-alltop.txt",
# 		inx="../snplist-alltop.txt",
# 		GENO_DATA_FILE
# 	output:
# 		"/mnt/HARVEST/top1-moba30k-dosage.csv.gz",
# 		"/mnt/HARVEST/top1x-moba30k-dosage.csv.gz"
# 	shell:
# 		"./extract-snps-moba30k.sh {inauto} top1 ; ./extract-snps-moba30k.sh {inx} top1x"
