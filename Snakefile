# In R, use snakemake@input[[1]] etc to get the variables
# on harvest:
# create snplist-alltop.txt, snplist-xtop.txt from meta results
# create betas_all.txt with PGS weights from meta

rule all:
	"Create all target files."
	input:
		"/home/julius/Documents/results/tv/table_main.tsv",
		"/home/julius/Documents/results/tv/plot_allmain.png",
		"/home/julius/Documents/results/tv/plot_pgs5.png",
		"/home/julius/Documents/results/tv/plot_explor.png",
		"/home/julius/Documents/results/tv/plot_suppcov.png",
		"/home/julius/Documents/results/tv/plot_supphaz.png",
		"/home/julius/Documents/results/tv/table_diag.tsv",
		"/home/julius/Documents/results/tv/plot_coxdiag.png",
		"/home/julius/Documents/results/tv/plot_pgs3.png"

rule preliminary:
	"Create files for initial analyses. No real need to run this."
	input:
		expand("/home/julius/Documents/results/tv/report_tvmodels{i}.pdf", i=range(1,23)),
		expand("/home/julius/Documents/results/tv/report_tvmodelsX{i}.pdf", i=range(1,3)),
		"/home/julius/Documents/results/tv/null-pvals.RData"

rule harvest:
	"Things to be produced directly in HARVEST."
	input:
		"/mnt/work2/jjuod/tmp/ga_cleaned.csv",
		"/mnt/work2/jjuod/tmp/ga_cleaned_f.csv",
		"/mnt/work2/jjuod/tmp/top1-moba30k-dosage.csv.gz",
		"/mnt/work2/jjuod/tmp/top1x-moba30k-dosage.csv.gz",
		"/mnt/work2/jjuod/tmp/top1f-moba30k-dosage.csv.gz",
		"/mnt/work2/jjuod/tmp/plinktests/res_top_regions.sscore",
		"/mnt/work2/jjuod/tmp/plinktests/res_burden_below0.001_rare.sscore",
		"/mnt/work2/jjuod/tmp/plinktests/res_burden_below0.001_rec.sscore",
		"/mnt/work2/jjuod/tmp/plinktests/pca.eigenvec",
		"/mnt/work2/jjuod/tmp/plinktests/chr6-maf1-M.bim"

# -----------------------
# harvest jobs:
rule prep_pheno:
	"Clean MFR, create a censoring indicator, export selected columns."
	input:
		"/mnt/archive2/p1724/v12/PDB1724_MFR_541_v12.csv",
		"/mnt/archive2/p1724/v12/parental_ID_to_PREG_ID.csv",
		"/mnt/archive2/p1724/v12/linkage_Mother_PDB1724.csv",
		"/mnt/archive2/p1724/v12/linkage_Child_PDB1724.csv",
		"/mnt/archive/MOBAGENETICS/genotypes-base/aux/flaglist-merged/mobagen-flaglist-n99259.txt",
		"/mnt/archive2/p1724/v12/PDB1724_Skjema1_v12.csv"  # could pre-cut first 100 columns from this if slow
	output:
		"/mnt/work2/jjuod/tmp/ga_cleaned.csv",
		"/mnt/work2/jjuod/tmp/ga_cleaned_f.csv"
	script:
		"clean-pheno.R"

# after this, the GWAS and PGS parts can be run.

rule extract_gt:
	"Extract the selected SNPs from complete genotyping data."
	input:
		inauto="snplists/snplist-alltop.txt",
		inx="snplists/snplist-xtop.txt",
		inf="snplists/snplist-ftop.txt",
		chrs=expand("/mnt/archive/MOBAGENETICS/genotypes-base/imputed/all/vcf/{chr}.vcf.gz", chr=range(1,22)),
		chrx="/mnt/archive/MOBAGENETICS/genotypes-base/imputed/all/vcf/X.vcf.gz"
	output:
		"/mnt/work2/jjuod/tmp/top1-moba30k-dosage.csv.gz",
		"/mnt/work2/jjuod/tmp/top1x-moba30k-dosage.csv.gz",
		"/mnt/work2/jjuod/tmp/top1f-moba30k-dosage.csv.gz"
	shell:
		"./extract-snps-moba30k.sh {inauto} top1 ; ./extract-snps-moba30k.sh {inx} top1x ; ./extract-snps-moba30k.sh {inf} top1f"

rule explore_pgs_plink:
	"Prepare the stats needed for the PGS analysis R script."
	input:
		expand("/mnt/archive/MOBAGENETICS/genotypes-base/imputed/all/plink/{chr}.bed", chr=range(1,22)),
		"/mnt/archive/MOBAGENETICS/genotypes-base/imputed/all/plink/X.bed",
		"/mnt/work2/jjuod/tmp/ga_cleaned.csv",
		"/mnt/work2/jjuod/tmp/betas_all.txt"
	output:
		"/mnt/work2/jjuod/tmp/plinktests/res_top_regions.sscore",
		"/mnt/work2/jjuod/tmp/plinktests/res_burden_below0.001_rare.sscore",
		"/mnt/work2/jjuod/tmp/plinktests/res_burden_below0.001_rec.sscore",
		"/mnt/work2/jjuod/tmp/plinktests/pca.eigenvec",
		"/mnt/work2/jjuod/tmp/plinktests/chr6-maf1-M.bim"
	shell:
		"bash pgs-bash.sh"

# -----------------------

rule analyze_all_main:
	output:
		maintable="/home/julius/Documents/results/tv/table_main.tsv",
		mainplot="/home/julius/Documents/results/tv/plot_allmain.png",
		suppplot="/home/julius/Documents/results/tv/plot_suppcov.png",
		supphaz="/home/julius/Documents/results/tv/plot_supphaz.png"
	input:
		mfr="/mnt/HARVEST/ga_cleaned.csv",
		mfrF="/mnt/HARVEST/ga_cleaned_f.csv",
		gt="/mnt/HARVEST/top1-moba30k-dosage.csv.gz",
		gtX="/mnt/HARVEST/top1x-moba30k-dosage.csv.gz",
		gtF="/mnt/HARVEST/top1f-moba30k-dosage.csv.gz",
		pgs="/mnt/HARVEST/PGS.txt",
		mobares="snplists/topsnps_meta_summaries.txt"
	script:
		"run-tvmodels-all.R"

rule analyze_pgs:
	output:
		mainplotpgs="/home/julius/Documents/results/tv/plot_pgs5.png",
		suppplotpgs="/home/julius/Documents/results/tv/plot_pgs3.png",
		mainplotpca="/home/julius/Documents/results/tv/plot_explor.png"
	input:
		mfr="/mnt/HARVEST/ga_cleaned.csv",
		pgs="/mnt/HARVEST/PGS.txt",
		pgstop="/mnt/HARVEST/plinktests/res_top_regions.sscore",
		rarescore="/mnt/HARVEST/plinktests/res_burden_below0.001_rare.sscore",
		recscore="/mnt/HARVEST/plinktests/res_burden_below0.001_rec.sscore",
		pca="/mnt/HARVEST/plinktests/pca.eigenvec",
		chr6pos="/mnt/HARVEST/plinktests/chr6-maf1-M.bim"
	script:
		"analyse-pgs.R"
		
rule analyze_diagnostics:
	output:
		diagtable="/home/julius/Documents/results/tv/table_diag.tsv",
		coxdiagplot="/home/julius/Documents/results/tv/plot_coxdiag.png"
	input:
		mfr="/mnt/HARVEST/ga_cleaned.csv",
		mfrF="/mnt/HARVEST/ga_cleaned_f.csv",
		gt="/mnt/HARVEST/top1-moba30k-dosage.csv.gz",
		gtX="/mnt/HARVEST/top1x-moba30k-dosage.csv.gz",
		gtF="/mnt/HARVEST/top1f-moba30k-dosage.csv.gz",
		maintable="/home/julius/Documents/results/tv/table_main.tsv"
	script:
		"run-tvmodels-diag.R"

# ----------------------------S
# unused:
		
rule prel_analyze_tv:
	"Run the preliminary reports for top SNPs in TV models (autosomal SNPs only)."
	input:
		mfr="/mnt/HARVEST/ga_cleaned.csv",
		gt="/mnt/HARVEST/top1-moba30k-dosage.csv.gz",
		moba="snplists/topsnps_meta_summaries.txt"
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
		moba="snplists/topsnps_meta_summaries.txt"
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

