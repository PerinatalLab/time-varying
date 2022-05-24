# In R, use snakemake@input[[1]] etc to get the variables

rule prep_pheno:
	"Clean MFR, create a censoring indicator, export selected columns."
	input:
		"/mnt/HARVEST/PDB1724_MFR_541_v12.csv",
		"/mnt/HARVEST/parental_ID_to_PREG_ID.csv",
		"/mnt/HARVEST/linkage_Mother_PDB1724.csv"
	output:
		"/mnt/HARVEST/ga_cleaned.csv"
	script:
		"clean-pheno.R"
		
