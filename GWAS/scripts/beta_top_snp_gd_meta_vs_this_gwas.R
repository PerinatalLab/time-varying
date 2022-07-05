library(dplyr)
library(data.table)
library(ggplot2)

snp_meta = fread(snakemake@input[[1]])
dat = fread(snakemake@input[[2]])

a = dat[dat$RSID %in% snp_meta$RSID & dat$POS %in% snp_meta$POS,]

b = left_join(a, snp_meta, by="POS" )
b = b %>% mutate(BETA = ifelse(a1_plink != EFF, BETA*-1,BETA))

p = ggplot(b, aes(x=BETA, y=Betadays)) + 
	geom_point(pch=19,color="blue",size=3) + 
	xlab("Beta (this GWAS)") + 
	ylab("Beta (GWAS meta-analysis)") + 
	ggtitle("Beta of top loci of gestational duration in meta-analysis against MoBa GWAS") + 
	geom_abline(intercept = 0, slope = 1, alpha = .5) 

ggsave(snakemake@output[[1]])























