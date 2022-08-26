library(dplyr)
library(data.table)
library(ggplot2)

dat = fread(snakemake@input[[1]])
dat2 = fread(snakemake@input[[2]])

dat$CHR = as.numeric(ifelse(dat$"#CHROM"=="X",23,dat$"#CHROM"))
dat2$CHR = as.numeric(ifelse(dat2$"#CHROM"=="X",23,dat2$"#CHROM"))

## PTD
# Sig level
dat_top = dat %>% filter(P<1*10**-5) %>% arrange(P)

final_pos_list = data.frame()
chrom_list = unique(dat_top$CHR)

# Extract sig snp, independent regions for top loci with a radius of 1.5 Mb, based on chromosome
for (chrom in chrom_list) {
  chromosome = filter(dat_top, CHR==chrom)
  positions = chromosome$POS
  min_P = positions[which.min(dat_top$P)]
  final_pos_list = bind_rows(final_pos_list, data.frame("POS"=min_P,"CHR"=chrom))
  for (pos in positions) {
    if ( !(pos %in% final_pos_list[,1] & chrom %in% final_pos_list[,2])){
      same_chrom = final_pos_list$POS[final_pos_list$CHR==chrom]
      if(min(abs(pos - same_chrom)) > 1.5e6){
        final_pos_list = bind_rows(final_pos_list, data.frame("POS"=pos,"CHR"=chrom))  
      }
    }
  }
}

#SNP to plot
sig_snp = inner_join(dat, final_pos_list, by = c('POS','CHR'))
sig_snp2 = inner_join(dat2, final_pos_list, by = c('POS','CHR'))

sig_snp$BETA_PTD = log(sig_snp$OR) # convert OR for sig_snp

dfb = inner_join(sig_snp, sig_snp2, by = c("POS", "CHR")) # 36 snps


## GD

# Sig level
dat2_top = dat2 %>% filter(P<1*10**-5) %>% arrange(P)

final_pos_list = data.frame()
chrom_list = unique(dat2_top$CHR)

# Extract sig snp, independent regions for top loci with a radius of 1.5 Mb, based on chromosome
for (chrom in chrom_list) {
  chromosome = filter(dat2_top, CHR==chrom)
  positions = chromosome$POS
  min_P = positions[which.min(dat2_top$P)]
  final_pos_list = bind_rows(final_pos_list, data.frame("POS"=min_P,"CHR"=chrom))
  for (pos in positions) {
    if ( !(pos %in% final_pos_list[,1] & chrom %in% final_pos_list[,2])){
      same_chrom = final_pos_list$POS[final_pos_list$CHR==chrom]
      if(min(abs(pos - same_chrom)) > 1.5e6){
        final_pos_list = bind_rows(final_pos_list, data.frame("POS"=pos,"CHR"=chrom))  
      }
    }
  }
}

#SNP to plot
sig_snp = inner_join(dat, final_pos_list, by = c('POS','CHR'))
sig_snp2 = inner_join(dat2, final_pos_list, by = c('POS','CHR'))

sig_snp$BETA_PTD = log(sig_snp$OR) # convert OR for sig_snp

dfc = inner_join(sig_snp, sig_snp2, by = c("POS", "CHR")) #227 snps


## Plot
# Join top variants from ptd and gd
a = bind_rows("PTD"=dfb, "GD"=dfc, .id="outcome")
write.table(a, snakemake@output[[1]], quote=F, row.names=F)

#plot
# p = ggplot(a, aes(x=Beta, y=BETA)) + 
# 	geom_point(pch=19,color="blue",size=3) +
#        	ylab("Beta (Gestational duration)") + 
# 	xlab("log OR(Preterm delivery)") + 
# 	ggtitle("Effect direction of top loci for preterm delivery and gestational duration") + 
# 	geom_abline(intercept = 0, slope = -1, alpha = .5) + 
# 	geom_hline(yintercept = 0, alpha = .5) +
#         geom_vline(xintercept = 0, alpha = .5)
# 
# ggsave(snakemake@output[[1]],p)
# 
# 
