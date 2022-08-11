library(dplyr)
library(data.table)
library(ggplot2)

dat = fread(snakemake@input[[1]])
dat2 = fread(snakemake@input[[2]])

## PTD
# Sig level
dat = dat %>% filter(P<1*10**-5) %>% arrange(P) %>% group_by('#CHROM')

final_pos_list = c()
chrom_list = unique(ifelse(dat$"#CHROM"=="X",23,dat$"#CHROM") )

# Extract sig snp, independent regions for top loci with a radius of 1.5 Mb, based on chromosome
for (chrom in chrom_list) {
        
        chromosome = dat[dat$"#CHROM" == ifelse(chrom == 23, "X",chrom),]
        positions = chromosome$POS
        min_P = chromosome %>% arrange(P) %>% filter(row_number()==1) %>% pull(POS)
        final_pos_list = rbind(final_pos_list, c(min_P,chrom) )
        for (pos in positions) {
		if ( ( (!pos %in% final_pos_list[,1] | (pos %in% final_pos_list[,1] & !chrom %in% as.numeric(final_pos_list[,2])) )  & (pos < as.numeric(final_pos_list[,1]) - (1.5*10**6) | pos > as.numeric(final_pos_list[,1]) + (1.5*10**6)) ) ) {
			final_pos_list = rbind(final_pos_list, c(pos,chrom))
                }
        }
}

#SNP to plot
colnames(final_pos_list) = c("POS","#CHROM")
final_pos_list = as.data.frame(final_pos_list)
final_pos_list$"#CHROM" = ifelse(final_pos_list$"#CHROM"==23,"X",final_pos_list$"#CHROM")

final_pos_list$POS = as.numeric(final_pos_list$POS)
dat$POS = as.numeric(dat$POS)
dat2$POS = as.numeric(dat2$POS)

sig_snp = inner_join(dat, final_pos_list, by = c('POS','#CHROM'))
sig_snp2 = inner_join(dat2, final_pos_list, by = c('POS','#CHROM'))

sig_snp$Beta = log(sig_snp$OR) # convert OR for sig_snp

b = inner_join(sig_snp, sig_snp2, by = "POS" ) # 165 snps


## GD
dat = fread(snakemake@input[[1]])
dat2 = fread(snakemake@input[[2]])

# Sig level
dat2 = dat2 %>% filter(P<1*10**-5) %>% arrange(P) %>% group_by('#CHROM')

final_pos_list = c()
chrom_list = unique(ifelse(dat2$"#CHROM"=="X",23,dat2$"#CHROM") )

# Extract sig snp, independent regions for top loci with a radius of 1.5 Mb, based on chromosome
for (chrom in chrom_list) {
        chromosome = dat2[dat2$"#CHROM" == ifelse(chrom == 23, "X",chrom),]
        positions = chromosome$POS
        min_POS = chromosome %>% arrange(P) %>% filter(row_number()==1) %>% pull(POS)
        final_pos_list = rbind(final_pos_list, c(min_POS,chrom) )

        for (pos in positions) {
                if ( ( (!pos %in% final_pos_list[,1] | (pos %in% final_pos_list[,1] & !chrom %in% as.numeric(final_pos_list[,2])) )  & (pos < as.numeric(final_pos_list[,1]) - (1.5*10**6) | pos > as.numeric(final_pos_list[,1]) + (1.5*10**6)) ) ) {
                       final_pos_list = rbind(final_pos_list, c(pos,chrom))

                }
        }
}

#SNP to plot
colnames(final_pos_list) = c("POS","#CHROM")
final_pos_list = as.data.frame(final_pos_list)
final_pos_list$"#CHROM" = ifelse(final_pos_list$"#CHROM"==23,"X",final_pos_list$"#CHROM")

final_pos_list$POS = as.numeric(final_pos_list$POS)
dat$POS = as.numeric(dat$POS)
dat2$POS = as.numeric(dat2$POS)


sig_snp = inner_join(dat, final_pos_list, by = c('POS','#CHROM'))
sig_snp2 = inner_join(dat2, final_pos_list, by = c('POS','#CHROM'))


# convert OR for sig_snp
sig_snp$Beta = log(sig_snp$OR)

c = inner_join(sig_snp, sig_snp2, by = "POS" ) #227 snps


## Plot
# Join top variants from ptd and gd
a = rbind(b, c)

#plot
p = ggplot(a, aes(x=Beta, y=BETA)) + 
	geom_point(pch=19,color="blue",size=3) +
       	ylab("Beta (Gestational duration)") + 
	xlab("log OR(Preterm delivery)") + 
	ggtitle("Effect direction of top loci for preterm delivery and gestational duration") + 
	geom_abline(intercept = 0, slope = -1, alpha = .5) + 
	geom_hline(yintercept = 0, alpha = .5) +
        geom_vline(xintercept = 0, alpha = .5)

ggsave(snakemake@output[[1]],p)



