library(dplyr)
library(data.table)
library(ggplot2)

dat = fread(snakemake@input[[1]])
topsnp = fread(snakemake@input[[2]])

colnames(dat)[1] = "CHROM"
dat$CHROM = as.numeric(ifelse(dat$CHROM=="X",23,dat$CHROM))

don <- dat %>%
    group_by(CHROM)      %>%
    summarise(chr_len= max(POS)) %>%
    mutate(tot= cumsum(as.numeric(chr_len))-chr_len) %>% # Calculate cumulative position of each chromosome
    select(-chr_len) %>%
    left_join(dat, ., by= 'CHROM') %>%
    arrange(CHROM, POS) %>% # Add a cumulative position of each SNP
    mutate(BPcum=POS+tot) %>%
         ungroup()

axisdf = don %>% group_by(CHROM) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  names(axisdf)= c('CHR', 'center')
HC= -log10(5*10**-8)

don = don %>% mutate(plotcolor = ifelse(CHROM %% 2 != 0 , "a","b")) %>% mutate(plotcolor = ifelse(POS %in% topsnp$POS & RSID %in% topsnp$RSID, "c",plotcolor))

p1= ggplot(don  %>% arrange(plotcolor), aes(x= BPcum, y= -log10(P), color = plotcolor)) +
  geom_point(size= 0.07) +   # Show all points
  theme_bw( ) +
  scale_colour_manual(values= c("#2c7fb8","#7fcdbb","#f60000"), guide= F) +
  scale_x_continuous(expand= c(0, 0.3),label = c(1:19, '', 21,'', 'X'), breaks= axisdf$center) + # label = ifelse(axisdf$CHR== 23, 'X', axisdf$CHR)
  scale_y_continuous(expand= c(0, 0.05), limits= c(0, max(-log10(don$P)) + 2), breaks= seq(0, 10, 5), labels= c(abs(seq(0, 10, 5)))) + # , sec.axis = sec_axis(~ ., name = derive())) +
  ylab('-log10(pvalue)') +
  xlab('Chromosome') +
  geom_hline(yintercept= 0, size= 0.25, colour= 'black') +
  geom_hline(yintercept= c(HC), size= 0.2, linetype= 2, colour= '#878787') +
  coord_cartesian(clip = "off") +
  theme(legend.position= 'none',
	plot.margin = unit(c(t= 0, r=0, b= 0, l=0), 'cm'),
        text= element_text(family="arial", size= 9),
	axis.line= element_line(size= 0.1),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank())

ggsave(snakemake@output[[1]], p1, width = 15, height = 7, units="cm")


