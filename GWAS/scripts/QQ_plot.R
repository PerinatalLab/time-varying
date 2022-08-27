library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)


gc = data.frame("lambda"=c(1.0375, 1.0075,  1.0557, 1.0195),
                "LDint"=c(1.0197, 0.999,  1.0294, 1.0005))

out = expand.grid("outcome"=c("GD", "PTD"), "geno"=c("Maternal", "Foetal"))

ps = list()

inputs = snakemake@input
# inputs = c("/mnt/HARVEST/GWAS/ga_cleaned_GD_filterd.txt",
#            "/mnt/HARVEST/GWAS/ga_cleaned_PTD_filterd.txt",
#            "/mnt/HARVEST/GWAS/ga_cleaned_f_GD_filterd.txt",
#            "/mnt/HARVEST/GWAS/ga_cleaned_f_PTD_filterd.txt")

for(i in 1:4){
  dat = fread(inputs[[i]])
  dat = arrange(dat, P) %>% mutate(exp1= -log10(1:length(P)/length(P)))
  
  ps[[i]] = ggplot(dat, aes(exp1, -log10(P))) +
    geom_point(size= 0.4, color= "#E69F00") +
    geom_abline(intercept = 0, slope = 1, alpha = .5) +
    labs(colour="") +
    annotate("text", x=0, y=6, hjust="left",
             label=sprintf("GC lambda=%.3f\nLD intercept=%.3f", gc$lambda[i], gc$LDint[i])) +
    annotate("text", x=0, y=7, hjust="left",
             label=sprintf("%s GWAS of %s", out$geno[i], out$outcome[i])) +
    theme_bw() + 
    xlab('Expected (-log10(p-value))') +
    ylab('Observed (-log10(p-value))') +
    theme(legend.position= 'bottom')
}

pall = plot_grid(plotlist = ps, ncol=2, labels="AUTO")
ggsave(snakemake@output[[1]],pall, width=20, height=20, units="cm")
# ggsave("~/Documents/results/tv/plot_qq_supp.png",pall, width=20, height=20, units="cm")

