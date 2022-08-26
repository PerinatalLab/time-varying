library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)

gc = data.frame("lambda"=c(1.0375, 1.0075,  1.0557, 1.0195),
                "LDint"=c(1.0197, 0.999,  1.0294, 1.0005))

out = expand.grid(c("GD", "PTD"), c("Maternal", "Foetal"))
out = out[,c("Var2", "Var1")]

ps = list()
for(i in 1:4){
  dat = fread(snakemake@input[[1]])
  dat = arrange(dat, P) %>% mutate(exp1= -log10(1:length(P)/length(P)))
  
  ps[[i]] = ggplot(dat, aes(exp1, -log10(P))) +
    geom_point(size= 0.4, color= "#E69F00") +
    geom_abline(intercept = 0, slope = 1, alpha = .5) +
    labs(colour="") +
    annotate("text", x=0, y=6, label=sprintf("GC lambda=%.3f\nLD intercept=%.3f", gc[i,])) +
    annotate("text", x=0, y=7, label=sprintf("%s GWAS of %s", out[i,]))
    theme_bw( ) +
    xlab('Expected (-log10(p-value))') +
    ylab('Observed (-log10(p-value))') +
    theme(legend.position= 'bottom')
}

pall = plot_grid(plotlist = ps, ncol=2, labels="AUTO")
ggsave(snakemake@output[[1]],p1)


