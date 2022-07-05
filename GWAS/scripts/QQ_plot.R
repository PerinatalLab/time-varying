library(dplyr)
library(data.table)
library(ggplot2)

dat = fread(snakemake@input[[1]])

dat= arrange(dat, P) %>% mutate(exp1= -log10(1:length(P)/length(P)))

p1= ggplot(dat, aes(exp1, -log10(P))) +
  geom_point(size= 0.4, color= "#E69F00") +
  geom_abline(intercept = 0, slope = 1, alpha = .5) +
labs(colour="") +
theme_bw( ) +
xlab('Expected (-log10(p-value))') +
ylab('Observed (-log10(p-value))') +
theme(legend.position= 'bottom')

ggsave(snakemake@output[[1]],p1)


