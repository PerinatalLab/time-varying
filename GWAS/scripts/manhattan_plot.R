library(dplyr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(cowplot)

HC= -log10(5e-8)

manhattan_prep = function(dat){
  colnames(dat)[1] = "CHROM"
  dat$CHROM = as.numeric(ifelse(dat$CHROM=="X",23,dat$CHROM))
  dat = dat[,c("CHROM", "POS", "RSID", "P")]
  
  don <- dat %>%
    group_by(CHROM) %>%
    summarize(chr_len= max(POS)) %>%
    mutate(tot = cumsum(as.numeric(chr_len))-chr_len) %>% # Calculate cumulative position of each chromosome
    select(-chr_len) %>%
    right_join(dat, by="CHROM") %>%
    mutate(BPcum=POS+tot) %>% # Add a cumulative position of each SNP
    ungroup()
  
  don = don %>% mutate(plotcolor = ifelse(CHROM %% 2 != 0 , "a","b"))
  
  return(don)
}

manhattan_plot = function(dat){
  axisdf = dat %>% group_by(CHROM) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 ) %>%
    mutate(labels=ifelse(CHROM %in% c(20,22), '', ifelse(CHROM==23, 'X', CHROM)))
  
  p = ggplot(dat, aes(x= BPcum, y= -log10(P), color = plotcolor)) +
    geom_point(size= 0.07) +   # Show all points
    geom_hline(yintercept= 0, size= 0.25, colour= 'black') +
    geom_hline(yintercept= HC, size= 0.2, linetype= 2, colour= '#878787') +
    theme_bw() +
    scale_colour_manual(values= c("#2c7fb8","#7fcdbb"), guide= "none") +
    scale_x_continuous(expand= c(0, 0.3),label= axisdf$labels, breaks= axisdf$center) +
    scale_y_continuous(expand= c(0, 0.05), limits= c(0, max(-log10(don$P)) + 2), breaks= seq(0, 10, 5)) +
    ylab('-log10(pvalue)') +
    xlab('Chromosome') +
    coord_cartesian(clip = "off") +
    theme(plot.margin = unit(c(t= 0, r=0, b= 1, l=0), 'cm'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  return(p)
}

# FETAL manhattans (Fig S1)

### ---- plot GA
topsnp = fread("GWAS/top_snp/top_snp_ga_cleaned_f.csv")
dat = fread(snakemake@input[[1]])
# dat = fread("/mnt/HARVEST/GWAS/ga_cleaned_f_GD_filterd.txt")

don = manhattan_prep(dat)

# attach Manhattan x coordinates to top snps
topsnpGA = left_join(topsnp, don, by=c("RSID", "POS", "CHR"="CHROM"))

p1 = manhattan_plot(don)

p1.final = p1 +
  geom_point(data=topsnpGA, color="#f60000", size=0.14) +
  geom_label_repel(data=mutate(topsnpGA, P=pmin(1, P*1.2)), aes(label=nearestGene),
                   ylim=c(-Inf,-0.9), fontface="italic",
                   arrow=arrow(angle=20, length=unit(0.02, "npc")), 
                   color="black", segment.color="#f60000")

### ---- same for PTD

dat = fread(snakemake@input[[2]])
# dat = fread("/mnt/HARVEST/GWAS/ga_cleaned_f_PTD_filterd.txt")

don = manhattan_prep(dat)
topsnpP = left_join(topsnp, don, by=c("RSID", "POS", "CHR"="CHROM"))
p2 = manhattan_plot(don)
p2.final = p2 +
  geom_point(data=topsnpP, color="#f60000", size=0.14) +
  geom_label_repel(data=mutate(topsnpP, P=pmin(1, P*1.2)), aes(label=nearestGene),
                   ylim=c(-Inf,-0.9), fontface="italic",
                   arrow=arrow(angle=20, length=unit(0.02, "npc")), 
                   color="black", segment.color="#f60000")
pall = plot_grid(p1.final, p2.final, nrow=2, labels="AUTO")

# ggsave("~/Documents/results/tv/plot_manh_s1.png", pall, width=15, height=14, units="cm")
ggsave(snakemake@output[[1]], pall, width = 15, height = 14, units="cm")

rm(pall, p1.final, p2.final, p1, p2, dat, don)
gc()


# MATERNAL manhattans (Fig 1)

# --- GA
topsnpGA = fread("GWAS/top_snp/top_snp_ga_cleaned_GD.csv")
dat = fread(snakemake@input[[3]])
# dat = fread("/mnt/HARVEST/GWAS/ga_cleaned_GD_filterd.txt")

don = manhattan_prep(dat)

# attach Manhattan x coordinates to top snps
topsnpGA$CHR = as.numeric(ifelse(topsnpGA$CHR=="X",23,topsnpGA$CHR))
topsnpGA = left_join(topsnpGA, don, by=c("RSID", "POS", "CHR"="CHROM"))
topsnpGA.bot = top_n(topsnpGA, 12, P)
topsnpGA.right = anti_join(topsnpGA, topsnpGA.bot) %>% filter(BPcum>620e6)
topsnpGA.left = anti_join(topsnpGA, topsnpGA.bot) %>% filter(BPcum<620e6)

p1 = manhattan_plot(don)

p1.final = p1 +
  geom_point(data=topsnpGA, color="#f60000", size=0.2) +
  geom_label_repel(data=mutate(topsnpGA.bot, P=pmin(1, P*1.2)), aes(label=nearestGene),
                   ylim=c(-Inf,-1.6), fontface="italic", force=2,
                   segment.size=0.3, arrow=arrow(angle=20, length=unit(0.02, "npc")), 
                   color="black", segment.color="#f60000") +
  geom_label_repel(data=mutate(topsnpGA.left, P=pmin(1, P*0.8)), aes(label=nearestGene),
                   ylim=c(9,13), xlim=c(NA, 600e6), fontface="italic",
                   segment.size=0.3, arrow=arrow(angle=20, length=unit(0.02, "npc")), 
                   color="black", segment.color="#f60000") +
  geom_label_repel(data=mutate(topsnpGA.right, P=pmin(1, P*0.8)), aes(label=nearestGene),
                   ylim=c(8,13), xlim=c(660e6, NA), fontface="italic",
                   segment.size=0.3, arrow=arrow(angle=20, length=unit(0.02, "npc")), 
                   min.segment.length=0.1,
                   color="black", segment.color="#f60000") +
  theme(plot.margin = unit(c(t=0, r=0, b=1.8, l=0), 'cm'))

p1.final

# --- PTD
topsnpPTD = fread("GWAS/top_snp/top_snp_ga_cleaned_PTD.csv")
dat = fread(snakemake@input[[4]])
# dat = fread("/mnt/HARVEST/GWAS/ga_cleaned_PTD_filterd.txt")

# Adjusted alpha values, precalculated in run-tvmodels-all.R:
alphas = data.frame(y=c(0.0558, 0.00749), power=c("80 %", "50 %"))

don = manhattan_prep(dat)

p2 = manhattan_plot(don)

topsnpPTD = left_join(topsnpPTD, don, by=c("RSID", "POS", "CHR"="CHROM"))
topsnpPTD.bot = top_n(topsnpPTD, 6, P)
topsnpPTD.top = anti_join(topsnpPTD, topsnpPTD.bot)

p2.final = p2 +
  geom_hline(data=alphas, aes(yintercept=-log10(y)), linetype="dotted", col="cornsilk1", size=0.5) +
  geom_text(data=alphas, aes(y=-log10(y), x=0, label=power),
            vjust="bottom", hjust="left", col="cornsilk1", size=3) +
  geom_point(data=topsnpPTD, color="#f60000", size=0.2) +
  geom_label_repel(data=mutate(topsnpPTD.bot, P=pmin(1, P*1.2)), aes(label=nearestGene),
                   ylim=c(-Inf,-1.1), fontface="italic", force=2,
                   segment.size=0.3, arrow=arrow(angle=20, length=unit(0.02, "npc")), 
                   color="black", segment.color="#f60000") +
  geom_label_repel(data=mutate(topsnpPTD.top, P=pmin(1, P*0.8)), aes(label=nearestGene),
                   ylim=c(7.7,13), xlim=c(660e6, NA), fontface="italic",
                   segment.size=0.3, arrow=arrow(angle=20, length=unit(0.02, "npc")), 
                   min.segment.length=0.1,
                   color="black", segment.color="#f60000") +
  theme(plot.margin = unit(c(t=0, r=0, b=0.8, l=0), 'cm'))

p2.final

pall = plot_grid(p1.final, p2.final, nrow=2, labels="AUTO", rel_heights=c(0.54,0.46))

ggsave("~/Documents/results/tv/plot_manh_main.png", pall, width=20, height=18, units="cm")
# ggsave(snakemake@output[[2]], pall, width = 20, height = 18, units="cm")

library(png)
library(magick)
plz = readPNG("~/Documents/results/tv/locuszoom_sm.png")

pbetas = 

plot_grid(ggdraw() + draw_image(plz), pbetas, labels=c("C", "D"))
