# PGS analyses

library(dplyr)
library(ggplot2)
library(mgcv)
library(pammtools)
library(cowplot)
library(tidyr)

# -------------------------------------------
# SETUP

# mfrfile="/mnt/HARVEST/ga_cleaned.csv"
# pgsfile = "/mnt/HARVEST/PGS.txt"
# pgstopfile = "/mnt/HARVEST/plinktests/res_top_regions.sscore"
# rarecountfile = "/mnt/HARVEST/plinktests/res_burden_below0.001_rare.sscore"
# raregenofile = "/mnt/HARVEST/plinktests/res_burden_below0.001_rec.sscore"
# pcafile = "/mnt/HARVEST/plinktests/pca.eigenvec"
# chr6posfile = "/mnt/HARVEST/plinktests/chr6-maf1-M.bim"
mfrfile = snakemake@input$mfr
pgsfile = snakemake@input$pgs
pgstopfile = snakemake@input$pgstop
rarecountfile = snakemake@input$rarescore
raregenofile = snakemake@input$recscore
pcafile = snakemake@input$pca
chr6posfile = snakemake@input$chr6pos

pcawtfile = paste0(pcafile, ".allele")

# some constants:
GA_START_TIME = 169

# make a palette with darker central color
pal5 = RColorBrewer::brewer.pal(5, "RdYlBu")
pal5[1] = "#FF410D"
pal5[2] = "#F7C530"
pal5[3] = "#ABB084"
# desaturated pal5 if need a bit calmer plotting
# pal5_desat = c("#FF9F85", "#FBE297", pal5[3], pal5[4], "#95BCDA")


# -------------------------------------------
# read in filtered phenotypes, MATERNAL
mfr_mid = read.table(mfrfile, h=T, sep=";")
nrow(mfr_mid)  # 26875

# Prepare covariates:
mfr_mid$PARITET_5[mfr_mid$PARITET_5==4] = 3  # parity honestly doesn't need 5 levels
mfr_mid$PARITET_5 = factor(mfr_mid$PARITET_5, levels=c(1,0,2,3))
# NOTE: using height will remove extra 1k missing values, so imputing:
mfr_mid$AA87[is.na(mfr_mid$AA87)] = mean(mfr_mid$AA87, na.rm=T)
# relevel to ensure the oldest and one of the largest batches is ref
mfr_mid$BATCH = factor(mfr_mid$BATCH, levels=c("M12A","FEB18","JAN15","JUN15",
                                               "M12B","M24","MAY16",
                                               "ROT1","ROT2"))

mfr_mid$MISD = factor(mfr_mid$MISD) # to prevent gg_slice complaining
mfr_mid$KJONN = factor(mfr_mid$KJONN)
mfr_mid$FAAR = mfr_mid$FAAR-2000  # center, for neater plot scales

# remove irrelevant timeframe for stability
mfr_mid$GAc = mfr_mid$SVLEN_DG-GA_START_TIME

# just some checks
range(mfr_mid$GAc)  # 1 139


# -------------------------------------------
# PGS. MAIN: 1) basic statistics

# merge maternal PGS w/ pheno data
pgs = read.table(pgsfile, h=T)
merged = inner_join(mfr_mid, pgs, by=c("SENTRIX_ID"="PREG_ID_1724"))
nrow(merged)  # 26875

# PGS categorized for visualization
PGSlabels = c("1st (shortest)", "2nd", "3rd", "4th", "5th (longest)")
merged$PGScat_rel = cut(merged$PGS, breaks=quantile(merged$PGS, c(0, 0.2, 0.4, 0.6, 0.8, 1)),
                    labels=PGSlabels, include.lowest=T)
merged$PGScat = factor(merged$PGScat_rel, levels=PGSlabels[c(3, 1,2,4,5)])  # set middle ref level


# variance explained: (about 0.0197)
summary(lm(GAc ~ PGS, data=merged[merged$hadevent,]))

# (slightly smaller r squared if all births included)
# summary(lm(GAc ~ PGS, data=merged))

# w/ covariates:
# (also checked w/ regressing out covariates and it doesn't change much)
summary(lm(GAc ~ PGS + BATCH + poly(MAGE,2) + FAAR + AA87 + KJONN + MISD + PARITET_5,
           data=merged[merged$hadevent,]))

# variance explained w/ a mini-PGS from top regions:
pgs_top = read.table(pgstopfile, h=T, comment.char = "")
merged2 = inner_join(merged, pgs_top[,c("IID", "SCORE1_AVG")], by=c("SENTRIX_ID"="IID"))
merged2$SCORE1_AVG = merged2$SCORE1_AVG * 2000 # just for my sanity
cor.test(merged2$PGS, merged2$SCORE1_AVG)
# variance explained: (about 0.006 & 0.014)
summary(lm(GAc ~ SCORE1_AVG, data=merged2[merged2$hadevent,]))
summary(lm(GAc ~ SCORE1_AVG + BATCH + poly(MAGE,2) + FAAR + AA87 + KJONN + MISD + PARITET_5,
           data=merged2[merged2$hadevent,]))


# 2) show overall time-varying effects of PGS groups
# create ped form for PAMs
ped = as_ped(merged, Surv(GAc, hadevent) ~ PGS + PGScat +
               BATCH + MAGE + FAAR + AA87 + KJONN + MISD + PARITET_5,
             id = "id", cut=c(0,seq(20, 137, by=7)))
nrow(ped)  # 389543

# run PAMM
mod.tv.pgs = bam(ped_status ~ ti(tend,bs='cr',k=11) + PGScat + ti(tend, by=as.ordered(PGScat),bs='cr') +
                   BATCH + poly(MAGE,2) + FAAR + AA87 + KJONN + MISD + PARITET_5,
                 data=ped, offset=offset, family=poisson())
summary(mod.tv.pgs)

# Plot partial effects
pred_df = make_newdata(ped, tend=unique(tend), PGScat=factor(PGSlabels)) %>%
  add_term(mod.tv.pgs, term="PGScat", se_mult=1.96) %>%
  mutate(tmid = 0.5*tstart+0.5*tend)

pmain = ggplot(pred_df, aes(x=(tmid+GA_START_TIME)/7, y=fit, col=PGScat, fill=PGScat)) +
  geom_vline(xintercept = 37, col="grey80") +
  geom_line(lwd=0.6) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),alpha=0.02,lwd=0.2,lty="dashed") +
  scale_color_manual(values=pal5, name="GA PGS quintile") +
  scale_fill_manual(values=pal5, name="GA PGS quintile") +
  scale_x_continuous(breaks=seq(23, 43, by=2), expand=c(0,0)) +
  ylim(c(-1.2, 2.2)) + 
  theme_bw() + xlab("gestational age, weeks") + ylab("log hazard ratio") +
  theme(legend.position=c(0.88, 0.75), legend.box.background= element_rect(colour="black"),
        panel.grid.minor.x=element_blank(), legend.title = element_text(size=10))
pmain

# Show the same with a basic, unadjusted binary test:
# summary(glm(SVLEN_DG<259 ~ PGScat + BATCH,
#             data=merged[merged$hadevent,], family="binomial"))  # only 1-2 are sign. vs 3
# summary(glm(SVLEN_DG<224 ~ PGScat + BATCH,
#             data=merged[merged$hadevent,], family="binomial"))  # only 1 is sign. vs 3

vecPropTest = Vectorize(prop.test, SIMPLIFY = F)
prisk = group_by(merged, PGScat_rel) %>%
  summarize(f_vptd=mean(SVLEN_DG<7*32), f_ptd=mean(SVLEN_DG<7*37), n=n()) %>%
  gather(key="pheno", value="f", f_vptd:f_ptd) %>%
  mutate(ci_lo = sapply(vecPropTest(f*n, n), function(x) x$conf.int[1]),
         ci_upp = sapply(vecPropTest(f*n, n), function(x) x$conf.int[2])) %>%
  mutate(pheno=factor(pheno, labels=c("PTD", "vPTD"))) %>%
  ggplot(aes(x=PGScat_rel, y=f*100)) +
  geom_col(aes(fill=PGScat_rel), alpha=0.5) + 
  geom_linerange(aes(ymin=ci_lo*100, ymax=ci_upp*100)) +
  facet_wrap(~pheno, scales="free_y") +
  xlab(NULL) + ylab("risk, %") +
  scale_fill_manual(values=pal5, name="GA PGS quintile", guide="none") +
  scale_y_continuous(expand=expansion(add=c(0, 0), mult=c(0,0.05))) +
  theme_bw() +
  theme(strip.background = element_rect(fill="#E7E8D3"),
        panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
prisk

plot_grid(pmain, prisk, labels="AUTO", nrow=2, rel_heights=c(3,2))

ggsave(snakemake@output$mainplotpgs, width=5.3, height=6, units="in")


# ----------------

# SUPPLEMENT: different cuts

## 7 quantiles
ped$PGScat = cut(ped$PGS, breaks=quantile(merged$PGS, (0:7)/7), labels=1:7)
ped$PGScat = factor(ped$PGScat, levels=c("4", "1", "2", "3", "5", "6", "7"))

mod.tv.pgs = bam(ped_status ~ ti(tend,bs='cr',k=11) + PGScat + ti(tend, by=as.ordered(PGScat),bs='cr') +
                   BATCH + poly(MAGE,2) + FAAR + AA87 + KJONN + MISD + PARITET_5,
                 data=ped, offset=offset, family=poisson())
summary(mod.tv.pgs)

pred_df = make_newdata(ped, tend=unique(tend), PGScat=factor(1:7)) %>%
  add_term(mod.tv.pgs, term="PGScat", se_mult=1.96) %>%
  mutate(tmid = 0.5*tstart+0.5*tend)

pal7 = RColorBrewer::brewer.pal(7, "RdBu")
pal7[4] = "grey60"

p1 = ggplot(pred_df, aes(x=(tmid+GA_START_TIME)/7, y=fit, col=PGScat, fill=PGScat)) +
  geom_vline(xintercept = 37, col="grey80") +
  geom_line(lwd=0.6) +
  # geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),alpha=0.04,lty=0) +
  scale_color_manual(values=pal7, name="PGS septile") +
  scale_fill_manual(values=pal7, name="PGS septile") +
  scale_x_continuous(breaks=seq(23, 43, by=2), expand=c(0,0)) +
  theme_bw() + xlab("gestational age, weeks") + ylab("log hazard ratio") +
  theme(legend.position=c(0.88, 0.75), legend.box.background= element_rect(colour="black"),
        legend.title = element_text(size=10))
p1

## 3 quantiles
ped$PGScat = cut(ped$PGS, breaks=quantile(merged$PGS, c(0, 0.333, 0.666, 1)),
                 labels=c("1st (shortest)", "2nd", "3rd (longest)"))
ped$PGScat = factor(ped$PGScat, levels=c("2nd", "1st (shortest)", "3rd (longest)"))

mod.tv.pgs = bam(ped_status ~ ti(tend,bs='cr',k=11) + PGScat + ti(tend, by=as.ordered(PGScat),bs='cr') +
                   BATCH + poly(MAGE,2) + FAAR + AA87 + KJONN + MISD + PARITET_5,
                 data=ped, offset=offset, family=poisson())
summary(mod.tv.pgs)

pred_df = make_newdata(ped, tend=unique(tend),
                       PGScat=factor(c("1st (shortest)", "2nd", "3rd (longest)"))) %>%
  add_term(mod.tv.pgs, term="PGScat", se_mult=1.96) %>%
  mutate(tmid = 0.5*tstart+0.5*tend)

p2 = ggplot(pred_df, aes(x=(tmid+GA_START_TIME)/7, y=fit, col=PGScat, fill=PGScat)) +
  geom_vline(xintercept = 37, col="grey80") +
  geom_line(lwd=0.6) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),alpha=0.05,lwd=0.2,lty="dashed") +
  scale_color_manual(values=pal5[c(1,3,5)], name="PGS tertile") +
  scale_fill_manual(values=pal5[c(1,3,5)], name="PGS tertile") +
  scale_x_continuous(breaks=seq(23, 43, by=2), expand=c(0,0)) +
  theme_bw() + xlab("gestational age, weeks") + ylab("log hazard ratio") +
  theme(legend.position=c(0.82, 0.83), legend.box.background= element_rect(colour="black"),
        legend.title = element_text(size=10))
p2

plot_grid(p2, p1+ylab(NULL), labels = "AUTO")

ggsave(snakemake@output$suppplotpgs, width=8, height=5, units="in")


# ----------------
# Exploratory part.
# To show:
# 1. rare alleles are not enriched in group 1
# 2. recessive genotypes are not enriched in group 1
# 3. PC3 is correlated w/ vPTD
# 4. PC3 contains lots of immunity

# 1. rare alleles are not enriched in group 1
rare = read.table(rarecountfile, h=T, comment.char = "")
rare$CNT = rare$NAMED_ALLELE_DOSAGE_SUM

rare = inner_join(merged, rare[,c("IID", "CNT")], by=c("SENTRIX_ID"="IID"))

# it is linearly related to PGS
# it is not helping to predict vPTD
summary(lm(SVLEN_DG ~ PGS + CNT + BATCH, data=rare[rare$hadevent,]))
summary(glm(SVLEN_DG<259 ~ PGS + CNT + BATCH,
            data=rare[rare$hadevent,], family="binomial"))
vptd_pval = summary(glm(SVLEN_DG<224 ~ PGS + CNT + BATCH,
                        data=rare[rare$hadevent,], family="binomial")) # note: 0/1 probs w/ many covars so simplified models here
vptd_pval
vptd_pval = vptd_pval$coefficients["CNT",4]
pgs_r = cor(rare$PGS, rare$CNT)

p_rare = mutate(rare, PGScat_rel=factor(PGScat_rel, labels=c("1st", "2nd", "3rd", "4th", "5th"))) %>%
  ggplot(aes(x=PGScat_rel, y=CNT)) + geom_boxplot() +
  ylab("minor allele count") + xlab("GA PGS quintile") +
  annotate("text", x="2nd", y=700, hjust=-0.0, size=3,
           label=sprintf("r=%.2g, vPTD p=%.2g", pgs_r, vptd_pval)) +
  theme_bw() + theme(panel.grid.major.x=element_blank(), 
                     axis.title = element_text(size=9))
p_rare


# 2. recessive genotypes are not enriched in group 1
rare2 = read.table(raregenofile, h=T, comment.char = "")
rare2$CNT_GT = rare2$NAMED_ALLELE_DOSAGE_SUM

rare = inner_join(rare, rare2[,c("IID", "CNT_GT")], by=c("SENTRIX_ID"="IID"))

# it is linearly related to PGS
# it is not helping to predict vPTD
summary(lm(SVLEN_DG ~ PGS + CNT_GT + BATCH, data=rare[rare$hadevent,]))
summary(glm(SVLEN_DG<259 ~ PGS + CNT_GT + BATCH,
            data=rare[rare$hadevent,], family="binomial"))
summary(glm(SVLEN_DG<224 ~ PGS + CNT_GT + BATCH,
            data=rare[rare$hadevent,], family="binomial"))

vptd_pval_gt = summary(glm(SVLEN_DG<224 ~ PGS + CNT_GT + BATCH,
                        data=rare[rare$hadevent,], family="binomial")) # note: 0/1 probs w/ many covars so simplified models here
vptd_pval_gt
vptd_pval_gt = vptd_pval_gt$coefficients["CNT_GT",4]
pgs_r_gt = cor(rare$PGS, rare$CNT_GT)

p_gt = mutate(rare, PGScat_rel=factor(PGScat_rel, labels=c("1st", "2nd", "3rd", "4th", "5th"))) %>%
  ggplot(aes(x=PGScat_rel, y=CNT_GT)) + geom_boxplot() +
  ylab("minor hom. genotype count") + xlab("GA PGS quintile") +
  annotate("text", x="2nd", y=830, hjust=-0.0, size=3,
           label=sprintf("r=%.2g, vPTD p=%.2g", pgs_r_gt, vptd_pval_gt)) +
  theme_bw() + theme(panel.grid.major.x=element_blank(),
                     axis.title = element_text(size=9))
p_gt

# we have also done this w/ different GT counting methods,
# positive alleles, different strength cutoffs, MAF limits.
# All were smoothly decreasing w/ PGS.

# a proxy for possible order 2 GxG interaction count
rare$SQCNT = rare$CNT*(rare$CNT-1)/2
ggplot(rare, aes(x=PGScat_rel, y=SQCNT)) + geom_boxplot()  # still really linear


# 3. PC3 is correlated w/ vPTD
pcs = read.table(pcafile, h=T, comment.char = "")
dim(pcs)
pcs = inner_join(merged, pcs, by=c("SENTRIX_ID"="IID"))

summary(lm(GAc ~ PGS + BATCH + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=pcs[pcs$hadevent,]))
summary(glm(SVLEN_DG<259 ~ PGS + PC1 + PC2 + PC3 + PC4 + PC5 + BATCH,
            data=pcs[pcs$hadevent,], family="binomial"))
summary(glm(SVLEN_DG<224 ~ PGS + PC1 + PC2 + PC3 + PC4 + PC5 + BATCH,
            data=pcs[pcs$hadevent,], family="binomial"))
summary(glm(SVLEN_DG<224 ~ PGS + PC3 + BATCH,
            data=pcs[pcs$hadevent,], family="binomial"))
summary(glm(SVLEN_DG<259 ~ PGS + PC3 + BATCH,
            data=pcs[pcs$hadevent,], family="binomial"))

# PC3 is still sign. even w/ non-linear PGS and clinical covars (except parity to allow fitting)
summary(glm(SVLEN_DG<224 ~ PGS + PGScat + PC3 + BATCH +
              poly(MAGE,2) + FAAR + AA87 + KJONN + MISD,
            data=pcs[pcs$hadevent,], family="binomial"))

ggplot(pcs, aes(x=PGScat_rel, y=PC3)) + geom_boxplot()
cor.test(pcs$PC3, pcs$PGS)  # only weakly related to PGS

ggplot(pcs, aes(x=PC2, y=PC3, col=BATCH)) + geom_point() # not a batch effect

# We also checked that PC3 is not related to reported birth country
# (also none of the birth countries related to vPTD)
# and adding the PC to the PAM doesn't change anything.


vptd_pval = summary(glm(SVLEN_DG<224 ~ PGS + PC3 + BATCH,
                        data=pcs[pcs$hadevent,], family="binomial"))
vptd_pval
vptd_pval = vptd_pval$coefficients["PC3",4]

# mutate(pcs, PC3cat=cut(-PC3, quantile(-pcs$PC3, (0:5)/5), include.lowest=T, labels=1:5)) %>%
#   group_by(PC3cat) %>%
#   summarize(f_vptd=mean(SVLEN_DG<7*32), f_ptd=mean(SVLEN_DG<7*37), n=n()) %>%
#   gather(key="pheno", value="f", f_vptd:f_ptd) %>%
#   mutate(ci_lo = sapply(vecPropTest(f*n, n), function(x) x$conf.int[1]),
#          ci_upp = sapply(vecPropTest(f*n, n), function(x) x$conf.int[2])) %>%
#   mutate(pheno=factor(pheno, labels=c("PTD", "vPTD"))) %>%
#   ggplot(aes(x=PC3cat, y=f)) +
#   geom_col(aes(fill=PC3cat), alpha=0.5) + 
#   geom_linerange(aes(ymin=ci_lo, ymax=ci_upp)) +
#   facet_wrap(~pheno, scales="free_y") +
#   xlab("GA PGS quintile") + ylab("risk") +
#   scale_fill_manual(values=pal5, name="GA PGS quintile", guide="none") +
#   theme_bw() + theme(strip.background = element_rect(fill="#E7E8D3"))

# let's see what goes into PC3
# loadings:
loads = read.table(pcawtfile, h=T, comment.char = "")

# positions:
chr6pos = read.table(chr6posfile)
chr6pos = inner_join(loads, chr6pos, by=c("ID"="V2", "A1"="V5", "X.CHROM"="V1"))

# high loadings in chr 6, in particular in HLA region!
p_allchr = filter(loads, abs(PC3)>1.0) %>%
  ggplot() + geom_boxplot(aes(x=X.CHROM, group=X.CHROM, y=abs(PC3), col=X.CHROM==6)) +
  xlab("chromosome") + ylab("loading magnitude") +
  scale_y_continuous(expand = c(0,0,0,0.7)) + 
  scale_color_manual(values=c("black", "#B54A4A"), guide="none") + 
  theme_bw() + theme(panel.grid.major.x=element_blank(),
                     axis.title=element_text(size=10))

p_chr6zoom = filter(chr6pos, X.CHROM==6) %>%
  mutate(hla=V4>28477797 & V4<33448354) %>%
  ggplot() + geom_point(aes(x=V4/1e6, y=abs(PC3), col=hla), size=0.6) +
  scale_y_continuous(expand = c(0,0,0,0.1)) + scale_x_continuous(expand=c(0,1,0,1)) +
  xlab("position, Mbp") + ylab(NULL) + # ylab("loading magnitude") +
  scale_color_manual(values=c("black", "#B54A4A"), guide="none") + 
  theme_bw() + theme(panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank(),
                     axis.title=element_text(size=9))
# around 6p21.3 marked, positions from https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37

p_pc3 = ggdraw(p_allchr) +
  draw_plot(p_chr6zoom, x=0.43, y=0.57, width=0.55, height=0.45)
  # draw_line(x=c(0.345, 0.475), y=c(0.8, 0.96), col="grey60") +  # maybe connect w/ lines?
  # draw_line(x=c(0.33, 0.36), y=c(0.8, 0.8), col="grey60")

p_rare_gt = plot_grid(p_rare, p_gt)
plot_grid(p_rare_gt, p_pc3, labels="AUTO", nrow=2, rel_heights=c(2.2,2))

ggsave(snakemake@output$mainplotpca, width=4.5, height=5.8, units="in")



#### experimenting
# pcs?...
# 
# # check if range-limited "PC3" can be used to predict vPTD
# # ("scores" based on chr6 alleles w/ weights from pc3)
# tmp1 = read.table("/mnt/HARVEST/plinktests/score-pc3-nochr6.profile", h=T)
# colnames(tmp1)[6] = "SCOREno6"
# pcs2 = inner_join(pcs, tmp1[,c("IID", "SCOREno6")], by=c("SENTRIX_ID"="IID"))
# 
# tmp1 = read.table("/mnt/HARVEST/plinktests/score-pc3-chr6.profile", h=T)
# colnames(tmp1)[6] = "SCOREonly6"
# pcs2 = inner_join(pcs2, tmp1[,c("IID", "SCOREonly6")], by=c("SENTRIX_ID"="IID"))
# 
# tmp1 = read.table("/mnt/HARVEST/plinktests/score-pc3-hla.profile", h=T)
# colnames(tmp1)[6] = "SCOREonlyHLA"
# pcs2 = inner_join(pcs2, tmp1[,c("IID", "SCOREonlyHLA")], by=c("SENTRIX_ID"="IID"))
# 
# head(pcs2)
# 
# cor.test(pcs2$PC3, pcs2$SCOREno6) # although still r=0.99
# cor.test(pcs2$PC3, pcs2$SCOREonly6) # v strong corr
# cor.test(pcs2$PC3, pcs2$SCOREonlyHLA) # v strong corr
# cor.test(pcs2$SCOREno6, pcs2$SCOREonly6)  # hmm
# 
# summary(glm(SVLEN_DG<224 ~ PGS + SCOREno6 + BATCH,
#             data=pcs2[pcs2$hadevent,], family="binomial"))  # not predictive of vPTD though
# 
# 
# cor.test(pgsnohla$SCORESUM, pgsnohla$PC3) # not related to chr1-5,7-22-X PGS
# cor.test(pgsnohla$PGS-pgsnohla$SCORESUM, pgsnohla$PC3) # but related to chr6 PGS!
# 
# 
# 
# # "scores" based on chr6 alleles w/ weights from PGS_GA
# pgsnohla = read.table("/mnt/HARVEST/plinktests/pgs-nochr6.profile", h=T)
# pgsnohla = inner_join(pcs, pgsnohla[,c("IID", "SCORESUM")], by=c("SENTRIX_ID"="IID"))
# head(pgsnohla)
# cor.test(pgsnohla$SCORESUM, pgsnohla$PGS) # PGS mostly preserved
# summary(glm(SVLEN_DG<224 ~ SCORESUM + PC3 + BATCH,  # and still predictive of vPTD
#             data=pgsnohla[pgsnohla$hadevent,], family="binomial"))
# cor.test(pgsnohla$SCORESUM, pgsnohla$PC3) # not related to chr1-5,7-22-X PGS
# 
# pgsonlyhla = read.table("/mnt/HARVEST/plinktests/pgs-onlychr6.profile", h=T)
# pgsonlyhla = inner_join(pcs, pgsonlyhla[,c("IID", "SCORESUM")], by=c("SENTRIX_ID"="IID"))
# cor.test(pgsonlyhla$SCORESUM, pgsonlyhla$PC3) # but related to chr6 PGS!
# 
# 
# # "score" based on chr6 alleles w/ weights from pc3:
# hlascore = read.table("/mnt/HARVEST/plinktests/score-pc3-chr6.profile", h=T)
# colnames(hlascore)[6] = "HLASCORE"
# pgsnohla = inner_join(pgsnohla, hlascore[,c("IID", "HLASCORE")], by=c("SENTRIX_ID"="IID"))
# head(pgsnohla)
# cor.test(pgsnohla$PC3, pgsnohla$HLASCORE) # v strong corr, chr6 explains a lot of PC3
# cor.test(pgsnohla$PGS, pgsnohla$PC3) # related to PGS?... r=0.03 though
# pgsnohla$HLASCORE = pgsnohla$HLASCORE/100
# 
# summary(glm(SVLEN_DG<224 ~ PGS + HLASCORE + BATCH,
#             data=pgsnohla[pgsnohla$hadevent,], family="binomial"))  # not predictive of vPTD though
# 
