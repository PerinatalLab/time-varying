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
mfrfile = snakemake@input$mfr
pgsfile = "/mnt/HARVEST/PGS.txt"

# some constants:
GA_START_TIME = 169

# make a palette with darker central color
pal5 = RColorBrewer::brewer.pal(5, "RdYlBu")
pal5[1] = "#FF410D"
pal5[2] = "#F7C530"
pal5[3] = "#ABB084"


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
# PGS. 1) basic statistics

# merge maternal PGS w/ pheno data
pgs = read.table(pgsfile, h=T)
merged = inner_join(mfr_mid, pgs, by=c("SENTRIX_ID"))
nrow(merged)  # 26875

# PGS categorized for visualization
PGSlabels = c("1st (shortest)", "2nd", "3rd", "4th", "5th (longest)")
merged$PGScat_rel = cut(merged$PGS, breaks=quantile(merged$PGS, c(0, 0.2, 0.4, 0.6, 0.8, 1)),
                    labels=PGSlabels, include.lowest=T)
merged$PGScat = factor(merged$PGScat_rel, levels=PGSlabels[c(3, 1,2,4,5)])  # set middle ref level


# variance explained: (about 0.018)
summary(lm(GAc ~ PGS, data=merged[merged$hadevent,]))

# (slightly smaller r squared if all births included)
summary(lm(GAc ~ PGS, data=merged))

# w/ covariates:
# (also checked w/ regressing out covariates and it doesn't change much)
summary(lm(GAc ~ PGS + BATCH + poly(MAGE,2) + FAAR + AA87 + KJONN + MISD + PARITET_5,
           data=merged[merged$hadevent,]))


# 2) show overall time-varying effects of PGS groups
# create ped form for PAMs
ped = as_ped(merged, Surv(GAc, hadevent) ~ PGS + PGScat +
               BATCH + MAGE + FAAR + AA87 + KJONN + MISD + PARITET_5,
             id = "id", cut=c(0,seq(20, 130, by=7)))
nrow(ped)  # 368280

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
  theme_bw() + xlab("gestational age, weeks") + ylab("hazard ratio") +
  theme(legend.position=c(0.85, 0.8), legend.box.background= element_rect(colour="black"),
        panel.grid.minor.x=element_blank(), legend.title = element_text(size=10))
pmain

ggsave(snakemake@output$mainplotpgs, width=5.5, height=4.5, units="in")



#### experimenting

# predicting PTD or vPTD: polynomials only support linear term,
# in categories only the risk ones (1-2 or 1) are significant
summary(lm(GAc ~ poly(PGS, 3) + BATCH, data=merged[merged$hadevent,]))
summary(glm(SVLEN_DG<259 ~ poly(PGS,3) + BATCH,
            data=merged[merged$hadevent,], family="binomial"))
summary(glm(SVLEN_DG<259 ~ PGScat + BATCH,
            data=merged[merged$hadevent,], family="binomial"))

summary(glm(SVLEN_DG<224 ~ poly(PGS,3) + BATCH,
            data=merged[merged$hadevent,], family="binomial"))
summary(glm(SVLEN_DG<224 ~ PGScat + BATCH,
            data=merged[merged$hadevent,], family="binomial"))

vecPropTest = Vectorize(prop.test, SIMPLIFY = F)
group_by(merged, PGScat_rel) %>%
  summarize(f_vptd=mean(SVLEN_DG<7*32), f_ptd=mean(SVLEN_DG<7*37), n=n()) %>%
  gather(key="pheno", value="f", f_vptd:f_ptd) %>%
  mutate(ci_lo = sapply(vecPropTest(f*n, n), function(x) x$conf.int[1]),
         ci_upp = sapply(vecPropTest(f*n, n), function(x) x$conf.int[2])) %>%
  mutate(pheno=factor(pheno, labels=c("PTD", "vPTD"))) %>%
  ggplot(aes(x=PGScat_rel, y=f)) +
  geom_col(aes(fill=PGScat_rel), alpha=0.5) + 
  geom_linerange(aes(ymin=ci_lo, ymax=ci_upp)) +
  facet_wrap(~pheno, scales="free_y") +
  xlab("GA PGS quintile") + ylab("risk") +
  scale_fill_manual(values=pal5, name="GA PGS quintile", guide="none") +
  theme_bw() + theme(strip.background = element_rect(fill="#E7E8D3"))

# hmmm. 5group predictions are wobbly, but 3 group are decently linear


summary(glm(!is.na(ICTERUS) ~ SVLEN_DG.x + PGScat_rel, data=merged2))

rare = read.table("/mnt/HARVEST/plinktests/res_burdend_rare.profile", h=T)
head(rare)
qplot(rare$SCORESUM)
rare = inner_join(rare[,c("IID", "SCORESUM")], merged, by=c("IID"="SENTRIX_ID"))
cor.test(rare$SCORESUM, rare$PGS)  # count not related to pgs, weighted obv related
ggplot(rare, aes(x=PGScat_rel, y=SCORESUM)) + geom_boxplot()  # weighted is smooth w/ PGS groups

summary(lm(SVLEN_DG ~ PGS + SCORESUM + BATCH, data=rare[rare$hadevent,]))
summary(glm(SVLEN_DG<259 ~ PGS + SCORESUM + BATCH,
            data=rare[rare$hadevent,], family="binomial"))
summary(glm(SVLEN_DG<224 ~ PGS + SCORESUM + BATCH,
            data=rare[rare$hadevent,], family="binomial"))
# rare-only scores do not add anything to pgs w/ either weighting
# could be interesting to show the plot w/ deleterious rare burden
# (that it is smoothly decr w/ PGS groups)

# just to see that it is not predicting vPTD:
rare = mutate(rare, BURcat=cut(SCORESUM, breaks=quantile(rare$SCORESUM, (0:5)/5),
                               labels=1:5, include.lowest=T))
group_by(rare, BURcat) %>%
  summarize(f_vptd=mean(SVLEN_DG<7*32), f_ptd=mean(SVLEN_DG<7*37), n=n()) %>%
  gather(key="pheno", value="f", f_vptd:f_ptd) %>%
  mutate(ci_lo = sapply(vecPropTest(f*n, n), function(x) x$conf.int[1]),
         ci_upp = sapply(vecPropTest(f*n, n), function(x) x$conf.int[2])) %>%
  mutate(pheno=factor(pheno, labels=c("PTD", "vPTD"))) %>%
  ggplot(aes(x=BURcat, y=f)) +
  geom_col(aes(fill=BURcat), alpha=0.5) + 
  geom_linerange(aes(ymin=ci_lo, ymax=ci_upp)) +
  facet_wrap(~pheno, scales="free_y") +
  xlab("BURDEN_del quintile") + ylab("risk") +
  scale_fill_manual(values=pal5, guide="none") +
  theme_bw() + theme(strip.background = element_rect(fill="#E7E8D3"))

pgs1 = merged$SENTRIX_ID[merged$PGScat=="1st (shortest)"]
write.table(data.frame(pgs1, pgs1), col.names = F, sep="\t", row.names = F, quote=F, file="/mnt/HARVEST/plinktests/pgs1.fam")
pgs1 = merged$SENTRIX_ID[merged$PGScat=="2nd"]
write.table(data.frame(pgs1, pgs1), col.names = F, sep="\t", row.names = F, quote=F, file="/mnt/HARVEST/plinktests/pgs2.fam")
pgs1 = merged$SENTRIX_ID[merged$PGScat=="3rd"]
write.table(data.frame(pgs1, pgs1), col.names = F, sep="\t", row.names = F, quote=F, file="/mnt/HARVEST/plinktests/pgs3.fam")
pgs1 = merged$SENTRIX_ID[merged$PGScat=="4th"]
write.table(data.frame(pgs1, pgs1), col.names = F, sep="\t", row.names = F, quote=F, file="/mnt/HARVEST/plinktests/pgs4.fam")
pgs1 = merged$SENTRIX_ID[merged$PGScat=="5th (longest)"]
write.table(data.frame(pgs1, pgs1), col.names = F, sep="\t", row.names = F, quote=F, file="/mnt/HARVEST/plinktests/pgs5.fam")

pgs1wt = read.table("/mnt/HARVEST/plinktests/pgs1counts.frqx", h=T, sep="\t")
pgs3wt = read.table("/mnt/HARVEST/plinktests/pgs3counts.frqx", h=T, sep="\t")
pgs5wt = read.table("/mnt/HARVEST/plinktests/pgs5counts.frqx", h=T, sep="\t")
betas = read.table("/mnt/HARVEST/plinktests/betas_all.txt", h=T)
pgs1wt = inner_join(pgs1wt, betas, by=c("SNP"="rsid", "A1"="a1"))
pgs3wt = inner_join(pgs3wt, betas, by=c("SNP"="rsid", "A1"="a1"))
pgs5wt = inner_join(pgs5wt, betas, by=c("SNP"="rsid", "A1"="a1"))
mean(pgs1wt$C.HOM.A1.); mean(pgs3wt$C.HOM.A1.); mean(pgs5wt$C.HOM.A1.)  # minor hom gt counts are basically same
mean(pgs1wt$C.HOM.A2.); mean(pgs3wt$C.HOM.A2.); mean(pgs5wt$C.HOM.A2.)  # major hom gt counts are basically same

gtcnt = bind_rows("1"=pgs1wt, "3"=pgs3wt, "5"=pgs5wt, .id="PGS")
group_by(gtcnt, PGS) %>%
  summarize(mean(C.HOM.A1.), mean(C.HOM.A2.))

group_by(gtcnt, PGS) %>%
  filter(W_Beta<0) %>%
  summarize(mean(C.HOM.A1.), mean(C.HOM.A2.)) # these counts again are smooth w/ PGS

group_by(gtcnt, PGS) %>%
  filter(W_Beta < -0.01) %>%
  summarize(mean(C.HOM.A1.), mean(C.HOM.A2.)) # even more so for stronger SNPs

group_by(gtcnt, PGS) %>%
  filter(W_Beta < -0.0, C.HOM.A1.<500) %>%
  summarize(mean(C.HOM.A1.), mean(C.HOM.A2.), n())  # hmmm

filter(gtcnt, C.HOM.A1.<300, W_Beta < -0.01) %>%
  ggplot(aes(x=PGS, y=C.HOM.A1.)) + geom_boxplot()  # visually - smooth w/ PGS again

# TODO counts of rare positive alleles???



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
  theme_bw() + xlab("gestational age, weeks") + ylab("hazard ratio") +
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
  theme_bw() + xlab("gestational age, weeks") + ylab("hazard ratio") +
  theme(legend.position=c(0.82, 0.83), legend.box.background= element_rect(colour="black"),
        legend.title = element_text(size=10))
p2

plot_grid(p2, p1+ylab(NULL), labels = "AUTO")

ggsave(snakemake@output$suppplotpgs, width=8, height=5, units="in")




# -----------------------------

# read in plink's heterozygosity output
het_all = list()
for(chr in c(8:12,16:22)){
  het1 = read.table(paste0("/mnt/HARVEST/plinktests/het",chr,"-all.het"), h=T)
  het1$CHR = chr
  het1$SNPS = "all"
  het2 = read.table(paste0("/mnt/HARVEST/plinktests/het",chr,"-pgs.het"), h=T)
  het2$CHR = chr
  het2$SNPS = "pgs"
  het_all[[length(het_all)+1]] = bind_rows(het1, het2)
}
het_all = bind_rows(het_all)


# recalculate F statistic over all chromosomes:
# (summing plink's fstats gives basically the same)
het = het_all %>%
  group_by(IID, SNPS) %>%
  summarize(Frecalc=(sum(O.HOM.)-sum(E.HOM.))/(sum(N.NM.)-sum(O.HOM.)),
            Fsum=max(F))

het = inner_join(merged, het, by=c("SENTRIX_ID"="IID"))

head(het)
# PGS and F totally uncorrelated
cor.test(het$PGS[het$SNPS=="all"], het$Fsum[het$SNPS=="all"])
cor.test(het$PGS[het$SNPS=="pgs"], het$Fsum[het$SNPS=="pgs"])

het %>%
  group_by(PGScat_rel, SNPS) %>%
  summarize(Fmean=mean(Fsum)) %>%
  ggplot(aes(x=PGScat_rel, y=Fmean, col=SNPS)) + geom_point()

het = filter(het, SNPS=="pgs")
summary(lm(GAc ~ PGS + Fsum + BATCH, data=het[het$hadevent,]))
summary(glm(SVLEN_DG<259 ~ PGS + Fsum + BATCH,
            data=het[het$hadevent,], family="binomial"))
summary(glm(SVLEN_DG<224 ~ PGS + Fsum + BATCH,
            data=het[het$hadevent,], family="binomial"))

summary(lm(GAc ~ PGS + Fsum, data=het[het$hadevent & het$Fsum<0.07,]))
summary(glm(SVLEN_DG<224 ~ PGS + Fsum + BATCH,
            data=het[het$hadevent & het$Fsum<0.07,], family="binomial"))

# --------- experiments ------------
# 
# miniPGS = as.matrix(gt[,1:22]) %*% metares$BETA[1:22] # this one doesn't show the same effect
# # or much effect overall
# miniPGS = data.frame(miniPGS, gt$SENTRIX_ID)
# merged = inner_join(merged, miniPGS, by=c("SENTRIX_ID"="gt.SENTRIX_ID"))
# 
# 
# merged = inner_join(mfr_mid, gt, by=c("SENTRIX_ID"))
# miniPGS = lm(GAc ~ X1 + X2 + X3 + X4 + X5 + X6 +
#                X7 + X8 + X9 + X10 + X11 + X12 +
#                X13 + X14 + X15 + X16 + X17 + X18 +         
#                X19 + X20 + X21 + X22 + X23 +
#                X24 + X25, data=merged)$fitted.values # for now, including all
# merged$miniPGS = miniPGS
# 
# ped = as_ped(merged, Surv(GAc, hadevent)~ #PGS +
#                miniPGS +
#                BATCH + MAGE + FAAR + AA87 + KJONN + MISD + PARITET_5,
#              id = "id", cut=c(0,seq(20, 130, by=7)))
# # ped$PGScat = cut(ped$PGS, breaks=quantile(merged$PGS, c(0, 0.333, 0.666, 1)),
# #                  labels=c("1st (shortest)", "2nd", "3rd (longest)"))
# # ped$PGScat = factor(ped$PGScat, levels=c("2nd", "1st (shortest)", "3rd (longest)"))
# ped$miniPGScat = cut(ped$miniPGS, breaks=quantile(merged$miniPGS, c(0, 0.333, 0.666, 1)),
#                  labels=c("1st (shortest)", "2nd", "3rd (longest)"))
# ped$miniPGScat = factor(ped$miniPGScat, levels=c("2nd", "1st (shortest)", "3rd (longest)"))
# 
# mod.tv.pgs = bam(ped_status ~ ti(tend,bs='cr',k=11) + miniPGScat + ti(tend, by=as.ordered(miniPGScat),bs='cr') +
#                    BATCH + poly(MAGE,2) + FAAR + AA87 + KJONN + MISD + PARITET_5,
#                  data=ped, offset=offset, family=poisson())
# summary(mod.tv.pgs)
# gg_slice(ped, mod.tv.pgs, "miniPGScat", tend=unique(tend),
#          miniPGScat=factor(c("1st (shortest)", "2nd", "3rd (longest)"),
#                        levels=c("1st (shortest)", "2nd", "3rd (longest)"))) +
#   scale_color_brewer(type="div", palette="RdBu") +
#   scale_fill_brewer(type="div", palette="RdBu") +
#   theme_bw()
# 
# 
# # BW
# mfr = read.table("/mnt/HARVEST/PDB1724_MFR_541_v12.csv", sep=";", h=T)
# merged = inner_join(mfr_mid, mfr, by="PREG_ID_1724")
# merged = inner_join(merged, pgs, by="SENTRIX_ID")
# nrow(merged)
# merged$PGScat = cut(merged$PGS, breaks=quantile(merged$PGS, c(0, 0.333, 0.666, 1)),
#                     labels=c("1st (shortest)", "2nd", "3rd (longest)"), include.lowest = T)
# ggplot(merged, aes(x=ZSCORE_BW_GA)) + geom_density(aes(col=PGScat, fill=PGScat), alpha=0.1) + 
#   theme_bw()
# 
# merged = mutate(merged, GAgr=cut(SVLEN_DG.x, c(169, 230, 260, 280, 292, 310)))
# table(merged$GAgr, merged$PGScat)
# merged %>%
#   ggplot(aes(x=ZSCORE_BW_GA)) + geom_density(aes(col=GAgr, fill=GAgr), alpha=0.1) + 
#   facet_wrap(~PGScat) + 
#   theme_bw()
# merged %>%
#   ggplot(aes(x=ZSCORE_BW_GA)) + geom_density(aes(col=PGScat, fill=PGScat), alpha=0.1) + 
#   facet_wrap(~GAgr) + 
#   theme_bw()
# 
# group_by(merged, GAgr, PGScat) %>%
#   filter(!is.na(ZSCORE_BW_GA)) %>%
#   summarize(mean(ZSCORE_BW_GA), sd(ZSCORE_BW_GA), n())

# burden score doesn't seem to be correlated to anything
# rare = read.table("/mnt/HARVEST/plinktests/res_burden_rare-22.profile", h=T)[,c("FID", "CNT2")]
# rare = inner_join(merged, rare, by=c("SENTRIX_ID"="FID"))


# just to check that PGS actually predicts PTD as well
summary(glm(SVLEN_DG<259 ~ PGS + BATCH + poly(MAGE,2) + FAAR + AA87 + KJONN + MISD + PARITET_5,
            data=merged[merged$hadevent,], family="binomial"))

