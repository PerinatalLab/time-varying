library(dplyr)
library(ggplot2)

setwd("/mnt/HARVEST")

# --------------------
# read in phenotypes
mfr = read.table("PDB1724_MFR_541_v12.csv", h=T, sep=";")

# clean a bit
mfr = filter(mfr, is.na(ART), is.na(DAAR) | DAAR!=FAAR,
             SVLEN_DG<295, is.na(FLERFODSEL),
             is.na(ZSCORE_BW_GA) | abs(ZSCORE_BW_GA)<10)
# 0 or NA APGARs almost always indicate big problems and aren't genotyped
mfr = filter(mfr, APGAR1>0, APGAR5>0)
# TODO Malformations skipped - might have a better variable later.

# create censoring:
# Might considering using SVLEN_DG>=295 here instead.
mfr = mutate(mfr, hadevent=!is.na(FSTART) & FSTART==1)
mfr = mfr[,c("PREG_ID_1724", "SVLEN_DG", "hadevent")]
nrow(mfr)  # 99765


# -------------------
# ID conversion
link2 = read.table("parental_ID_to_PREG_ID.csv", h=T, sep=";")
link2 = link2[,c("M_ID_1724", "PREG_ID_1724")]
link2$M_ID_1724 = trimws(link2$M_ID_1724)

# some n lost, probably not genotyped parents:
mfr_mid = inner_join(mfr, link2, by="PREG_ID_1724")
nrow(mfr_mid)

# sentrix-pregid converter:
link = read.table("linkage_Mother_PDB1724.csv", h=T, sep=";")
link = link[,c("M_ID_1724", "SENTRIX_ID")]
mfr_mid = inner_join(mfr_mid, link, by="M_ID_1724")

# Remove repeated pregnancies & genotyping duplicates:
mfr_mid = mfr_mid[!duplicated(mfr_mid$M_ID_1724),]
nrow(mfr_mid)  # 79194


# --------------------
# read in genotypes
gt = read.table("ebf1-moms-dosage.csv.gz", h=T, sep=" ")
gtinfo = gt[,1:6]
gt = t(gt[,7:ncol(gt)])
gt = data.frame(gt)
gt$SENTRIX_ID = substring(rownames(gt), 2)
rownames(gt) = NULL
gt = filter(gt, !is.na(X1)) # tends to add an empty line in the end


# merge w/ pheno data
merged = inner_join(gt, mfr_mid, by=c("SENTRIX_ID"))
nrow(merged)  # 9500

# just some checks
sum(duplicated(merged$SENTRIX_ID))
sum(duplicated(merged$PREG_ID_1724))

# remove irrelevant timeframe for stability
merged$GAc = merged$SVLEN_DG-170
# TODO might need to tweak this later based on genotyped data


# TODO temp
merged$GT = merged$X1

# NOTE: flipping alleles to effect=minor alignment,
# because the last analysis is sensitive to that.
if(mean(merged$GT)>1){
  merged$GT = 2-merged$GT
  print("NOTE: allele flipped")
  # TODO might want to flip the gtinfo row accordingly
}

merged$GTcat = factor(round(merged$GT))


GTpalette =  c("#A6D96A", "#FD8D3C", "#D7191C")

# make a plot & a table to check the distributions visually
ggplot(merged[merged$hadevent,], aes(x=SVLEN_DG, group=GTcat)) +
  geom_vline(xintercept=259, col="grey60") +
  geom_density(aes(col=GTcat, fill=GTcat), alpha=0.05) +
  scale_fill_manual(values=GTpalette) + scale_color_manual(values=GTpalette) +
  theme_bw()
merged %>%
  group_by(GTcat) %>%
  summarize(n(), sum(SVLEN_DG<259), sum(SVLEN_DG<240), mean(SVLEN_DG<259))

# ----------------
# models

# test to confirm effect size & direction
summary(lm(SVLEN_DG ~ GT, data=merged[merged$hadevent,]))
summary(lm(SVLEN_DG ~ GTcat, data=merged[merged$hadevent,]))
summary(glm(SVLEN_DG<259 ~ GT, family="binomial", data=merged[merged$hadevent,]))
summary(glm(SVLEN_DG<259 ~ GTcat, family="binomial", data=merged[merged$hadevent,]))

library(survival)

# 1. Cox PH
mod.ph = coxph(Surv(SVLEN_DG, hadevent) ~ GT, data=merged)
summary(mod.ph)

km.ph = survfit(coxph(Surv(GAc, hadevent) ~ GTcat, data=merged), newdata=data.frame(GTcat=factor(0:2)))
plot(km.ph, conf.int=F, col=4:2)
legend(5, 0.4, levels(merged$GTcat), col=4:2, lty=1)
plot(km.ph, conf.int=F, col=4:2, cumhaz=T, log=T)
legend(5, 0.4, levels(merged$GTcat), col=4:2, lty=1)

# Default transform supposedly downweighs outliers.
# Could use the identity one, but that is identical to the interaction w/ t below.
# zp.ph = cox.zph(mod.ph, transform='identity')
zp.ph = cox.zph(mod.ph)
zp.ph
plot(zp.ph, resid=F)


# 2. Time-varying

# KM curve check
km.tv = survfit(Surv(GAc, hadevent) ~ GTcat, data=merged)
plot(km.tv, col=4:2)
legend(5, 0.4, levels(merged$GTcat), col=4:2, lty=1)
plot(km.tv, conf.int=F, col=4:2, cumhaz=T, log=T)
legend(5, 0.4, levels(merged$GTcat), col=4:2, lty=1)

# via specified forms in survival
mod.tv.1 = coxph(Surv(GAc, hadevent) ~ GT + tt(GT), data=merged,
                 tt=function(x, t, ...) x*(t/100))
summary(mod.tv.1)


# contrast for PTD
mod.tv.2 = coxph(Surv(GAc, hadevent) ~ GT + tt(GT), data=merged,
                 tt=function(x, t, ...) x*(t<89))
summary(mod.tv.2)


# via pammtools
library(mgcv)
library(pammtools)
# NOTE: using cubic regr splines for consistency betw. ti() and s()
# convert to piecewise format - i.e. 1 row per each period per indiv.:
ped = as_ped(merged, Surv(GAc, hadevent)~GT + GTcat, id = "id", cut=c(0,seq(20, 130, by=7)))
nrow(ped)  # 135052

# fit a piecewise PH model for testing
mod.pc = gam(ped_status ~ s(tend,bs='cr',k=11) + GT, data=ped, offset=offset, family=poisson())
summary(mod.pc)
plot(mod.pc, pages=1, all.terms=T, scale=0)


# time-varying by prespecified shape
mod.tv.pc.i = gam(ped_status ~ s(tend,bs='cr',k=11) + GT + GT:tend,
                  data=ped, offset=offset, family=poisson())
summary(mod.tv.pc.i)
plot(mod.tv.pc.i, pages=1, all.terms=T, scale=0)


# time-varying smoothly (interaction with spline)
# explanation: ti(tend) is the baseline h(x), and ti(tend, by=GT)
# is a h_1(x)*GT term with potentially different shape.
mod.tv.pc.sm = gam(ped_status ~ ti(tend,bs='cr',k=11) + GT + ti(tend, by=GT,bs='cr'),
                     data=ped, offset=offset, family=poisson())
summary(mod.tv.pc.sm)
plot(mod.tv.pc.sm, pages=1, all.terms=T, scale=0)
# TODO: might need to change GAM to BAM for bigger data, but be aware
# that there are some diffs.


# NOTE: this one is sensitive to allele choice, so should try both ways
# non-linear in time AND non-linear in GT
mod.tv.pc.sm2 = gam(ped_status ~ ti(tend,bs='cr',k=11) + GTcat + ti(tend, by=as.ordered(GTcat),bs='cr'),
                   data=ped, offset=offset, family=poisson())
summary(mod.tv.pc.sm2)
plot(mod.tv.pc.sm2, pages=1, all.terms=T, scale=0)

ped$GTcatRev = factor(round(2-ped$GT))
mod.tv.pc.sm2.r = gam(ped_status ~ ti(tend,bs='cr',k=11) + GTcatRev + ti(tend, by=as.ordered(GTcatRev),bs='cr'),
                    data=ped, offset=offset, family=poisson())
summary(mod.tv.pc.sm2.r)
plot(mod.tv.pc.sm2.r, pages=1, all.terms=T, scale=0)

AIC(mod.pc)
AIC(mod.tv.pc.i); AIC(mod.tv.pc.sm)
AIC(mod.tv.pc.sm2); AIC(mod.tv.pc.sm2.r)

# ---------------
# Useful for printing the summaries:
# library(itsadug)
# gamtabs(mod.tv.pc.sm2.r, type="HTML")
# report_stats(mod.tv.pc.sm2.r)

# ---------------
# VIZ
vis.gam(mod.tv.pc.i, theta=50)
vis.gam(mod.tv.pc.sm, theta=50)
vis.gam(mod.tv.pc.sm2, theta=50)
vis.gam(mod.tv.pc.sm2.r, theta=50)

# library(gratia)
# smooths(mod.pc)
# smooths(mod.tv.pc.i)
# smooths(mod.tv.pc.sm)
# smooths(mod.tv.pc.sm2)
# smooths(mod.tv.pc.sm2.r)
# esthaz.ph = smooth_estimates(mod.pc, smooth=c("s(tend)"))
# esthaz.i = smooth_estimates(mod.tv.pc.i, smooth=c("s(tend)"))
# esthaz.sm = smooth_estimates(mod.tv.pc.sm, smooth=c("ti(tend):GT"))
# esthaz.sm2 = smooth_estimates(mod.tv.pc.sm2,
#                               smooth=c("ti(tend):as.ordered(GTcat)1",
#                                        "ti(tend):as.ordered(GTcat)2"))
# esthaz.sm2.r = smooth_estimates(mod.tv.pc.sm2.r,
#                               smooth=c("ti(tend):as.ordered(GTcatRev)1",
#                                        "ti(tend):as.ordered(GTcatRev)2"))

# will plot 2SE intervals by default
gg_slice(ped, mod.pc, "GT", tend=unique(tend), GT=0:2) + theme_bw()
gg_slice(ped, mod.tv.pc.i, "GT", tend=unique(tend), GT=0:2) + theme_bw()
gg_slice(ped, mod.tv.pc.sm, "GT", tend=unique(tend), GT=0:2) + theme_bw()
gg_slice(ped, mod.tv.pc.sm2, "GTcat", tend=unique(tend), GTcat=factor(0:2)) + theme_bw()
gg_slice(ped, mod.tv.pc.sm2, "GTcat", tend=unique(tend), GTcat=factor(0:2)) + theme_bw() +
  coord_cartesian(ylim=c(-2,2))
gg_slice(ped, mod.tv.pc.sm2.r, "GTcatRev", tend=unique(tend), GTcatRev=factor(0:2)) + theme_bw()

# works like this:
# retrieves predict.gam(type="terms") values
# ped %>%
#   make_newdata(tend=unique(tend), GTcat=factor(0:2)) %>%
#   add_term(mod.tv.pc.sm2, term = "GTcat")
# tmpdata = expand.grid(tend=unique(ped$tend), GT=0:2)
# cbind(tmpdata, predict(mod.pc, tmpdata, se.fit=T, type="iterms")) %>%
#   ggplot(aes(x=tend,col=factor(GT),group=GT)) + geom_line(aes(y=fit.GT)) +
#   geom_ribbon(aes(ymin=fit.GT-2*se.fit.GT, ymax=fit.GT+2*se.fit.GT), alpha=0.1) +
#   theme_bw()


# TODO figure out what of this needs to be stored as the output