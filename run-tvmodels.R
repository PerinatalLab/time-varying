library(dplyr)
library(ggplot2)

setwd("/mnt/HARVEST")

# --------------------
# read in phenotypes
mfr = read.table("harvest_mfr.csv", h=T, sep=";")

# clean a bit
mfr = filter(mfr, is.na(IVF), is.na(DAAR) | DAAR!=FAAR, SVLEN_DG<295, FLERFODSEL==0, C00_MALF_ALL==0)
# Malformations - might have a better variable later.

# create censoring:
# Might considering using SVLEN_DG>=295 here instead.
mfr = mutate(mfr, hadevent=!is.na(FSTART) & FSTART==1)
mfr = mfr[,c("PREG_ID_1724", "SVLEN_DG", "hadevent")]

# sentrix-pregid converter:
link = read.table("harvest_linkage.csv", h=T, sep=";")
link = link[,c("PREG_ID_1724", "SentrixID_1")]



# --------------------
# read in genotypes
gt = read.table("ebf1-moms-dosage.csv.gz", h=T, sep=" ")
gtinfo = gt[,1:6]
gt = t(gt[,7:ncol(gt)])
gt = data.frame(PID=substring(rownames(gt), 2), GT=as.numeric(gt))

# merge w/ pheno data
gt = inner_join(link, gt, by=c("SentrixID_1"="PID"))

# NOTE: some genotypes will be duplicated (have multiple pregs).
# currently deleting the duplicates
gt = filter(gt, !duplicated(gt$SentrixID_1))
nrow(gt)

merged = inner_join(gt, mfr, by=c("PREG_ID_1724"))
nrow(merged)


# TODO look here
merged$GAc = merged$SVLEN_DG-170
merged$GTeff = 2-merged$GT  # recode to have 0 as the reference allele
merged$GTcat = factor(round(merged$GTeff))

# make a plot & a table to check the distributions visually
ggplot(merged[merged$hadevent,], aes(x=SVLEN_DG, group=GTcat)) +
  geom_density(aes(col=GTcat, fill=GTcat), alpha=0.1)
merged %>%
  group_by(GTcat) %>%
  summarize(n(), sum(SVLEN_DG<255), sum(SVLEN_DG<240), mean(SVLEN_DG<259))

# ----------------
# models

# test to confirm effect size & direction
summary(lm(SVLEN_DG ~ GT, data=merged[merged$hadevent,]))
summary(lm(SVLEN_DG ~ GTcat, data=merged[merged$hadevent,]))
summary(glm(SVLEN_DG<259 ~ GT, family="binomial", data=merged[merged$hadevent,]))

library(survival)

# 1. Cox PH
mod.ph = coxph(Surv(SVLEN_DG, hadevent) ~ GT, data=merged)
summary(mod.ph)

km.ph = survfit(coxph(Surv(GAc, hadevent) ~ GTcat, data=merged), newdata=data.frame(GTcat=factor(0:2)))
plot(km.ph, conf.int=F, col=2:4)
legend(5, 0.4, levels(merged$GTcat), col=2:4, lty=1)
plot(km.ph, conf.int=F, col=2:4, cumhaz=T, log=T)
legend(5, 0.4, levels(merged$GTcat), col=2:4, lty=1)

# using identity b/c the default tr removes outliers (possibly the early births of interest)
zp.ph = cox.zph(mod.ph, transform='identity')
zp.ph
plot(zp.ph, resid=F)


# 2. Time-varying

# KM curve check
km.tv = survfit(Surv(GAc, hadevent) ~ GTcat, data=merged)
plot(km.tv, col=2:4)
legend(5, 0.4, levels(merged$GTcat), col=2:4, lty=1)
plot(km.tv, conf.int=F, col=2:4, cumhaz=T, log=T)
legend(5, 0.4, levels(merged$GTcat), col=2:4, lty=1)

# via specified forms in survival
mod.tv.1 = coxph(Surv(GAc, hadevent) ~ GT + tt(GT), data=merged,
                 tt=function(x, t, ...) x*(t/100))
summary(mod.tv.1)


# hand-tuned form for the suspected location
mod.tv.2 = coxph(Surv(GAc, hadevent) ~ GT + tt(GT), data=merged,
                 tt=function(x, t, ...) x*(t<89))
summary(mod.tv.2)


# via pammtools
library(mgcv)
library(pammtools)
# convert to piecewise format - i.e. 1 row per each period per indiv.:
ped = as_ped(merged, Surv(GAc, hadevent)~GT, id = "id", cut=c(0,seq(20, 130, by=7)))
nrow(ped)
ped$GTcat = factor(round(ped$GT))

# fit a piecewise PH model for testing
mod.pc = gam(ped_status ~ s(tend) + GT, data=ped, offset=offset, family=poisson())
summary(mod.pc)
plot(mod.pc, pages=1, all.terms=T, scale=0)


# time-varying by prespecified shape
mod.tv.pc.i = gam(ped_status ~ s(tend) + GT + GT:tend, data=ped, offset=offset, family=poisson())
summary(mod.tv.pc.i)
plot(mod.tv.pc.i, pages=1, all.terms=T, scale=0)

ped_times = int_info(ped)$tend
ped_effs = coef(mod.tv.pc.i)["GT"] + coef(mod.tv.pc.i)["GT:tend"] * ped_times
plot(x = ped_times, y = ped_effs, type = "s")


# time-varying smoothly (interaction with spline)
# explanation: ti(tend) is the baseline h(x), and ti(tend, by=GT)
# is a h_1(x)*GT term with potentially different shape.
mod.tv.pc.sm = gam(ped_status ~ ti(tend) + GT + ti(tend, by=GT),
                     data=ped, offset=offset, family=poisson())
summary(mod.tv.pc.sm)
plot(mod.tv.pc.sm, pages=1, all.terms=T, scale=0)
# TODO: might need to change GAM to BAM for bigger data, but be aware
# that there are some diffs.


# via rstpm2
library(rstpm2)
mod.ph.stpm = stpm2(Surv(GAc, hadevent)~GT, data=merged, df=21)
summary(mod.ph.stpm)

mod.tv.stpm = stpm2(Surv(GAc, hadevent)~GT, data=merged, tvc=list(GT=11), df=21)
summary(mod.tv.stpm)
plot(mod.tv.stpm, newdata=data.frame(GT=0), line.col=2)
plot(mod.tv.stpm, newdata=data.frame(GT=1), line.col=3, add=T, ci=F)
plot(mod.tv.stpm, newdata=data.frame(GT=2), line.col=4, add=T, ci=F)

# Somewhat questionable way of getting one p-value
anova(mod.ph.stpm, mod.tv.stpm)

# The high dfs currently are selected based on test controls.
# Ideally this should pick them w/ AIC or such.

# TODO check one more SNP
# p=0.14 for non-prop test in Cox, 0.143 for the linear int term,
# 0.297 or 0.493 for timereg non-prop process (empirical) via two tests
# 0.107 for linear int term in GAM, 0.15 & 0.48 for strata in GAM,
# 0.247 for interaction spline non-linearity in GAM
# around 0.5 at best for spline terms in stmp2, 0.108 for ANOVA vs stmp2 PH

