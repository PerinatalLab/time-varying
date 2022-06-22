library(dplyr)
library(ggplot2)

library(survival)
library(flexsurv)

library(mgcv)
library(pammtools)

mfr_mid = read.table("/mnt/HARVEST/ga_cleaned.csv", h=T, sep=";")
nrow(mfr_mid)  # 26875

# fit a distr
mod = flexsurvreg(Surv(GAc, hadevent) ~ KJONN, data=mfr_mid, dist="gompertz")
summary(mod)
sh = coef(mod)[1]
rt = exp(coef(mod)[2])
sims = rgompertz(nrow(mfr_mid), shape=sh, rate=rt)
qplot(sims)
qplot(mfr_mid$GAc)

mean(mfr_mid$GAc); mean(sims)
var(mfr_mid$GAc); var(sims)
summary(mfr_mid$GAc); summary(sims)

# -----------------
# simulate PH effects
pgs = rnorm(nrow(mfr_mid))
BETA = -0.2
y.noeff = rgompertz(nrow(mfr_mid), shape=sh, rate=rt)
y.ph = rgompertz(nrow(mfr_mid), shape=sh, rate=rt*exp(BETA*pgs))

df = data.frame(y.noeff, y.ph, pgs, ev=TRUE)

summary(coxph(Surv(y.noeff, ev) ~ pgs, data=df))  # preserves alpha
summary(coxph(Surv(y.ph, ev) ~ pgs, data=df))  # actually good estimates of the beta
summary(lm(y.ph ~ pgs, data=df))  # realistic rsq ~2 %

df$pgscat = cut(df$pgs, breaks=quantile(df$pgs, c(0, 0.333, 0.666, 1)), include.lowest=T,
                labels=c("1st (shortest)", "2nd", "3rd (longest)"))
ggplot(df, aes(x=y.ph, col=pgscat, fill=pgscat)) + geom_density(alpha=0.2) + 
  theme_bw()

# try fitting a PAM
ped = as_ped(df, Surv(y.ph, ev) ~ pgscat,
             id = "id", cut=c(0,seq(20, 130, by=7)))
ped$pgscat = factor(ped$pgscat, levels=c("2nd", "1st (shortest)", "3rd (longest)"))

mod.sim = bam(ped_status ~ ti(tend,bs='cr',k=11) + pgscat + ti(tend, by=as.ordered(pgscat),bs='cr'),
                 data=ped, offset=offset, family=poisson())
summary(mod.sim)
gg_slice(ped, mod.sim, "pgscat", tend=unique(tend),
         pgscat=factor(c("1st (shortest)", "2nd", "3rd (longest)"),
                       levels=c("1st (shortest)", "2nd", "3rd (longest)"))) +
  scale_color_brewer(type="div", palette="RdBu") +
  scale_fill_brewer(type="div", palette="RdBu") +
  theme_bw()
# can get totally horizontal lines actually!


# -----------------
# simulate FRAILTY + PH
pgs = rnorm(nrow(mfr_mid))
fra = rnorm(nrow(mfr_mid))
BETA = -0.3
BETAfra = 1.0
y.ph = rgompertz(nrow(mfr_mid), shape=sh, rate=rt*exp(BETA*pgs + BETAfra*fra))

df = data.frame(y.ph, pgs, ev=TRUE)

summary(coxph(Surv(y.ph, ev) ~ pgs, data=df))  # estimates of the beta slightly towards 0
summary(lm(y.ph ~ pgs, data=df))  # realistic rsq ~2 %

df$pgscat = cut(df$pgs, breaks=quantile(df$pgs, c(0, 0.333, 0.666, 1)), include.lowest=T,
                labels=c("1st (shortest)", "2nd", "3rd (longest)"))
ggplot(df, aes(x=y.ph, col=pgscat, fill=pgscat)) + geom_density(alpha=0.2) + 
  theme_bw()

# try fitting a PAM
ped = as_ped(df, Surv(y.ph, ev) ~ pgscat,
             id = "id", cut=c(0,seq(20, 130, by=7)))
ped$pgscat = factor(ped$pgscat, levels=c("2nd", "1st (shortest)", "3rd (longest)"))

mod.sim = bam(ped_status ~ ti(tend,bs='cr',k=11) + pgscat + ti(tend, by=as.ordered(pgscat),bs='cr'),
              data=ped, offset=offset, family=poisson())
summary(mod.sim)
gg_slice(ped, mod.sim, "pgscat", tend=unique(tend),
         pgscat=factor(c("1st (shortest)", "2nd", "3rd (longest)"),
                       levels=c("1st (shortest)", "2nd", "3rd (longest)"))) +
  scale_color_brewer(type="div", palette="RdBu") +
  scale_fill_brewer(type="div", palette="RdBu") +
  theme_bw()


# -----------------
# simulate SCM w/ bad PGS and another cause
pgs = rnorm(nrow(mfr_mid))
c2 = rnorm(nrow(mfr_mid))
BETA = -1.0
hascause = pgs< -1.5 #| c2< -1.5
hascause = pgs< -1.0 & c2>0
hascause = c2>1.0

y.ph = rgompertz(nrow(mfr_mid), shape=sh, rate=rt*exp(BETA*hascause*pgs-0.1*pgs))

df = data.frame(y.ph, pgs, ev=TRUE)

summary(coxph(Surv(y.ph, ev) ~ pgs, data=df))  # estimates of the beta slightly towards 0
summary(lm(y.ph ~ pgs, data=df))  # realistic rsq ~2 %

df$pgscat = cut(df$pgs, breaks=quantile(df$pgs, c(0, 0.333, 0.666, 1)), include.lowest=T,
                labels=c("1st (shortest)", "2nd", "3rd (longest)"))

ggplot(df, aes(x=y.ph, col=pgscat, fill=pgscat)) + geom_density(alpha=0.2) + 
  theme_bw()

# try fitting a PAM
ped = as_ped(df, Surv(y.ph, ev) ~ pgscat,
             id = "id", cut=c(0,seq(20, 130, by=7)))
ped$pgscat = factor(ped$pgscat, levels=c("2nd", "1st (shortest)", "3rd (longest)"))

mod.sim = bam(ped_status ~ ti(tend,bs='cr',k=11) + pgscat + ti(tend, by=as.ordered(pgscat),bs='cr'),
              data=ped, offset=offset, family=poisson())
summary(mod.sim)
gg_slice(ped, mod.sim, "pgscat", tend=unique(tend),
         pgscat=factor(c("1st (shortest)", "2nd", "3rd (longest)"),
                       levels=c("1st (shortest)", "2nd", "3rd (longest)"))) +
  scale_color_brewer(type="div", palette="RdBu") +
  scale_fill_brewer(type="div", palette="RdBu") +
  theme_bw()

