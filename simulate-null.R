library(dplyr)
library(ggplot2)
library(broom)
library(mgcv)
library(pammtools)

mfrfile = snakemake@input[[1]] # "/mnt/HARVEST/ga_cleaned.csv"
outfile = snakemake@output[[1]] # "~/Documents/results/tv/null-pvals.RData"

## Run simulations to show that alpha is maintained under the null
## for the main pammtools model that we use
mfr_true = read.table(mfrfile, h=T, sep=";")
nrow(mfr_true)  # 26908
mfr_true$GAc = mfr_true$SVLEN_DG - 169
mfr_true = mfr_true[,c("GAc", "hadevent")]

# simulate genotypes
MAF = 0.3
gt_sim = factor(rbinom(nrow(mfr_true), 2, MAF))

NITER = 400 # 1 h
pvals1 = pvals2 = rep(NA, NITER)

for(i in 1:NITER){
  print(i)
  # bootstrap the outcomes (to get a bit more varied distr than by permuting)
  mfr_boot = sample_frac(mfr_true, 1, replace=T)
  mfr_boot$GTcat = gt_sim
  
  # run the model
  ped = as_ped(mfr_boot, Surv(GAc, hadevent) ~ GTcat, id = "id", cut=c(0,seq(20, 130, by=7)))
  mod = bam(ped_status ~ ti(tend,bs='cr',k=11) + GTcat + ti(tend, by=as.ordered(GTcat),bs='cr'),
            data=ped, offset=offset, family=poisson())
  
  ptable = tidy(mod)
  pvals1[i] = ptable$p.value[2]  # for GTcat=1
  pvals2[i] = ptable$p.value[3]  # for GTcat=2
}

save(pvals1, pvals2, file=outfile)
