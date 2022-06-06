# Main script for generating the publication outputs

library(dplyr)
library(ggplot2)
library(survival)
library(mgcv)
library(pammtools)
library(cowplot)
library(kableExtra)
library(tidyr)
library(broom)

# palette for plotting. colors based on gg_sci::pal_tron
GTpalette = c("#55AEC4", "#F7C530", "#FF410D")


# mfrfile="/mnt/HARVEST/ga_cleaned.csv"
# gtfile="/mnt/HARVEST/top1-moba30k-dosage.csv.gz"
# gtfileX="/mnt/HARVEST/top1x-moba30k-dosage.csv.gz"
# mobaresfile="/mnt/HARVEST/topsnps_moba_summaries.txt"
mfrfile = snakemake@input[[1]]
gtfile = snakemake@input[[2]]
gtfileX = snakemake@input[[3]]
mobaresfile = snakemake@input[[4]]



# -------------------------------------------
# read in filtered phenotypes
mfr_mid = read.table(mfrfile, h=T, sep=";")
nrow(mfr_mid)  # 25648

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
GA_START_TIME = 169
mfr_mid$GAc = mfr_mid$SVLEN_DG-GA_START_TIME

# just some checks
range(mfr_mid$GAc)  # 1 125

# metaanalysis result table
metares = read.table(mobaresfile, h=T)


# -------------------------------------------
# read in genotypes
gt = data.table::fread(gtfile, sep=" ", h=T)
gtinfo = gt[,1:6]
gt = t(gt[,7:ncol(gt)])
gt = data.frame(gt)
gt$SENTRIX_ID = rownames(gt)
rownames(gt) = NULL
gt = filter(gt, !is.na(X1)) # tends to add an empty line in the end
dim(gt)

# merge w/ pheno data
merged = inner_join(gt, mfr_mid, by=c("SENTRIX_ID"))
nrow(merged)  # all 25648 for autosomes

# create ped form for PAMs
ped = as_ped(merged, Surv(GAc, hadevent)~ X1 + X2 + X3 + X4 + X5 + X6 +
               X7 + X8 + X9 + X10 + X11 + X12 +
               X13 + X14 + X15 + X16 + X17 + X18 +         
               X19 + X20 + X21 + X22 + X23 +  
               BATCH + MAGE + FAAR + AA87 + KJONN + MISD + PARITET_5,
             id = "id", cut=c(0,seq(20, 130, by=7)))
nrow(ped)  # 367824


# -------------------------------------------
# Analysis loop

out = data.frame(snpnum=1:23, rsid=NA, ref=NA, eff=NA, locus=NA,
                 lin.beta=NA, lin.se=NA, lin.p=NA)
out_plots = vector("list", length(out$snpnum))

for(snpnum in out$snpnum){
  print(paste("Working on SNP", snpnum))
  # assign the right SNP to GT
  merged$GT = merged[,paste0("X",snpnum)]
  ped$GT = ped[,paste0("X",snpnum)]
  
  # store marker info
  out[snpnum, c("rsid", "chr", "pos")] = gtinfo[snpnum, c("ID", "CHROM", "POS")]
  out[snpnum,"locus"] = metares$gene[metares$POS==gtinfo$POS[snpnum]]
  
  # NOTE: flipping alleles to effect=minor alignment,
  # because the last (factor) analysis is sensitive to that.
  if(mean(merged$GT)>1){
    merged$GT = 2-merged$GT
    ped$GT = 2-ped$GT
    # store alleles
    out[snpnum,c("ref", "eff")] = gtinfo[snpnum, c("ALT", "REF")]
  } else {
    out[snpnum,c("ref", "eff")] = gtinfo[snpnum, c("REF", "ALT")]  
  }
  out[snpnum, "EAF"] = mean(merged$GT)/2
  
  # run the linear model
  print("linear models...")
  mod.lin = lm(SVLEN_DG ~ GT + BATCH + poly(MAGE,2) + FAAR + AA87 + KJONN + MISD + PARITET_5,
             data=merged[merged$hadevent,])
  mod.lin = tidy(mod.lin)
  # report the beta, se and pval for the GT effect
  out[snpnum,c("lin.beta", "lin.se", "lin.p")] = mod.lin[mod.lin$term=="GT", c("estimate", "std.error", "p.value")]
  
  # run the PAM w/ smooth interaction on GT dosage
  print("PAM models...")
  mod.tv.pc.sm = bam(ped_status ~ ti(tend,bs='cr',k=11) + GT + ti(tend, by=GT,bs='cr') +
                       BATCH + poly(MAGE,2) + FAAR + AA87 + KJONN + MISD + PARITET_5,
                     data=ped, offset=offset, family=poisson())
  
  # report the p-value and edf for the smooth
  out.tmp = tidy(mod.tv.pc.sm)
  out[snpnum,c("pc.sm.edf", "pc.sm.p")] = out.tmp[out.tmp$term=="ti(tend):GT",c("edf", "p.value")]
  
  # create & store the predictions and CIs
  print("predicting...")
  plottitle = sprintf("%s / %s, p=%0.2g", out[snpnum,"rsid"], out[snpnum,"locus"], out[snpnum,"pc.sm.p"])
  pred_df = make_newdata(ped, tend=unique(tend), GT=0:2) %>%
    add_term(mod.tv.pc.sm, term="GT", se_mult=1.96) %>%
    mutate(tmid = 0.5*tstart+0.5*tend, locus=plottitle)
  out_plots[[snpnum]] = pred_df[,c("tmid", "locus", "GT", "fit", "ci_lower", "ci_upper")]
  
}

# -------------------------------------------
# read in genotypes for X snps
gt = data.table::fread(gtfileX, sep=" ", h=T)
gtinfo = gt[,1:6]
gt = t(gt[,7:ncol(gt)])
gt = data.frame(gt)
gt$SENTRIX_ID = rownames(gt)
rownames(gt) = NULL
gt = filter(gt, !is.na(X1)) # tends to add an empty line in the end
dim(gt)

# merge w/ pheno data
merged = inner_join(gt, mfr_mid, by=c("SENTRIX_ID"))
nrow(merged)  # all 25648

# create ped form for PAMs
ped = as_ped(merged, Surv(GAc, hadevent)~ X1 + X2 +
               BATCH + MAGE + FAAR + AA87 + KJONN + MISD + PARITET_5,
             id = "id", cut=c(0,seq(20, 130, by=7)))
nrow(ped)  # 367824


# -------------------------------------------
# Analysis for X snps

outX = data.frame(snpnum=1:2, rsid=NA, ref=NA, eff=NA, locus=NA,
                 lin.beta=NA, lin.se=NA, lin.p=NA)
out_plotsX = vector("list", length(outX$snpnum))

for(snpnum in outX$snpnum){
  print(paste("Working on X chr SNP", snpnum))
  # assign the right SNP to GT
  merged$GT = merged[,paste0("X",snpnum)]
  ped$GT = ped[,paste0("X",snpnum)]
  
  # store marker info
  outX[snpnum, c("rsid", "chr", "pos")] = gtinfo[snpnum, c("ID", "CHROM", "POS")]
  outX[snpnum,"locus"] = metares$gene[metares$POS==gtinfo$POS[snpnum]]
  
  # NOTE: flipping alleles to effect=minor alignment,
  # because the last (factor) analysis is sensitive to that.
  if(mean(merged$GT)>1){
    merged$GT = 2-merged$GT
    ped$GT = 2-ped$GT
    # store alleles
    outX[snpnum,c("ref", "eff")] = gtinfo[snpnum, c("ALT", "REF")]
  } else {
    outX[snpnum,c("ref", "eff")] = gtinfo[snpnum, c("REF", "ALT")]  
  }
  outX[snpnum, "EAF"] = mean(merged$GT)/2
  
  # run the linear model
  print("linear models...")
  mod.lin = lm(SVLEN_DG ~ GT + BATCH + poly(MAGE,2) + FAAR + AA87 + KJONN + MISD + PARITET_5,
               data=merged[merged$hadevent,])
  mod.lin = tidy(mod.lin)
  # report the beta, se and pval for the GT effect
  outX[snpnum,c("lin.beta", "lin.se", "lin.p")] = mod.lin[mod.lin$term=="GT", c("estimate", "std.error", "p.value")]
  
  # run the PAM w/ smooth interaction on GT dosage
  print("PAM models...")
  mod.tv.pc.sm = bam(ped_status ~ ti(tend,bs='cr',k=11) + GT + ti(tend, by=GT,bs='cr') +
                       BATCH + poly(MAGE,2) + FAAR + AA87 + KJONN + MISD + PARITET_5,
                     data=ped, offset=offset, family=poisson())
  
  # report the p-value and edf for the smooth
  out.tmp = tidy(mod.tv.pc.sm)
  outX[snpnum,c("pc.sm.edf", "pc.sm.p")] = out.tmp[out.tmp$term=="ti(tend):GT",c("edf", "p.value")]
  
  # create & store the predictions and CIs
  print("predicting...")
  plottitle = sprintf("%s / %s, p=%0.2g", outX[snpnum,"rsid"], outX[snpnum,"locus"], outX[snpnum,"pc.sm.p"])
  pred_df = make_newdata(ped, tend=unique(tend), GT=0:2) %>%
    add_term(mod.tv.pc.sm, term="GT", se_mult=1.96) %>%
    mutate(tmid = 0.5*tstart+0.5*tend, locus=plottitle)
  out_plotsX[[snpnum]] = pred_df[,c("tmid", "locus", "GT", "fit", "ci_lower", "ci_upper")]
}

# attach to autosome results
out = bind_rows(out, outX)
out_plots = c(out_plots, out_plotsX)


# TODO later will probably want to separate them by pvals into two plots.

out
write.table(out, snakemake@output$mainplot, sep="\t")

out_plots = bind_rows(out_plots)
ggplot(out_plots, aes(x=(tmid+GA_START_TIME)/7, y=fit, col=factor(GT), fill=factor(GT))) +
  geom_vline(xintercept = 37, col="grey80") +
  facet_wrap(~locus, scales="free_y") +
  geom_line(lwd=0.6) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),alpha=0.05,lwd=0.2,lty="dashed") +
  scale_color_manual(values=GTpalette, name="allele count") +
  scale_fill_manual(values=GTpalette, name="allele count") +
  scale_x_continuous(breaks=seq(23, 43, by=2)) +
  theme_bw() + xlab("gestational age, weeks") + ylab("hazard ratio") +
  theme(legend.position=c(0.9, 0.2), legend.box.background= element_rect(colour="black"),
        strip.background = element_rect(fill="#E7E8D3"))
ggsave(snakemake@output$maintable, width=7, height=8, units="in")



# later, for geno effects
# merged$GTcat = factor(round(merged$GT))

