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
if(ncol(gt)!=24) stop("Currently hardcoded for 23 SNPs, make sure you adapt all the code appropriately if this changes!")

# read in genotypes for X snps
gtX = data.table::fread(gtfileX, sep=" ", h=T)
gtXinfo = gtX[,1:6]
gtX = t(gtX[,7:ncol(gtX)])
gtX = data.frame(gtX)
gtX$SENTRIX_ID = rownames(gtX)
rownames(gtX) = NULL
gtX = filter(gtX, !is.na(X1)) # tends to add an empty line in the end
dim(gtX)
colnames(gtX) = c("X24", "X25", "SENTRIX_ID")

gt = inner_join(gt, gtX, by=c("SENTRIX_ID"))
gtinfo = rbind(gtinfo, gtXinfo)

# merge w/ pheno data
merged = inner_join(gt, mfr_mid, by=c("SENTRIX_ID"))
nrow(merged)  # all 25648 for autosomes

# create ped form for PAMs
ped = as_ped(merged, Surv(GAc, hadevent)~ X1 + X2 + X3 + X4 + X5 + X6 +
               X7 + X8 + X9 + X10 + X11 + X12 +
               X13 + X14 + X15 + X16 + X17 + X18 +         
               X19 + X20 + X21 + X22 + X23 +
               X24 + X25 + PGS +
               BATCH + MAGE + FAAR + AA87 + KJONN + MISD + PARITET_5,
             id = "id", cut=c(0,seq(20, 130, by=7)))
nrow(ped)  # 367824

# read and attach PGS
pgs = read.table("/mnt/HARVEST/PGS.txt", h=T)
merged = inner_join(merged, pgs, by=c("SENTRIX_ID"))
nrow(merged)  # 25648


# -------------------------------------------
# Analysis loop

out = data.frame(snpnum=1:25, rsid=NA, chr=NA, pos=NA, ref=NA, eff=NA, EAF=NA, locus=NA,
                 lin.beta=NA, lin.se=NA, lin.p=NA, pc.sm.edf=NA, pc.sm.p=NA)
out_plots = vector("list", length(out$snpnum))

for(snpnum in out$snpnum){
  print(paste("Working on SNP", snpnum))
  # assign the right SNP to GT
  merged$GT = merged[,paste0("X",snpnum)]
  ped$GT = ped[,paste0("X",snpnum)]
  
  # store marker info
  out[snpnum, c("rsid", "pos")] = gtinfo[snpnum, c("ID", "POS")]
  out[snpnum, "chr"] = as.character(gtinfo[snpnum, "CHROM"])
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
  
  # SUPPLEMENTAL: check that adding PGS or other covariates doesn't change the effects
  mod.tv.allcov = bam(ped_status ~ ti(tend,bs='cr',k=11) + GT + ti(tend, by=GT,bs='cr') +
                        PGS + BATCH + poly(MAGE,2) + FAAR + AA87 + KJONN + MISD + PARITET_5,
                      data=ped, offset=offset, family=poisson())
  summary(mod.tv.allcov)
  
  mod.tv.pgscov = bam(ped_status ~ ti(tend,bs='cr',k=11) + GT + ti(tend, by=GT,bs='cr') + PGS,
                      data=ped, offset=offset, family=poisson())
  summary(mod.tv.pgscov)
  
  mod.tv.nocov = bam(ped_status ~ ti(tend,bs='cr',k=11) + GT + ti(tend, by=GT,bs='cr'),
                     data=ped, offset=offset, family=poisson())
  summary(mod.tv.nocov)
  
  locusname = sprintf("%s / %s", out[snpnum,"rsid"], out[snpnum,"locus"])
  # only using GT=1 now: baselines don't need plotting,
  # and as GT is continuous here, effect for GT=2 is just 2beta_{GT=1}.
  pred_data = make_newdata(ped, tend=unique(tend), GT=1) %>%
    mutate(tmid = 0.5*tstart+0.5*tend, locus=locusname)
  out_pred_all[[snpnum]] = add_term(pred_data, mod.tv.allcov, term="GT",
                         se_mult=1.96)[,c("tmid", "locus", "GT", "fit", "ci_lower", "ci_upper")]
  out_pred_pgs[[snpnum]] = add_term(pred_data, mod.tv.pgscov, term="GT",
                         se_mult=1.96)[,c("tmid", "locus", "GT", "fit", "ci_lower", "ci_upper")]
  out_pred_nocov[[snpnum]] = add_term(pred_data, mod.tv.nocov, term="GT",
                           se_mult=1.96)[,c("tmid", "locus", "GT", "fit", "ci_lower", "ci_upper")]
  
}


# TODO later will probably want to separate them by pvals into two plots.

out
write.table(out, snakemake@output$maintable, sep="\t", quote=F, row.names = F)

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
  theme(legend.position=c(0.9, 0.1), legend.box.background= element_rect(colour="black"),
        strip.background = element_rect(fill="#E7E8D3"))
ggsave(snakemake@output$mainplot, width=8, height=10, units="in")


# later, for geno effects
# merged$GTcat = factor(round(merged$GT))



# SUPPLEMENTAL covariate effect plots
out_pred_all = bind_rows(out_pred_all)
out_pred_pgs = bind_rows(out_pred_pgs)
out_pred_nocov = bind_rows(out_pred_nocov)
  
bind_rows("all"=out_pred_all, "clinical"=out_plots,
          "PGS"=out_pred_pgs, "none"=out_pred_nocov, .id="covariates") %>%
  mutate(covariates=factor(covariates, levels=c("none", "clinical", "PGS", "all"))) %>%
  ggplot(aes(x=(tmid+GA_START_TIME)/7, y=fit, col=covariates, fill=covariates)) +
  geom_vline(xintercept = 37, col="grey80") +
  facet_wrap(~locus, scales="free_y") +
  geom_line(lwd=0.6) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),alpha=0.02,lwd=0.2,lty="dashed") +
  scale_color_brewer(type="qual", palette="Set2") +
  scale_fill_brewer(type="qual", palette="Set2") +
  scale_x_continuous(breaks=seq(23, 43, by=2)) +
  theme_bw() + xlab("gestational age, weeks") + ylab("hazard ratio") +
  theme(legend.position=c(0.9, 0.1), legend.box.background= element_rect(colour="black"),
        strip.background = element_rect(fill="#E7E8D3"))





# -------------------------------------------
# PGS

# create ped form for PAMs
ped = as_ped(merged, Surv(GAc, hadevent)~ X1 + X2 + X3 + X4 + X5 + X6 +
               X7 + X8 + X9 + X10 + X11 + X12 +
               X13 + X14 + X15 + X16 + X17 + X18 +         
               X19 + X20 + X21 + X22 + X23 + PGS +  
               BATCH + MAGE + FAAR + AA87 + KJONN + MISD + PARITET_5,
             id = "id", cut=c(0,seq(20, 130, by=7)))
nrow(ped)  # 367824

# PGS-only analyses
summary(lm(GAc ~ PGS, data=merged[merged$hadevent,]))
# r squared is 0.017. Checked w/ regressing out covariates and it doesn't change much.
summary(lm(GAc ~ PGS + BATCH + poly(MAGE,2) + FAAR + AA87 + KJONN + MISD + PARITET_5,
           data=merged[merged$hadevent,]))
summary(glm(SVLEN_DG<259 ~ PGS + BATCH + poly(MAGE,2) + FAAR + AA87 + KJONN + MISD + PARITET_5,
            data=merged[merged$hadevent,], family="binomial"))

mod.tv.pgs = bam(ped_status ~ ti(tend,bs='cr',k=11) + PGS + ti(tend, by=PGS,bs='cr') +
                   BATCH + poly(MAGE,2) + FAAR + AA87 + KJONN + MISD + PARITET_5,
                 data=ped, offset=offset, family=poisson())
summary(mod.tv.pgs)
gg_slice(ped, mod.tv.pgs, "PGS", tend=unique(tend), PGS=c(-1,0,1)) + theme_bw()

# can verify that categorizing PGS also gives the same patterns
ped$PGScat = ifelse(ped$PGS< -1, "lo", ifelse(ped$PGS>1, "hi", "mid"))
ped$PGScat = factor(ped$PGScat, levels=c("mid", "hi", "lo"))
mod.tv.pgs = bam(ped_status ~ ti(tend,bs='cr',k=11) + PGScat + ti(tend, by=as.ordered(PGScat),bs='cr') +
                   BATCH + poly(MAGE,2) + FAAR + AA87 + KJONN + MISD + PARITET_5,
                 data=ped, offset=offset, family=poisson())
summary(mod.tv.pgs)
gg_slice(ped, mod.tv.pgs, "PGScat", tend=unique(tend),
         PGScat=factor(c("hi", "mid", "lo"), levels=c("hi","mid","lo"))) +
  scale_color_manual(values=GTpalette) +
  scale_fill_manual(values=GTpalette) +
  theme_bw()

