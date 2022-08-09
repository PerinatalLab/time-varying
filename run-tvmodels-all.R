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

# -------------------------------------------
# SETUP

# mfrfile="/mnt/HARVEST/ga_cleaned.csv"
# mfrfileF="/mnt/HARVEST/ga_cleaned_f.csv"
# gtfile="/mnt/HARVEST/top1-moba30k-dosage.csv.gz"
# gtfileX="/mnt/HARVEST/top1x-moba30k-dosage.csv.gz"
# gtfileF="/mnt/HARVEST/top1f-moba30k-dosage.csv.gz"
# mobaresfile="snplists/topsnps_meta_summaries.txt"
mfrfile = snakemake@input$mfr
mfrfileF = snakemake@input$mfrF
gtfile = snakemake@input$gt
gtfileX = snakemake@input$gtX
gtfileF = snakemake@input$gtF
mobaresfile = snakemake@input$mobares

# some constants:
# palette for plotting. colors based on gg_sci::pal_tron
GTpalette = c("#55AEC4", "#F7C530", "#FF410D")
GA_START_TIME = 169

# metaanalysis result table
metares = read.table(mobaresfile, h=T)


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
# read in filtered phenotypes, FETAL
mfr_fid = read.table(mfrfileF, h=T, sep=";")
nrow(mfr_fid)  # 25515

# Prepare covariates:
mfr_fid$PARITET_5[mfr_fid$PARITET_5==4] = 3  # parity honestly doesn't need 5 levels
mfr_fid$PARITET_5 = factor(mfr_fid$PARITET_5, levels=c(1,0,2,3))
# NOTE: using height will remove extra 1k missing values, so imputing:
mfr_fid$AA87[is.na(mfr_fid$AA87)] = mean(mfr_fid$AA87, na.rm=T)
# relevel to ensure the oldest and one of the largest batches is ref
mfr_fid$BATCH = factor(mfr_fid$BATCH, levels=c("M12A","FEB18","JAN15","JUN15",
                                               "M12B","M24","MAY16",
                                               "ROT1","ROT2"))

mfr_fid$MISD = factor(mfr_fid$MISD) # to prevent gg_slice complaining
mfr_fid$KJONN = factor(mfr_fid$KJONN)
mfr_fid$FAAR = mfr_fid$FAAR-2000  # center, for neater plot scales

# remove irrelevant timeframe for stability
mfr_fid$GAc = mfr_fid$SVLEN_DG-GA_START_TIME

# just some checks
range(mfr_fid$GAc)  # 13 139


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

# read in genotypes for fetal snps
gtF = data.table::fread(gtfileF, sep=" ", h=T)
gtFinfo = gtF[,1:6]
gtF = t(gtF[,7:ncol(gtF)])
gtF = data.frame(gtF)
gtF$SENTRIX_ID = rownames(gtF)
rownames(gtF) = NULL
gtF = filter(gtF, !is.na(X1)) # tends to add an empty line in the end
dim(gtF)
colnames(gtF) = c("X26", "X27", "X28", "X29", "SENTRIX_ID")

gt = full_join(gt, gtX, by=c("SENTRIX_ID"))
gt = full_join(gt, gtF, by=c("SENTRIX_ID"))
gtinfo = rbind(gtinfo, gtXinfo, gtFinfo)


# merge w/ pheno data
merged = inner_join(gt, mfr_mid, by=c("SENTRIX_ID"))
nrow(merged)  # all 26875 for autosomes

# read and attach PGS
pgs = read.table("/mnt/HARVEST/PGS.txt", h=T)
merged = inner_join(merged, pgs, by=c("SENTRIX_ID"))
nrow(merged)  # 26875

# create ped form for PAMs
ped = as_ped(merged, Surv(GAc, hadevent)~ X1 + X2 + X3 + X4 + X5 + X6 +
               X7 + X8 + X9 + X10 + X11 + X12 +
               X13 + X14 + X15 + X16 + X17 + X18 +         
               X19 + X20 + X21 + X22 + X23 +
               X24 + X25 + PGS +
               BATCH + MAGE + FAAR + AA87 + KJONN + MISD + PARITET_5,
             id = "id", cut=c(0,seq(20, 130, by=7)))
nrow(ped)  # 368280


# -------------------------------------------
# Analysis loop MATERNAL

out = data.frame(snpnum=1:29, rsid=NA, chr=NA, pos=NA, ref=NA, eff=NA, EAF=NA, locus=NA,
                 lin.beta=NA, lin.se=NA, lin.p=NA, pc.sm.edf=NA, pc.sm.p=NA)
out_plots = vector("list", length(out$snpnum))
out_pred_all = out_pred_pgs = out_pred_nocov = out_plots

for(snpnum in 1:25){
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
  print("supplemental checks...")
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
  pred_data = make_newdata(ped, tend=unique(tend), GT=c(1)) %>%
    mutate(tmid = 0.5*tstart+0.5*tend, locus=locusname)
  out_pred_all[[snpnum]] = add_term(pred_data, mod.tv.allcov, term="GT",
                         se_mult=1.96)[,c("tmid", "locus", "GT", "fit", "ci_lower", "ci_upper")]
  out_pred_pgs[[snpnum]] = add_term(pred_data, mod.tv.pgscov, term="GT",
                         se_mult=1.96)[,c("tmid", "locus", "GT", "fit", "ci_lower", "ci_upper")]
  out_pred_nocov[[snpnum]] = add_term(pred_data, mod.tv.nocov, term="GT",
                           se_mult=1.96)[,c("tmid", "locus", "GT", "fit", "ci_lower", "ci_upper")]
  
}


# merge w/ pheno data
merged = inner_join(gt, mfr_fid, by=c("SENTRIX_ID"))
nrow(merged)  # all 25515 for autosomes

# create ped form for PAMs
ped = as_ped(merged, Surv(GAc, hadevent)~ X26 + X27 + X28 + X29 +
               BATCH + MAGE + FAAR + AA87 + KJONN + MISD + PARITET_5,
             id = "id", cut=c(0,seq(20, 130, by=7)))
nrow(ped)  # 369140


# -------------------------------------------
# Analysis loop FETAL
for(snpnum in 26:29){
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
  # mod.tv.allcov = bam(ped_status ~ ti(tend,bs='cr',k=11) + GT + ti(tend, by=GT,bs='cr') +
  #                       PGS + BATCH + poly(MAGE,2) + FAAR + AA87 + KJONN + MISD + PARITET_5,
  #                     data=ped, offset=offset, family=poisson())
  # summary(mod.tv.allcov)
  # 
  # mod.tv.pgscov = bam(ped_status ~ ti(tend,bs='cr',k=11) + GT + ti(tend, by=GT,bs='cr') + PGS,
  #                     data=ped, offset=offset, family=poisson())
  # summary(mod.tv.pgscov)
  
  mod.tv.nocov = bam(ped_status ~ ti(tend,bs='cr',k=11) + GT + ti(tend, by=GT,bs='cr'),
                     data=ped, offset=offset, family=poisson())
  summary(mod.tv.nocov)
  
  # only using GT=1 now: baselines don't need plotting,
  # and as GT is continuous here, effect for GT=2 is just 2beta_{GT=1}.
  pred_data = make_newdata(ped, tend=unique(tend), GT=c(1)) %>%
    mutate(tmid = 0.5*tstart+0.5*tend, locus=plottitle)
  # out_pred_all[[snpnum]] = add_term(pred_data, mod.tv.allcov, term="GT",
  #                                   se_mult=1.96)[,c("tmid", "locus", "GT", "fit", "ci_lower", "ci_upper")]
  # out_pred_pgs[[snpnum]] = add_term(pred_data, mod.tv.pgscov, term="GT",
  #                                   se_mult=1.96)[,c("tmid", "locus", "GT", "fit", "ci_lower", "ci_upper")]
  out_pred_nocov[[snpnum]] = add_term(pred_data, mod.tv.nocov, term="GT",
                                      se_mult=1.96)[,c("tmid", "locus", "GT", "fit", "ci_lower", "ci_upper")]
  
}


out
write.table(out, snakemake@output$maintable, sep="\t", quote=F, row.names = F)

out_plots = bind_rows(out_plots)

# SUPP PLOT for insignificant hazards (Fig S2)
out_plots %>%
  mutate(rsid=unlist(lapply(strsplit(locus, " "), "[[", 1))) %>%
  filter(rsid %in% out$rsid[out$lin.p>=0.05]) %>%
  ggplot(aes(x=(tmid+GA_START_TIME)/7, y=fit, col=factor(GT), fill=factor(GT))) +
  geom_vline(xintercept = 37, col="grey80") +
  facet_wrap(~locus, scales="free_y") +
  geom_line(lwd=0.6) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),alpha=0.05,lwd=0.2,lty="dashed") +
  scale_color_manual(values=GTpalette, name="minor allele count") +
  scale_fill_manual(values=GTpalette, name="minor allele count") +
  scale_x_continuous(breaks=seq(23, 43, by=2)) +
  theme_bw() + xlab("gestational age, weeks") + ylab("log hazard ratio") +
  theme(legend.position=c(0.85, 0.1), legend.box.background= element_rect(colour="black"),
        strip.background = element_rect(fill="#E7E8D3"))

ggsave(snakemake@output$supphaz, width=8, height=7, units="in")
# ggsave("~/Documents/results/tv/plot_supphaz.png", width=8, height=6, units="in")


# MAIN PLOT (Fig 2)
facetorder = sort(unique(out_plots$locus))
facetp = as.numeric(unlist(lapply(strsplit(facetorder, "p="), "[[", 2)))
old_plot = out_plots %>%
  mutate(rsid=unlist(lapply(strsplit(locus, " "), "[[", 1))) %>%
  mutate(locus=factor(locus, levels=facetorder[order(facetp)])) %>%
  filter(rsid %in% out$rsid[out$lin.p<0.05]) %>%
  ggplot(aes(x=(tmid+GA_START_TIME)/7, y=fit, col=factor(GT), fill=factor(GT))) +
  geom_vline(xintercept = 37, col="grey80") +
  facet_wrap(~locus, scales="free_y", ncol=4) +
  geom_line(lwd=0.6) +
  # geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),alpha=0.05,lwd=0.2,lty="dashed") +
  scale_color_manual(values=GTpalette, name="minor allele count") +
  scale_fill_manual(values=GTpalette, name="minor allele count") +
  scale_x_continuous(breaks=seq(23, 43, by=2)) +
  theme_bw() + xlab("gestational age, weeks") + ylab("log hazard ratio") +
  scale_y_continuous(expand = expansion(mult=0.3)) +
  theme(legend.position=c(0.78, 0.06), legend.direction="horizontal",
        strip.background = element_rect(fill="#E7E8D3"))

# hax for setting tiny bit tighter y axes
new_plot = old_plot +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),alpha=0.05,lwd=0.2,lty="dashed")
old_plot_data = ggplot_build(old_plot)
new_plot_data = ggplot_build(new_plot)
new_plot_data$layout$panel_params = old_plot_data$layout$panel_params

plot(ggplot_gtable(new_plot_data))
ggsave(snakemake@output$mainplot, plot=ggplot_gtable(new_plot_data), width=9.5, height=9, units="in")


# SUPPLEMENTAL covariate effect plots
out_pred_all = bind_rows(out_pred_all)
out_pred_pgs = bind_rows(out_pred_pgs)
out_pred_nocov = bind_rows(out_pred_nocov)
  
bind_rows("all"=out_pred_all, "clinical"=out_plots,
          "PGS"=out_pred_pgs, "none"=out_pred_nocov, .id="covariates") %>%
  filter(GT==1) %>%
  mutate(covariates=factor(covariates, levels=c("none", "clinical", "PGS", "all")),
         locus=unlist(lapply(strsplit(locus, ","), "[[", 1))) %>%
  ggplot(aes(x=(tmid+GA_START_TIME)/7, y=fit, col=covariates, fill=covariates)) +
  geom_vline(xintercept = 37, col="grey80") +
  geom_hline(yintercept = 0, col="lightblue") +
  facet_wrap(~locus, scales="free_y", ncol=5) +
  geom_line(lwd=0.6) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),alpha=0.02,lwd=0.2,lty="dashed") +
  scale_color_brewer(type="qual", palette="Set2") +
  scale_fill_brewer(type="qual", palette="Set2") +
  scale_x_continuous(breaks=seq(23, 43, by=2)) +
  theme_bw() + xlab("gestational age, weeks") + ylab("log hazard ratio") +
  theme(legend.position=c(0.9, 0.07), legend.box.background= element_rect(colour="black"),
        strip.background = element_rect(fill="#E7E8D3"))

ggsave(snakemake@output$suppplot, width=9, height=9, units="in")


# print some outputs for the text
filter(out, lin.p>=0.05)
filter(out, lin.p<0.05) %>% mutate(pc.sm.p<0.05)


# "zoom in" on IL1A: where is it sign. diff. from 0?
filter(out_plots, grepl("IL1", locus), GT>0, sign(ci_lower)==sign(ci_upper)) %>%
  mutate(tmid=(tmid+GA_START_TIME)/7)


# TODO
# later, for geno effects
# merged$GTcat = factor(round(merged$GT))


