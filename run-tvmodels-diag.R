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

# some constants:
# palette for plotting. colors based on gg_sci::pal_tron
GTpalette = c("#55AEC4", "#F7C530", "#FF410D")
GA_START_TIME = 169

out_prev = read.table(snakemake@input$maintable, h=T)


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
colnames(gtF) = c("X26", "X27", "X28", "SENTRIX_ID")

gt = full_join(gt, gtX, by=c("SENTRIX_ID"))
gt = full_join(gt, gtF, by=c("SENTRIX_ID"))
gtinfo = rbind(gtinfo, gtXinfo, gtFinfo)


# merge w/ pheno data
merged = inner_join(gt, mfr_mid, by=c("SENTRIX_ID"))
nrow(merged)  # all 26875 for autosomes

# create ped form for PAMs
ped = as_ped(merged, Surv(GAc, hadevent)~ X1 + X2 + X3 + X4 + X5 + X6 +
               X7 + X8 + X9 + X10 + X11 + X12 +
               X13 + X14 + X15 + X16 + X17 + X18 +         
               X19 + X20 + X21 + X22 + X23 +
               X24 + X25 + 
               BATCH + MAGE + FAAR + AA87 + KJONN + MISD + PARITET_5,
             id = "id", cut=c(0,seq(20, 130, by=7)))
nrow(ped)  # 368280


# -------------------------------------------
# Analysis loop MATERNAL

out = data.frame(snpnum=1:28, rsid=NA, ref=NA, eff=NA,
                 pc.int.beta=NA, pc.int.p=NA, cox.sch.p=NA,
                 cox.beta=NA, cox.se=NA, pc.beta=NA, pc.se=NA)

for(snpnum in 1:25){
  print(paste("Working on SNP", snpnum))
  # assign the right SNP to GT
  merged$GT = merged[,paste0("X",snpnum)]
  ped$GT = ped[,paste0("X",snpnum)]
  
  # store marker info
  out[snpnum,"rsid"] = gtinfo[snpnum, "ID"]

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
  
  # run the PAM w/ linear interaction on GT dosage
  print("PAM model w/ linear interaction...")
  mod.tv.pc.int = bam(ped_status ~ ti(tend,bs='cr',k=11) + GT + GT:tend +
                       BATCH + poly(MAGE,2) + FAAR + AA87 + KJONN + MISD + PARITET_5,
                     data=ped, offset=offset, family=poisson())
  
  # report the p-value and edf for the smooth
  out.tmp = tidy(mod.tv.pc.int, parametric=T)
  out[snpnum,c("pc.int.beta", "pc.int.p")] = out.tmp[out.tmp$term=="GT:tend",c("estimate", "p.value")]
  
  
  # run PAM w/o smooth terms for comparison v Cox
  print("no-interaction models...")
  mod.pc = bam(ped_status ~ s(tend,bs='cr',k=11) + GT + BATCH +
                 poly(MAGE,2) + FAAR + AA87 + KJONN + MISD + PARITET_5,
               data=ped, offset=offset, family=poisson())
  mod.ph = coxph(Surv(GAc, hadevent) ~ GT + BATCH + poly(MAGE,2) + FAAR +
                   AA87 + KJONN + MISD + PARITET_5, data=merged)
  
  # confirm that doesn't differ much
  out.tmp = tidy(mod.pc,parametric=T)
  out[snpnum,c("pc.beta", "pc.se")] = out.tmp[out.tmp$term=="GT", c("estimate", "std.error")]
  
  out.tmp = tidy(mod.ph)
  out[snpnum,c("cox.beta", "cox.se")] = out.tmp[out.tmp$term=="GT", c("estimate", "std.error")]
  
  # do the schoenfeld residuals test with KM transform
  out.tmp = cox.zph(mod.ph)$table
  out[snpnum,"cox.sch.p"] = out.tmp[rownames(out.tmp)=="GT", "p"]
}


# merge w/ pheno data
merged = inner_join(gt, mfr_fid, by=c("SENTRIX_ID"))
nrow(merged)  # all 25515 for autosomes

# create ped form for PAMs
ped = as_ped(merged, Surv(GAc, hadevent)~ X26 + X27 + X28 + 
               BATCH + MAGE + FAAR + AA87 + KJONN + MISD + PARITET_5,
             id = "id", cut=c(0,seq(20, 130, by=7)))
nrow(ped)  # 369140


# -------------------------------------------
# Analysis loop FETAL
for(snpnum in 26:28){
  print(paste("Working on SNP", snpnum))
  # assign the right SNP to GT
  merged$GT = merged[,paste0("X",snpnum)]
  ped$GT = ped[,paste0("X",snpnum)]
  
  # store marker info
  out[snpnum,"rsid"] = gtinfo[snpnum, "ID"]

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
  
  # run the PAM w/ linear interaction on GT dosage
  print("PAM model w/ linear interaction...")
  mod.tv.pc.int = bam(ped_status ~ ti(tend,bs='cr',k=11) + GT + GT:tend +
                        BATCH + poly(MAGE,2) + FAAR + AA87 + KJONN + MISD + PARITET_5,
                      data=ped, offset=offset, family=poisson())
  
  # report the p-value and edf for the smooth
  out.tmp = tidy(mod.tv.pc.int, parametric=T)
  out[snpnum,c("pc.int.beta", "pc.int.p")] = out.tmp[out.tmp$term=="GT:tend",c("estimate", "p.value")]
  
  
  # run PAM w/o smooth terms for comparison v Cox
  print("no-interaction models...")
  mod.pc = bam(ped_status ~ s(tend,bs='cr',k=11) + GT + BATCH +
                 poly(MAGE,2) + FAAR + AA87 + KJONN + MISD + PARITET_5,
               data=ped, offset=offset, family=poisson())
  mod.ph = coxph(Surv(GAc, hadevent) ~ GT + BATCH + poly(MAGE,2) + FAAR +
                   AA87 + KJONN + MISD + PARITET_5, data=merged)
  
  # confirm that doesn't differ much
  out.tmp = tidy(mod.pc,parametric=T)
  out[snpnum,c("pc.beta", "pc.se")] = out.tmp[out.tmp$term=="GT", c("estimate", "std.error")]
  
  out.tmp = tidy(mod.ph)
  out[snpnum,c("cox.beta", "cox.se")] = out.tmp[out.tmp$term=="GT", c("estimate", "std.error")]
  
  # do the schoenfeld residuals test with KM transform
  out.tmp = cox.zph(mod.ph)$table
  out[snpnum,"cox.sch.p"] = out.tmp[rownames(out.tmp)=="GT", "p"]
}


# attach locus and pval from the maintext PAM model
out = left_join(out, out_prev[,c("rsid", "locus", "ref", "eff", "pc.sm.p")],
                by=c("rsid", "ref", "eff"))

out

write.table(out, snakemake@output$diagtable, sep="\t", quote=F, row.names = F)

ggplot(out, aes(x=cox.beta, y=pc.beta)) +
  geom_vline(xintercept = 0, col="grey80") + geom_hline(yintercept = 0, col="grey80") +
  geom_abline(slope = 1, col="grey70") +
  geom_errorbar(aes(ymin=pc.beta-pc.se, ymax=pc.beta+pc.se), col="#FF410D", alpha=0.3) +
  geom_errorbarh(aes(xmin=cox.beta-cox.se, xmax=cox.beta+cox.se), col="#FF410D", alpha=0.3) +
  geom_point(col="grey20", pch=18) +
  xlab("log hazard ratio, Cox regression") + ylab("log hazard ratio, PAM") +
  coord_fixed() + theme_bw()

ggsave(snakemake@output$coxdiagplot, width=6, height=5, units="in")
ggsave("~/Documents/results/tv/plot_coxdiag.png", width=4, height=4, units="in")

out[,c("rsid", "locus", "ref", "eff", "pc.sm.p", "pc.int.p", "cox.sch.p")] %>%
  arrange(pc.sm.p) %>%
  mutate_at(c("pc.sm.p", "pc.int.p", "cox.sch.p"), prettyNum, digits=3) %>%
  kable(format="simple")

  
# later, for geno effects
# merged$GTcat = factor(round(merged$GT))

