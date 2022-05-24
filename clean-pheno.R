library(dplyr)

# Cleans MFR, creates censoring indicator, merges with sentrix ids
# mfrfile = "PDB1724_MFR_541_v12.csv"
# parentalidfile = "parental_ID_to_PREG_ID.csv"
# sentrixlinkfile = "linkage_Mother_PDB1724.csv"
mfrfile = snakemake@input[[1]] # "PDB1724_MFR_541_v12.csv"
parentalidfile = snakemake@input[[2]] # "parental_ID_to_PREG_ID.csv"
sentrixlinkfile = snakemake@input[[3]] # "linkage_Mother_PDB1724.csv"

# read in phenotypes
mfr = read.table(mfrfile, h=T, sep=";")

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

# ID conversion
link2 = read.table(parentalidfile, h=T, sep=";")
link2 = link2[,c("M_ID_1724", "PREG_ID_1724")]
link2$M_ID_1724 = trimws(link2$M_ID_1724)

# some n lost, probably not genotyped parents:
mfr_mid = inner_join(mfr, link2, by="PREG_ID_1724")
nrow(mfr_mid)  # 99707

# sentrix-pregid converter:
link = read.table(sentrixlinkfile, h=T, sep=";")
link = link[,c("M_ID_1724", "SENTRIX_ID")]
mfr_mid = inner_join(mfr_mid, link, by="M_ID_1724")
nrow(mfr_mid)  # 95388

# genotyping QC flags for the 30k data:
flags = read.table("mobagen-flaglist-n99259.txt", h=T)

# keep only genotyped samples
mfr_mid = inner_join(mfr_mid, flags, by=c("SENTRIX_ID"="IID"))
nrow(mfr_mid)  # 40240

# possibly some sample mixup
mfr_mid = filter(mfr_mid, ROLE=="FOUNDER")
mfr_mid = filter(mfr_mid, genotypesOK, phenoOK)
nrow(mfr_mid)  # 39673

# TODO sort out coreUNRELATED, i.e. relatedness and ancestry problems.
# using that flag is the easiest but costs 3k samples (10%).
# Using PCs instead may save that, plus account for batch effects?..
# seems to be fine when looking at e.g. EBF1 in HARVEST though.
mfr_mid = filter(mfr_mid, coreUNRELATED)
nrow(mfr_mid)  # 35710

# Remove repeated pregnancies & genotyping duplicates:
mfr_mid = mfr_mid[!duplicated(mfr_mid$M_ID_1724),]
nrow(mfr_mid)  # 26908 (or 29763 before core flag)

mfr_mid = mfr_mid[,c("PREG_ID_1724", "SVLEN_DG", "hadevent", "BATCH", "SENTRIX_ID")]

write.csv(mfr_mid, snakemake@output[[1]], quote=F, sep=";", row.names=F)
