library(dplyr)

# Cleans MFR, creates censoring indicator, merges with sentrix ids
# mfrfile = "PDB1724_MFR_541_v12.csv"
# parentalidfile = "parental_ID_to_PREG_ID.csv"
# sentrixlinkfile = "linkage_Mother_PDB1724.csv"
mfrfile = snakemake@input[[1]] # "PDB1724_MFR_541_v12.csv"
parentalidfile = snakemake@input[[2]] # "parental_ID_to_PREG_ID.csv"
sentrixlinkfile = snakemake@input[[3]] # "linkage_Mother_PDB1724.csv"
genoflagfile = snakemake@input[[4]] # "mobagen-flaglist-n99259.txt"
q1file = snakemake@input[[5]] # "/mnt/HARVEST/PDB1724_Q1_v12.csv"

# read in phenotypes
mfr = read.table(mfrfile, h=T, sep=";")

# clean a bit
mfr = filter(mfr, is.na(ART), is.na(DAAR) | DAAR!=FAAR,
             SVLEN_DG<295, is.na(FLERFODSEL),
             is.na(DIABETES_MELLITUS),
             is.na(ZSCORE_BW_GA) | abs(ZSCORE_BW_GA)<10)
# 0 or NA APGARs almost always indicate big problems and aren't genotyped
mfr = filter(mfr, APGAR1>0, APGAR5>0)
# For malformations, most are removed via APGAR,
# the remaining ones are mostly irrelevant or unspecified.

# create censoring:
# Might considering using SVLEN_DG>=295 here instead.
mfr = mutate(mfr, hadevent=!is.na(FSTART) & FSTART==1)
nrow(mfr)  # 98247

# Attach maternal height
q1 = read.table(q1file, h=T, sep=";")
q1 = filter(q1, AA87>140)
mfr = left_join(mfr, q1[,c("PREG_ID_1724", "AA87")], by="PREG_ID_1724")

# ID conversion
link2 = read.table(parentalidfile, h=T, sep=";")
link2 = link2[,c("M_ID_1724", "PREG_ID_1724")]
link2$M_ID_1724 = trimws(link2$M_ID_1724)

# some n lost, probably not genotyped parents:
mfr_mid = inner_join(mfr, link2, by="PREG_ID_1724")
nrow(mfr_mid)  # 98190

# sentrix-pregid converter:
link = read.table(sentrixlinkfile, h=T, sep=";")
link = link[,c("M_ID_1724", "SENTRIX_ID")]
mfr_mid = inner_join(mfr_mid, link, by="M_ID_1724")
nrow(mfr_mid)  # 93948

# genotyping QC flags for the 30k data:
flags = read.table(genoflagfile, h=T)

# keep only genotyped samples
mfr_mid = inner_join(mfr_mid, flags, by=c("SENTRIX_ID"="IID"))
nrow(mfr_mid)  # 39670

# possibly some sample mixup
mfr_mid = filter(mfr_mid, ROLE=="FOUNDER")
mfr_mid = filter(mfr_mid, genotypesOK, phenoOK)
nrow(mfr_mid)  # 39110

# coreUNRELATED filters for relatedness and ancestry problems.
# using that flag is the easiest but costs 3k samples (10%).
# seems to be fine when looking at e.g. EBF1 in HARVEST though.
mfr_mid = filter(mfr_mid, coreUNRELATED)
nrow(mfr_mid)  # 35213

# remove the case-control batch
mfr_mid = filter(mfr_mid, BATCH!="TED")
nrow(mfr_mid)  # 33648

# prepare covariates for modelling
mfr_mid$MISD = !is.na(mfr_mid$MISD)
mfr_mid$MAGE = mfr_mid$FAAR - mfr_mid$MOR_FAAR

# Remove repeated pregnancies & genotyping duplicates:
# (currently sorted by the meaningless pregid so removing a random one)
mfr_mid = mfr_mid[!duplicated(mfr_mid$M_ID_1724),]
nrow(mfr_mid)  # 25964 (or ~28700 before core flag)

mfr_mid = mfr_mid[,c("PREG_ID_1724", "SVLEN_DG", "hadevent", "BATCH", "MAGE", "FAAR",
                     "AA87", "KJONN", "MISD", "PARITET_5", "SENTRIX_ID")]

write.table(mfr_mid, snakemake@output[[1]], quote=F, sep=";", row.names=F)
