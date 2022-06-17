library(dplyr)

# Cleans MFR, creates censoring indicator, merges with sentrix ids
# mfrfile = "/mnt/HARVEST/PDB1724_MFR_541_v12.csv"
# parentalidfile = "/mnt/HARVEST/parental_ID_to_PREG_ID.csv"
# sentrixlinkfile = "/mnt/HARVEST/linkage_Mother_PDB1724.csv"
# sentrixlinkfileF = "/mnt/HARVEST/linkage_Child_PDB1724.csv"
# genoflagfile = "/mnt/HARVEST/mobagen-flaglist-n99259.txt"
# q1file = "/mnt/HARVEST/PDB1724_Q1_v12.csv"
mfrfile = snakemake@input[[1]]
parentalidfile = snakemake@input[[2]]
sentrixlinkfile = snakemake@input[[3]]
sentrixlinkfileF = snakemake@input[[4]]
genoflagfile = snakemake@input[[5]]
q1file = snakemake@input[[6]]

# read in phenotypes
mfr = read.table(mfrfile, h=T, sep=";")

# clean a bit
mfr = filter(mfr, is.na(ART), is.na(DAAR) | DAAR!=FAAR,
             SVLEN_DG<309,
             is.na(FLERFODSEL),
             is.na(DIABETES_MELLITUS),
             is.na(ZSCORE_BW_GA) | abs(ZSCORE_BW_GA)<10)
# 0 or NA APGARs almost always indicate big problems and aren't genotyped
mfr = filter(mfr, APGAR1>0, APGAR5>0)
# For malformations, most are removed via APGAR,
# the remaining ones are mostly irrelevant or unspecified.

# create censoring:
# Might considering using SVLEN_DG>=295 here instead.
mfr = mutate(mfr, hadevent=!is.na(FSTART) & FSTART==1)
nrow(mfr)  # 104454

# Attach maternal height
q1 = read.table(q1file, h=T, sep=";")
q1 = filter(q1, AA87>140)
mfr = left_join(mfr, q1[,c("PREG_ID_1724", "AA87")], by="PREG_ID_1724")

# prepare covariates for modelling
mfr$MISD = !is.na(mfr$MISD)
mfr$MAGE = mfr$FAAR - mfr$MOR_FAAR


# genotyping QC flags for the 30k data:
flags = read.table(genoflagfile, h=T)
nrow(flags)  # 99259
table(flags$ROLE)

# bad/mixed up samples
flags = filter(flags, genotypesOK, phenoOK)
nrow(flags)  # 97679

# coreUNRELATED filters for relatedness and ancestry problems.
# using that flag is the easiest but costs 3k samples in the end (10%).
# seems to be fine when looking at e.g. EBF1 in HARVEST though.
flags = filter(flags, coreUNRELATED)
nrow(flags)  # 88658

# remove the case-control batch
flags = filter(flags, BATCH!="TED")
nrow(flags)  # 83886


# extra ID conversion for mothers
link2 = read.table(parentalidfile, h=T, sep=";")
link2 = link2[,c("M_ID_1724", "PREG_ID_1724")]
link2$M_ID_1724 = trimws(link2$M_ID_1724)

# some n lost, probably not genotyped parents:
mfr_mid = inner_join(mfr, link2, by="PREG_ID_1724")
nrow(mfr_mid)  # 104393

# sentrix-pregid converter, MATERNAL:
linkM = read.table(sentrixlinkfile, h=T, sep=";")
linkM = linkM[,c("M_ID_1724", "SENTRIX_ID")]
mfr_mid = inner_join(mfr_mid, linkM, by="M_ID_1724")
nrow(mfr_mid)  # 99913


# sentrix-pregid converter, FETAL:
linkF = read.table(sentrixlinkfileF, h=T, sep=";")
linkF = linkF[,c("PREG_ID_1724", "SENTRIX_ID")]
mfr_fid = inner_join(mfr, linkF, by="PREG_ID_1724")
nrow(mfr_fid)  # 79113


# keep only genotyped samples:
mfr_mid = inner_join(mfr_mid, flags, by=c("SENTRIX_ID"="IID"))
nrow(mfr_mid)  # 35785

mfr_fid = inner_join(mfr_fid, flags, by=c("SENTRIX_ID"="IID"))
nrow(mfr_fid)  # 25515


# Remove repeated pregnancies & genotyping duplicates:
# (currently sorted by the meaningless pregid so removing a random one)
mfr_mid = mfr_mid[!duplicated(mfr_mid$M_ID_1724),]
nrow(mfr_mid)  # 26875

# none of the kids appear duplicated
mfr_fid = mfr_fid[!duplicated(mfr_fid$PREG_ID_1724),]
nrow(mfr_fid)  # 25515


mfr_mid = mfr_mid[,c("PREG_ID_1724", "SVLEN_DG", "hadevent", "BATCH", "MAGE", "FAAR",
                     "AA87", "KJONN", "MISD", "PARITET_5", "SENTRIX_ID")]
mfr_fid = mfr_fid[,c("PREG_ID_1724", "SVLEN_DG", "hadevent", "BATCH", "MAGE", "FAAR",
                     "AA87", "KJONN", "MISD", "PARITET_5", "SENTRIX_ID")]

write.table(mfr_mid, snakemake@output[[1]], quote=F, sep=";", row.names=F)
write.table(mfr_fid, snakemake@output[[2]], quote=F, sep=";", row.names=F)
