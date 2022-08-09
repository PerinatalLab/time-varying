### Re-analyses of public expression data

# read in the RNAseq counts (as directly downloaded from GEO)
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE150830&format=file&file=GSE150830_allCounts.xlsx",
              "~/Documents/results/tv/GSE150830_allCounts.xlsx")
ex = readxl::read_excel("~/Documents/results/tv/GSE150830_allCounts.xlsx")

# need to prepare the ensembl-hugo gene id conversion table:
# BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
geneIDs = ensembldb::select(EnsDb.Hsapiens.v86, keys=ex$EnsemblId, keytype="GENEID", columns=c("SYMBOL","GENEID", "GENEBIOTYPE"))
nrow(geneIDs)  # 52699
dim(ex)  # 57905 probes x 125 samples (+1 column)

# read in the sample data
# BiocManager::install("GEOquery")
library(GEOquery)
gse = getGEO("GSE150830", GSEMatrix = TRUE)
gse = gse[[1]]
phenos = pData(gse)
phenos = phenos[,c("description.1", "gestational age:ch1")]
colnames(phenos)[1] = "ID"
colnames(phenos)[2] = "ga"
phenos$batch = !grepl("PAC", phenos$ID)  # data was created in two batches, mark them
phenos$ga = as.numeric(phenos$ga)

# clean up
unloadNamespace("EnsDb.Hsapiens.v86")
unloadNamespace("ensembldb")
unloadNamespace("GEOquery")

library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)

our_genes = read.table("~/Documents/gitrep/time-varying/snplists/topsnps_meta_summaries.txt", h=T)
our_genes = our_genes$gene
our_genes = our_genes[our_genes!="EBF1.2"]  # second locus for this gene
our_genes[our_genes=="LRATD2"] = "FAM84B"  # alternative name

table(our_genes %in% geneIDs$SYMBOL)  # all are found now

# first, calculate the median counts of each transcript
rowmeds = apply(ex[,2:ncol(ex)], 1, median)
df = data.frame(id=ex$EnsemblId, med=rowmeds)
our_genes_ensembl = geneIDs$GENEID[geneIDs$SYMBOL %in% our_genes]
pseudogenes = geneIDs$GENEID[grep("pseudogene", geneIDs$GENEBIOTYPE)]
proteingenes = geneIDs$GENEID[grep("protein_coding", geneIDs$GENEBIOTYPE)]
df = mutate(df, ours=ifelse(id %in% our_genes_ensembl, "GWAS hits",
                    ifelse(id %in% pseudogenes, "pseudogenes",
                    ifelse(id %in% proteingenes, "protein coding", "other genes"))))
group_by(df, ours) %>%
  summarize(median(med))

# plot, excluding pseudogenes
filter(df, ours!="pseudogenes", ours!="other genes") %>%
  ggplot(aes(x=log10(med+1), group=ours, col=ours, fill=ours)) +
  geom_density(alpha=0.1) + xlab("log10(median count + 1)") +
  scale_color_manual(name=NULL, values=c("#FF410D", "#2C7BB6")) + 
  scale_fill_manual(name=NULL, values=c("#FF410D", "#2C7BB6")) + 
  theme_bw() + theme(legend.position=c(0.82, 0.83))
ggsave("~/Documents/results/tv/plot_suppexpr.png", width=4, height=3, units="in")  

# extracting our genes now
ex2 = filter(ex, EnsemblId %in% our_genes_ensembl)

ex2 = inner_join(ex2, geneIDs, by=c("EnsemblId"="GENEID"))
ex2 = gather(ex2, key="sample", value="cnt", -EnsemblId, -SYMBOL)
ex2$cnt = as.numeric(ex2$cnt)
ex2 = inner_join(ex2, phenos, by=c("sample"="ID"))

nrow(ex2)/28  # 96+29=125 samples x 28 genes. amazing

range(ex2$ga)

# all our transcripts were detected in at least 15 samples:
# or: all were present in >=120 samples, except AGTR2 which had low counts
group_by(ex2, SYMBOL) %>%
  summarize(min(cnt), mean(cnt), max(cnt), sum(cnt>0), sum(cnt>=5)) %>%
  arrange(`mean(cnt)`) %>%
  print.data.frame


#### experimenting

# not sure if these slopes are telling anything?..
mods = group_by(ex2, SYMBOL) %>%
  do(fit = glm(cnt ~ ga + batch, data=., family="poisson")) %>%
  summarize(gene=max(SYMBOL), tidy(fit)) %>%
  filter(term!="(Intercept)")
mods %>% View

# plot a few
group_by(ex2, SYMBOL) %>%
  summarize(min(cnt), mean(cnt), max(cnt)) 
mods = group_by(ex2, SYMBOL) %>%
  do(fit = glm(cnt ~ ga, data=., family="poisson")) %>%
  summarize(gene=max(SYMBOL), tidy(fit)) %>%
  filter(term!="(Intercept)")
print.data.frame(mods)

not_interesting = c("LPP", "SPATA6", "LRP5", "MRPS22", "HAND2", "LSM3")
filter(ex2, SYMBOL %in% not_interesting) %>%
  ggplot(aes(x=ga, y=log(cnt+1))) + geom_point() +
  facet_wrap(~SYMBOL) + geom_smooth(method="lm")
filter(ex2, !SYMBOL %in% not_interesting) %>%
  ggplot(aes(x=ga, y=log(cnt+1))) + geom_point() +
  facet_wrap(~SYMBOL) + geom_smooth(method="lm")

filter(ex2, SYMBOL=="WNT4") %>%
  ggplot(aes(x=ga, y=log(cnt+1), col=batch)) + geom_point() +
  facet_wrap(~SYMBOL) + geom_smooth(method="lm")


# btw, these are well correlated with the Danish/UPenn Science study counts:
download.file("https://github.com/miramou/pregnancy_cfRNA/tree/master/raw_data/rnaseq_counts.csv",
              "~/Documents/results/tv/rnaseq_counts.csv")
rnaseq = read.table("~/Documents/results/tv/rnaseq_counts.csv", h=T,sep=",")
danmeds = apply(rnaseq[,2:ncol(rnaseq)], 1, median)
danmeds = data.frame("ID"=rnaseq$external_gene_name, danmeds)
df2 = inner_join(df, geneIDs, by=c("id"="GENEID"))
df2 = filter(df2, GENEBIOTYPE=="protein_coding")
df2 = inner_join(danmeds, df2, by=c("ID"="SYMBOL"))
cor.test(log(df2$danmeds+1), log(df2$med+1))  # r=0.65