
library(data.table)
library(dplyr)

df <- fread(file = snakemake@input[[1]])

df$FID <- df$SENTRIX_ID

df$IID <- df$SENTRIX_ID

df$GA <- df$SVLEN_DG

df <- df %>% select(FID, IID, GA, BATCH)

print("pheno done")

write.table(df, snakemake@output[[1]], 
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)