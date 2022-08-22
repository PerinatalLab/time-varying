
library(data.table)
library(dplyr)

df <- fread(file = snakemake@input[[1]])

df <- df %>% select(FID, SCORESUM) %>% rename (PREG_ID_1724 = FID, PGS = SCORESUM)

write.table(df, snakemake@output[[1]], 
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)