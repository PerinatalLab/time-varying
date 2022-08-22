
library(data.table)
library(dplyr)

fam <- fread(file = snakemake@input[[1]])
grid <- fread(file = snakemake@input[[2]])
df <- fread(file = snakemake@input[[3]])


prs <- cbind(fam, grid) %>% select(FID, V1) %>% rename(PRS = V1)

df1 <- full_join(df, prs, by=c('SENTRIX_ID' = 'FID'))


write.table(df, snakemake@output[[1]], 
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)