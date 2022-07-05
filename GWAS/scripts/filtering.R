library(data.table)
library(dplyr)

dat = fread(snakemake@input[[1]])

dat = dat %>% filter(A1_FREQ !=0) %>% mutate(MAF = ifelse(A1_FREQ < 0.5 ,A1_FREQ,1-A1_FREQ)) %>% filter(MAF>0.01)
dat = dat %>% rename(a1_plink = A1) %>% mutate(A1 = ifelse(A1_FREQ <0.5,a1_plink,REF)) %>% mutate(A2 = ifelse(A1_FREQ >0.5,a1_plink,REF)) #A1 (effect allele) A2 (non-effect allele)

dat = dat %>% rename (RSID = ID,N = OBS_CT)

# change A1 and A2 to lower case
dat$A1 <- tolower(dat$A1)
dat$A2 <- tolower(dat$A2)

fwrite(as.data.frame(dat),snakemake@output[[1]], sep=" ")

