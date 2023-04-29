# install_packages("gprofiler2")
library(gprofiler2)
library(tidyr)
library(dplyr)
loc.all = c("WNT4",
            "HIVEP3",
            "FAF1", 
            "TET3",
            "LSM3",
            "ADCY5",
            "EEFSEC",
            "MRPS22",
            "ZBTB38",
            "KCNAB1",
            "LEKR1",
            "KDR",
            "HAND2",
            "EBF1",
            "EBF1",   # duplicate - gprofiler sorts it out
            "HGNC:4942", # HLA-DQA1 resolved
            "GDAP1",
            "FBXO32",
            "COL27A1",
            "TFAP4",
            "MYOCD",
            "TCEA2",
            "AGTR2",
            "RAP2C",
            "LRP5",
            "LRATD2",
            "IL1A",
            "LPP",
            "SPATA6")
res.all = gost(loc.all, organism="hsapiens")
# check that all were mapped
res.all$meta$genes_metadata

res.all$result

# loci w/ T-V effects in this cohort
loc.varying = c("WNT4",
              "HIVEP3",
              "TET3",
              "ADCY5",
              "EEFSEC",
              "MRPS22",
              "ZBTB38",
              "KCNAB1",
              "LEKR1",
              "KDR",
              "EBF1",
              "GDAP1",
              "FBXO32",
              "COL27A1",
              "TFAP4",
              "MYOCD",
              "TCEA2",
              "AGTR2",
              "RAP2C",
              "LRATD2",
              "IL1A")
res.varying = gost(loc.varying, organism="hsapiens")
res.varying$result

full_join(res.all$result[,c("term_name", "p_value")],
          res.varying$result[,c("term_name", "p_value")],
          by="term_name", suffix=c("all", "tv")) %>%
    View

# removed loci w/ later effect peaks
loc.early = c("WNT4",
                "HIVEP3",
                "ADCY5",
                "EEFSEC",
                "MRPS22",
                "ZBTB38",
                "KCNAB1",
                "KDR",
                "GDAP1",
                "FBXO32",
                "TFAP4",
                "MYOCD",
                "TCEA2",
                "LEKR1",
                "AGTR2",
                "RAP2C",
                "LRATD2")
res.early = gost(loc.early, organism="hsapiens")
res.early$result

res.wide = bind_rows("all"=res.all$result[,c("term_name", "p_value")],
          "tv"=res.varying$result[,c("term_name", "p_value")],
          "early"=res.early$result[,c("term_name", "p_value")],
          .id="gr") %>%
    mutate(p_value = signif(p_value, digits=2)) %>%
    spread(key="gr", value="p_value")

View(res.wide)
arrange(res.wide[,c("term_name", "all", "tv", "early")], early) %>%
    print
