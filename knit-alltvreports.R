#!/usr/bin/Rscript

for(snpnum in 1:22){
	rmarkdown::render('run-tvmodels.Rmd', params=list(i=snpnum), output_format="pdf_document", output_file=paste0('/home/julius/Documents/results/tv/report_tvmodels', snpnum))
}
