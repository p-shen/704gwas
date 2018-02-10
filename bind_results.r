library(tidyverse)
gwaResult <- data.frame(`rsid`=character(), `MinorAllele`=character(),  `MajorAllele`=character(), `DiseaseFrequency`=numeric(), `ControlFrequency`=numeric(), `OR`=numeric(), `pvalue`=numeric(), `hwpvalue`=numeric(), stringsAsFactors=F)
for(i in 1:22) {
  gwaResult <- bind_rows(gwaResult, read.delim(paste0("./results/gwa", i, ".csv"), sep="\t", header=T))
}

write.table(gwaResult, file=paste0("gwa_results.csv"), sep="\t")