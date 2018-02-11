library(tidyverse)

gwaResult <- do.call(rbind, lapply(1:22, function(i) {
  read.delim(paste0("./results/gwa", i, ".csv"), sep="\t", header=T)
}))

write.table(gwaResult, file=paste0("gwa_results.csv"), sep="\t")