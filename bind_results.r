library(tidyverse)

diseases <- c("BD", "CAD", "CD", "HT", "RA", "T1D", "T2D")

gwaResult.All <- do.call(rbind, lapply(diseases, function(disease) {
    gwaResult <- do.call(rbind, lapply(1:22, function(i) {
      read.delim(paste0(disease,"/gwa", i, ".csv"), sep="\t", header=T)
    }))
    gwaResult$Disease <- disease
    return(gwaResult)
  })
)

write.table(gwaResult.All, file=paste0("gwa_results_all.csv"), sep="\t")