#!/usr/bin/env Rscript
# load a library for option parsing
library(optparse)

# define the options that the script should support
opt.list <- list(
  make_option(c("-a","--a"),type='character',help="file for the 1st control sample"),
  make_option(c("-b","--b"),type='character',help="file for the 2nd control sample"),
  make_option(c("-d", "--disease"), type='character',help="file for the diseased sample"),
  make_option(c("-n", "--chromosome"), type='integer',help="chromosome number"),
  make_option(c("-s", "--snp"), type='character',help="SNP rsid file")
)

# do the argument parsing
pa <- parse_args(OptionParser(option_list=opt.list),positional_arguments=T)
opt <- pa$options;

control1 <- opt$a
control2 <- opt$b
disease <- opt$disease
chrom <- opt$chromosome
snps <- opt$snp

print(paste(control1, control2, disease, chrom, snps, sep = ", "))

# check to see if the arguments are both provided
if(is.null(control1) || is.null(control2) || is.null(disease) || is.null(chrom) || is.null(snps)) {
  stop("Control files or disease file, chromosome number, or snp file is not provided. Run 'Rscript execute_gwas.r -h' to see the help menu.")
}

# does it exist?
if (!file.exists(control1)) {
  print(paste(control1, " does not exist"))
} else if (!file.exists(control2)) {
  print(paste(control2, " does not exist"))
}else if (!file.exists(disease)) {
  print(paste(disease, " does not exist"))
} else if (!file.exists(snps)) {
  print(paste(snps, " does not exist"))
}

#----------------------------
# Start processing GWAS data
library(tidyverse)

control1 <- read.delim(control1, sep="\t", header=F)
control2 <- read.delim(control2, sep="\t", header=F)
disease <- read.delim(disease, sep="\t", header=F)
snps <- read.delim(snps, header=F)

control <- control1 %>% bind_cols(control2[,5:ncol(control2)])

GWA <- function(i) {
  
  csnp <- control[i,]
  dsnp <- disease[i,]
  
  # build contingency table
  snpTable <- table(t(csnp[1,5:ncol(csnp)]), useNA="no") %>% bind_rows(table(t(dsnp[1,5:ncol(dsnp)]), useNA="no"))
  snpTable[is.na(snpTable)] <- 0
  
  # check for edge cases with alleles
  if (ncol(snpTable)==1){
    # if there is only 1 SNP, there's not going to be any contributional effect from this snp so just skip it
    return(NA)
  } else if(ncol(snpTable)==2) {
    # if there are 2 SNPS found in the samples, then we just add an empty count for the third possible configuration of the SNP
    snpTable <- bind_cols(snpTable, data.frame(`N_N`=c(0,0)))
  }
  
  snpTable.colNames <- unlist(strsplit(colnames(snpTable[,2]), " "))
  snpTable.additive <- snpTable %>% transmute(a1=.[[1]]*2+.[[2]], a2=.[[2]]+.[[3]]*2)
  colnames(snpTable.additive) <- snpTable.colNames
  
  # determine major and minor alleles
  majorAllele <- max.col(snpTable.additive[1,])
  minorAllele <- colnames(snpTable.additive[1,-majorAllele])
  majorAllele <- colnames(snpTable.additive[1,majorAllele])
  
  return(c(`WTCCC`=as.character(csnp[1,2]), 
           `MinorAllele`=minorAllele,  
           `MajorAllele`=majorAllele,
           `ControlMajorAlleleCount`=snpTable.additive[[1,majorAllele]],
           `ControlMinorAlleleCount`=snpTable.additive[[1,minorAllele]],
           `DiseaseMajorAlleleCount`=snpTable.additive[[2,majorAllele]],
           `DiseaseMinorAlleleCount`=snpTable.additive[[2,minorAllele]],
           `AACount`=snpTable[[1,1]],
           `ABCount`=snpTable[[1,2]],
           `BBCount`=snpTable[[1,3]]))
  
}

gwaResult <- data.frame(t(sapply(seq_len(nrow(disease)), GWA)), stringsAsFactors = F)

# Class conversions
gwaResult <- gwaResult %>% mutate(ControlMajorAlleleCount=as.numeric(ControlMajorAlleleCount),
                                  ControlMinorAlleleCount=as.numeric(ControlMinorAlleleCount),
                                  DiseaseMajorAlleleCount=as.numeric(DiseaseMajorAlleleCount),
                                  DiseaseMinorAlleleCount=as.numeric(DiseaseMinorAlleleCount),
                                  AACount=as.numeric(AACount),
                                  ABCount=as.numeric(ABCount),
                                  BBCount=as.numeric(BBCount)) %>%
  mutate(rsid=snps[which(snps[,4] %in% WTCCC), 5]) %>% ## get rsid
  ## frequency calculations
  mutate(ControlMinAlleleFreq=ControlMinorAlleleCount/(ControlMinorAlleleCount+ControlMajorAlleleCount)) %>% 
  mutate(DiseaseMinAlleleFreq=DiseaseMinorAlleleCount/(DiseaseMinorAlleleCount+DiseaseMajorAlleleCount)) %>%
  ## Odds Ratio Calculation
  mutate(OR=(DiseaseMinorAlleleCount*ControlMajorAlleleCount)/(DiseaseMajorAlleleCount*ControlMinorAlleleCount)) %>%
  ## generate HW p and q values
  mutate(p=(2*AACount+ABCount)/(2*AACount+2*ABCount+2*BBCount), total = AACount + ABCount + BBCount) %>% mutate(q=1-p) %>% mutate(HWAA=total*p^2, HWAB=total*2*p*q, HWBB=total*q^2)


# Chi-sq tests for SNP and HW
gwaResult$PValue <- apply(gwaResult, 1, function(x) {
  chisq.test(
    matrix(
      as.numeric(c(x['ControlMajorAlleleCount'], x['ControlMinorAlleleCount'],
                   x['DiseaseMajorAlleleCount'], x['DiseaseMinorAlleleCount'])),
      nrow=2, ncol=2, byrow = T)
    , correct = F)$p.value
})

# HW chi-sq test
gwaResult$HWPValue <- apply(gwaResult, 1, function(x) {
  chisq.test(
    matrix(
      as.numeric(c(x['AACount'], x['ABCount'], x['BBCount'],
                   x['HWAA'], x['HWAB'], x['HWBB'])),
      nrow=2, ncol=3, byrow = T)
    , correct = F)$p.value
})

gwaResult <- gwaResult %>% select(rsid, MinorAllele, MajorAllele, DiseaseMinAlleleFreq, ControlMinAlleleFreq, OR, PValue, HWPValue)
gwaResult$Chrom <- chrom

write.table(gwaResult, file=paste0("./results/gwa", chrom, ".csv"), sep="\t")
