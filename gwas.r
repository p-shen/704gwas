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

print(opt)

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

GWA <- function(csnp, dsnp, snpMap) {
  
  # WTCCC snp to rsid
  rsid <- snps[which(snps[,4]==csnp[1,2]), 5]
  
  # build contingency table
  snpTable <- table(t(csnp[1,5:ncol(csnp)]), useNA="no") %>% bind_rows(table(t(dsnp[1,5:ncol(dsnp)]), useNA="no"))
  snpTable[is.na(snpTable)] <- 0
  
  if(ncol(snpTable)<3) {snpTable <- bind_cols(snpTable, data.frame(`N_N`=c(0,0)))}
  
  snpTable.colNames <- unlist(strsplit(colnames(snpTable[,2]), " "))
  snpTable.additive <- snpTable %>% transmute(a1=.[[1]]*2+.[[2]], a2=.[[2]]+.[[3]]*2)
  colnames(snpTable.additive) <- snpTable.colNames
  
  # determine major and minor alleles
  majorAllele <- max.col(snpTable.additive[1,])
  minorAllele <- colnames(snpTable.additive[,-majorAllele])
  majorAllele <- colnames(snpTable.additive[,majorAllele])
  
  # minor allele frequency in controls
  minAlleleFreqControl <- snpTable.additive[1,minorAllele]/sum(snpTable.additive[1])
  minAlleleFreqDisease <- snpTable.additive[2,minorAllele]/sum(snpTable.additive[2])
  
  # oddsRatio
  OR <- (snpTable.additive[2,minorAllele]*snpTable.additive[1,majorAllele])/(snpTable.additive[2,majorAllele]*snpTable.additive[1,minorAllele])
  
  # chisq test for the snp
  snpTable.additive
  pvalue.snp <- chisq.test(snpTable.additive, correct = F)$p.value
  
  # HW test
  hw.control <- snpTable[1,]
  total <- sum(hw.control)
  colnames(hw.control) <- c("AA", "AB", "BB")
  hw <- hw.control %>% transmute(p=(2*AA+AB)/(2*AA+2*AB+2*BB)) %>% mutate(q=1-p) %>% transmute(AA=total*p^2,   AB=total*2*p*q, BB=total*q^2)
  hw <- hw.control %>% bind_rows(hw)
  pvalue.hw <- chisq.test(hw, correct = F)$p.value
  
  df <- data.frame(`rsid`=as.character(rsid), `MinorAllele`=as.character(minorAllele),  `MajorAllele`=as.character(majorAllele), `DiseaseFrequency`=as.numeric(minAlleleFreqDisease), `ControlFrequency`=as.numeric(minAlleleFreqControl), `OR`=as.numeric(OR), `pvalue`=as.numeric(pvalue.snp), `hwpvalue`=as.numeric(pvalue.hw), stringsAsFactors=F)
  
  return(df)
}

gwaResult <- data.frame(`rsid`=character(), `MinorAllele`=character(),  `MajorAllele`=character(), `DiseaseFrequency`=numeric(), `ControlFrequency`=numeric(), `OR`=numeric(), `pvalue`=numeric(), `hwpvalue`=numeric(), stringsAsFactors=F)

for (i in min(nrow(control), nrow(disease))) {
  csnp <- control[i,]
  dsnp <- t2d22[i,]
  print(i)
  
  res <- GWA(csnp, dsnp, snpMap = snps)
  gwaResult <- bind_rows(gwaResult, res)
}

gwaResult <- gwaResult %>% mutate(`chrom`=chrom)

write.table(gwaResult, file=paste0("gwa", chrom, ".csv"), sep="\t")
