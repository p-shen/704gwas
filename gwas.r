#!/usr/bin/env Rscript
# load a library for option parsing
library(optparse)

# define the options that the script should support
opt.list <- list(
  make_option(c("-a","--a"),type='character',help="file for the 1st control sample"),
  make_option(c("-b","--b"),type='character',help="file for the 2nd control sample"),
  make_option(c("-d", "--disease"), type='character',help="file for the diseased sample"),
  make_option(c("-n", "--chromosome"), type='integer',help="chromosome number"),
  make_option(c("-s", "--snp"), type='character',help="SNP rsid file"),
  make_option(c("-g", "--disease_name"), type='character',help="Disease name (eg. T2D)")
)

# do the argument parsing
pa <- parse_args(OptionParser(option_list=opt.list),positional_arguments=T)
opt <- pa$options;

control1 <- opt$a
control2 <- opt$b
disease <- opt$disease
chrom <- opt$chromosome
snps <- opt$snp
diseaseName <- opt$disease_name

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

control1 <- read.delim(control1, sep="\t", header=F, stringsAsFactors = F)
control2 <- read.delim(control2, sep="\t", header=F, stringsAsFactors = F)
disease <- read.delim(disease, sep="\t", header=F, stringsAsFactors = F)
snps <- read.delim(snps, header=F, stringsAsFactors = F)

control <- control1 %>% bind_cols(control2[,5:ncol(control2)])

GWA <- function(csnp, dsnp) {
  
  # build contingency table
  snpTable <- table(t(csnp[5:ncol(csnp)]), useNA="no") %>% bind_rows(table(t(dsnp[5:ncol(dsnp)]), useNA="no"))
  snpTable[is.na(snpTable)] <- 0
  
  # check for edge cases with alleles
  if (ncol(snpTable)==1){
    # if there is only 1 SNP, there's not going to be any contributional effect from this snp so just skip it
    return(NA)
  } else if(ncol(snpTable)==2) {
    # if there are 2 SNPS found in the samples, then we just add an empty count for the third possible configuration of the SNP
    snpTable <- bind_cols(snpTable, data.frame(`N N`=c(0,0)))
  }
  
  snpTable.colNames <- unlist(strsplit(colnames(snpTable[,2]), " "))
  snpTable.additive <- snpTable %>% transmute(a1=.[[1]]*2+.[[2]], a2=.[[2]]+.[[3]]*2)
  colnames(snpTable.additive) <- snpTable.colNames
  
  # determine major and minor alleles
  majorAllele <- max.col(snpTable.additive[1,])
  minorAllele <- colnames(snpTable.additive[1,-majorAllele])
  majorAllele <- colnames(snpTable.additive[1,majorAllele])
  
  return(c (`Chrom`=as.numeric(csnp[1]),
            `WTCCC`=as.character(csnp[[2]]), 
            `MinorAllele`=minorAllele,  
            `MajorAllele`=majorAllele,
            `ControlMajorAlleleCount`=snpTable.additive[[1,majorAllele]],
            `ControlMinorAlleleCount`=snpTable.additive[[1,minorAllele]],
            `DiseaseMajorAlleleCount`=snpTable.additive[[2,majorAllele]],
            `DiseaseMinorAlleleCount`=snpTable.additive[[2,minorAllele]],
            `AACount`=snpTable[[1,1]],
            `ABCount`=snpTable[[1,2]],
            `BBCount`=snpTable[[1,3]]) )
  
}

# QC for files
print(paste("Performing GWAS for ", diseaseName, " on chromosome ", chrom))
print(paste("Performing QC for chromosome", chrom))
control.misread <- which(control$V1==chrom)
control <- control[control.misread,]
disease <- disease[control.misread,]

disease.misread <- which(disease$V1==chrom)
control <- control[disease.misread,]
disease <- disease[disease.misread,]

print(paste("Performing GWA for chromosome", chrom))
gwaResult <- do.call(rbind, lapply(seq_len(nrow(control)), function(i){GWA(control[i,], disease[i,])}))
gwaResult <- data.frame(gwaResult, stringsAsFactors = F)

print("Finished GWAS, writing out intermediate data")
write.table(gwaResult, file=paste0("./imm/", diseaseName, "/", "gwa", chrom, ".csv"), sep="\t")

# Remove any rows with NA values
print("Remove rows with NA")
gwaResult <- gwaResult[complete.cases(gwaResult),]

# Class conversions
print("Class conversion")
gwaResult <- gwaResult %>% mutate(Chrom=as.numeric(Chrom),
                                  ControlMajorAlleleCount=as.numeric(ControlMajorAlleleCount),
                                  ControlMinorAlleleCount=as.numeric(ControlMinorAlleleCount),
                                  DiseaseMajorAlleleCount=as.numeric(DiseaseMajorAlleleCount),
                                  DiseaseMinorAlleleCount=as.numeric(DiseaseMinorAlleleCount),
                                  AACount=as.numeric(AACount),
                                  ABCount=as.numeric(ABCount),
                                  BBCount=as.numeric(BBCount))

# get RSID
print("Get RSID")
snps <- snps %>% select(V2, V4, V5) %>% transmute(WTCCC=V4, rsid=V5, pos=V2)
gwaResult <- gwaResult %>% left_join(snps, by='WTCCC')

print("Frequency Calculations")
gwaResult <- gwaResult %>% 
  mutate(ControlMinAlleleFreq=ControlMinorAlleleCount/(ControlMinorAlleleCount+ControlMajorAlleleCount)) %>%
  mutate(DiseaseMinAlleleFreq=DiseaseMinorAlleleCount/(DiseaseMinorAlleleCount+DiseaseMajorAlleleCount))

print("Calculate OR Ratios")
gwaResult <- gwaResult %>% mutate(OR=(DiseaseMajorAlleleCount*ControlMinorAlleleCount)/(DiseaseMinorAlleleCount*ControlMajorAlleleCount))

print("generate HW p and q values")
gwaResult <- gwaResult %>% 
  mutate(p=(2*AACount+ABCount)/(2*AACount+2*ABCount+2*BBCount), total = AACount + ABCount + BBCount) %>% mutate(q=1-p) %>% mutate(HWAA=total*p^2, HWAB=total*2*p*q, HWBB=total*q^2)


# Chi-sq tests for SNP and HW
print("Chi-sq Test for SNPs")
gwaResult$PValue <- apply(gwaResult, 1, function(x) {
  chisq.test(matrix(
    as.numeric(c(x['ControlMajorAlleleCount'], x['ControlMinorAlleleCount'],
                 x['DiseaseMajorAlleleCount'], x['DiseaseMinorAlleleCount'])),
    nrow=2, ncol=2, byrow = T), correct = F)$p.value
})

# HW chi-sq test
print("Chi-sq Test for HW")
gwaResult$HWPValue <- apply(gwaResult, 1, function(x) {
  m <-  matrix(
    as.numeric(c(x['AACount'], x['ABCount'], x['BBCount'],
                 x['HWAA'], x['HWAB'], x['HWBB'])),
    nrow=2, ncol=3, byrow = T)
  hw.test.q <- sum((m[1,]-m[2,])^2 / m[2,])
  return(pchisq(hw.test.q, df=1, lower.tail = F))
})

# concatenate to the results we want
print("Select results table")
gwaResult <- gwaResult %>% select(Chrom, rsid, pos, MinorAllele, MajorAllele, DiseaseMinAlleleFreq, ControlMinAlleleFreq, OR, PValue, HWPValue)

print("Write out results to file")
write.table(gwaResult, file=paste0("./results/",diseaseName, "/gwa",  chrom, ".csv"), sep="\t")