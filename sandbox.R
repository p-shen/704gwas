
gwa5.imm <- read.delim('./imm_bk/gwa5.csv', sep="\t", header=T)
gwa5.res <- read.delim('./results/gwa5.csv', sep="\t", header=T)
snps <- read.delim('imm_bk/snps/snps_22', sep="\t", header=F)

gwa22.imm <- read.delim('./imm_bk/gwa22.csv', sep="\t", header=T)
gwa22.res <- read.delim('./results/gwa22.csv', sep="\t", header=T)


test.table <- bind_rows(apply(control[1:5,-1:-4], 1, function(x){
  as.list(table(t(x)))
}))

test.table.counts <- bind_cols(test.table, bind_rows(apply(test.table, 1, function(x){
  cols <- names(x)[!is.na(x)]
  data <- x[cols]
  uniqchars <- unique(strsplit(paste(cols, collapse = ''), "")[[1]])
  uniqchars <- uniqchars[uniqchars!=" "]
  counts <- sapply(uniqchars, function(y) sum(data*str_count(cols, y)))
  names(counts) <- uniqchars
  
  # determine major and minor alleles
  majorAllele <- names(which(counts==max(counts)))
  minorAllele <- names(which(counts==min(counts)))
  majorAlleleCount <- counts[[majorAllele]]
  minorAlleleCount <- counts[[minorAllele]]
  
  as.list(c(`MajorAllele`=majorAllele, `MinorAllele`=minorAllele, `MajorAlleleCount`=majorAlleleCount, `MinorAlleleCount`=minorAlleleCount))
})))

print(append(uniqchars, ret))

snpTable.colNames <- unlist(strsplit(colnames(snpTable[,2]), " "))
snpTable.additive <- snpTable %>% transmute(a1=.[[1]]*2+.[[2]], a2=.[[2]]+.[[3]]*2)
colnames(snpTable.additive) <- snpTable.colNames

# determine major and minor alleles
majorAllele <- max.col(snpTable.additive[1,])
minorAllele <- colnames(snpTable.additive[1,-majorAllele])
majorAllele <- colnames(snpTable.additive[1,majorAllele])
