
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
