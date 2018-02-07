#!/usr/bin/env Rscript
# load a library for option parsing
library(optparse)

# define the options that the script should support
opt.list <- list(
  make_option(c("-c1","--control1"),type='character',help="file name for the 1st control sample"),
  make_option(c("-c2","--control2"),type='character',help="file name for the 2nd control sample"),
  make_option(c("-d", "--disease"), type='character',help="file name for the diseased sample"),
  make_option(c("-n", "--chromosome"), type='integer',help="chromosome number")
)

# do the argument parsing
pa <- parse_args(OptionParser(option_list=opt.list),positional_arguments=T)
opt <- pa$options;

control1 <- opt$control1
control2 <- opt$control2
disease <- opt$disease
chrom <- opt$chromosome

# check to see if the arguments are both provided
if(is.null(control1) || is.null(control2) || is.null(disease) || is.null(chrom)) {
  stop("Control files or disease file or chromosome number is not provided. Run 'Rscript execute_gwas.r -h' to see the help menu.")
}

# does it exist?
if (!file.exists(control1)) {
  print(paste(control1, " does not exist"))
} else if (!file.exists(control2)) {
  print(paste(control2, " does not exist"))
}else if (!file.exists(disease)) {
  print(paste(disease, " does not exist"))
}

#----------------------------
# Start processing GWAS data
library(tidyverse)


quit()