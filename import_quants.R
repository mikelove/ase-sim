cmd_args <- commandArgs(TRUE)
nsamp <- as.numeric(cmd_args[1])
library(readr)
library(fishpond)
library(SummarizedExperiment)
samps <- list.files("quants")
stopifnot(length(samps) == nsamp) # check that we have same number
coldata <- data.frame(files=file.path("quants",samps,"quant.sf"),
                      names=samps)
# txp-level data
se <- importAllelicCounts(coldata, a1="P", a2="M", format="wide")                      
save(se, file="txp_allelic_se.rda")
