cmd_args <- commandArgs(TRUE)
tracefile <- cmd_args[1]
mfile <- cmd_args[2]
pfile <- cmd_args[3]
library(readr)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)
trace <- read_delim(tracefile, col_names=FALSE, skip=1)
# the last column is extra
trace <- trace[,-ncol(trace)]
names <- scan(tracefile, what="char", n=ncol(trace))
for (i in 1:2) {
  outfile <- sub(".gz","",c(mfile, pfile)[i])
  idx <- grep(c("M","P")[i], names)
  tracenames <- names[idx]
  tracenames <- sub("_.","",tracenames)
  write(tracenames, file=outfile, ncolumns=length(idx))
  write.table(trace[,idx], file=outfile,
              quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", append=TRUE)
  R.utils::gzip(outfile)
}
# also the identical trace file
identfile <- sub(".trace_gibbs.gz", ".identical.trace_gibbs.gz", tracefile)
itrace <- read_delim(identfile, col_names=FALSE, skip=1)
itrace <- itrace[,-ncol(itrace)]
inames <- scan(identfile, what="char", n=ncol(itrace))
for (i in 1:2) {
  ioutfile <- sub(".trace_gibbs.gz",".identical.trace_gibbs",c(mfile, pfile)[i])
  idx <- grep(c("M","P")[i], inames)
  itracenames <- inames[idx]
  itracenames <- gsub("_.","",itracenames)
  write(itracenames, file=ioutfile, ncolumns=length(idx))
  write.table(itrace[,idx], file=ioutfile,
              quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", append=TRUE)
  R.utils::gzip(ioutfile)
}
