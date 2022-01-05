cmd_args <- commandArgs(TRUE)
mmseqfile <- cmd_args[1]
mfile <- cmd_args[2]
pfile <- cmd_args[3]
library(readr)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)
trace <- read_delim(mmseqfile, col_names=FALSE, skip=1)
# the last column is extra
trace <- trace[,-ncol(trace)]
names <- scan(mmseqfile, what="char", n=ncol(trace))
for (i in 1:2) {
  outfile <- sub(".gz","",c(mfile, pfile)[i])
  idx <- grep(c("M","P")[i], names)
  write(names[idx], file=outfile, ncolumns=length(idx))
  write.table(trace[,idx], file=outfile,
              quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", append=TRUE)
  R.utils::gzip(outfile)
}
