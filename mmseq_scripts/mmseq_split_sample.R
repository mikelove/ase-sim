cmd_args <- commandArgs(TRUE)
mmseqfile <- cmd_args[1]
mfile <- cmd_args[2]
pfile <- cmd_args[3]
header <- scan(mmseqfile, what="char", n=2, sep="\n")
mmseq <- read.table(mmseqfile, skip=2)
midx <- grepl("M", mmseq$V1)
pidx <- grepl("P", mmseq$V1)
mmseq$V1 <- sub("_.","",mmseq$V1)
for (i in 1:2) {
  outfile <- c(mfile, pfile)[i]
  idx <- list(midx, pidx)[[i]]
  write(header, outfile)
  write.table(mmseq[idx,], file=outfile,
              append=TRUE, quote=FALSE, sep="\t",
              row.names=FALSE, col.names=FALSE)
}
