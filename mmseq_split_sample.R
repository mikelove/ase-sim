cmd_args <- commandArgs(TRUE)
mmseqfile <- cmd_args[1]
mfile <- cmd_args[2]
pfile <- cmd_args[3]
header <- scan(mmseqfile, what="char", n=2, sep="\n")
mmseq <- read.table(mmseqfile, skip=2)
n <- nrow(mmseq)/2
stopifnot(all(grepl("M", mmseq$V1[1:n])))
stopifnot(all(grepl("P", mmseq$V1[(n+1):(2*n)])))
mmseq$V1 <- sub("_.","",mmseq$V1)
write(header, mfile)
write.table(mmseq[1:n,], file=mfile,
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)
write(header, pfile)
write.table(mmseq[(n+1):(2*n),], file=pfile,
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)
