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
# also the identical file
identfile <- sub(".mmseq", ".identical.mmseq", mmseqfile)
iheader <- scan(identfile, what="char", n=2, sep="\n")
ident <- read.table(identfile, skip=2)
midx <- grepl("_M", ident$V1)
ident$V1 <- sub("_.","",ident$V1)
imfile <- sub(".mmseq", ".identical.mmseq", mfile)
ipfile <- sub(".mmseq", ".identical.mmseq", pfile)
write(iheader, imfile)
write.table(ident[midx,], file=imfile,
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)
write(iheader, ipfile)
write.table(ident[!midx,], file=ipfile,
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)
