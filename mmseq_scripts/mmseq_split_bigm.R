cmd_args <- commandArgs(TRUE)
bigmfile <- cmd_args[1]
mfile <- cmd_args[2]
pfile <- cmd_args[3]
header <- strsplit(scan(bigmfile, what="char", n=1, sep="\n"), "\t")[[1]][-1]
sparse <- read.table(bigmfile, skip=1)
library(Matrix)
incidence <- sparseMatrix(i=sparse$V1, j=sparse$V2, index1=FALSE)
midx <- grepl("_M", header)
mtxp <- gsub("_.","",header[midx])
ptxp <- gsub("_.","",header[!midx])
mincidence <- incidence[,midx]
pincidence <- incidence[,!midx]
common <- intersect(mtxp, ptxp)
mincidence <- mincidence[,match(common, mtxp)]
pincidence <- pincidence[,match(common, ptxp)]
together <- mincidence + pincidence
together <- as(together, "dgTMatrix")
dat <- data.frame(i=together@i, j=together@j)
dat <- dat[order(dat$i),]
for (i in 1:2) {
  outfile <- c(mfile, pfile)[i]
  write(paste0("#\t",paste(mtxp, collapse="\t")), file=outfile)
  write.table(dat, file=outfile, append=TRUE, quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE)
}
