cmd_args=commandArgs(TRUE)
fastafile <- cmd_args[1]
grangesfile <- cmd_args[2]

library(AnnotationHub)
library(ensembldb)
library(dplyr)
library(ggplot2)
library(GenomicFeatures)
library(Biostrings)

# set random seed
set.seed(1)

# use AHub to obtain Drosophila genes and genome
ah <- AnnotationHub()

# Ensembl database
release <- 100
q <- query(ah, c("EnsDb",release,"Drosophila"))
edb <- q[[1]]

# genome file for the Ensembl database
q <- query(ah, c("BDGP6.28",release,"TwoBit","dna_sm"))
genome <- q[[1]]

# STEP 1 - flip one basepair within each gene to make alt chromosomes

# exons by transcript and gene, transcripts by gene
ebt <- exonsBy(edb, by="tx")
ebg <- exonsBy(edb, by="gene")
tbg <- transcriptsBy(edb, by="gene")

# just standard chromosomes -- excluding Y
chroms <- c("2L","2R","3L","3R","4","X")
ebt <- keepSeqlevels(ebt, value=chroms, pruning.mode = "coarse")
ebg <- keepSeqlevels(ebg, value=chroms, pruning.mode = "coarse")
tbg <- keepSeqlevels(tbg, value=chroms, pruning.mode = "coarse")

# order 'ebt' by 'tbg'
ebt <- ebt[ mcols(unlist(tbg))$tx_name ]

# pick 'nexons' exon per gene, then take the middle positions
nexons <- 3
exons <- unlist(ebg)
spl_ids <- split( mcols(exons)$exon_id, names(exons) )
which_id <- unname(unlist(lapply(spl_ids, function(ids) {
  if (length(ids) <= nexons) {
    ids
  } else {
    sample(ids, nexons, replace=FALSE)
  }
})))
nexons * length(ebg)
length(which_id) # less than nexons x ngenes, bc replace=FALSE

# grab the exons
which_exons <- exons[match(which_id, mcols(exons)$exon_id)]
# grab the middle position
at <- mid(which_exons)

# plot the positions
## dat <- data.frame(at=at, chrom=seqnames(which_exons))
## png(file="positions.png")
## ggplot(dat, aes(at)) + geom_histogram() + facet_wrap(~chrom)
## dev.off()

# split the positions by chromosome
spl_at <- split(at, seqnames(which_exons))

# make DNAStringSet of the genome (only standard chromosomes)
dna <- import(genome)[c("2L","2R","3L","3R","4","X")]

# flip a letter within each gene
dna_alt <- lapply(names(dna), function(c) {
  x <- spl_at[[c]]
  # Views() is a lightweight way to pull out positions of a DNAString
  flip <- as( complement(Views(dna[[c]], x, x)), "character" )
  replaceLetterAt(dna[[c]], x, letter=flip)
})
names(dna_alt) <- names(dna)
dna_alt <- DNAStringSet(dna_alt)

# write a 2bit file for the alt chromosomes
rtracklayer::export(dna_alt, con="drosophila_alt.2bit")

# point to a 2bit file
genome_alt <- TwoBitFile("drosophila_alt.2bit")

# extract transcript sequence for both alleles
cdna <- extractTranscriptSeqs(genome, ebt)
cdna_alt <- extractTranscriptSeqs(genome_alt, ebt)

# STEP 2 - create allelic read counts

# group transcripts by TSS
tss_pos <- start(resize(tbg, width=1))

ntxp <- lengths(tbg)
ntss <- sapply(tss_pos, function(x) length(unique(x)))

## dat <- data.frame(ntxp, ntss)
## png(file="ntss_over_ntxp.png")
## dat %>% filter(ntxp < 10 & ntxp > 1) %>%
##   ggplot(aes(ntxp, ntss)) + geom_jitter(alpha=.1) + geom_abline()
## dev.off()

# how many genes are we working with
txp_idx <- ntxp %in% 3:6 & ntss < ntxp & ntss > 1
sum(txp_idx)

# put tss position in order of 'ebt'
tss_pos_vector <- unlist(tss_pos)
names(tss_pos_vector) <- unlist(tbg)$tx_name
tss_pos_vector <- tss_pos_vector[names(ebt)]

# this will be the altered abundance of alt transcripts
abundance <- rep(2, length(ebt))
names(abundance) <- names(ebt)

# alter expression of 1,000 genes
genes_to_alter <- sample(names(tbg)[txp_idx], 1000)
for (i in seq_along(genes_to_alter)) {
  if (i %% 100 == 0) print(i)
  g <- genes_to_alter[i]
  promoter <- sample(unique(tss_pos[[g]]), 1)
  idx <- mcols(tbg[[g]])$tx_name[tss_pos[[g]] == promoter]
  neg_idx <- mcols(tbg[[g]])$tx_name[tss_pos[[g]] != promoter]
  abundance[idx] <- abundance[idx] + 1/length(idx)
  abundance[neg_idx] <- abundance[neg_idx] - 1/length(neg_idx)
}

# STEP 3 write out files

# combine the original and altered transcripts
cdna_both <- c(cdna, cdna_alt)

# get gene id
tbg_reorder <- unlist(tbg)
names(tbg_reorder) <- tbg_reorder$tx_name
tbg_reorder <- tbg_reorder[ names(ebt) ]
gene <- mcols(tbg_reorder)$gene_id

# get the mutated basepair per gene
names(at) <- names(which_exons)
at_along_txps <- at[ gene ]

# build a suffix for the transcript names
suffix <- paste0(tss_pos_vector, "_", round(abundance,2), "_", gene, "_", at_along_txps)

# transcript headers are: 'name_allele|suffix'
names(cdna_both) <- paste0(names(cdna_both), "_",
                           rep(c("M","P"),each=length(suffix)),
                           "|", suffix)

# save abundance on unlisted 'tbg'
txps <- unlist(tbg)
names(txps) <- mcols(txps)$tx_name
mcols(txps)$abundance <- abundance
mcols(txps)$tss <- tss_pos_vector
mcols(txps)$width <- width(cdna)
mcols(txps)$snp_loc <- at_along_txps
# FASTA and GRanges (with abundance)
writeXStringSet(cdna_both, file=fastafile)
save(ebt, tbg, txps, file=grangesfile)
