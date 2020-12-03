cmd_args=commandArgs(TRUE)
outfile <- cmd_args[1]

library(AnnotationHub)
library(ensembldb)
library(dplyr)
library(ggplot2)
library(GenomicFeatures)
library(Biostrings)

# use AHub to obtain Drosophila genes and genome
ah <- AnnotationHub()

# Ensembl database 
q <- query(ah, c("EnsDb",101,"Drosophila"))
edb <- q[[1]]

# genome file for the Ensembl database
q <- query(ah, c("BDGP6.28",101,"TwoBit","dna_sm"))
genome <- q[[1]]

# STEP 1 - flip one basepair within each gene to make alt chromosomes

# exons by transcript and gene
ebt <- exonsBy(edb, by="tx")
ebg <- exonsBy(edb, by="gene")

# just standard chromosomes
chroms <- c("2L","2R","3L","3R","4","X","Y")
ebt <- keepSeqlevels(ebt, value=chroms, pruning.mode = "coarse")
ebg <- keepSeqlevels(ebg, value=chroms, pruning.mode = "coarse")

# pick an exon per gene, then take its middle position
exons <- unlist(ebg)
spl_ids <- split( mcols(exons)$exon_id, names(exons) )
which_id <- unname(sapply(spl_ids, function(ids) sample(ids, 1)))
which_exons <- exons[match(which_id, mcols(exons)$exon_id)]
at <- mid(which_exons)

# plot the positions
dat <- data.frame(at=at, chrom=seqnames(which_exons))
png(file="positions.png")
ggplot(dat, aes(at)) + geom_histogram() + facet_wrap(~chrom)
dev.off()

# split the positions by chromosome
spl_at <- split(at, seqnames(which_exons))

# make DNAStringSet of the genome (only standard chromosomes)
dna <- import(genome)[1:7]

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

# STEP 2 - come up with allelic read counts

# group transcripts by TSS
tbg <- transcriptsBy(edb, by="gene")
tss_pos <- start(resize(tbg, width=1))

ntxp <- lengths(tbg)
ntss <- sapply(tss_pos, function(x) length(unique(x)))

dat <- data.frame(ntxp, ntss)
png(file="ntss_over_ntxp.png")
dat %>% filter(ntxp < 10 & ntxp > 1) %>%
  ggplot(aes(ntxp, ntss)) + geom_jitter(alpha=.1) + geom_abline()
dev.off()

# how many genes are we working with
sum(ntxp %in% 2:5 & ntss < ntxp)

# put tss position in order of 'ebt'
tss_pos_vector <- unlist(tss_pos)
names(tss_pos_vector) <- unlist(tbg)$tx_name
tss_pos_vector <- tss_pos_vector[names(ebt)]

# start to build a suffix for the transcripts
suffix <- as(tss_pos_vector, "character")
names(suffix) <- names(ebt)

# this will be the altered abundance of alt transcripts
abundance <- rep(2, length(ebt))
names(abundance) <- names(ebt)

# alter expression of 1,000 genes
genes_to_alter <- sample(names(tbg)[ntxp %in% 2:5 & ntss < ntxp], 1000)
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

txps <- c(cdna, cdna_alt)

tbg_reorder <- unlist(tbg)
names(tbg_reorder) <- tbg_reorder$tx_name
tbg_reorder <- tbg_reorder[names(ebt)]
gene <- mcols(tbg_reorder)$gene_id
suffix <- paste0(suffix, "_", abundance, "_", gene)
suffix <- c(paste0(suffix, "_ref"), paste0(suffix, "_alt"))
names(txps) <- suffix

writeXStringSet(txps, file=outfile)
