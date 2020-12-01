# use AHub to obtain Drosophila genes and genome
library(AnnotationHub)
ah <- AnnotationHub()

# Ensembl database 
q <- query(ah, c("EnsDb",101,"Drosophila"))
edb <- q[[1]]

# genome file for the Ensembl database
q <- query(ah, c("BDGP6.28",101,"TwoBit","dna_sm"))
genome <- q[[1]]

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
library(ggplot2)
dat <- data.frame(at=at, chrom=seqnames(which_exons))
ggplot(dat, aes(at)) + geom_histogram() + facet_wrap(~chrom)

# split the positions by chromosome
spl_at <- split(at, seqnames(which_exons))

# DNAStringSet of the genome (only standard chromosomes)
dna <- import(genome)[1:7]
# flip a letter:
library(Biostrings)
dna_alt <- lapply(names(dna), function(c) {
  x <- spl_at[[c]]
  flip <- as(complement(Views(dna[[c]], x, x)), "character")
  replaceLetterAt(dna[[c]], x, letter=flip)
})
names(dna_alt) <- names(dna)
dna_alt <- DNAStringSet(dna_alt)
# write a 2bit file
rtracklayer::export(dna, con="~/Desktop/test.2bit")
# point to a 2bit file
genome_alt <- TwoBitFile("~/Desktop/test.2bit")

# extract transcript sequence for both alleles
library(GenomicFeatures)
cdna <- extractTranscriptSeqs(genome, ebt)
cdna <- extractTranscriptSeqs(genome_alt, ebt)

library(polyester)
simulate_experiment(fasta=fasta, outdir=".",
                    num_reps=c(1,1),
                    reads_per_transcript=300,
                    size=100,
                    fold_changes=fold_changes,
                    paired=TRUE)
