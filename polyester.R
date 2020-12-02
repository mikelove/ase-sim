# use AHub to obtain Drosophila genes and genome
library(AnnotationHub)
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
library(ggplot2)
dat <- data.frame(at=at, chrom=seqnames(which_exons))
ggplot(dat, aes(at)) + geom_histogram() + facet_wrap(~chrom)

# split the positions by chromosome
spl_at <- split(at, seqnames(which_exons))

# make DNAStringSet of the genome (only standard chromosomes)
dna <- import(genome)[1:7]

# flip a letter within each gene
library(Biostrings)
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
library(GenomicFeatures)
cdna <- extractTranscriptSeqs(genome, ebt)
cdna_alt <- extractTranscriptSeqs(genome_alt, ebt)

# STEP 2 - come up with allelic read counts

# group transcripts by TSS
tbg <- transcriptsBy(edb, by="gene")
tss_pos <- start(resize(tbg, width=1))

ntxp <- lengths(tbg)
ntss <- sapply(tss_pos, function(x) length(unique(x)))

library(dplyr)
library(ggplot2)
dat <- data.frame(ntxp, ntss)
dat %>% filter(ntxp < 10 & ntxp > 1) %>%
  ggplot(aes(ntxp, ntss)) + geom_jitter(alpha=.1) + geom_abline()

# how many genes are we working with
sum(ntxp %in% 2:5 & ntss < ntxp)

# this will be the altered abundance of alt transcripts
abundance <- rep(2, length(ebt))
names(abundance) <- names(ebt)

# put tss position in order of 'ebt'
tss_pos_vector <- unlist(tss_pos)
names(tss_pos_vector) <- unlist(tbg)$tx_name
tss_pos_vector <- tss_pos_vector[names(ebt)]

# start to build a suffix for the transcripts
suffix <- as(tss_pos_vector, "character")
names(suffix) <- names(ebt)

genes_to_alter <- names(tbg)[ntxp %in% 2:5 & ntss < ntxp]
for (g in genes_to_alter) {
  promoter <- sample(unique(tss_pos[[g]]), 1)
  idx <- mcols(tbg[[g]])$tx_name[tss_pos[[g]] == promoter]
  neg_idx <- mcols(tbg[[g]])$tx_name[tss_pos[[g]] != promoter]
  abundance[idx] <- abundance[idx] + 1/length(idx)
  abundance[neg_idx] <- abundance[neg_idx] - 1/length(neg_idx)
}

library(polyester)
fold_changes <- matrix(1, nrow=length(reads), ncol=2)
simulate_experiment(fasta=fasta, outdir=".",
                    num_reps=c(1,1),
                    reads_per_transcript=300, # read counts
                    size=100, # dispersion
                    fold_changes=fold_changes,
                    paired=TRUE)
