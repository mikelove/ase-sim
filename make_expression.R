cmd_args <- commandArgs(TRUE)
fastafile <- cmd_args[1]
grangesfile <- cmd_args[2]
chrfile <- cmd_args[3]
vcffile <- cmd_args[4]
testtargetfile <- cmd_args[5]
mmseqfile <- cmd_args[6]
t2gfile <- cmd_args[7]
a2tfile <- cmd_args[8]

library(AnnotationHub)
library(ensembldb)
library(dplyr)
library(ggplot2)
library(GenomicFeatures)
library(Biostrings)
library(rtracklayer)

# use AHub to obtain Drosophila genes and genome
ah <- AnnotationHub()

# Ensembl database
release <- 100
q <- query(ah, c("EnsDb",release,"Drosophila"))
edb <- q[[1]]

# genome file for the Ensembl database
q <- query(ah, c("BDGP6.28",release,"TwoBit","dna_sm"))
genome <- q[[1]]

# define 'standard chromosomes', exluding Y
chroms <- c("2L","2R","3L","3R","4","X")

# STEP 1 - flip one basepair within each gene to make alt chromosomes

# exons by transcript and gene, transcripts by gene
g <- genes(edb)
ebt <- exonsBy(edb, by="tx")
ebg <- exonsBy(edb, by="gene")
tbg <- transcriptsBy(edb, by="gene")

# just standard chromosomes -- excluding Y
g <- keepSeqlevels(g, value=chroms, pruning.mode = "coarse")
ebt <- keepSeqlevels(ebt, value=chroms, pruning.mode = "coarse")
ebg <- keepSeqlevels(ebg, value=chroms, pruning.mode = "coarse")
tbg <- keepSeqlevels(tbg, value=chroms, pruning.mode = "coarse")

# subset to protein_coding ~14k and ncRNA ~2.4k
g <- g[sort(names(g))]
table(g$gene_biotype)
g <- g[g$gene_biotype %in% c("protein_coding","ncRNA")]

# filtering the 'by-gene' objects
ebg <- ebg[names(g)]
tbg <- tbg[names(g)]

# subset to genes not overlapping large simple tandem repeats (keeping ~14.8k genes)
rep <- import("data/simpleRepeat_dm6_Aug2014.bed.gz")
rep <- keepSeqlevels(rep, paste0("chr", chroms), pruning.mode="coarse")
seqlevels(rep) <- sub("chr","",seqlevels(rep))
big_rep <- rep[width(rep) >= 50]
over_rep <- names(ebg)[overlapsAny(ebg, big_rep, minoverlap=50)]

g <- g[!names(g) %in% over_rep]
ebg <- ebg[!names(ebg) %in% over_rep]
tbg <- tbg[!names(tbg) %in% over_rep]

# filter, and order 'ebt' by 'tbg'
ebt <- ebt[ mcols(unlist(tbg))$tx_name ]
stopifnot(length(ebt) == sum(lengths(tbg)))

# set random seed for the sampling below
set.seed(1)

# pick 'nexons' exon per gene, then take the middle positions
nexons <- 5
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

# split the positions per gene
at_per_gene <- as(split(at, names(which_exons)), "IntegerList")

# plot the positions
## dat <- data.frame(at=at, chrom=seqnames(which_exons))
## png(file="positions.png")
## ggplot(dat, aes(at)) + geom_histogram() + facet_wrap(~chrom)
## dev.off()

# split the positions by chromosome
spl_at <- split(at, seqnames(which_exons))

# make DNAStringSet of the genome (only standard chromosomes)
dna <- import(genome)[chroms]

# flip a letter within each gene
dna_alt <- lapply(names(dna), function(chr) {
  x <- spl_at[[chr]]
  # Views() is a lightweight way to pull out positions of a DNAString
  flip <- as( complement(Views(dna[[chr]], x, x)), "character" )
  replaceLetterAt(dna[[chr]], x, letter=flip)
})
names(dna_alt) <- names(dna)
dna_alt <- DNAStringSet(dna_alt)

# repeat to record ref and  alt alleles
ref_allele <- lapply(names(dna), function(chr) {
  x <- spl_at[[chr]]
  as( Views(dna[[chr]], x, x), "character" )
})
alt_allele <- lapply(names(dna), function(chr) {
  x <- spl_at[[chr]]
  as( complement(Views(dna[[chr]], x, x)), "character" )
})
names(ref_allele) <- names(alt_allele) <- names(dna)

# write a 2bit file for the alt chromosomes
rtracklayer::export(dna_alt, con="data/drosophila_alt.2bit")

# point to a 2bit file
genome_alt <- TwoBitFile("data/drosophila_alt.2bit")

# extract transcript sequence for both alleles
cdna <- extractTranscriptSeqs(genome, ebt)
cdna_alt <- extractTranscriptSeqs(genome_alt, ebt)

# STEP 2 - write out various genome-based files for downstream tools

# write out individual chroms for WASP
for (chr in names(dna)) {
  rtracklayer::export(dna[chr], con=sub("2L",chr,chrfile))
}

# alt alleles table
stopifnot(all(lengths(spl_at) == lengths(alt_allele)))
alt_chr <- rep(names(spl_at), lengths(spl_at))
alt_table <- data.frame(chr=alt_chr,
                        pos=unlist(spl_at),
                        ref=unlist(ref_allele),
                        alt=unlist(alt_allele))
alt_table <- alt_table[order(alt_table$chr, alt_table$pos),]
rownames(alt_table) <- paste0("rs",seq_along(alt_chr))

# VCF files, one per chrom, for WASP
vcf_table <- data.frame(CHROM=alt_table$chr, POS=alt_table$pos, ID=rownames(alt_table),
                        REF=alt_table$ref, ALT=alt_table$alt,
                        QUAL=".", FILTER="PASS", INFO=".", FORMAT="GT:GL",
                        sample="0|1:-100,0,-100")
for (chr in names(dna)) {
  write.table(vcf_table[vcf_table$CHROM == chr,], file=sub("2L",chr,vcffile),
              quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
  system(paste0("sed -i -e 's/CHROM/#CHROM/' ",sub("2L",chr,vcffile)))
}

# VCF file (just one) for WASP2
write.table(vcf_table, file=sub("chr_2L","wg",vcffile),
            quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
system(paste0("sed -i -e 's/CHROM/#CHROM/' ",sub("chr_2L","wg",vcffile)))

# the test SNP and target regions for WASP
red_ebg <- reduce(ebg)
left_pos <- sapply(at_per_gene, min)
wasp_test_target <- data.frame(
  chr=seqnames(g),
  pos=left_pos,
  endpos=left_pos,
  strand="+",
  name=paste0(seqnames(g),".",left_pos),
  tstart=sapply(start(red_ebg), paste0, collapse=";"),
  tend=sapply(end(red_ebg), paste0, collapse=";")
)
wasp_test_target <- merge(wasp_test_target, alt_table)
wasp_test_target <- wasp_test_target[,c("chr","pos","endpos","ref","alt",
                                        "strand","name","tstart","tend")]
write.table(wasp_test_target, file=testtargetfile, quote=FALSE,
            sep=" ", row.names=FALSE, col.names=FALSE)

# STEP 3 - create allelic read counts

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
# - half will be isoform level AI
# - half will be gene level AI
# set random seed for the sampling below
set.seed(1)
genes_to_alter <- sample(names(tbg)[txp_idx], 1000)
for (i in seq_along(genes_to_alter)) {
  if (i %% 100 == 0) print(i)
  g2a <- genes_to_alter[i]
  if (i <= 500) {
    # isoform AI - balanced within gene
    promoter <- sample(unique(tss_pos[[g2a]]), 1)
    idx <- mcols(tbg[[g2a]])$tx_name[tss_pos[[g2a]] == promoter]
    neg_idx <- mcols(tbg[[g2a]])$tx_name[tss_pos[[g2a]] != promoter]
    abundance[idx] <- abundance[idx] + 1/length(idx)
    abundance[neg_idx] <- abundance[neg_idx] - 1/length(neg_idx)
  } else {
    # gene AI
    idx <- mcols(tbg[[g2a]])$tx_name
    up_or_down <- c(2.5, 1.6)
    abundance[idx] <- rep(sample(up_or_down, 1), length(idx))
  }
}

# paternal allele is close to maternal
sum(abundance) / (2 * length(abundance))

# STEP 4 write out simulation files

# get gene id
tbg_reorder <- unlist(tbg)
names(tbg_reorder) <- tbg_reorder$tx_name
tbg_reorder <- tbg_reorder[ names(ebt) ]
gene <- mcols(tbg_reorder)$gene_id

# build a suffix for the transcript names
suffix <- paste0(tss_pos_vector, "_", round(abundance,2), "_", gene)

# save abundance on unlisted 'tbg'
txps <- unlist(tbg)
names(txps) <- mcols(txps)$tx_name
mcols(txps)$abundance <- abundance
mcols(txps)$tss <- tss_pos_vector
mcols(txps)$snp_loc <- at_per_gene[ mcols(txps)$gene_id ]

# save the AI status at txp level
mcols(txps)$isoAI <- FALSE
mcols(txps)$isoAI[ mcols(txps)$gene_id %in% genes_to_alter[1:500] ] <- TRUE
mcols(txps)$geneAI <- FALSE
mcols(txps)$geneAI[ mcols(txps)$gene_id %in% genes_to_alter[501:1000] ] <- TRUE

# save the AI status at gene level
mcols(g)$isoAI <- FALSE
mcols(g)[ genes_to_alter[1:500], "isoAI"] <- TRUE
mcols(g)$geneAI <- FALSE
mcols(g)[ genes_to_alter[501:1000], "geneAI"] <- TRUE

# width of cDNA
mcols(txps)$width <- width(cdna)

# write out GRanges (with abundance)
save(g, ebg, ebt, tbg, txps, genes_to_alter, file=grangesfile)

# combine the original and altered transcripts
cdna_both <- c(cdna, cdna_alt)

# transcript headers are: 'name_allele|suffix'
names(cdna_both) <- paste0(names(cdna_both), "_",
                           rep(c("M","P"),each=length(suffix)),
                           "|", suffix)

# write out FASTA
writeXStringSet(cdna_both, file=fastafile)

# write out the t2g and a2t files for terminus
terminus_t2g <- data.frame(txp=names(cdna_both),
                           gene=sub(".*(FBgn.*)","\\1",names(cdna_both)))
terminus_a2t <- data.frame(allele=names(cdna_both),
                           txp=sub("(FBtr.*?)_.*","\\1",names(cdna_both)))
write.table(terminus_t2g, file=t2gfile,
            quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(terminus_a2t, file=a2tfile,
            quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

# write out FASTA reference for mmseq
cdna_mmseq <- cdna_both
names(cdna_mmseq) <- paste0(names(cdna_mmseq), "_",
                           rep(c("M","P"),each=length(gene)),
                           " gene:", gene, "_",
                           rep(c("M","P"),each=length(gene)))
writeXStringSet(cdna_mmseq, file=mmseqfile)
writeXStringSet(cdna_mmseq[!duplicated(cdna_mmseq)],
                file=sub(".fa","_nodup.fa",mmseqfile))
