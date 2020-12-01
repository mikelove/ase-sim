library(AnnotationHub)
ah <- AnnotationHub()

q <- query(ah, c("EnsDb",101,"Drosophila"))
edb <- q[[1]] # Ensembl database 

q <- query(ah, c("BDGP6.28",101,"TwoBit","dna_rm"))
genome <- q[[1]] # genome file

dna <- import(genome)[1:7]
# flip a letter:
# replaceLetterAt(dna, at, letter=as.character(complement(subseq(dna, at, at))))
rtracklayer::export(dna, con="~/Desktop/test.2bit")
genome_alt <- TwoBitFile("~/Desktop/test.2bit")

ebt <- exonsBy(edb, by="tx")
chroms <- c("2L","2R","3L","3R","4","X","Y")
ebt <- keepSeqlevels(ebt, value=chroms, pruning.mode = "coarse")

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
