cmd_args <- commandArgs(TRUE)
fasta <- cmd_args[1] # fasta file
load(cmd_args[2]) # granges file
libsize <- as.numeric(cmd_args[3]) # number of reads
pair <- cmd_args[4] # "A", "B", etc.

library(GenomicRanges)
library(polyester)

# set random seed, "A" -> 1, "B" -> 2, etc.
set.seed(match(pair, LETTERS))

# from abundance to read counts
abundance <- mcols(txps)$abundance

# combine the ref (constant 2) with the alt allelic abundance
abundance <- c(rep(2,length(abundance)), abundance)
nucfrac <- abundance * mcols(txps)$width
nucfrac <- nucfrac / sum(nucfrac)
counts <- nucfrac * libsize
counts[counts < 1] <- 0 # zero out when expected count is < 1

# no fold changes across sample
fold_changes <- matrix(1, nrow=2 * length(txps), ncol=2)

# just 1 per group, for Snakemake
n_per_group <- 1

simulate_experiment(fasta=fasta, outdir=file.path("reads",pair),
                    num_reps=c(n_per_group, n_per_group), # 1 per group
                    fraglen=400, fragsd=25, # long fragments
                    readlen=150, # 2x 150bp reads
                    reads_per_transcript=counts, # read counts
                    size=100, # 1/dispersion
                    fold_changes=fold_changes, # no DE -> FC all are 1
                    paired=TRUE)
