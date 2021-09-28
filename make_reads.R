cmd_args=commandArgs(TRUE)
fasta <- cmd_args[1] # fasta file
load(cmd_args[2]) # granges file
n <- as.numeric(cmd_args[3]) # total number of samples
n_per_group <- n/2
libsize <- as.numeric(cmd_args[4]) # number of reads

library(GenomicRanges)
library(polyester)

# set random seed
set.seed(1)

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

simulate_experiment(fasta=fasta, outdir="reads",
                    num_reps=c(n_per_group, n_per_group),
                    fraglen=400, fragsd=25, readlen=150,
                    reads_per_transcript=counts, # read counts
                    size=100, # dispersion
                    fold_changes=fold_changes,
                    paired=TRUE)
