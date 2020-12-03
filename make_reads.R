cmd_args=commandArgs(TRUE)
fasta <- cmd_args[1] # fasta file
load(cmd_args[2]) # granges file
n <- as.numeric(cmd_args[3]) # per group
libsize <- as.numeric(cmd_args[4]) # number of reads

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

library(polyester)
simulate_experiment(fasta=fasta, outdir=".",
                    num_reps=c(n,n),
                    reads_per_transcript=counts, # read counts
                    size=100, # dispersion
                    fold_changes=fold_changes,
                    paired=TRUE)

# writes out:
# sample_01_1.fasta
# sample_01_2.fasta
# sample_02_1.fasta
# sample_02_2.fasta
