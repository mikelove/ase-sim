cmd_args=commandArgs(TRUE)
n <- as.numeric(cmd_args[1]) # per group

library(polyester)
fold_changes <- matrix(1, nrow=length(reads), ncol=2)
simulate_experiment(fasta=fasta, outdir=".",
                    num_reps=c(n,n),
                    reads_per_transcript=300, # read counts
                    size=100, # dispersion
                    fold_changes=fold_changes,
                    paired=TRUE)
