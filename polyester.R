library(polyester)
simulate_experiment(fasta=fasta, outdir=".",
                    num_reps=c(1,1),
                    reads_per_transcript=300,
                    size=100,
                    fold_changes=fold_changes,
                    paired=TRUE)
