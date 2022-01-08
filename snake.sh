#!/bin/bash
#
#SBATCH --job-name=snake
#SBATCH --time=720
#SBATCH --mem=1000
#SBATCH --mail-user=milove@email.unc.edu
#SBATCH --mail-type=ALL

module load python
module load r/4.1.0
#module load pigz
#module load hisat2/2.2.1
#module load wasp
#module load samtools
#module load bowtie
snakemake -j 10 --latency-wait 30 --cluster "sbatch --mem=5000 -N 1 -n 1 --time=30"
