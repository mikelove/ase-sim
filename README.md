# ase-sim

## Simulation for creating allelic expression from fly genome

Author: Michael Love
Last modified: Nov 12 2021
Version: 0.0.3

## Software versions used:

* `anaconda/2019.10` (for snakemake)
* `r/4.1.0`
* `hisat2/2.2.1`
* `samtools/1.13`

The following was used to build the splice site file that is included in the repo.

```
wget ftp://ftp.ensembl.org/pub/release-100/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.28.100.chr.gtf.gz
gunzip Drosophila_melanogaster.BDGP6.28.100.chr.gtf.gz
hisat2_extract_splice_sites.py Drosophila_melanogaster.BDGP6.28.100.chr.gtf > Drosophila_melanogaster.BDGP6.28.100.ss
```
