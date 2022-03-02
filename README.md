# ase-sim

## Simulation for creating allelic expression from fly genome

Author: Michael Love and Euphy Wu
Last modified: Feb 26 2022
Version: 0.0.6

## Software versions used:

* `r/4.1.0`
* `hisat2/2.2.1`
* `samtools/1.13`
* `wasp/2019-12`
* `bowtie/1.3.1`
* `mmseq/1.0.10a`
* `anaconda/2019.10`

## Directories

The following directory structure is needed:

```
anno
data
ht2_align
mmseq
quants
reads
wasp_cht
wasp_mapping
wasp2_counts
wasp2_tests
```

## Provenance of included files

`simpleRepeat_dm6_Aug2014.bed` is from UCSC Genome Browser:

<https://genome.ucsc.edu/cgi-bin/hgTables?db=dm6&hgta_group=varRep&hgta_track=simpleRepeat&hgta_table=simpleRepeat&hgta_doSchema=describe+table+schema>

The Drosophila GTF file is downloaded from:

<http://ftp.ensembl.org/pub/release-100/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.28.100.chr.gtf.gz>
