# ase-sim

## Simulation for creating allelic expression from fly genome

Author: Michael Love
Last modified: Nov 25 2021
Version: 0.0.3

## Software versions used:

* `r/4.1.0`
* `bowtie2/2.4.1`
* `samtools/1.13`
* `wasp/2019-12`

## Directories

The following directory structure is needed:

```
anno
bt2_align
data
quants
reads
wasp
```

## Provenance of included files

`simpleRepeat_dm6_Aug2014.bed` is from UCSC Genome Browser:

<https://genome.ucsc.edu/cgi-bin/hgTables?db=dm6&hgta_group=varRep&hgta_track=simpleRepeat&hgta_table=simpleRepeat&hgta_doSchema=describe+table+schema>

Bowtie2 index for BDGP6 was downloaded directly from the software website.
