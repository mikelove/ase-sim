configfile: "config.json"

SALMON = "/proj/milovelab/bin/salmon-1.5.2_linux_x86_64/bin/salmon"

rule all:
    input: 
        summarized_experiment = "txp_allelic_se.rda",
        alignments = expand("align/{sample}.sorted.bam", sample=config["samples"]),
        seq = "data/drosophila_seq.h5"

rule make_expression:
    output:
        txps_fa = "transcripts.fa",
	granges = "granges.rda",
        ref = "data/drosophila_ref.fasta",
        chr = "data/drosophila_chr_2L.fasta",
        alt = "data/drosophila_alt_zero-based.tsv",
        vcf = "data/drosophila_chr_2L.vcf",
        haps = "data/drosophila_alt.haps"
    shell:
        "R CMD BATCH --no-save --no-restore '--args {output.txps_fa} {output.granges} "
        "{output.ref} {output.chr} {output.alt} {output.vcf} {output.haps}' make_expression.R"

rule make_reads:
    input:
        txps_fa = "transcripts.fa",
	granges = "granges.rda"
    params:
        n = len(config["samples"]),
	libsize = "50e6"
    output:
        expand("reads/{sample}_{read}.fasta", sample=config["samples"], read=config["reads"])
    shell:
        "R CMD BATCH --no-save --no-restore '--args {input.txps_fa} {input.granges} "
	"{params.n} {params.libsize}' make_reads.R"

rule shuffle:
    input:
        r1 = "reads/{sample}_1.fasta",
	r2 = "reads/{sample}_2.fasta"
    output:
        r1 = "reads/{sample}_1.shuffled.fa",
	r2 = "reads/{sample}_2.shuffled.fa"
    shell:
        "./shuffle.sh -l {input.r1} -r {input.r2}"

rule salmon_index:
    input: "transcripts.fa"
    output: directory("anno/salmon-index-1.5.2")
    params:
        threads = "12"
    shell: "{SALMON} index --keepDuplicates -p {params.threads} -t {input} -i {output}"

rule salmon_quant:
    input:
        r1 = "reads/{sample}_1.shuffled.fa",
        r2 = "reads/{sample}_2.shuffled.fa",
        index = "anno/salmon-index-1.5.2"
    output:
        "quants/{sample}/quant.sf"
    params:
        dir = "quants/{sample}",
        threads = "12",
	nboot = "30"
    shell:
        "{SALMON} quant -i {input.index} -l IU -p {params.threads} "
	"--numBootstraps {params.nboot} "
        "-o {params.dir} -1 {input.r1} -2 {input.r2}"

rule import_quants:
    input: expand("quants/{sample}/quant.sf", sample=config["samples"])
    params:
        nsamp = len(config["samples"])
    output: "txp_allelic_se.rda"
    shell:
        "R CMD BATCH --no-save --no-restore '--args {params.nsamp}' import_quants.R"

rule hisat_ss:
    input: "data/Drosophila_melanogaster.BDGP6.28.100.ss.gz"
    output: "data/Drosophila_melanogaster.BDGP6.28.100.ss"
    shell:
        "gunzip -c {input} > {output}"

rule hisat_index:
    input: 
        ref = "data/drosophila_ref.fasta",
        alt = "data/drosophila_alt_zero-based.tsv",
        haps = "data/drosophila_alt.haps",
        ss = "data/Drosophila_melanogaster.BDGP6.28.100.ss"
    output: "anno/bdgp6_sim/genome.1.ht2"
    params:
        threads = "12"
    shell:
        "hisat2-build -p {params.threads} -f --snp {input.alt} --haplotype {input.haps} --ss {input.ss} {input.ref} anno/bdgp6_sim/genome"

rule hisat_align:
    input:
        index = "anno/bdgp6_sim/genome.1.ht2",
        r1 = "reads/{sample}_1.shuffled.fa",
        r2 = "reads/{sample}_2.shuffled.fa"
    output: "align/{sample}.sam"
    params:
        threads = "12"
    shell:
        "hisat2 -p {params.threads} -f -x bdgp6_sim/genome -1 {input.r1} -2 {input.r2} -S {output}"

rule sort_alignments:
    input: "align/{sample}.sam"
    output: "align/{sample}.sorted.bam"
    params:
        unsorted = "align/{sample}.bam",
        threads = "12"
    shell:
        """
        samtools view -@ {params.threads} -bS {input} > {params.unsorted}
        samtools sort -@ {params.threads} {params.unsorted} -o {output}
        samtools index {output}
        rm {input}
        """

rule snp2h5:
    input: "data/drosophila_chr_2L.vcf"
    output:
        hap = "data/drosophila_hap.h5",
        index = "data/drosophila_snp_index.h5",
        tab = "data/drosophila_snp_tab.h5"
    shell:
        "snp2h5 --chrom data/drosophila_chromInfo.txt "
        "--format vcf --haplotype {output.hap} --snp_index {output.index} --snp_tab {output.tab} "
        "data/drosophila_*.vcf "

rule fasta_h5:
    input: "data/drosophila_chr_2L.fasta"
    output: "data/drosophila_seq.h5"
    shell: 
        "fasta2h5 --chrom data/drosophila_chromInfo.txt "
        "--seq {output} data/drosophila_chr_*.fasta"
