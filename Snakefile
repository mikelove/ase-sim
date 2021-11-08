configfile: "config.json"

SALMON = "/proj/milovelab/bin/salmon-1.5.2_linux_x86_64/bin/salmon"

rule all:
    input: 
        summarized_experiment = "txp_allelic_se.rda",
        alignments = expand("align/{sample}.sorted.bam", sample=config["samples"])

rule make_expression:
    output:
        txps_fa = "transcripts.fa",
	granges = "granges.rda"
    shell:
        "R CMD BATCH --no-save --no-restore '--args {output.txps_fa} {output.granges}' make_expression.R"

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

rule hisat_index:
    output: "bdgp6/genome.1.ht2"
    shell:
        """
        wget https://genome-idx.s3.amazonaws.com/hisat/bdgp6.tar.gz
        tar -xvf bdgp6.tar.gz
        rm -f bdgp6.tar.gz
        """

rule hisat_align:
    input:
        index = "bdgp6/genome.1.ht2",
        r1 = "reads/{sample}_1.shuffled.fa",
        r2 = "reads/{sample}_2.shuffled.fa"
    output: "align/{sample}.sam"
    params:
        threads = "12"
    shell:
        "hisat2 -p {params.threads} -f -x bdgp6/genome -1 {input.r1} -2 {input.r2} -S {output}"

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
        """
