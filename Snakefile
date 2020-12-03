configfile: "config.json"

rule all:
    input:
        expand("{sample}_{read}.shuffled.fa", sample=config["samples"], read=config["reads"])

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
        n = "1",
	libsize = "1e5"
    output:
        s1r1 = "sample_01_1.fasta",
	s1r2 = "sample_01_2.fasta",
        s2r1 = "sample_02_1.fasta",
	s2r2 = "sample_02_2.fasta"
    shell:
        "R CMD BATCH --no-save --no-restore '--args {input.txps_fa} {input.granges} "
	"{params.n} {params.libsize}' make_reads.R"

rule shuffle:
    input:
        r1 = "{sample}_1.fasta",
	r2 = "{sample}_2.fasta"
    output:
        r1 = "{sample}_1.shuffled.fa",
	r2 = "{sample}_2.shuffled.fa"
    shell:
        "./shuffle.sh -l {input.r1} -r {input.r2}"
