configfile: "config.json"

SALMON = "~/bin/salmon-1.5.2_linux_x86_64/bin/salmon"

rule all:
    input:
        expand("quants/{sample}/quant.sf", sample=config["samples"])

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
	libsize = "1e5"
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
    shell: "{SALMON} index -p 2 -t {input} -i {output}"

rule salmon_quant:
    input:
        r1 = "reads/{sample}_1.shuffled.fa",
        r2 = "reads/{sample}_2.shuffled.fa",
        index = "anno/salmon-index-1.5.2"
    output:
        "quants/{sample}/quant.sf"
    params:
        dir = "quants/{sample}"
    shell:
        "{SALMON} quant -i {input.index} -l UI -p 2 "
        "-o {params.dir} -1 {input.r1} -2 {input.r2}"
