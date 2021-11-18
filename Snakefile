configfile: "config.json"

SALMON = "/proj/milovelab/bin/salmon-1.5.2_linux_x86_64/bin/salmon"

BAM2H5 = "python3.5 /nas/longleaf/apps/wasp/2019-12/WASP/CHT/bam2h5.py"
EXTRACTHAP = "python3.5 /nas/longleaf/apps/wasp/2019-12/WASP/CHT/extract_haplotype_read_counts.py"
UPDATEDEP = "python3.5 /nas/longleaf/apps/wasp/2019-12/WASP/CHT/update_total_depth.py"
UPDATEHET = "python3.5 /nas/longleaf/apps/wasp/2019-12/WASP/CHT/update_het_probs.py"

rule all:
    input: 
        summarized_experiment = "txp_allelic_se.rda",
        alignments = expand("align/{sample}.bam", sample=config["samples"]),
        wasp = expand("wasp/hap_read_counts.{sample}.hetp", sample=config["samples"]),

rule make_expression:
    output:
        txps_fa = "transcripts.fa",
	granges = "granges.rda",
        ref = "data/drosophila_ref.fasta",
        chr = "data/drosophila_chr_2L.fasta",
        alt = "data/drosophila_alt_zero-based.tsv",
        vcf = "data/drosophila_chr_2L.vcf",
        haps = "data/drosophila_alt.haps",
        tt = "data/drosophila_test_target.txt"
    shell:
        "R CMD BATCH --no-save --no-restore '--args {output.txps_fa} {output.granges} "
        "{output.ref} {output.chr} {output.alt} {output.vcf} {output.haps} {output.tt}' "
        "make_expression.R"

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
        "hisat2-build-s -p {params.threads} -f --snp {input.alt} --haplotype {input.haps} "
        "--ss {input.ss} {input.ref} anno/bdgp6_sim/genome"

rule hisat_align:
    input:
        index = "anno/bdgp6_sim/genome.1.ht2",
        r1 = "reads/{sample}_1.shuffled.fa",
        r2 = "reads/{sample}_2.shuffled.fa"
    output: "align/{sample}.bam"
    params:
        threads = "12",
        mem_per_thread = "1G"
    shell:
        """
        hisat2 -p {params.threads} -f -x anno/bdgp6_sim/genome \
          -1 {input.r1} -2 {input.r2} | samtools sort -@ {params.threads} \
          -m {params.mem_per_thread} -o {output}
        samtools index {output}
        """

rule filter_and_tally_alignments:
    input: "align/{sample}.bam"
    output: 
        lowqual = "align/{sample}.lowqual",
        filter = "align/{sample}.filt.bam"
    params:
        threads = "12",
        mem_per_thread = "1G"
    shell:
        """
        samtools view {input} | awk '$5 < 60 {{print $1}}' | grep -o '_.|' | \
          sed 's/[_|]//g' | sort | uniq -c > {output.lowqual}
        samtools view -q 60 -@ {params.threads} -m {params.mem_per_thread} -b {input} > {output.filter}
        samtools index {output.filter}
        """

rule wasp_snp2h5:
    input: "data/drosophila_chr_2L.vcf"
    output:
        genoprob = "data/drosophila_genoprob.h5",
        index = "data/drosophila_snp_index.h5",
        tab = "data/drosophila_snp_tab.h5",
        hap = "data/drosophila_haps.h5"
    shell:
        "snp2h5 --chrom data/drosophila_chromInfo.txt --format vcf "
        "--geno_prob {output.genoprob} --snp_index {output.index} "
        "--snp_tab {output.tab} --haplotype {output.hap} "
        "data/drosophila_*.vcf "

rule wasp_fasta_h5:
    input: "data/drosophila_chr_2L.fasta"
    output: "data/drosophila_seq.h5"
    shell: 
        "fasta2h5 --chrom data/drosophila_chromInfo.txt "
        "--seq {output} data/drosophila_chr_*.fasta"

rule wasp_read_count:
    input:
        index = "data/drosophila_snp_index.h5",
        tab = "data/drosophila_snp_tab.h5",
        bam = "align/{sample}.filt.bam"
    output:
        ref = "wasp/ref_as_counts.{sample}.h5",
        alt = "wasp/alt_as_counts.{sample}.h5",
        other = "wasp/other_as_counts.{sample}.h5",
        count = "wasp/read_counts.{sample}.h5"
    shell:
        """
        {BAM2H5} --chrom data/drosophila_chromInfo.txt \
         --snp_index {input.index} --snp_tab {input.tab} \
         --ref_as_counts {output.ref} --alt_as_counts {output.alt} \
         --other_as_counts {output.other} --read_counts {output.count} \
         {input.bam} 2> /dev/null
         """

rule wasp_extract:
    input:
        genoprob = "data/drosophila_genoprob.h5",
        hap = "data/drosophila_haps.h5",
        index = "data/drosophila_snp_index.h5",
        tab = "data/drosophila_snp_tab.h5",
        ref = "wasp/ref_as_counts.{sample}.h5",
        alt = "wasp/alt_as_counts.{sample}.h5",
        other = "wasp/other_as_counts.{sample}.h5",
        count = "wasp/read_counts.{sample}.h5",
        tt = "data/drosophila_test_target.txt"
    output: "wasp/hap_read_counts.{sample}.txt"
    shell:
        """
        {EXTRACTHAP}  --chrom data/drosophila_chromInfo.txt \
         --geno_prob {input.genoprob} --haplotype {input.hap} \
         --snp_index {input.index} --snp_tab {input.tab} \
         --samples data/samples --individual sample \
         --ref_as_counts {input.ref} --alt_as_counts {input.alt} \
         --other_as_counts {input.other} --read_counts {input.count} \
         {input.tt} > {output}
        """

rule wasp_adjust_read_count:
    input:
        seq = "data/drosophila_seq.h5",
        hap_counts = expand("wasp/hap_read_counts.{sample}.txt", sample=config["samples"])
    output:
        expand("wasp/hap_read_counts.{sample}.adj", sample=config["samples"])
    shell:
        """
        ls -1 {input.hap_counts} > wasp/preadj
        sed 's/.txt/.adj/' wasp/preadj > wasp/postadj
        {UPDATEDEP} --seq {input.seq} wasp/preadj wasp/postadj
        """

rule wasp_adjust_het_prob:
    input:
        ref = "wasp/ref_as_counts.{sample}.h5",
        alt = "wasp/alt_as_counts.{sample}.h5",
        adj = "wasp/hap_read_counts.{sample}.adj"
    output: "wasp/hap_read_counts.{sample}.hetp"
    shell:
        "{UPDATEHET} --ref_as_counts {input.ref} --alt_as_counts {input.alt} "
        "{input.adj} {output}"
