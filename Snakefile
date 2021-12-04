configfile: "config.json"

SALMON = "/proj/milovelab/bin/salmon-1.5.2_linux_x86_64/bin/salmon"

CHT = "python3.5 /nas/longleaf/apps/wasp/2019-12/WASP/CHT"

MAPPING = "python3.5 /nas/longleaf/apps/wasp/2019-12/WASP/mapping"

rule all:
    input: 
        reads = expand("reads/sample_{pair}_{sample}_{read}.fasta", 
                       pair=config["pairs"], sample=config["samples"], 
                       read=config["reads"]),
        summarized_experiment = "txp_allelic_se.rda",
        hisat = expand("ht2_align/sample_{pair}_{sample}.bam",
                       pair=config["pairs"], sample=config["samples"]),
        wasp = expand("wasp_cht/ref_as_counts.sample_{pair}_{sample}.h5",
                      pair=config["pairs"], sample=config["samples"]),
        wasp_result = "wasp_cht/cht_results.txt"

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
        libsize = "50e6",
        pair = "{pair}"
    output:
        r11 = "reads/sample_{pair}_1_1.fasta",
        r12 = "reads/sample_{pair}_1_2.fasta",
        r21 = "reads/sample_{pair}_2_1.fasta",
        r22 = "reads/sample_{pair}_2_2.fasta"
    shell:
        """
        R CMD BATCH --no-save --no-restore \
          '--args {input.txps_fa} {input.granges} {params.libsize} {params.pair}' make_reads.R
        mv reads/{params.pair}/sample_01_1.fasta {output.r11}
        mv reads/{params.pair}/sample_01_2.fasta {output.r12}
        mv reads/{params.pair}/sample_02_1.fasta {output.r21}
        mv reads/{params.pair}/sample_02_2.fasta {output.r22}
        """

rule compress_reads:
    input: "{sample}.fasta"
    output: "{sample}.fasta.gz"
    params:
        threads = "12"
    shell: "pigz -p {params.threads} {input}"

rule salmon_index:
    input: "transcripts.fa"
    output: directory("anno/salmon-index-1.5.2")
    params:
        threads = "12"
    shell: "{SALMON} index --keepDuplicates -p {params.threads} -t {input} -i {output}"

rule salmon_quant:
    input:
        r1 = "reads/{sample}_1.fasta.gz",
        r2 = "reads/{sample}_2.fasta.gz",
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
    input: 
        expand("quants/sample_{pair}_{sample}/quant.sf", 
               pair=config["pairs"], sample=config["samples"])
    params:
        nsamp = len(config["pairs"]) * len(config["samples"])
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
        r1 = "reads/{sample}_1.fasta.gz",
        r2 = "reads/{sample}_2.fasta.gz"
    output: "ht2_align/{sample}.bam"
    params:
        xdir = "anno/bdgp6_sim/genome",
        temp_sam = "ht2_align/{sample}.unfilt.sam",
        threads = "12",
        mem = "1G" # per thread
    shell:
        """
        hisat2 -p {params.threads} -f -x {params.xdir} -1 {input.r1} -2 {input.r2} > {params.temp_sam}
        samtools view -b -q 60 {params.temp_sam} | samtools sort -@ {params.threads} -m {params.mem} -o {output}
        samtools index {output}
        rm -f {params.temp_sam}
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

rule find_intersecting_snps:
    input: 
        bam = "ht2_align/{sample}.bam",
        index = "data/drosophila_snp_index.h5",
        tab = "data/drosophila_snp_tab.h5",
        hap = "data/drosophila_haps.h5"
    output: 
        keep = "wasp_mapping/{sample}.keep.bam",
        remap_bam = "wasp_mapping/{sample}.to.remap.bam",
        remap_fq1 = "wasp_mapping/{sample}.remap.fq1.gz",
        remap_fq2 = "wasp_mapping/{sample}.remap.fq2.gz"
    shell:
        """
        {MAPPING}/find_intersecting_snps.py \
          --is_paired_end \
          --is_sorted \
          --output_dir wasp_mapping \
          --snp_tab {input.tab} \
          --snp_index {input.index} \
          --haplotype {input.hap} \
          --samples data/samples \
          {input.bam}
        """

rule hisat_align_flip:
    input:
        remap_fq1 = "wasp_mapping/{sample}.remap.fq1.gz",
        remap_fq2 = "wasp_mapping/{sample}.remap.fq2.gz"
    output: "wasp_mapping/{sample}.map2.bam"
    params:
        xdir = "anno/bdgp6_sim/genome",
        temp_sam = "wasp_mapping/{sample}.unfilt.sam",
        threads = "12",
        mem = "1G"
    shell:
        """
        hisat2 -p {params.threads} -x {params.xdir} \
        -1 {input.remap_fq1} -2 {input.remap_fq2} > {params.temp_sam}
        samtools view -b -q 60 {params.temp_sam} | samtools sort -@ {params.threads} -m {params.mem} -o {output}
        samtools index {output}
        rm -f {params.temp_sam}
        """

rule filter_remapped:
    input:
        remap_bam = "wasp_mapping/{sample}.to.remap.bam",
        map2_bam = "wasp_mapping/{sample}.map2.bam"
    output: "wasp_mapping/{sample}.keep2.bam"
    shell:
        "{MAPPING}/filter_remapped_reads.py {input.remap_bam} {input.map2_bam} {output}"

rule wasp_merge:
    input:
        keep = "wasp_mapping/{sample}.keep.bam",        
        keep2 = "wasp_mapping/{sample}.keep2.bam"
    output: "wasp_mapping/{sample}.merge.bam"
    params:
        presort = "wasp_mapping/{sample}.presort.bam",
        threads = "12",
        mem = "1G"
    shell:
        """
        samtools merge -@ {params.threads} -o {params.presort} {input.keep} {input.keep2}
        samtools sort -@ {params.threads} -m {params.mem} -o {output} {params.presort}
        samtools index {output}
        rm -f {params.presort}
        """

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
        bam = "wasp_mapping/{sample}.merge.bam"
    output:
        ref = "wasp_cht/ref_as_counts.{sample}.h5",
        alt = "wasp_cht/alt_as_counts.{sample}.h5",
        other = "wasp_cht/other_as_counts.{sample}.h5",
        count = "wasp_cht/read_counts.{sample}.h5"
    shell:
        """
        {CHT}/bam2h5.py --chrom data/drosophila_chromInfo.txt \
          --snp_index {input.index} --snp_tab {input.tab} \
          --ref_as_counts {output.ref} --alt_as_counts {output.alt} \
          --other_as_counts {output.other} --read_counts {output.count} \
          {input.bam} 2>&1 | grep -v "partially"
        """

rule wasp_extract:
    input:
        genoprob = "data/drosophila_genoprob.h5",
        hap = "data/drosophila_haps.h5",
        index = "data/drosophila_snp_index.h5",
        tab = "data/drosophila_snp_tab.h5",
        ref = "wasp_cht/ref_as_counts.{sample}.h5",
        alt = "wasp_cht/alt_as_counts.{sample}.h5",
        other = "wasp_cht/other_as_counts.{sample}.h5",
        count = "wasp_cht/read_counts.{sample}.h5",
        tt = "data/drosophila_test_target.txt"
    output: "wasp_cht/hap_read_counts.{sample}.txt"
    shell:
        """
        {CHT}/extract_haplotype_read_counts.py  --chrom data/drosophila_chromInfo.txt \
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
        hap_counts = expand("wasp_cht/hap_read_counts.sample_{pair}_{sample}.txt", 
                            pair=config["pairs"], sample=config["samples"])
    output:
        expand("wasp_cht/hap_read_counts.sample_{pair}_{sample}.adj", 
               pair=config["pairs"], sample=config["samples"])
    shell:
        """
        ls -1 {input.hap_counts} > wasp_cht/preadj
        sed 's/.txt/.adj/' wasp_cht/preadj > wasp_cht/postadj
        {CHT}/update_total_depth.py --seq {input.seq} wasp_cht/preadj wasp_cht/postadj
        """

rule wasp_adjust_het_prob:
    input:
        ref = "wasp_cht/ref_as_counts.{sample}.h5",
        alt = "wasp_cht/alt_as_counts.{sample}.h5",
        adj = "wasp_cht/hap_read_counts.{sample}.adj"
    output: "wasp_cht/hap_read_counts.{sample}.hetp"
    shell:
        "{CHT}/update_het_probs.py --ref_as_counts {input.ref} --alt_as_counts {input.alt} "
        "{input.adj} {output}"

rule CHT:
    input: 
        expand("wasp_cht/hap_read_counts.sample_{pair}_{sample}.hetp", 
               pair=config["pairs"], sample=config["samples"])
    output: "wasp_cht/cht_results.txt"
    shell:
        """
        ls -1 {input} > wasp_cht/cht_input_files.txt
        {CHT}/fit_as_coefficients.py wasp_cht/cht_input_files.txt wasp_cht/cht_as_coef.txt
        {CHT}/fit_bnb_coefficients.py --min_counts 50 --min_as_counts 10 \
        wasp_cht/cht_input_files.txt wasp_cht/cht_bnb_coef.txt
        {CHT}/combined_test.py --min_as_counts 10 \
        --bnb_disp wasp_cht/cht_bnb_coef.txt --as_disp wasp_cht/cht_as_coef.txt \
        wasp_cht/cht_input_files.txt {output}
        """

# rule shuffle:
#     input:
#         r1 = "reads/{sample}_1.fasta",
# 	r2 = "reads/{sample}_2.fasta"
#     output:
#         r1 = "reads/{sample}_1.shuffled.fa",
# 	r2 = "reads/{sample}_2.shuffled.fa"
#     shell:
#         "./shuffle.sh -l {input.r1} -r {input.r2}"
