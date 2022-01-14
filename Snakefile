configfile: "config.json"

SALMON = "/proj/milovelab/bin/salmon-1.5.2_linux_x86_64/bin/salmon"
CHT = "python3.5 /nas/longleaf/apps/wasp/2019-12/WASP/CHT"
MAPPING = "python3.5 /nas/longleaf/apps/wasp/2019-12/WASP/mapping"
MMSEQ = "/proj/milovelab/bin/mmseq-1.0.10a/bin"
TERMINUS = "/proj/milovelab/bin/terminus/target/release/terminus"

rule all:
    input: 
        # gr = "granges.rda",
        # reads = expand("reads/sample_{pair}_{sample}_{read}.shuffled.fasta.gz", 
        #                pair=config["pairs"], sample=config["samples"], 
        #                read=config["reads"]),
        # summarized_experiment = "txp_allelic_se.rda"
        # hisat = expand("ht2_align/sample_{pair}_{sample}.bam",
        #                pair=config["pairs"], sample=config["samples"]),
        # wasp_counts = expand("wasp_cht/ref_as_counts.sample_{pair}_{sample}.h5",
        #                      pair=config["pairs"], sample=config["samples"]),
        # wasp_result = "wasp_cht/cht_results.txt"
        # mmseq = "mmseq/mmdiff_results.txt"
        # mmseq = "mmseq/mmdiff_gene_results.txt"
        # mmseq = expand("mmseq_nodup/sample_{pair}_{sample}_{allele}.collapsed.mmseq",
        #               pair=config["pairs"], sample=config["samples"], allele=config["alleles"])
        terminus = expand("terminus/sample_{pair}_{sample}/groups.txt",
                          pair=config["pairs"], sample=config["samples"])

rule make_expression:
    output:
        txps_fa = "transcripts.fa",
	granges = "granges.rda",
        chr = "data/drosophila_chr_2L.fasta",
        vcf = "data/drosophila_chr_2L.vcf",
        tt = "data/drosophila_test_target.txt",
        mmseq = "data/transcripts_mmseq.fa",
        t2g = "data/t2g.tsv",
        a2t = "data/a2t.tsv"
    shell:
        "R CMD BATCH --no-save --no-restore '--args {output.txps_fa} {output.granges} "
        "{output.chr} {output.vcf} {output.tt} {output.mmseq} {output.t2g}' make_expression.R"

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

rule shuffle_reads:
    input:
        r1 = "{sample}_1.fasta",
	r2 = "{sample}_2.fasta"
    output:
        r1 = "{sample}_1.shuffled.fasta",
	r2 = "{sample}_2.shuffled.fasta"
    shell:
        """
        ./shuffle.sh -l {input.r1} -r {input.r2}
        rm -f {input.r1}
        rm -f {input.r2}
        """

rule compress_reads:
    input: "{sample}.shuffled.fasta"
    output: "{sample}.shuffled.fasta.gz"
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
        r1 = "reads/{sample}_1.shuffled.fasta.gz",
        r2 = "reads/{sample}_2.shuffled.fasta.gz",
        index = "anno/salmon-index-1.5.2"
    output:
        "quants/{sample}/quant.sf"
    params:
        dir = "quants/{sample}",
        threads = "12",
	nboot = "30"
    shell:
        "{SALMON} quant -i {input.index} -l IU -p {params.threads} "
	"--numBootstraps {params.nboot} -d "
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

rule terminus_group:
    input: "quants/{sample}/quant.sf"
    output: "terminus/{sample}/groups.txt"
    params:
        t2g = "data/t2g.tsv",
        a2t = "data/a2t.tsv",
        qdir = "quants/{sample}",
        tdir = "terminus"
    shell: "{TERMINUS} group --t2g --a2t --dir {params.qdir} --out {params.tdir}"

rule hisat_align:
    input:
        index = "anno/hisat2_bdgp6/genome.1.ht2",
        r1 = "reads/{sample}_1.shuffled.fasta.gz",
        r2 = "reads/{sample}_2.shuffled.fasta.gz"
    output: "ht2_align/{sample}.bam"
    params:
        xdir = "anno/hisat2_bdgp6/genome",
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

rule wasp_find_intersecting_snps:
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

rule wasp_hisat_align_flip:
    input:
        remap_fq1 = "wasp_mapping/{sample}.remap.fq1.gz",
        remap_fq2 = "wasp_mapping/{sample}.remap.fq2.gz"
    output: "wasp_mapping/{sample}.map2.bam"
    params:
        xdir = "anno/hisat2_bdgp6/genome",
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

rule wasp_bam_the_fake_bams:
    input:
        keep = "wasp_mapping/{sample}.keep.bam",
        to_remap = "wasp_mapping/{sample}.to.remap.bam",
    output:
        keep = "wasp_mapping/{sample}.keep.actual.bam",
        to_remap = "wasp_mapping/{sample}.to.remap.actual.bam",
    shell:
        """
        samtools view -b {input.keep} > {output.keep}
        samtools view -b {input.to_remap} > {output.to_remap}
        rm -f {input.keep}
        rm -f {input.to_remap}
        """

rule wasp_filter_remapped:
    input:
        remap_bam = "wasp_mapping/{sample}.to.remap.actual.bam",
        map2_bam = "wasp_mapping/{sample}.map2.bam"
    output: "wasp_mapping/{sample}.keep2.bam"
    shell:
        "{MAPPING}/filter_remapped_reads.py {input.remap_bam} {input.map2_bam} {output}"

rule wasp_merge:
    input:
        keep = "wasp_mapping/{sample}.keep.actual.bam",
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

rule wasp_CHT:
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

rule mmseq_index:
    input: "data/transcripts_mmseq.fa"
    output: "anno/bt_index/index.1.ebwt"
    params:
        dir = "anno/bt_index/index",
        threads = "12"
    shell:
        "bowtie-build --offrate 3 --threads {params.threads} {input} {params.dir}"

rule mmseq_index_nodup:
    input: "data/transcripts_mmseq_nodup.fa"
    output: "anno/bt_index_nodup/index.1.ebwt"
    params:
        dir = "anno/bt_index_nodup/index",
        threads = "12"
    shell:
        "bowtie-build --offrate 3 --threads {params.threads} {input} {params.dir}"

rule mmseq_align:
    input:
        r1 = "reads/{sample}_1.shuffled.fasta.gz",
        r2 = "reads/{sample}_2.shuffled.fasta.gz",
        index = "anno/bt_index/index.1.ebwt"
    output: "mmseq/{sample}.bam"
    params:
        dir = "anno/bt_index/index",
        threads = "12",
        mem = "1G",
        m = "100",
        X = "500",
        chunkmbs = "256"
    shell:
        """
        bowtie -a --best --strata -S -m {params.m} -X {params.X} --chunkmbs {params.chunkmbs} \
        -p {params.threads} -f -x {params.dir} \
        -1 <(gzip -dc {input.r1}) -2 <(gzip -dc {input.r2}) | \
        samtools view -F 0xC -bS - | samtools sort -@ {params.threads} -m {params.mem} -n - -o {output}
        """

rule mmseq_align_nodup:
    input:
        r1 = "reads/{sample}_1.shuffled.fasta.gz",
        r2 = "reads/{sample}_2.shuffled.fasta.gz",
        index = "anno/bt_index_nodup/index.1.ebwt"
    output: "mmseq_nodup/{sample}.bam"
    params:
        dir = "anno/bt_index_nodup/index",
        threads = "12",
        mem = "1G",
        m = "100",
        X = "500",
        chunkmbs = "256"
    shell:
        """
        bowtie -a --best --strata -S -m {params.m} -X {params.X} --chunkmbs {params.chunkmbs} \
        -p {params.threads} -f -x {params.dir} \
        -1 <(gzip -dc {input.r1}) -2 <(gzip -dc {input.r2}) | \
        samtools view -F 0xC -bS - | samtools sort -@ {params.threads} -m {params.mem} -n - -o {output}
        """

rule mmseq_bam2hits:
    input: 
        ref = "data/transcripts_mmseq.fa",
        bam = "mmseq/{sample}.bam"
    output: "mmseq/{sample}.hits"
    params:
        threads = "12"
    shell: "OMP_NUM_THREADS={params.threads} {MMSEQ}/bam2hits-linux {input.ref} {input.bam} > {output}"

rule mmseq_bam2hits_nodup:
    input: 
        ref = "data/transcripts_mmseq_nodup.fa",
        bam = "mmseq_nodup/{sample}.bam"
    output: "mmseq_nodup/{sample}.hits"
    params:
        threads = "12"
    shell: "OMP_NUM_THREADS={params.threads} {MMSEQ}/bam2hits-linux {input.ref} {input.bam} > {output}"

rule mmseq_quant:
    input: "{sample}.hits"
    output: "{sample}.mmseq"
    params:
        name = "{sample}",
        threads = "12"
    shell: "OMP_NUM_THREADS={params.threads} {MMSEQ}/mmseq-linux {input} {params.name}"

rule mmseq_split:
    input: 
        txp = "mmseq/{sample}.mmseq",
        gene = "mmseq/{sample}.gene.mmseq"
    output: 
        m = "mmseq/{sample}_M.mmseq",
        p = "mmseq/{sample}_P.mmseq",
        m_gene = "mmseq/{sample}_M.gene.mmseq",
        p_gene = "mmseq/{sample}_P.gene.mmseq"
    shell:
        """
        R CMD BATCH --no-save --no-restore \
          '--args {input.txp} {output.m} {output.p}' mmseq_scripts/mmseq_split_sample.R
        R CMD BATCH --no-save --no-restore \
          '--args {input.gene} {output.m_gene} {output.p_gene}' mmseq_scripts/mmseq_split_sample.R
        """

rule mmseq_split_nodup:
    input: "mmseq_nodup/{sample}.mmseq"
    output: 
        m = "mmseq_nodup/{sample}_M.mmseq",
        p = "mmseq_nodup/{sample}_P.mmseq"
    params:
        base = "mmseq_nodup/{sample}"
    shell:
        """
        R CMD BATCH --no-save --no-restore \
          '--args {input} {output.m} {output.p}' mmseq_scripts/mmseq_split_sample.R
        cp {params.base}.identical.mmseq {params.base}_M.identical.mmseq
        cp {params.base}.identical.mmseq {params.base}_P.identical.mmseq
        cp {params.base}.identical.trace_gibbs.gz {params.base}_M.identical.trace_gibbs.gz
        cp {params.base}.identical.trace_gibbs.gz {params.base}_P.identical.trace_gibbs.gz
        """

rule mmseq_split_trace:
    input: "mmseq_nodup/{sample}.trace_gibbs.gz"
    output: 
        m = "mmseq_nodup/{sample}_M.trace_gibbs.gz",
        p = "mmseq_nodup/{sample}_P.trace_gibbs.gz"
    shell:
        """
        R CMD BATCH --no-save --no-restore \
          '--args {input} {output.m} {output.p}' mmseq_scripts/mmseq_split_trace.R
        """

rule mmseq_split_bigm_and_k:
    input: 
        bigm = "mmseq_nodup/{sample}.M",
        k = "mmseq_nodup/{sample}.k"
    output: 
        bigm_m = "mmseq_nodup/{sample}_M.M",
        bigm_p = "mmseq_nodup/{sample}_P.M",
        k_m = "mmseq_nodup/{sample}_M.k",
        k_p = "mmseq_nodup/{sample}_P.k"
    shell:
        """
        R CMD BATCH --no-save --no-restore \
          '--args {input.bigm} {output.bigm_m} {output.bigm_p}' mmseq_scripts/mmseq_split_bigm.R
        cp {input.k} {output.k_m}
        cp {input.k} {output.k_p}
        """

# rule mmcollapse:
#     input: 
#         mmseq = expand("mmseq_nodup/sample_{pair}_{sample}_{allele}.mmseq", pair=config["pairs"], sample=config["samples"], allele=config["alleles"]),
#         trace = expand("mmseq_nodup/sample_{pair}_{sample}_{allele}.trace_gibbs.gz", pair=config["pairs"], sample=config["samples"], allele=config["alleles"]),
#         bigm = expand("mmseq_nodup/sample_{pair}_{sample}_{allele}.M", pair=config["pairs"], sample=config["samples"], allele=config["alleles"])
#     output: expand("mmseq_nodup/sample_{pair}_{sample}_{allele}.collapsed.mmseq", pair=config["pairs"], sample=config["samples"], allele=config["alleles"])
#     params:
#         base = expand("mmseq_nodup/sample_{pair}_{sample}_{allele}", pair=config["pairs"], sample=config["samples"], allele=config["alleles"]),
#         threads = "12" # requires 1.5 Gb per thread
#     shell: "OMP_NUM_THREADS={params.threads} {MMSEQ}/mmcollapse-linux {params.base}"

rule mmdiff:
    input: 
        m = expand("mmseq/sample_{pair}_{sample}_M.mmseq",
                   pair=config["pairs"], sample=config["samples"]),
        p = expand("mmseq/sample_{pair}_{sample}_P.mmseq",
                   pair=config["pairs"], sample=config["samples"])
    output: "mmseq/mmdiff_results.txt"
    params:
        n = 2 * len(config["pairs"])
    shell: "{MMSEQ}/mmdiff-linux -de {params.n} {params.n} {input.m} {input.p} > {output}"

rule mmdiff_gene:
    input: 
        m = expand("mmseq/sample_{pair}_{sample}_M.gene.mmseq",
                   pair=config["pairs"], sample=config["samples"]),
        p = expand("mmseq/sample_{pair}_{sample}_P.gene.mmseq",
                   pair=config["pairs"], sample=config["samples"])
    output: "mmseq/mmdiff_gene_results.txt"
    params:
        n = 2 * len(config["pairs"])
    shell: "{MMSEQ}/mmdiff-linux -de {params.n} {params.n} {input.m} {input.p} > {output}"
