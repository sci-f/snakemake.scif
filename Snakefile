configfile: "config.yaml"
workdir: "/scif/data"

rule bwa_index:
    input:
        "{ref}.fa"
    output:
        "{ref}.fa.bwt"
    shell:
        "scif run bwa index {input}"

rule bwa_map:
    input:
        index="genome.fa.bwt",
        ref="genome.fa",
        reads="samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.sam"
    shell:
        "scif run bwa mem -o {output} {input.ref} {input.reads}"


rule samtools_sam_to_bam:
    input:
        "mapped_reads/{sample}.sam"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "scif --quiet run samtools 'view -bS $SCIF_DATA/{input} > $SCIF_DATA/{output}'"


rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "scif --quiet run samtools 'sort -T $SCIF_DATA/sorted_reads/{wildcards.sample} -O bam $SCIF_DATA/{input} > $SCIF_DATA/{output}'"


rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "scif run samtools 'index $SCIF_DATA/{input}'"


rule bcftools_call:
    input:
        fa="genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
    output:
        "calls/all.vcf"
    shell:
        "scif --quiet run samtools 'mpileup -g -f $SCIF_DATA/{input.fa} $SCIF_DATA/{input.bam} | bcftools call -mv - > $SCIF_DATA/{output}'"


rule report:
    input:
        "calls/all.vcf"
    output:
        "report.html"
    run:
        from snakemake.utils import report
        with open(input[0]) as vcf:
            n_calls = sum(1 for l in vcf if not l.startswith("#"))

        report("""
        An example variant calling workflow
        ===================================

        Reads were mapped to the Yeast
        reference genome and variants were called jointly with
        SAMtools/BCFtools.

        This resulted in {n_calls} variants (see Table T1_).
        """, output[0], T1=input[0])


rule all:
    input:
        "report.html"
