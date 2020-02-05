# cutAnd snakemake pipeline for Cut&Run or Cut&Tag data
import glob
import os

# import configuration
configfile: "src/config.yml"

# map samples to fastqs
def get_samples(dir):
    '''Matches samples to their fastq files.'''
    samples = config["SAMPLES"]
    sdict = {s: glob.glob(f"{dir}/*{s}*") for s in samples}
    return sdict

sampdict = get_samples("data/raw/test")
print(sampdict)

# define target output
rule all:
    input:
        expand(["data/aligned/{sample}.bam",
                "data/sort/{sample}.sorted.bam",
                "data/mrkdup/{sample}.sorted.mrkdup.bam",
                "data/mrkdup/{sample}.sorted.mrkdup.bam.bai"],
                 sample=sampdict.keys())

# perform trimming and fastqc
rule trim_galore:
    input:
        lambda wildcards: sampdict[wildcards.sample]
    output:
        "data/trim_galore/{sample}_val_1.fastq.gz", "data/trim_galore/{sample}_val_2.fastq.gz"
    params:
        extra="--illumina -q 20",
        outdir="data/trim_galore"
    conda:
        "envs/trim.yml"
    log:
        "data/logs/trim_{sample}.log"
    shell:
        "trim_galore --suppress_warn --no_report_file --basename "
        "{wildcards.sample} --fastqc --paired -o {params.outdir} " 
        "{input} > {log} 2>&1"

# align samples to genome
rule bowtie2:
    input:
        "data/trim_galore/{sample}_val_1.fastq.gz", "data/trim_galore/{sample}_val_2.fastq.gz"
    output:
        "data/aligned/{sample}.bam"
    log:
        "data/logs/bowtie2_{sample}.log"
    conda:
        "envs/align.yml"
    threads: 8
    shell:
        "bowtie2 --local --very-sensitive-local "
        "--no-unal --no-mixed --threads {threads} "
        "--no-discordant --phred33 "
        "-I 10 -X 700 -x {config[GENOME_IDX]} "
        "-1 {input[0]} -2 {input[1]} 2>&1 | tee {log} | samtools view -Sbh > {output}"

# sort
rule sort:
    input:
        rules.bowtie2.output 
    output:
        "data/sort/{sample}.sorted.bam"
    conda:
        "envs/sam.yml"
    log: 
        "data/logs/sort_{sample}.log"
    threads: 2
    shell:
        "sambamba sort {input} -t {threads} -o {output} > {log} 2>&1"

# mark pcr duplicates
rule mrkdup:
    input:
        rules.sort.output
    output:
        "data/mrkdup/{sample}.sorted.mrkdup.bam"
    conda:
        "envs/sam.yml"
    log:
        "data/logs/mrkdup_{sample}.log"
    threads: 2
    shell:
        "sambamba mrkdup -t {threads} {input} {output} > {log} 2>&1"

# generate index file
rule index:
    input:
       rules.mrkdup.output 
    output:
        "data/mrkdup/{sample}.sorted.mrkdup.bam.bai"
    conda:
        "envs/sam.yml"
    log: 
        "data/logs/index_{sample}.log"
    threads: 2
    shell:
        "sambamba index -t 2 {input} > {log} 2>&1"

