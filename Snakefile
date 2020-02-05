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

sampdict = get_samples("data/raw")
readpair=[[k+"_R1", k+"_R2"] for k in sampdict.keys()]
reads=[v for s in readpair for v in s]

# define target output
rule all:
    input:
         "data/multiqc/multiqc_report.html",
         expand("data/fastqc/{read}.html", read=reads),
         expand([
                "data/aligned/{sample}.bam",
                "data/mrkdup/{sample}.sorted.mrkdup.bam",
                "data/mrkdup/{sample}.sorted.mrkdup.bam.bai"],
                 sample=sampdict.keys())

rule fastqc:
    input:
        "data/raw/{read}.fastq.gz"
    output:
        html="data/fastqc/{read}.html",
	zip="data/fastqc/{read}_fastqc.zip" 
    log:
        "data/logs/fastqc_{read}.log"
    threads: 4
    wrapper:
        "0.49.0/bio/fastqc"
        

# align samples to genome
rule bowtie2:
    input:
        lambda wildcards: sampdict[wildcards.sample]
    output:
        "data/aligned/{sample}.bam"
    log:
        out="data/logs/bowtie2_{sample}.log",
        err="data/logs/bowtie2_{sample}.err"
    conda:
        "envs/align.yml"
    threads: 8
    shell:
        "bowtie2 --local --very-sensitive-local "
        "--no-unal --no-mixed --threads {threads} "
        "--no-discordant --phred33 "
        "-I 10 -X 700 -x {config[GENOME_IDX]} "
        "-1 {input[0]} -2 {input[1]} 2>{log.err} | tee {log.out} | samtools view -Sbh - > {output}"

# sort
rule sort:
    input:
        rules.bowtie2.output 
    output:
        temp("data/sort/{sample}.sorted.bam")
    conda:
        "envs/sam.yml"
    log: 
        "data/logs/sort_{sample}.log"
    threads: 4
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
    threads: 4
    shell:
        "sambamba markdup -t {threads} {input} {output} > {log} 2>&1"

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
    threads: 4
    shell:
        "sambamba index -t 2 {input} > {log} 2>&1"

rule multiqc:
    input:
        expand("data/mrkdup/{sample}.sorted.mrkdup.bam.bai", sample=sampdict.keys())
    output:
        "data/multiqc/multiqc_report.html"
    conda:
        "envs/multiqc.yml"
    log:
        "data/logs/multiqc.log"
    shell:
        "multiqc --force -o data/multiqc data/ > {log} 2>&1"
    

