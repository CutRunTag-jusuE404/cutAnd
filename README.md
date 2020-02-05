# cutAnd
Snakemake pipeline for cut and tag QC and peak calling with mark specific sensitivity

This pipeline will do the inital QC steps for cut&Tag or cut&Run data. 

The configuration all happens in config.yml 

Steps:
    Fastqc
    Trimmomatic
    Bowtie2
    Sambamba dedup
    Macs2 callpeak
