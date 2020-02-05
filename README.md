# cutAnd
Snakemake pipeline for cut and tag QC and peak calling with mark specific sensitivity.

This pipeline will do the inital QC steps for cut&Tag or cut&Run data. 

# SETUP 

Clone this repository into your project directory:

```
git@github.com:maxsonBraunLab/cutAnd.git

# create a directory to hold your data
cd cutAnd

# create directory for your fastq files
mkdir -p data/raw

# link your fastqs to here
ln -s /path/to/fastq/files/* data/raw

```

The configuration happens in config.yml 

 - Specify the bowtie index for the genome you are aligning to. 

 - A list of sample names, the names must already be part of the filename


# Execution

Setup snakemake profile to run on compute cluster:

SLURM: follow [these instructions](https://github.com/Snakemake-Profiles/slurm)

Example of how to run pipeline

```
snakemake --use-conda --profile slurm -j 60 --latency-wait 60
```

Steps:
    Fastqc
    Bowtie2
    Sort
    Markdups
    Index
    MultiQC

To Compare duplication rates, alignment statistics, and fastqc results across samples, open then file `data/multiqc/multiqc_report.html` in your web browser.
