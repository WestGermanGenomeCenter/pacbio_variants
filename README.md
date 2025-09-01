# pb_variants 
A snakemake 9 based Pipeline for hifi snp, sv, cnv calling, phasing and more

Only PacBio data for now

__!!THIS PIPLINE IS IN-DEVELOPMENT AND EXPERIMENTAL, USE AT YOUR OWN RISK!!__

## what this tools aims to deliver:
    - newest and best tools suited for HiFi data (only for now)
    - singletons and trio analysis (trio is coming sometime...)
    - human-first (hg38 for now), others should be possible (untested...)

## included tools:
- deepvariant or bcftools for snp calling
- snps get used for phasing with whatshap and longphase
- paraphase 
- trgt
- hificnv 
- mitorsaw (just hg38)
- sniffles for sv calls that get phased with longphase
- sawfish for svs and cnvs (results phased by sawfish)
- mosdepth, multiqc
- pbmm2 for mapping
- kraken2 for contamination detection (downsamples massively, needs kraken2 database)
- demultiplexing of input as option, will not split the files per barcode.
- for now one .bam per sample
- NanoCaller for phased snp/indel calls

## how to run
- make sure you have a conda environment active with snakemake9+ (called smk9 in the runPipeline_local.sh)
    - this can also be achieved by running the included setupPipeline_hpc.sh
        - that script uses conda to create the env smk9 - with snakemake 9 installed already (check file smk9_minimal.txt)
- cp/mv/ln your unmapped .bam file into the root folder of this directory (pb_variants/.)
- edit samplesheet.csv with your filename 
    - one sample per line, do not delete the header line
- edit config.yaml to your liking/ folder structure
- make sure you are in an interactive terminal session inside a screen / tmux or similar
- bash runPipeline_local.sh for local installment on single-server setups, 
- bash runPipeline.sh on HPC 
- non-hpc users need to edit the command line deepvariant in rules/snakefile.smk: the rule is currently only working on the HPC, not on other hardware. Please adjust.


# DAG
This DAG was made:
- with snakemake --rulegraph option
- demultiplexing disabled through config.yaml setting
- deepvariant enabled through config.yaml setting
- kraken enabled through config.yaml setting

![alt text](<Bildschirmfoto vom 2025-08-28 16-36-02.png>)




## output files
- the first step of the pipeline is to strip the kinetics data out of the .bam input file, but keep the methylation data. This makes all following processes much faster without any real data loss. 
- for each input sample:
    - mosdepth and kraken (optional) report that get summarized with multiqc
    - mapped .bam file haplotaged with whatshap and longphase
    - .vcf(.gz) file for:
        - snps / indels from deepvariant or bcftools, phased with whatshap and longphase
        - trgt
        - paraphase
        - mitorsaw
        - hificnv
        - sawfish sv / cnv 
        - svs from sniffles phased with longphase
        - snp / svs / indels from nanocaller
- snakemake report, rulegraph, copy of samplesheet and config.yaml with timestamp


## roadmap:
- trio calling : deeptrio, glnexus
- methylation tools (only if .bam is input format)
- testing cuteSV, svim Clair3, delly 
- snp /sv annotation: annotsv,


## why this work is being done:
- nf-core/pacvar: https://nf-co.re/pacvar/1.0.1/
    - does not run without sudo for us
    - seems not mature enough (imho)
    - not newest tools included
    - not all wanted tools included

- pacbios wdl-based workflow: https://github.com/PacificBiosciences/HiFi-somatic-WDL
    - doesnt run on our hardware

- other, locally developed snakemake-based workflows: (eg: https://github.com/core-unit-bioinformatics/workflow-smk-longread-variant-calling)
    - not all wanted tools included

- Radboud's Valentine workflow:
    - not available to us
    - not all wanted tools included

- smrtlinks internal pipeline:
    - singularity not working, limited tool options
    - not all wanted tools included
