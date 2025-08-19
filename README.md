# pb_variants 
A snakemake 9 based Pipeline for hifi snp, sv, cnv calling, phasing and more

Only PacBio data for now, only single samples for now


## why this work is being done:
- nf-core/pacvar: https://nf-co.re/pacvar/1.0.1/
    - does not run without sudo for us
    - seems not mature enough
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


!!THIS PIPLINE IS IN-DEVELOPMENT AND EXPERIMENTAL, USE AT YOUR OWN RISK!!

## what this tools aims to deliver:
    - newest and best tools suited for HiFi data (only for now)
    - singletons and trio analysis (trio is coming sometime...)
    - human-first (hg38 for now), others should be possible (untested...)

## included tools:
- deepvariant or bcftools for snp calling
- snps get used for phasing with whatshap and hiphase 
- paraphase 
- trgt
- hificnv 
- sniffles for sv calls
- sawfish for svs and cnvs (results phased)
- mosdepth, multiqc
- pbmm2 for mapping
- .fastq.gz files or .bam (unmapped) as input file
- demultiplexing of input as option, will not split the files per barcode.
- for now one .bam per sample
- NanoCaller for phased snp/sv calls

## how to run
- make sure you have a conda environment active with snakemake9+ (called smk9 in the runPipeline_local.sh)
- cp/mv/ln your unmapped .bam or .fastq file into the root folder of this directory
- edit samplesheet.csv with your filename minus the ending, so for sample.bam you should enter sample 
    - one sample per line, do not delete the header line
- edit config.yaml to your liking/ folder structure
- make sure you are in an interactive terminal session inside a screen / tmux or similar
- bash runPipeline_local.sh for local installment on single-server setups, 
- bash runPipeline.sh on HPC 


## roadmap:
- trio calling : deeptrio, glnexus
- methylation tools (only if .bam is input format)
- testing cuteSV, svim Clair3, delly 
- kraken2