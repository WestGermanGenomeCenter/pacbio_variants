#!/usr/bin/env bash
source ~/conda/etc/profile.d/conda.sh
conda activate smk9
# get the output dir from the config.yaml, and create the folder
out="$(grep output_dir config.yaml | tail -n 1 | awk '{print $2}' | sed 's/"//g')"
mkdir -p $out
# for reproducibility, copy used config and samplesheet with timestamp into output dir
start_time="`date +"%Y_%m_%d_%I_%M_%p"`"
cp config.yaml $out/config_$start_time.yaml
cp samplesheet.csv $out/samplesheet_$start_time.csv

# create a rulegraph before executing the actual pipeline
snakemake -s rules/snakefile.smk --forceall --rulegraph | dot -Tpdf > $out/pb_variants_rulegraph.$start_time.pdf
nice snakemake -c 96 --jobs 4 -s rules/snakefile.smk --use-conda --rerun-incomplete
# after the run a report is created with task runtime and more info
snakemake -s rules/snakefile.smk --report $out/pb_variants_report.$start_time.html

