#!/usr/bin/env bash
source ~/conda/etc/profile.d/conda.sh
conda activate smk9
# get the output dir from the config.yaml, and create the folder
out="$(grep output_dir config.yaml | tail -n 1 | awk '{print $2}' | sed 's/"//g')"
mkdir -p $out
# for reproducibility, copy used config and samplesheet with timestamp into output dir
start_time="`date +"%Y_%m_%d_%I_%M_%p"`"



# print before execution what is on and what is off
echo "Pre-run: options enabled: "
echo "(config.yaml option set to True)"
echo "==================================="
grep True config.yaml | awk '{print $2":",$1}' | sed 's/use_//g'
echo ""
echo ""
echo "Pre-run: options disabled: "
echo "(config.yaml option set to False)"
echo "==================================="
grep False config.yaml | awk '{print $2":",$1}' | sed 's/use_//g'
echo ""
echo ""



# if the pipeline is not executed for the first time with that outputfolder, then move previous report/config/samplesheet/rulegraph into a new folder: outputfolder/logs/previous_executions/.

if ls $out/config*.yaml 1> /dev/null 2>&1; then
    echo "Found files from previous execution, moving them to $out/logs/previous_executions"
    mkdir -p $out/logs/previous_executions
    mv -f $out/pb_variants_* $out/logs/previous_executions/. 2>/dev/null
    mv -f $out/config*.yaml $out/logs/previous_executions/. 2>/dev/null
    mv -f $out/samplesheet*.csv $out/logs/previous_executions/. 2>/dev/null
    mv -f $out/*.sha256 $out/logs/previous_executions/. 2>/dev/null
    mv -f $out/*.filelist $out/logs/previous_executions/. 2>/dev/null

    echo "Files from old execution moved."
else
    echo "No files from a previous execution found. starting..."
fi 

# for reproducibility, copy used config and samplesheet with timestamp into output dir
cp config.yaml $out/config_$start_time.yaml
cp samplesheet.csv $out/samplesheet_$start_time.csv

# create a rulegraph before executing the actual pipeline
snakemake -s rules/snakefile.smk --forceall --rulegraph | dot -Tpdf > $out/pb_variants_rulegraph.$start_time.pdf
nice snakemake -c 120 --jobs 12 -s rules/snakefile.smk --use-conda --rerun-incomplete
# after the run a report is created with task runtime and more info
snakemake -s rules/snakefile.smk --report $out/pb_variants_report.$start_time.html



# create listfile and checksum file of complete output, can be switched off.

if [[ "$1" == "--no-checksums" ]]; then
    echo "Skipping checksum file creation."
else
    echo "Creating checksum and filelist. You can skip this step by providing the --no-checksums parameter: bash runPipeline --no-checksums ."
    echo "Creating Filelist of $out"
    find $out -type f -exec ls  -alth --time-style=long-iso {} \; | sort > $out/filelist_project_$start_time.filelist
    echo "Last Task: creating checksums:"
    echo "Creating checksumfile $out/checksums_$start_time.sha256 ..."
    find $out -type f -exec sha256sum {} \; | sort > $out/checksums_$start_time.sha256
    echo "Done."

fi




echo "Completed the run."
