import pandas as pd
configfile: "config.yaml"
samples_df = pd.read_csv("samplesheet.csv")

# Define the samples
samples = samples_df['sample'].tolist()
filenames_without_extension = [file.rstrip('.bam') for file in samples]
#print("filenames without extension:")
#print(filenames_without_extension)
def get_folder_name(samplename_with_bam): # only one sample name for each run allowed#
        only_folder_name=samplename_with_bam.split(".bam")[0]
        return (only_folder_name)

output_dir=config["output_dir"]


rule index_reference:
    input:
        reference=config["reference"]
    output:
        index="{output_dir}/reference.mmi"
    conda:
        "../envs/pbmm2.yaml"
    log:
        "{output_dir}/logs/index_reference.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 12,
        time_hrs=lambda wildcards, attempt: attempt * 1,
        mem_gb=lambda wildcards, attempt: 36 + (attempt * 12)
    message:
        "Indexing reference genome..."
    shell:
        """
        pbmm2 index {input.reference} {output} --num-threads {resources.threads} >> {log} 2>&1
        samtools faidx {input.reference} >> {log} 2>&1
        """

rule demultiplex:
    input:
        full_bam="{sample}.bam"
    output:
        demux_fastq="{output_dir}/bams/{sample}_demux.bam"
    conda:
        "../envs/lima.yaml"
    params:
        barcodes=config["barcodes_lima_fasta"]
    log:
        "{output_dir}/logs/demultiplex_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 12,
        time_hrs=lambda wildcards, attempt: attempt * 1,
        mem_gb=lambda wildcards, attempt: 12 + (attempt * 12)

    message:
        "Demultiplexing  {sample}..."
    shell:
        """
        lima {input.full_bam} {params.barcodes} {output} --num-threads {resources.threads} >{log} 2>&1
        pbindex {output} --num-threads {resources.threads} >{log} 2>&1
        samtools index {output} -@ {resources.threads} >{log} 2>&1
        """

rule kinetics_removal:
    input:
        bam = "{sample}.bam" if not config["demultiplex"] else "{output_dir}/bams/{sample}_demux.bam",
    output:
        bam_no_kinetics="{output_dir}/bams/{sample}_no_kinetics.bam",
    conda:
        "../envs/samtools.yaml"
    log:
        "{output_dir}/logs/remove_kinetics_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: 6 + (attempt * 10),
        time_hrs=lambda wildcards, attempt: attempt * 2,
        mem_gb=lambda wildcards, attempt: 2 + (attempt * 24)
    message:
        "removing kinetics of {input.bam} ..."
    shell:
        """
        samtools view -@ {resources.threads} --bam --remove-tag fi,fp,fn,ri,rp,rn --output {params.bam_no_kinetics} {input.bam} >{log} 2>&1
        """




rule map:
    input:
        bam_no_kinetics="{output_dir}/bams/{sample}_no_kinetics.bam",
        index ="{output_dir}/reference.mmi",
    output:
        bam ="{output_dir}/bams/{sample}_aligned.bam",
        bai ="{output_dir}/bams/{sample}_aligned.bam.bai",
    conda:
        "../envs/pbmm2.yaml"
    log:
        "{output_dir}/logs/preprocess_{sample}.log"

    resources:
        threads=lambda wildcards, attempt: 6 + (attempt * 10),
        time_hrs=lambda wildcards, attempt: attempt * 4,
        mem_gb=lambda wildcards, attempt: 72 + (attempt * 24)
    params:
        picard_tmp="{output_dir}/bams/",
        bam_no_kinetics="{output_dir}/bams/{sample}_no_kinetics.bam",
        bai="{output_dir}/bams/{sample}_aligned.bam.bai",
        bamstats="{output_dir}/bams/{sample}_bamstats.txt",
        pref="{output_dir}/bams/{sample}_bamstats",

    message:
        "Aligning reads for {input.bam} ..."
    shell:
        """
        pbmm2 align {input.index} {input.bam_no_kinetics} --num-threads {resources.threads} --sort {output.bam} >{log} 2>&1
        pbindex {output.bam} --num-threads {resources.threads} >{log} 2>&1
        #picard BuildBamIndex -I {output.bam} -TMP_DIR {params.picard_tmp} -O {params.bai} >{log} 2>&1
        samtools index {output.bam} -@ {resources.threads} >{log} 2>&1
        samtools stats {output.bam} > {params.bamstats} 2>{log}
        plot-bamstats {params.bamstats} -p {params.pref} 2>{log}
        """

#samtools stats {input.bams} >{output.samtools_stats_file} 2>{log}
#        samtools coverage {input.bams} >{params.samtools_coverage_file} 2>{log}
#        plot-bamstats {output.samtools_stats_file} -p {params.samtools_plot_prefix}   2>{log}   
#



rule mosdepth:
    input:
        bam="{output_dir}/bams/{sample}_aligned.bam",
    output:
        depth_file="{output_dir}/mosdepth/{sample}.mosdepth.summary.txt"
    params:
        prefix="{output_dir}/mosdepth/{sample}"
    conda:
        "../envs/mosdepth.yaml"
    log:
        "{output_dir}/logs/mosdepth_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 8,
        time_hrs=lambda wildcards, attempt: attempt * 1,
        mem_gb=lambda wildcards, attempt: 4 + (attempt * 12)

    message:
        "reporting mapping depth of  {input.bam}..."
    shell:
        """
        mosdepth --threads {resources.threads} {params.prefix} {input.bam} >{log} 2>&1
        """


# adding fastqc or fastp



rule fastqc:
    input:
        bam = "{sample}.bam" if not config["demultiplex"] else "{output_dir}/bams/{sample}_demux.bam",
    output:
        zip = temp("{output_dir}/fastqc/{sample}_fastqc.zip"), 
        html = "{output_dir}/fastqc/{sample}_fastqc.html",
    params:
        subsetted_fastq=temp("{output_dir}/fastqc/{sample}_subsetted.fastq"),
        fastp_report = "{output_dir}/fastqc/{sample}_fastp.html",
        fastp_json = "{output_dir}/fastqc/{sample}_fastp.json",
        output_dir="{output_dir}/fastqc/",
    conda:
        "../envs/fastqc.yaml"
    log:
        "{output_dir}/logs/fastqc_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 4,
        time_hrs=lambda wildcards, attempt: attempt * 3,
        mem_gb=lambda wildcards, attempt: 4 + (attempt * 8)

    message:
        "reporting read quality input for {input.bam}..."
    shell:
        """
        samtools fastq {input.bam} -@ {resources.threads} >{params.subsetted_fastq} 2>{log}
        fastqc -q --threads {resources.threads} {params.subsetted_fastq} -o {params.output_dir} >> {log} 2>&1
        fastp -i {params.subsetted_fastq} -h {params.fastp_report} -j {params.fastp_json} >> {log} 2>&1
        """



rule kraken2:
    input:
        bam = "{sample}.bam" if not config["demultiplex"] else "{output_dir}/bams/{sample}_demux.bam",
    output:
        kraken2_report="{output_dir}/kraken2/{sample}_kraken2.report",
        kraken2_outfile=temp("{output_dir}/kraken2/{sample}_kraken2.kraken2")
    params:
        subsetted_bam=temp("{output_dir}/kraken2/{sample}_subsetted.bam"),
        subsetted_fastq=temp("{output_dir}/kraken2/{sample}_subsetted.fastq"),
        kraken_db_folder=config["kraken2_db_folder"],
    conda:
        "../envs/kraken2.yaml"
    log:
        "{output_dir}/logs/kraken2_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 8,
        time_hrs=lambda wildcards, attempt: attempt * 1,
        mem_gb=lambda wildcards, attempt: 170 + (attempt * 16)

    message:
        "reporting sample species with kraken2 for {input.bam}..."
    shell:
        """
        samtools view -s 0.01 -b {input.bam} -@ {resources.threads} >{params.subsetted_bam} 2>{log}
        samtools fastq {params.subsetted_bam} -@ {resources.threads} >{params.subsetted_fastq} 2>{log}
        kraken2 --use-names --db {params.kraken_db_folder} --threads {resources.threads} --confidence 0.05 --report {output.kraken2_report} {params.subsetted_fastq} >{output.kraken2_outfile} 2>{log}
        """