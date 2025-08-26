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
        lima {input.full_bam} {params.barcodes} {output} --num-threads {resources.threads} >> {log} 2>&1
        pbindex {output} --num-threads {resources.threads} >> {log} 2>&1
        samtools index {output} -@ {resources.threads} >> {log} 2>&1
        """
#
rule map:
    input:
        bam = "{sample}.bam" if not config["demultiplex"] else "{output_dir}/bams/{sample}_demux.bam",
        index ="{output_dir}/reference.mmi",
    output:
        bam ="{output_dir}/bams/{sample}_aligned.bam"
    conda:
        "../envs/pbmm2.yaml"
    log:
        "{output_dir}/logs/preprocess_{sample}.log"

    resources:
        threads=lambda wildcards, attempt: attempt * 24,
        time_hrs=lambda wildcards, attempt: attempt * 1,
        mem_gb=lambda wildcards, attempt: 48 + (attempt * 12)

    message:
        "Aligning reads for {input.bam} ..."
    shell:
        """
        pbmm2 align {input.index} {input.bam} --num-threads {resources.threads} --sort {output} >> {log} 2>&1
        pbindex {output} --num-threads {resources.threads} >> {log} 2>&1
        samtools index {output} >> {log} 2>&1
        """


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
        mosdepth --threads {resources.threads} {params.prefix} {input.bam} >> {log} 2>&1
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
        mem_gb=lambda wildcards, attempt: 170 + (attempt * 12)

    message:
        "reporting sample species with kraken2  {input.bam}..."
    shell:
        """
        samtools view -s 0.01 -b {input.bam} -@ {resources.threads} >{params.subsetted_bam} 2>{log}
        samtools fastq {params.subsetted_bam} -@ {resources.threads} >{params.subsetted_fastq} 2>{log}
        kraken2 --use-names --db {params.kraken_db_folder} --threads {resources.threads} --confidence 0.05 --report {output.kraken2_report} {params.subsetted_fastq} >{output.kraken2_outfile} 2>{log}
        """