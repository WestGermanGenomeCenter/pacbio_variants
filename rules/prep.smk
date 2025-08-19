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
        bam = "{sample}.bam" if config["demultiplex"] else "{sample}_demux.bam",
        index ="{output_dir}/reference.mmi",
    output:
        bam ="{output_dir}/bams/{sample}_aligned.bam"
    conda:
        "../envs/pbmm2.yaml"
    log:
        "{output_dir}/logs/preprocess_{sample}.log"

    resources:
        threads=lambda wildcards, attempt: attempt * 96,
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

# using .fastq.gz files is also possible, but loses a lot of data included in the .bam files, especially the kinetics!
rule map_fq_gz:
    input:
        bam = "{sample}.fastq.gz" if config["demultiplex"] else "{sample}_demux.bam",
        index ="{output_dir}/reference.mmi",
    output:
        bam ="{output_dir}/bams/{sample}_aligned.bam"
    conda:
        "../envs/pbmm2.yaml"
    log:
        "{output_dir}/logs/preprocess_{sample}.log"

    resources:
        threads=lambda wildcards, attempt: attempt * 96,
        time_hrs=lambda wildcards, attempt: attempt * 1,
        mem_gb=lambda wildcards, attempt: 48 + (attempt * 12)

    message:
        "Aligning reads for {input.bam} ..."
    shell:
        """
        pbmm2 align {input.index} {input.bam} --num-threads {resources.threads} --sort {output} >> {log} 2>&1
        pbindex {output} --num-threads {resources.threads} >> {log} 2>&1
        samtools index {output} -@ {resources.threads} >> {log} 2>&1
        """

# using .fastq.gz files is also possible, but loses a lot of data included in the .bam files, especially the kinetics!
rule demultiplex_fq_gz:
    input:
        full_bam="{sample}.fastq.gz"
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

