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





rule hiphase: # phases snps, svs and more but is not compatible with eiter bcftools or sniffles output files
    input:
        reference=config["reference"], # must be fasta
        gz_file= "{output_dir}/variants/deepvariant_{sample}/{sample}_variants.vcf.gz" if config["use_deepvariant"] else "{output_dir}/variants/bcftools_{sample}/{sample}_bcft_snps.vcf.gz",
        svs="{output_dir}/variants/sniffles_{sample}/{sample}_svs.vcf.gz",
        bam="{output_dir}/bams/{sample}_aligned.bam"
    output:
        vcf_phased="{output_dir}/variants/hiphase_{sample}/{sample}_snps_phased.vcf.gz"
    conda:
        "../envs/hiphase.yaml"
    log:
        "{output_dir}/logs/hiphase_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 12,
        time_hrs=lambda wildcards, attempt: attempt * 1,
        mem_gb=lambda wildcards, attempt: 48 + (attempt * 12)
    message:
        "Phasing snps and svs with hiphase for {input.bam} ..."
    shell:
        """
        tabix -f {input.gz_file} 2>{log}
        tabix -f {input.svs} 2>{log}
        hiphase --bam {input.bam} --reference {input.reference} --threads {resources.threads} --vcf {input.gz_file} --output-vcf {output.vcf_phased}  --ignore-read-groups --min-vcf-qual 40 >> {log} 2>&1
        """    




rule nanocaller: # output snps are already haplotaged
    input:
        reference=config["reference"], # must be fasta
        bam="{output_dir}/bams/{sample}_aligned.bam"
    output:
        vcf_nano="{output_dir}/variants/nanocaller_{sample}/{sample}_nanocaller.vcf.gz"
    params:
        path_out="{output_dir}/variants/nanocaller_{sample}/",
        prefix="{sample}_nanocaller"
    conda:
        "../envs/nanocaller.yaml"
    log:
        "{output_dir}/logs/nanocaller_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 12,
        time_hrs=lambda wildcards, attempt: attempt * 2,
        mem_gb=lambda wildcards, attempt: 48 + (attempt * 12)
    message:
        "Calling SNPs and SVs for {input.bam} using NanoCaller..."
    shell:
        """
        NanoCaller --bam {input.bam} --ref {input.reference} --cpu {resources.threads} --mode all --preset ccs --output {params.path_out} --prefix {params.prefix} --phase >> {log} 2>&1
        """    


rule whatshap: # only able to haplotype snps, cannot use svs. for this longphase is used
    input: 
        vcf= "{output_dir}/variants/deepvariant_{sample}/{sample}_variants.vcf.gz" if config["use_deepvariant"] else "{output_dir}/variants/bcftools_{sample}/{sample}_bcft_snps.vcf.gz",
        bam="{output_dir}/bams/{sample}_aligned.bam",
        reference=config["reference"], # must be fasta

    output:
        phased_vcf="{output_dir}/variants/whatshap_{sample}/{sample}_phased.vcf.gz",
        haplotaged_bam="{output_dir}/variants/whatshap_{sample}/{sample}_haplotaged.bam",

    conda:
        "../envs/whatshap.yaml"
    log:
        "{output_dir}/logs/whatshap_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 12,
        time_hrs=lambda wildcards, attempt: attempt * 12,
        mem_gb=lambda wildcards, attempt: 52 + (attempt * 12)
    params:
        sorted="{output_dir}/variants/whatshap_{sample}/{sample}_phased_sorted.vcf.gz",
        folder="{output_dir}/variants/whatshap_{sample}",
        stats_file="{output_dir}/variants/whatshap_{sample}/whatshap_{sample}_stats.tsv",
        packed="{output_dir}/variants/whatshap_{sample}/{sample}_phased.vcf.gz",
    message:
        "Phasing haplotypes for {input.bam}..."
    shell:
        """
        whatshap phase -o {output.phased_vcf} --reference {input.reference} {input.vcf} {input.bam} --ignore-read-groups  2>{log}
        tabix -f {output.phased_vcf} 2>{log}
        whatshap haplotag {params.packed} {input.bam} --output {output.haplotaged_bam} --reference {input.reference} --output-threads {resources.threads} --ignore-read-groups 2>{log}
        whatshap stats {params.packed} --tsv > {params.stats_file} 2>{log}
        """

rule longphase: # phases snps, svs and more
    input:
        reference=config["reference"], # must be fasta
        gz_file= "{output_dir}/variants/deepvariant_{sample}/{sample}_variants.vcf.gz" if config["use_deepvariant"] else "{output_dir}/variants/bcftools_{sample}/{sample}_bcft_snps.vcf.gz",
        svs="{output_dir}/variants/sniffles_{sample}/{sample}_svs.vcf.gz",
        bam="{output_dir}/bams/{sample}_aligned.bam"
    output:
        vcf_phased="{output_dir}/variants/longphase_{sample}/{sample}_phased.vcf",
        svs_phased="{output_dir}/variants/longphase_{sample}/{sample}_phased_SV.vcf",
        long_hap_bam="{output_dir}/variants/longphase_{sample}/{sample}_aligned_haplotaged.bam",
    conda:
        "../envs/longphase.yaml"
    log:
        "{output_dir}/logs/longphase_{sample}.log"
    params:
        prefix="{output_dir}/variants/longphase_{sample}/{sample}_phased",
        prefix_bam="{output_dir}/variants/longphase_{sample}/{sample}_aligned_haplotaged",
    resources:
        threads=lambda wildcards, attempt: attempt * 12,
        time_hrs=lambda wildcards, attempt: attempt * 3,
        mem_gb=lambda wildcards, attempt: 48 + (attempt * 12)
    message:
        "Phasing snps and svs with longphase for {input.bam} ..."
    shell:
        """
        tabix -f {input.gz_file} 2>{log}
        tabix -f {input.svs} 2>{log}
        longphase phase -s {input.gz_file} -b {input.bam} -r {input.reference} --sv-file={input.svs} --pb --indels -t {resources.threads} -o {params.prefix} 2>{log}
        longphase haplotag -r {input.reference} -s {output.vcf_phased} --sv-file {output.svs_phased} -b {input.bam} -t {resources.threads} -o {params.prefix_bam} 2>{log}
        """