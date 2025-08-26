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

include:"prep.smk"
include:"phasers.smk"

# this is the same as the get_output_files except minus the multiqc report
def get_mqc_files():
    all=list()
    if config["use_deepvariant"]:
        all.extend(expand("{output_dir}/variants/deepvariant_{sample}/{sample}_variants.vcf.gz", sample=filenames_without_extension, output_dir=config["output_dir"])),

    if not config["use_deepvariant"]:
        all.extend(expand("{output_dir}/variants/bcftools_{sample}/{sample}_bcft_snps.vcf.gz", sample=filenames_without_extension, output_dir=config["output_dir"])),

    all.extend(expand("{output_dir}/variants/sniffles_{sample}/{sample}_svs.vcf", sample=filenames_without_extension, output_dir=config["output_dir"])),
    all.extend(expand("{output_dir}/variants/nanocaller_{sample}/{sample}_nanocaller.vcf.gz", sample=filenames_without_extension, output_dir=config["output_dir"])),
    all.extend(expand("{output_dir}/variants/sawfish_phased_{sample}/{sample}_genotyped.sv.vcf.gz", sample=filenames_without_extension, output_dir=config["output_dir"])),

    if config["use_kraken2"]:
        all.extend(expand("{output_dir}/kraken2/{sample}_kraken2.report", sample=filenames_without_extension, output_dir=config["output_dir"])),

    all.extend(expand("{output_dir}/variants/trgt_{sample}/{sample}.vcf.gz", sample=filenames_without_extension, output_dir=config["output_dir"])),
    all.extend(expand("{output_dir}/variants/paraphase_{sample}/{sample}_done.flag", sample=filenames_without_extension, output_dir=config["output_dir"])), 
    all.extend(expand("{output_dir}/variants/whatshap_{sample}/{sample}_phased.vcf.gz", sample=filenames_without_extension, output_dir=config["output_dir"])),
    all.extend(expand("{output_dir}/variants/mitorsaw_{sample}/{sample}_mitochondiral_variants.vcf.gz", sample=filenames_without_extension, output_dir=config["output_dir"])),
    all.extend(expand("{output_dir}/variants/hificnv_{sample}/{sample}_hificnv_done.flag", sample=filenames_without_extension, output_dir=config["output_dir"])),    
    all.extend(expand("{output_dir}/bams/{sample}_aligned.bam", sample=filenames_without_extension, output_dir=config["output_dir"])),    
    all.extend(expand("{output_dir}/mosdepth/{sample}.mosdepth.summary.txt", sample=filenames_without_extension, output_dir=config["output_dir"])),   
    return all


def get_output_files():
    all=list()
    if config["use_deepvariant"]:
        all.extend(expand("{output_dir}/variants/deepvariant_{sample}/{sample}_variants.vcf.gz", sample=filenames_without_extension, output_dir=config["output_dir"])),

    if not config["use_deepvariant"]:
        all.extend(expand("{output_dir}/variants/bcftools_{sample}/{sample}_bcft_snps.vcf.gz", sample=filenames_without_extension, output_dir=config["output_dir"])),

    all.extend(expand("{output_dir}/variants/sniffles_{sample}/{sample}_svs.vcf", sample=filenames_without_extension, output_dir=config["output_dir"])),
    all.extend(expand("{output_dir}/variants/nanocaller_{sample}/{sample}_nanocaller.vcf.gz", sample=filenames_without_extension, output_dir=config["output_dir"])),
    all.extend(expand("{output_dir}/variants/sawfish_phased_{sample}/{sample}_genotyped.sv.vcf.gz", sample=filenames_without_extension, output_dir=config["output_dir"])),

    if config["use_kraken2"]:
        all.extend(expand("{output_dir}/kraken2/{sample}_kraken2.report", sample=filenames_without_extension, output_dir=config["output_dir"])),

    all.extend(expand("{output_dir}/variants/longphase_{sample}/{sample}_phased.vcf", sample=filenames_without_extension, output_dir=config["output_dir"])),
    all.extend(expand("{output_dir}/variants/trgt_{sample}/{sample}.vcf.gz", sample=filenames_without_extension, output_dir=config["output_dir"])),
    all.extend(expand("{output_dir}/variants/paraphase_{sample}/{sample}_done.flag", sample=filenames_without_extension, output_dir=config["output_dir"])), 
    all.extend(expand("{output_dir}/variants/whatshap_{sample}/{sample}_phased.vcf.gz", sample=filenames_without_extension, output_dir=config["output_dir"])),
    all.extend(expand("{output_dir}/variants/hificnv_{sample}/{sample}_hificnv_done.flag", sample=filenames_without_extension, output_dir=config["output_dir"])),    
    all.extend(expand("{output_dir}/bams/{sample}_aligned.bam", sample=filenames_without_extension, output_dir=config["output_dir"])),    
    all.extend(expand("{output_dir}/mosdepth/{sample}.mosdepth.summary.txt", sample=filenames_without_extension, output_dir=config["output_dir"])),   
    all.extend(expand("{output_dir}/multiqc_report.html", output_dir=config["output_dir"])) # this line is the only difference to multiqc input, for now
    all.extend(expand("{output_dir}/variants/mitorsaw_{sample}/{sample}_mitochondiral_variants.vcf.gz", sample=filenames_without_extension, output_dir=config["output_dir"])),

    return all

rule all:
    input:
        get_output_files()
 


rule multiqc:
    input:
        get_mqc_files()
    output:
        mqc_report="{output_dir}/multiqc_report.html"
    conda:
        "../envs/multiqc.yaml"
    resources:
        threads=lambda wildcards, attempt: attempt * 2,
        time_hrs=lambda wildcards, attempt: attempt * 1,
        mem_gb=lambda wildcards, attempt: 8 + (attempt * 12)
    message:
        "Creating MultiQC Report..."
    params:
        dir="./{output_dir}"
    log:
        "{output_dir}/logs/Multiqc.log"
    shell:
        """
        multiqc {params.dir} --filename {output} --no-data-dir >> {log} 2>&1
        """


rule deepvariant:
    input:
        bam="{output_dir}/bams/{sample}_aligned.bam"
    output:
        vcf="{output_dir}/variants/deepvariant_{sample}/{sample}_variants.vcf.gz",
        gvcf="{output_dir}/variants/deepvariant_{sample}/{sample}_variants.gvcf.gz",
    conda:
        "../envs/deepvariant.yaml"
    log:
        "{output_dir}/logs/deepvariant_{sample}.log"
    params:
        intermediate_dir="{output_dir}/bams/{sample}_deepvariant_workdir",
        sif_dir=config["sif_image_deepvariant"],
        ref=config["reference"],
        checkpoint_dir=config["deepvariant_checkpoint_dir"],
    resources:
        threads=lambda wildcards, attempt: attempt * 20,
        time_hrs=lambda wildcards, attempt: attempt * 2,
        mem_gb=lambda wildcards, attempt: 36 + (attempt * 12)

    message:
        "Calling short variants for {input.bam} using DeepVariant..."
    shell:
        """
        # only for hpc 
        module load deepvariant/1.9.0 2>{log}
        run_deepvariant --model_type=PACBIO --ref={params.ref} --reads={input.bam} --vcf_stats_report=true --output_vcf={output.vcf} --output_gvcf={output.gvcf} --num_shards {resources.threads} --intermediate_results_dir {params.intermediate_dir} >> {log} 2>&1
        """

rule bcftools_snp:
    input:
        reference=config["reference"], # must be fasta
        bam="{output_dir}/bams/{sample}_aligned.bam"
    output:
        vcf_bcf="{output_dir}/variants/bcftools_{sample}/{sample}_bcft_snps_unfiltered.vcf.gz",
        gz_file="{output_dir}/variants/bcftools_{sample}/{sample}_bcft_snps.vcf.gz"
    params:
        min_qual_filter=config["bcftools_min_q_filter"],
        in_between_file="{output_dir}/variants/bcftools_{sample}/{sample}_bcft_mpileup.bcf",

    conda:
        "../envs/bcftools.yaml"
    log:
        "{output_dir}/logs/bcftools_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 12,
        time_hrs=lambda wildcards, attempt: attempt * 12,
        mem_gb=lambda wildcards, attempt: 48 + (attempt * 12)
    message:
        "Calling SNPs for {input.bam} using bcftools..."
    shell:
        """
        bcftools mpileup -O b -f {input.reference} {input.bam} --threads {resources.threads} -o {params.in_between_file} -Q {params.min_qual_filter} -q {params.min_qual_filter} >> {log} 2>&1
        bcftools call -mv --threads {resources.threads} -Oz -o {output.vcf_bcf} {params.in_between_file} >> {log} 2>&1
        bcftools view {output.vcf_bcf} -i 'QUAL>={params.min_qual_filter}' --threads {resources.threads} -Oz -o {output.gz_file} --write-index >> {log} 2>&1
        """    


# todo: mitorsaw: https://github.com/PacificBiosciences/mitorsaw


rule mitorsaw: # mitochondrial variants, only hg38 compatible
    input:
        reference=config["reference"], # must be fasta
        bam="{output_dir}/bams/{sample}_aligned.bam"
    output:
        mit_vcf="{output_dir}/variants/mitorsaw_{sample}/{sample}_mitochondiral_variants.vcf.gz"
    conda:
        "../envs/mitorsaw.yaml"
    params:
        stats_json="{output_dir}/variants/mitorsaw_{sample}/{sample}_hap_stats.json"
    log:
        "{output_dir}/logs/mitorsaw_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 24,
        time_hrs=lambda wildcards, attempt: attempt * 2,
        mem_gb=lambda wildcards, attempt: 48 + (attempt * 12)
    message:
        "Mitochondrial variant detection for {input.bam} using mitorsaw..."
    shell:
        """
        mitorsaw haplotype --reference {input.reference} --bam {input.bam} --output-vcf {output.mit_vcf} --output-hap-stats {params.stats_json} --output-debug {log}
        tabix -f {output.mit_vcf} 2>{log}
        """    




rule sawfish: # svs and cnv, instead of pbsv + more does only minimal phasing information 
    input:
        reference=config["reference"], # must be fasta
        bam="{output_dir}/bams/{sample}_aligned.bam"
    output:
        #stats="{output_dir}/variants/sawfish_{sample}/run.stats.json", # might now be missing
        phased_cnv_and_svs="{output_dir}/variants/sawfish_phased_{sample}/{sample}_genotyped.sv.vcf.gz" 
    params:
        output_dir="{output_dir}/variants/sawfish_{sample}",
        call_output="{output_dir}/variants/sawfish_phased_{sample}",
        file_to_rename="{output_dir}/variants/sawfish_phased_{sample}/genotyped.sv.vcf.gz" 
    conda:
        "../envs/sawfish.yaml"
    log:
        "{output_dir}/logs/sawfish_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 24,
        time_hrs=lambda wildcards, attempt: attempt * 2,
        mem_gb=lambda wildcards, attempt: 48 + (attempt * 12)
    message:
        "Calling SVs and CNVs for {input.bam} using sawfish..."
    shell:
        """
        sawfish discover --threads {resources.threads} --bam {input.bam} --ref {input.reference} --clobber --output-dir {params.output_dir} >> {log} 2>&1
        sawfish joint-call --threads {resources.threads} --sample {params.output_dir} --output-dir {params.call_output} >> {log} 2>&1
        rm -rf {params.output_dir} >> {log} 2>&1 # no need to store in-between files
        mv {params.file_to_rename} {output} >> {log} 2>&1
        """    

rule paraphase:
    input:
        reference=config["reference"], # must be fasta
        bam="{output_dir}/bams/{sample}_aligned.bam"
    output:
        done_flag= "{output_dir}/variants/paraphase_{sample}/{sample}_done.flag"
    params:
        prefix="{sample}_",
        dir="{output_dir}/variants/paraphase_{sample}/"
    conda:
        "../envs/paraphase.yaml"
    log:
        "{output_dir}/logs/paraphase_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 24,
        time_hrs=lambda wildcards, attempt: attempt * 2,
        mem_gb=lambda wildcards, attempt: 48 + (attempt * 12)
    message:
        "Paralog annotation of {input.bam} ..."
    shell:
        """
        paraphase -b {input.bam} -r {input.reference} -t {resources.threads} -p {params.prefix} -o {params.dir} >> {log} 2>&1 
        touch {output.done_flag}
        """    


rule trgt:
    input:
        reference=config["reference"],
        bam="{output_dir}/bams/{sample}_aligned.bam"
    output:
        vcf="{output_dir}/variants/trgt_{sample}/{sample}.vcf.gz"
        #done_flag= "{output_dir}/variants/trgt_{sample}/{sample}_done.flag"
    params:
        repeats_bed=config["repeats_file"],
        prefix="{output_dir}/variants/trgt_{sample}/{sample}"
    conda:
        "../envs/trgt.yaml"
    log:
        "{output_dir}/logs/trgt_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 24,
        time_hrs=lambda wildcards, attempt: attempt * 2,
        mem_gb=lambda wildcards, attempt: 48 + (attempt * 12)
    message:
        "Tandem repeats genotyping {input.bam} ..."
    shell:
        """
        trgt genotype --genome {input.reference} --reads {input.bam} --repeats {params.repeats_bed} --threads {resources.threads} --output-prefix {params.prefix} >> {log} 2>&1
        """    

rule sniffles:
    input:
        reference=config["reference"], # must be fasta
        bam="{output_dir}/bams/{sample}_aligned.bam"
    output:
        vcf="{output_dir}/variants/sniffles_{sample}/{sample}_svs.vcf",
        gziped_file="{output_dir}/variants/sniffles_{sample}/{sample}_svs.vcf.gz"

    params:
        plot_out="{output_dir}/sniffles_{sample}/{sample}_sniffles_plots/",
        snf="{output_dir}/variants/sniffles_{sample}/{sample}_sniffles.snf",
    conda:
        "../envs/sniffles.yaml"
    log:
        "{output_dir}/logs/sniffles_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 32,
        time_hrs=lambda wildcards, attempt: attempt * 2,
        mem_gb=lambda wildcards, attempt: 48 + (attempt * 12)
    message:
        "Calling structural variants for {input.bam} using Sniffles..."
    shell:
        """
        sniffles --input {input.bam} --vcf {output.vcf} --reference {input.reference} --snf {params.snf} --allow-overwrite >> {log} 2>&1
        # python3 -m sniffles2_plot -i {output.vcf} -o {params.plot_out} >> {log} 2>&1
        bgzip -c {output.vcf} > {output.gziped_file} 2>{log}
        """


rule hificnv: # todo: use haplotagged .bam file as input instead
    input:
        bam="{output_dir}/bams/{sample}_aligned.bam",
        reference=config["reference"], # must be fasta

    output:
        flag_done="{output_dir}/variants/hificnv_{sample}/{sample}_hificnv_done.flag"
    params:
        prefix="{output_dir}/variants/hificnv_{sample}/{sample}"
        
    conda:
        "../envs/hificnv.yaml"
    log:
        "{output_dir}/logs/hificnv_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 32,
        time_hrs=lambda wildcards, attempt: attempt * 2,
        mem_gb=lambda wildcards, attempt: 48 + (attempt * 12)

    message:
        "Detecting CNVs in  {input.bam}..."
    shell:
        """
        hificnv --bam {input.bam} --ref {input.reference} --threads {resources.threads} --output-prefix {params.prefix} >> {log} 2>&1
        touch {output}
        """
