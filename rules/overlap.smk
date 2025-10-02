import pandas as pd
configfile: "config.yaml"
samples_df = pd.read_csv("samplesheet.csv")

# Define the samples
samples = samples_df['sample'].tolist()
filenames_without_extension = [file.rstrip('.bam') for file in samples]
def get_folder_name(samplename_with_bam): # only one sample name for each run allowed#
        only_folder_name=samplename_with_bam.split(".bam")[0]
        return (only_folder_name)

output_dir=config["output_dir"]



rule truvari: # overlap unannotated svs
    input:
        phased_cnv_and_svs="{output_dir}/variants/sawfish_phased_{sample}/{sample}_genotyped.sv.vcf.gz",
        gziped_file="{output_dir}/variants/sniffles_{sample}/{sample}_svs.vcf.gz",
        reference=config["reference"],
 

    output:
        truvari_done="{output_dir}/overlaped_variants/svs_{sample}/{sample}_sniffles_vs_sawfish/summary.json"
    conda:
        "../envs/truvari.yaml"
    log:
        "{output_dir}/logs/truvari_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 2,
        time_hrs=lambda wildcards, attempt: attempt * 2,
        mem_gb=lambda wildcards, attempt: 2 + (attempt * 10)
    params:
        sawf="{output_dir}/overlaped_variants/svs_{sample}/{sample}_sawfish.vcf.gz",
        snfs="{output_dir}/overlaped_variants/svs_{sample}/{sample}_sniffles.vcf.gz",
        dir="{output_dir}/overlaped_variants/svs_{sample}/{sample}_sniffles_vs_sawfish"
    message:
        "Overlaping svs with truvari: {input.phased_cnv_and_svs} and {input.gziped_file}..."
    shell:
        """
        rm -rf {params.dir} 2>{log}
        cp {input.phased_cnv_and_svs} {params.sawf} 2>{log}
        cp {input.gziped_file} {params.snfs} 2>{log}
        tabix -f {params.sawf} 2>{log}
        tabix -f {params.snfs} 2>{log}

        truvari bench -b {params.snfs} -c {params.sawf} -f {input.reference} -o {params.dir} >{log} 2>&1
        """





rule create_rtg_ref:
    input:
        reference=config["reference"],
    output:
        ref_file="{output_dir}/overlaped_variants/ref/summary.txt",
    conda:
        "../envs/rtgtools.yaml"
    log:
        "{output_dir}/logs/rtgtools_reference.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 2,
        time_hrs=lambda wildcards, attempt: attempt * 2,
        mem_gb=lambda wildcards, attempt: 2 + (attempt * 10)
    params:
        ref_dir="{output_dir}/overlaped_variants/ref",
    message:
        "Building ref for snp overlap with rtg-tools: {input.reference}..."
    shell:
        """
        rm -rf {params.ref_dir} 2>{log}
        rtg format {input.reference} -o {params.ref_dir} >{log} 2>&1
        """


rule snp_overlap: # snps
    input:
        ref_file="{output_dir}/overlaped_variants/ref/summary.txt",
        gz_file= "{output_dir}/variants/deepvariant_{sample}/{sample}_variants.vcf.gz" if config["use_deepvariant"] else "{output_dir}/variants/bcftools_{sample}/{sample}_bcft_snps.vcf.gz",
        vcf_nano="{output_dir}/variants/nanocaller_{sample}/{sample}_nanocaller.vcf.gz" # same here
    output:
        summary="{output_dir}/overlaped_variants/snps_{sample}/overlap/summary.txt"
    conda:
        "../envs/rtgtools.yaml"
    log:
        "{output_dir}/logs/rtgtools_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 2,
        time_hrs=lambda wildcards, attempt: attempt * 2,
        mem_gb=lambda wildcards, attempt: 2 + (attempt * 22)
    params:
        out_dir="{output_dir}/overlaped_variants/snps_{sample}/overlap",
        ref_file="{output_dir}/overlaped_variants/ref",
        snps="{output_dir}/overlaped_variants/snps_{sample}/snps_{sample}.vcf.gz",
        nano="{output_dir}/overlaped_variants/snps_{sample}/nano_{sample}.vcf.gz",
    message:
        "comparing snps with rtg-tools: {input.gz_file} and {input.vcf_nano}..."
    shell:
        """
        rm -rf {params.out_dir} 2>{log}
        cp {input.gz_file} {params.snps} 2>{log}
        cp {input.vcf_nano} {params.nano} 2>{log}
        tabix -f {params.snps} 2>{log}
        tabix -f {params.nano} 2>{log}
        rtg vcfeval -t {params.ref_file} -b {params.snps} -c {params.nano} -o {params.out_dir} >{log} 2>&1
        """
 
