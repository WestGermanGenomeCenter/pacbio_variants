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


rule sansa: # annotate svs
    input:
        svs_phased="{output_dir}/variants/longphase_{sample}/{sample}_phased_SV.vcf",
        phased_cnv_and_svs="{output_dir}/variants/sawfish_phased_{sample}/{sample}_genotyped.sv.vcf.gz" 
    output:
        annotated_bcf="{output_dir}/annotated_variants/sansa_svs_sniffles_{sample}/{sample}_sniffles_longphase_annotated.bcf",
        tsv="{output_dir}/annotated_variants/sansa_svs_sniffles_{sample}/{sample}_sniffles_longphase_annotated.tsv.gz",
        sawfish_bcf="{output_dir}/annotated_variants/sansa_svs_cnvs_sawfish_{sample}/{sample}_sawfish_annotated.bcf",
        sawfish_tsv="{output_dir}/annotated_variants/sansa_svs_sniffles_{sample}/{sample}_sawfish_annotated.tsv.gz",
    conda:
        "../envs/sansa.yaml"
    log:
        "{output_dir}/logs/sansa_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 2,
        time_hrs=lambda wildcards, attempt: attempt * 2,
        mem_gb=lambda wildcards, attempt: 2 + (attempt * 10)
    params:
        annotation_sv_file=config["sv_annotation_file"] # maybe start here with the hgsvc data?
    message:
        "Annotating the svs from longphase (sniffles) and sawfish with sansa: {input.svs_phased} and {input.phased_cnv_and_svs}..."
    shell:
        """
        sansa annotate -d {params.annotation_sv_file} {input.svs_phased} -a {output.annotated_bcf} -o {output.tsv} >> {log} 2>&1
        sansa annotate -d {params.annotation_sv_file} {input.phased_cnv_and_svs} -a {output.sawfish_bcf} -o {output.sawfish_tsv} >> {log} 2>&1
        """


rule snpsift: # snps
    input:
        vcf_phased_longp="{output_dir}/variants/longphase_{sample}/{sample}_phased.vcf",
        phased_vcf_whatsh="{output_dir}/variants/whatshap_{sample}/{sample}_phased.vcf.gz", # need to unpack, maybe
        vcf_nano="{output_dir}/variants/nanocaller_{sample}/{sample}_nanocaller.vcf.gz" # same here
    output:
        longp_snp="{output_dir}/annotated_variants/snps_longphase_{sample}/{sample}_snp_longphase_annotated.vcf",
        whatsh_snp="{output_dir}/annotated_variants/snps_whatshap_{sample}/{sample}_snp_whatshap_annotated.vcf",
        nanoc_smp="{output_dir}/annotated_variants/snps_nanocaller_{sample}/{sample}_snp_nanocaller_annotated.vcf"
    conda:
        "../envs/snpeff.yaml"
    log:
        "{output_dir}/logs/snpeff_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 2,
        time_hrs=lambda wildcards, attempt: attempt * 2,
        mem_gb=lambda wildcards, attempt: 2 + (attempt * 10)
    params:
        annotation_snp_file=config["snp_annotation_file"] # dnsnp, hg38
    message:
        "Annotating the snps with snpsift: {input.vcf_phased_longp}, {input.phased_vcf_whatsh} and {input.vcf_nano}..."
    shell:
        """
        SnpSift annotate {params.annotation_snp_file} {input.vcf_phased_longp} >{output.longp_snp} 2>{log}
        SnpSift annotate {params.annotation_snp_file} {input.phased_vcf_whatsh} >{whatsh_snp} 2>{log}
        SnpSift annotate {params.annotation_snp_file} {input.vcf_nano} >{output.nanoc_smp} 2>{log}
        """


# maybe future sv annotation: https://strvctvre.berkeley.edu/