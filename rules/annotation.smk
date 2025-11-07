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
        sawfish_tsv="{output_dir}/annotated_variants/sansa_svs_cnvs_sawfish_{sample}/{sample}_sawfish_annotated.tsv.gz",
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
        sansa annotate -d {params.annotation_sv_file} {input.svs_phased} -a {output.annotated_bcf} -o {output.tsv} >{log} 2>&1
        sansa annotate -d {params.annotation_sv_file} {input.phased_cnv_and_svs} -a {output.sawfish_bcf} -o {output.sawfish_tsv} >{log} 2>&1
        """


rule snpsift: # snps
    input:
        vcf_phased_longp="{output_dir}/variants/longphase_{sample}/{sample}_phased.vcf",
        phased_vcf_whatsh="{output_dir}/variants/whatshap_{sample}/{sample}_phased.vcf.gz", # need to unpack, maybe
        vcf_nano="{output_dir}/variants/nanocaller_{sample}/{sample}_nanocaller.vcf.gz" # same here
    output:
        longp_snp="{output_dir}/annotated_variants/snpsift_longphase_{sample}/{sample}_snpsift_longphase_annotated.vcf",
        whatsh_snp="{output_dir}/annotated_variants/snpsift_whatshap_{sample}/{sample}_snpsift_whatshap_annotated.vcf",
        nanoc_smp="{output_dir}/annotated_variants/snpsift_nanocaller_{sample}/{sample}_snpsift_nanocaller_annotated.vcf"
    conda:
        "../envs/snpeff.yaml"
    log:
        "{output_dir}/logs/snpsift_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 2,
        time_hrs=lambda wildcards, attempt: attempt * 2,
        mem_gb=lambda wildcards, attempt: 2 + (attempt * 22)
    params:
        annotation_snp_file=config["snp_annotation_file"] # dnsnp, hg38
    message:
        "Annotating the snps with snpsift: {input.vcf_phased_longp}, {input.phased_vcf_whatsh} and {input.vcf_nano}..."
    shell:
        """
        export _JAVA_OPTIONS="-Xmx12g" # otherwise annotation will run out of java heapspace at least 12 
        export _JAVA_OPTIONS="-Xmx{resources.mem_gb}g" # otherwise annotation will run out of java heapspace dynamically as much as the job can give

        SnpSift annotate {params.annotation_snp_file} {input.vcf_phased_longp} >{output.longp_snp} 2>{log}
        SnpSift annotate {params.annotation_snp_file} {input.phased_vcf_whatsh} >{output.whatsh_snp} 2>{log}
        SnpSift annotate {params.annotation_snp_file} {input.vcf_nano} >{output.nanoc_smp} 2>{log}
        """


rule vep:
    input:
        vcf_phased_longp="{output_dir}/variants/longphase_{sample}/{sample}_phased.vcf",
        phased_vcf_whatsh="{output_dir}/variants/whatshap_{sample}/{sample}_phased.vcf.gz", # need to unpack, maybe
        vcf_nano="{output_dir}/variants/nanocaller_{sample}/{sample}_nanocaller.vcf.gz" # same here
    output:
        longp_snp="{output_dir}/annotated_variants/vep_longphase_{sample}/{sample}_vep_longphase_annotated.vcf",
        whatsh_snp="{output_dir}/annotated_variants/vep_whatshap_{sample}/{sample}_vep_whatshap_annotated.vcf",
        nanoc_smp="{output_dir}/annotated_variants/vep_nanocaller_{sample}/{sample}_vep_nanocaller_annotated.vcf"
    conda:
        "../envs/vep.yaml"
    log:
        "{output_dir}/logs/vep_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 12,
        time_hrs=lambda wildcards, attempt: attempt * 5,
        mem_gb=lambda wildcards, attempt: 2 + (attempt * 10)
    params:
        cache_dir=config["vep_cache_dir"],
        whatsh_unzp="{output_dir}/variants/whatshap_{sample}/{sample}_phased.vcf",
        nanocaller_unzipped="{output_dir}/variants/nanocaller_{sample}/{sample}_nanocaller.vcf",
        dir_long="{output_dir}/annotated_variants/vep_longphase_{sample}",
        dir_whats="{output_dir}/annotated_variants/vep_whatshap_{sample}",
        dir_nanoc="{output_dir}/annotated_variants/vep_nanocaller_{sample}",
    message:
        "Annotating the snps with vep: {input.vcf_phased_longp}, {input.phased_vcf_whatsh} and {input.vcf_nano}..."
    shell:
        """
        gunzip {input.phased_vcf_whatsh} -f -c >{params.whatsh_unzp}
        gunzip {input.vcf_nano} -f -c >{params.nanocaller_unzipped}
        rm -f {params.dir_long}/*
        rm -f {params.dir_whats}/*
        rm -f {params.dir_nanoc}/*
        vep --offline --dir_cache {params.cache_dir} -i {input.vcf_phased_longp} -o {output.longp_snp} --everything --force_overwrite --fork {resources.threads} --warning_file {log} >>{log} 2>&1
        vep --offline --dir_cache {params.cache_dir} -i {params.whatsh_unzp} -o {output.whatsh_snp} --everything --force_overwrite --fork {resources.threads} --warning_file {log} >>{log} 2>&1
        vep --offline --dir_cache {params.cache_dir} -i {params.nanocaller_unzipped} -o {output.nanoc_smp} --everything --force_overwrite --fork {resources.threads} --warning_file {log} >>{log} 2>&1
        """


rule annotsv:
    input:
        svs_phased="{output_dir}/variants/longphase_{sample}/{sample}_phased_SV.vcf",
        phased_cnv_and_svs="{output_dir}/variants/sawfish_phased_{sample}/{sample}_genotyped.sv.vcf.gz" 
    output:
        snfls="{output_dir}/annotated_variants/annotsv_sniffles_{sample}/{sample}_phased_SV.annotated.tsv",
        sawfs="{output_dir}/annotated_variants/annotsv_sawfish_{sample}/{sample}_genotyped.sv.annotated.tsv",
    conda:
        "../envs/annotsv.yaml"
    log:
        "{output_dir}/logs/annotsv_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 2,
        time_hrs=lambda wildcards, attempt: attempt * 2,
        mem_gb=lambda wildcards, attempt: 2 + (attempt * 10)
    params:
        annotsv_data=config["annotsv_data_dir"],
        dir_out_snfls="{output_dir}/annotated_variants/annotsv_sniffles_{sample}",
        dir_out_sawfs="{output_dir}/annotated_variants/annotsv_sawfish_{sample}",
        unziped_safw="{output_dir}/variants/sawfish_phased_{sample}/{sample}_genotyped.sv.vcf",
        parental_dir="{output_dir}",
        output_file1="{output_dir}/{sample}_snfls.tsv",
        output_file2="{output_dir}/{sample}_sawf.tsv",
    message:
        "Annotating the svs from longphase (sniffles) and sawfish with annotsv: {input.svs_phased} and {input.phased_cnv_and_svs}..."
    shell:
        """
        rm -rf {params.dir_out_snfls} >>{log} 2>&1
        rm -rf {params.dir_out_sawfs} >>{log} 2>&1
        mkdir -p {params.dir_out_snfls} >>{log} 2>&1 # annotsv needs these dirs to be present before it starts
        mkdir -p {params.dir_out_sawfs} >>{log} 2>&1
        gunzip {input.phased_cnv_and_svs} -f -c >{params.unziped_safw}
        AnnotSV -annotationsDir {params.annotsv_data} -SVinputFile {input.svs_phased} -outputDir {params.parental_dir} -outputFile {params.output_file1}  >>{log} 2>&1        
        AnnotSV -annotationsDir {params.annotsv_data} -SVinputFile {params.unziped_safw} -outputDir {params.parental_dir} -outputFile {params.output_file2}  >>{log} 2>&1
        mv {params.output_file1} {output.snfls} >>{log} 2>&1
        mv {params.output_file2} {output.sawfs} >>{log} 2>&1
        rm -f {params.parental_dir}/*unannotated.tsv >>{log} 2>&1
        rm -f {params.parental_dir}/*.bash >>{log} 2>&1
        rm -f {params.parental_dir}/*.bed >>{log} 2>&1
        """


