import pandas as pd
configfile: "config.yaml"
samples_df = pd.read_csv("samplesheet.csv")

# Define the samples
samples = samples_df['sample'].tolist()
filenames_without_extension = [file.rstrip('.bam') for file in samples]
#print("filenames without extension:")
#print(filenames_without_extension)
#filenames_without_extension = [os.path.splitext(samples)[0] for sample in samples]
#bamfiles = dict(zip(samples_df['sample_name'], samples_df['bam_file']))
def get_folder_name(samplename_with_bam): # only one sample name for each run allowed#
        only_folder_name=samplename_with_bam.split(".bam")[0]
        return (only_folder_name)
# use this function later to remove the .bam files from the output folder name

# how to run: snakemake -c 12 --jobs 2 -s rules/snakefile.smk --jobscript scripts/jobscript.sh --use-conda

#get_sample_short=
#wildcard_constraints:
#        sample = "[^/]+\\.bam$"  # Matches any string that ends with a single .bam

output_dir=config["output_dir"]
# get the correct files to be named and remove the .bam_aligned.bam -> aligned.bam, already created filenames_without_extension, but can do a function



def get_output_files():
    all=list()
    if config["use_deepvariant"]:
        all.extend(expand("{output_dir}/variants/deepvariant_{sample}/{sample}_variants.vcf.gz", sample=filenames_without_extension, output_dir=config["output_dir"])),

    if not config["use_deepvariant"]:
        all.extend(expand("{output_dir}/variants/bcftools_{sample}/{sample}_bcft_snps.vcf.gz", sample=filenames_without_extension, output_dir=config["output_dir"])),

    all.extend(expand("{output_dir}/variants/sniffles_{sample}/{sample}_svs.vcf", sample=filenames_without_extension, output_dir=config["output_dir"])),
    all.extend(expand("{output_dir}/variants/nanocaller_{sample}/{sample}_nanocaller.vcf.gz", sample=filenames_without_extension, output_dir=config["output_dir"])),
    all.extend(expand("{output_dir}/variants/sawfish_phased_{sample}/genotyped.sv.vcf.gz", sample=filenames_without_extension, output_dir=config["output_dir"])),
    #  phased_cnv_and_svs="{output_dir}/variants/sawfish_phased_{sample}/genotyped.sv.vcf.gz" 
 #   all.extend(expand("{output_dir}/variants/hiphase_{sample}/{sample}_svs_phased.vcf.gz", sample=filenames_without_extension, output_dir=config["output_dir"])),
    #fix that either bcftools or deepvariant output is used
    
    #         vcf="{output_dir}/variants/trgt_{sample}/{sample}.vcf.gz"

    all.extend(expand("{output_dir}/variants/trgt_{sample}/{sample}.vcf.gz", sample=filenames_without_extension, output_dir=config["output_dir"])),
    all.extend(expand("{output_dir}/variants/paraphase_{sample}/{sample}_done.flag", sample=filenames_without_extension, output_dir=config["output_dir"])), 
    all.extend(expand("{output_dir}/variants/whatshap_{sample}/{sample}_phased.vcf", sample=filenames_without_extension, output_dir=config["output_dir"])),
    all.extend(expand("{output_dir}/variants/hificnv_{sample}/{sample}_hificnv_done.flag", sample=filenames_without_extension, output_dir=config["output_dir"])),    
    all.extend(expand("{output_dir}/bams/{sample}_aligned.bam", sample=filenames_without_extension, output_dir=config["output_dir"])),    
    all.extend(expand("{output_dir}/mosdepth/{sample}.mosdepth.summary.txt", sample=filenames_without_extension, output_dir=config["output_dir"])),   
    all.extend(expand("{output_dir}/multiqc_report.html", output_dir=config["output_dir"])) # this line is the only difference to multiqc input, for now

    return all

rule all:
    input:
        get_output_files()
 


rule multiqc:
    input:
    #    expand("{output_dir}/variants/{sample}_variants.vcf", sample=filenames_without_extension, output_dir=config["output_dir"]), # deepvariant

    # the hiphase output should be the actually phased snps and svs, but hiphase does not work for now
      #  expand("{output_dir}/variants/hiphase_{sample}/{sample}_svs_phased.vcf.gz", sample=filenames_without_extension, output_dir=config["output_dir"]),
      # fix that either bcftools or deepvariant output is used



        expand("{output_dir}/variants/sniffles_{sample}/{sample}_svs.vcf", sample=filenames_without_extension, output_dir=config["output_dir"]),
        expand("{output_dir}/variants/nanocaller_{sample}/{sample}_nanocaller.vcf.gz", sample=filenames_without_extension, output_dir=config["output_dir"]),
        expand("{output_dir}/variants/whatshap_{sample}/{sample}_phased.vcf", sample=filenames_without_extension, output_dir=config["output_dir"]),
        expand("{output_dir}/variants/hificnv_{sample}/{sample}_hificnv_done.flag", sample=filenames_without_extension, output_dir=config["output_dir"]),    
        expand("{output_dir}/bams/{sample}_aligned.bam", sample=filenames_without_extension, output_dir=config["output_dir"]),
        expand("{output_dir}/mosdepth/{sample}.mosdepth.summary.txt", sample=filenames_without_extension, output_dir=config["output_dir"]),   
   
 

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
        ref=config["reference"]
    resources:
        threads=lambda wildcards, attempt: attempt * 48,
        time_hrs=lambda wildcards, attempt: attempt * 1,
        mem_gb=lambda wildcards, attempt: 64 + (attempt * 12)

    message:
        "Calling short variants for {input.bam} using DeepVariant..."
    shell:
        """
        # only for hpc 
        module load deepvariant 2>{log}
        deepvariant --model_type=PACBIO --ref={params.ref} --reads={input.bam} ---vcf_stats_report=true --output_vcf={output.vcf} --output_gvcf={output.gvcf} --num_shards {resources.threads} --intermediate_results_dir {params.intermediate_dir} >> {log} 2>&1
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
        threads=lambda wildcards, attempt: attempt * 24,
        time_hrs=lambda wildcards, attempt: attempt * 1,
        mem_gb=lambda wildcards, attempt: 48 + (attempt * 12)
    message:
        "Calling SNPs for {input.bam} using bcftools..."
    shell:
        """
        bcftools mpileup -O b -f {input.reference} {input.bam} --threads {resources.threads} -o {params.in_between_file} >> {log} 2>&1
        bcftools call -mv --threads {resources.threads} -Oz -o {output.vcf_bcf} {params.in_between_file} >> {log} 2>&1
        bcftools view {output.vcf_bcf} -i 'QUAL>={params.min_qual_filter}' --threads {resources.threads} -Oz -o {output.gz_file} -W >> {log} 2>&1
        """    

# command > out 2>error


rule sawfish: # svs and cnv, instead of pbsv + more
    input:
        reference=config["reference"], # must be fasta
        bam="{output_dir}/bams/{sample}_aligned.bam"
    output:
        #stats="{output_dir}/variants/sawfish_{sample}/run.stats.json", # might now be missing
        phased_cnv_and_svs="{output_dir}/variants/sawfish_phased_{sample}/genotyped.sv.vcf.gz" 
    params:
        output_dir="{output_dir}/variants/sawfish_{sample}",
        call_output="{output_dir}/variants/sawfish_phased_{sample}"
    conda:
        "../envs/sawfish.yaml"
    log:
        "{output_dir}/logs/sawfish_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 24,
        time_hrs=lambda wildcards, attempt: attempt * 1,
        mem_gb=lambda wildcards, attempt: 48 + (attempt * 12)
    message:
        "Calling SVs and CNVs for {input.bam} using sawfish..."
    shell:
        """
        sawfish discover --threads {resources.threads} --bam {input.bam} --ref {input.reference} --clobber --output-dir {params.output_dir} >> {log} 2>&1
        sawfish joint-call --threads {resources.threads} --sample {params.output_dir} --output-dir {params.call_output} >> {log} 2>&1
        rm -rf {params.output_dir} >> {log} 2>&1 # no need to store in-between files
        """    
    # new pacbio stuff

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
        time_hrs=lambda wildcards, attempt: attempt * 1,
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
        reference=config["reference"], # must be fasta
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
        time_hrs=lambda wildcards, attempt: attempt * 1,
        mem_gb=lambda wildcards, attempt: 48 + (attempt * 12)
    message:
        "Tandem repeats genotyping {input.bam} ..."
    shell:
        """
        trgt genotype --genome {input.reference} --reads {input.bam} --repeats {params.repeats_bed} --threads {resources.threads} --output-prefix {params.prefix} >> {log} 2>&1
        """    



rule hiphase: # phases snps, svs and more
    input:
        reference=config["reference"], # must be fasta
        gz_file= "{output_dir}/variants/deepvariant_{sample}/{sample}_variants.vcf.gz" if config["use_deepvariant"] else "{output_dir}/variants/bcftools_{sample}/{sample}_bcft_snps.vcf.gz",
# add svs 
        svs="{output_dir}/variants/sniffles_{sample}/{sample}_svs.vcf.gz",
        #gz_file="{output_dir}/variants/bcftools_{sample}/{sample}_bcft_snps.vcf.gz" if config["use_deepvariant"] == 'False' else "{output_dir}/variants/deepvariant_{sample}/{sample}_variants.vcf.gz",
        bam="{output_dir}/bams/{sample}_aligned.bam"
    output: # vcf_bcf="{output_dir}/variants/bcftools_{sample}/{sample}_bcft_snps.vcf"
        vcf_phased="{output_dir}/variants/hiphase_{sample}/{sample}_svs_phased.vcf.gz"
    conda:
        "../envs/hiphase.yaml"
    log:
        "{output_dir}/logs/hiphase_{sample}.log"
    params:
        merged_file= "{output_dir}/variants/hiphase_{sample}/{sample}_merged_svs_and_snps.vcf.gz"
    resources:
        threads=lambda wildcards, attempt: attempt * 24,
        time_hrs=lambda wildcards, attempt: attempt * 1,
        mem_gb=lambda wildcards, attempt: 48 + (attempt * 12)
    message:
        "Phasing snps and svs with hiphase for {input.bam} ..."
    shell:
        """
        tabix -f {input.gz_file} 2>{log}
        tabix -f {input.svs} 2>{log}
        # needs a singular input file
        hiphase --bam {input.bam} --reference {input.reference} --threads {resources.threads} --vcf {input.svs} --output-vcf {output.vcf_phased}  --ignore-read-groups --min-vcf-qual 40 >> {log} 2>&1
        """    




# only active once 
# phasing snps for 
# rework output: for each tool output one folder
# including the log file maybe
# in that folder all samples the same output files


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
        threads=lambda wildcards, attempt: attempt * 24,
        time_hrs=lambda wildcards, attempt: attempt * 1,
        mem_gb=lambda wildcards, attempt: 48 + (attempt * 12)
    message:
        "Calling SNPs and SVs for {input.bam} using NanoCaller..."
    shell:
        """
        NanoCaller --bam {input.bam} --ref {input.reference} --cpu {resources.threads} --mode all --preset ccs --output {params.path_out} --prefix {params.prefix} --phase >> {log} 2>&1
        """    




rule sniffles:
    input:
        reference=config["reference"], # must be fasta
        bam="{output_dir}/bams/{sample}_aligned.bam"
    output:
        vcf="{output_dir}/variants/sniffles_{sample}/{sample}_svs.vcf"
    params:
        plot_out="{output_dir}/sniffles_{sample}/{sample}_sniffles_plots/",
        snf="{output_dir}/variants/sniffles_{sample}/{sample}_sniffles.snf",
        gziped_file="{output_dir}/variants/sniffles_{sample}/{sample}_svs.vcf.gz"
    conda:
        "../envs/sniffles.yaml"
    log:
        "{output_dir}/logs/sniffles_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 32,
        time_hrs=lambda wildcards, attempt: attempt * 1,
        mem_gb=lambda wildcards, attempt: 48 + (attempt * 12)
    message:
        "Calling structural variants for {input.bam} using Sniffles..."
    shell:
        """
        sniffles --input {input.bam} --vcf {output} --reference {input.reference} --snf {params.snf} --allow-overwrite >> {log} 2>&1
        # python3 -m sniffles2_plot -i {output} -o {params.plot_out} >> {log} 2>&1
        bgzip -c {output.vcf} > {params.gziped_file} 2>{log}
        """
#         bcftools view {output.vcf} -Oz -o {params.gziped_file} >> {log} 2>&1



rule whatshap:
    input: 
        vcf= "{output_dir}/variants/deepvariant_{sample}/{sample}_variants.vcf.gz" if config["use_deepvariant"] else "{output_dir}/variants/bcftools_{sample}/{sample}_bcft_snps.vcf.gz",
        bam="{output_dir}/bams/{sample}_aligned.bam",
        reference=config["reference"], # must be fasta

    output:
        phased_vcf="{output_dir}/variants/whatshap_{sample}/{sample}_phased.vcf"
    conda:
        "../envs/whatshap.yaml"
    log:
        "{output_dir}/logs/whatshap_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 32,
        time_hrs=lambda wildcards, attempt: attempt * 1,
        mem_gb=lambda wildcards, attempt: 48 + (attempt * 12)
    params:
        haplotaged_bam="{output_dir}/variants/whatshap_{sample}/{sample}_haplotaged.bam",
        sorted="{output_dir}/variants/whatshap_{sample}/{sample}_phased_sorted.vcf.gz",
        folder="{output_dir}/variants/whatshap_{sample}"
    message:
        "Phasing haplotypes for {input.bam}..."
    shell:
        """
        whatshap phase --output {output} --reference {input.reference} {input.vcf} {input.bam} --ignore-read-groups >> {log} 2>&1
        bcftools sort {output} -Oz -o {params.sorted} >> {log} 2>&1
        bcftools index -t  {params.sorted} >> {log} 2>&1
        mkdir -p {params.folder} >> {log} 2>&1
        whatshap haplotag {params.sorted} {input.bam} --output {params.haplotaged_bam} --reference {input.reference} --output-threads {resources.threads} --ignore-read-groups >> {log} 2>&1
        whatshap stats {output.phased_vcf} >> {log} 2>&1
        """

rule hificnv:
    input:
        bam="{output_dir}/bams/{sample}_aligned.bam",
        reference=config["reference"], # must be fasta

    output:
        flag_done="{output_dir}/variants/hificnv_{sample}/{sample}_hificnv_done.flag"
        #cnv_track="{output_dir}/variants/{sample}.copynum.bedgraph"
    params:
        prefix="{output_dir}/variants/hificnv_{sample}/{sample}"
        
        # hificnv for some reason adds the sample number into the output filename. need to get it first with samtools
    conda:
        "../envs/hificnv.yaml"
    log:
        "{output_dir}/logs/hificnv_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 32,
        time_hrs=lambda wildcards, attempt: attempt * 1,
        mem_gb=lambda wildcards, attempt: 48 + (attempt * 12)

    message:
        "Detecting CNVs in  {input.bam}..."
    shell:
        """
        hificnv --bam {input.bam} --ref {input.reference} --threads {resources.threads} --output-prefix {params.prefix} >> {log} 2>&1
        touch {output}
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



########################
#
# Tools to add:
# cpg meth tools:
#Genotype tandem repeats - produce spanning bams and vcf (TRGT)
# add a whole set of rules for trio calling
# 
# 
# sv filters: svpack, silvar, bcftools, phrank
# volcanosv, get deepvariant on hpc working
