# tertiary tools like:
# svtopo
# kivvi
# sawshark 
# humanatee
# svpack


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


rule svtopo:
    input:
        bam="{output_dir}/bams/{sample}_aligned.bam",
        phased_cnv_and_svs="{output_dir}/variants/sawfish_phased_{sample}/{sample}_genotyped.sv.vcf.gz"
    output:
        svtopo_flag="{output_dir}/visualizations/svtopo_{sample}/{sample}_done.flag"
    params:
        exclude_regions_file=config["svtopo_exclude_regions"],
        svtopo_dir="{output_dir}/visualizations/svtopo_{sample}/",
        prefix="{sample}"
    conda:
        "../envs/svtopo.yaml"
    log:
        "{output_dir}/logs/svtopo_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 1,
        time_hrs=lambda wildcards, attempt: attempt * 2,
        mem_gb=lambda wildcards, attempt: 4 + (attempt * 16)
    message:
        "Visualizing Sawfish SVs with svtopo for sample {input.bam} ..."
    shell:
        """
        rm -rf {params.svtopo_dir} 2>{log}
        svtopo --bam {input.bam} --svtopo-dir {params.svtopo_dir} --prefix {params.prefix} --vcf {input.phased_cnv_and_svs} --exclude-regions {params.exclude_regions_file} >{log} 2>&1
        svtopovz --svtopo-dir {params.svtopo_dir} >{log} 2>&1
        touch {output} >{log} 2>&1
        """


rule paraviewer:
    input:
        done_flag= "{output_dir}/variants/paraphase_{sample}/{sample}_done.flag",
        reference=config["reference"], # must be fasta
        bam="{output_dir}/bams/{sample}_aligned.bam"
    output:
        viewer_done="{output_dir}/visualizations/paraviewer_{sample}/{sample}_done.flag"
    params:
        prefix="{sample}_",
        dir="{output_dir}/visualizations/paraviewer_{sample}/",
        para_dir="{output_dir}/variants/paraphase_{sample}/"
    conda:
        "../envs/paraviewer.yaml"
    log:
        "{output_dir}/logs/paraviewer_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 2,
        time_hrs=lambda wildcards, attempt: attempt * 1,
        mem_gb=lambda wildcards, attempt: 4 + (attempt * 12)
    message:
        "Paralog annotation / visualization for {input.bam} ..."
    shell:
        """
        paraviewer --outdir {params.dir} --paraphase-dir {params.para_dir} --genome hg38 >{log} 2>&1
        touch {output.viewer_done}
        """


# chromograph, this is test-only for now



rule mosdepth_windows:
    input:
        long_hap_bam="{output_dir}/variants/longphase_{sample}/{sample}_aligned_haplotaged.bam",
    output:
        bed = config["output_dir"] + "/qc/{sample}/{sample}.windows.regions.bed.gz"
    params:
        prefix = config["output_dir"] + "/qc/{sample}/{sample}.windows",
        step = config.get("chromograph_step", 5000)
    conda:
        "../envs/mosdepth.yaml"
    log:
        config["output_dir"] + "/logs/mosdepth_windows_{sample}.log"
    resources:
        threads = lambda wildcards, attempt: 4,
        time_hrs = lambda wildcards, attempt: attempt * 2,
        mem_gb = lambda wildcards, attempt: 4 * attempt
    message:
        "Windowed coverage for chromograph, sample {wildcards.sample}..."
    shell:
        "mosdepth -t {resources.threads} --by {params.step} --no-per-base --fast-mode {params.prefix} {input.bam} > {log} 2>&1"


rule chromograph_coverage:
    input:
        bed = config["output_dir"] + "/qc/{sample}/{sample}.windows.regions.bed.gz"
    output:
        directory(config["output_dir"] + "/qc/{sample}/plots")
    params:
        script = os.path.join(workflow.basedir, "scripts", "mosdepth_regions_to_wig.py"),
        step = config.get("chromograph_step", 5000),
        strip = "1" if config.get("chromograph_strip_chr", True) else "0",
        norm = "--norm" if config.get("chromograph_norm", True) else "",
        wig = config["output_dir"] + "/qc/{sample}/{sample}.coverage.wig"
    conda:
        "../envs/mosdepth.yaml"
    log:
        config["output_dir"] + "/logs/chromograph_{sample}.log"
    resources:
        threads = lambda wildcards, attempt: 1,
        time_hrs = lambda wildcards, attempt: attempt * 1,
        mem_gb = lambda wildcards, attempt: 4 * attempt
    message:
        "Chromograph coverage plots for sample {wildcards.sample}..."
    shell:
        """
        python {params.script} {input.bed} {params.wig} {params.step} {params.strip}
        mkdir -p {output}
        chromograph --coverage {params.wig} --step {params.step} {params.norm} --euploid --outd {output} > {log} 2>&1
        """