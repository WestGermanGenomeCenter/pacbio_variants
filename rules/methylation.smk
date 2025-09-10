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



rule pb_cpg_tools: # only able to haplotype snps, cannot use svs. for this longphase is used
    input: 
        haplotaged_bam="{output_dir}/variants/whatshap_{sample}/{sample}_haplotaged.bam",
    output:
        bed_track="{output_dir}/variants/cpg_tools_{sample}/{sample}.combined.bed.gz"
    conda:
        "../envs/pb_cpg_tools.yaml"
    log:
        "{output_dir}/logs/cpg_tools_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: attempt * 12,
        time_hrs=lambda wildcards, attempt: attempt * 2,
        mem_gb=lambda wildcards, attempt: 12 + (attempt * 12)
    params:
        prefix="{output_dir}/variants/cpg_tools_{sample}/{sample}"
    message:
        "Creating Methylation tracks for the phased whatshap output: {input.haplotaged_bam}..."
    shell:
        """
        tabix -f {input.haplotaged_bam} 2>{log}
        aligned_bam_to_cpg_scores --bam {input.haplotaged_bam} --output-prefix {params.prefix} --threads {resources.threads} >{log} 2>&1
        """
# next methbat, probably only https://github.com/PacificBiosciences/MethBat/blob/main/docs/profile_guide.md#rare-methylation-analysis