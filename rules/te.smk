import os
from snakemake.exceptions import WorkflowError

# --- which BAM sTELLeR reads: "longphase" (haplotagged + indexed), "whatshap" (haplotagged), "aligned" (plain) ---
def te_input_bam(wildcards):
    src = config.get("te_bam_source", "aligned")
    if src == "longphase":
        return f"{output_dir}/variants/longphase_{wildcards.sample}/{wildcards.sample}_aligned_haplotaged.bam"
    if src == "whatshap":
        return f"{output_dir}/variants/whatshap_{wildcards.sample}/{wildcards.sample}_haplotaged.bam"
    return f"{output_dir}/bams/{wildcards.sample}_aligned.bam"

def te_sv_vcf(wildcards):      # optional: Sniffles calls used to confirm the TE insertions
    return f"{output_dir}/variants/sniffles_{wildcards.sample}/{wildcards.sample}_svs.vcf.gz"

def te_sample_sex(wildcards):  # optional "sex" column (male/female/unknown) in samplesheet.csv
    if "sex" not in samples_df.columns:
        return "unknown"
    names = samples_df["sample"].map(get_folder_name)
    sex = dict(zip(names, samples_df["sex"].fillna("unknown").astype(str).str.lower()))
    s = sex.get(wildcards.sample, "unknown")
    return s if s in ("male", "female") else "unknown"
# ----------------------------------------------------------------------------------------------------

if config.get("run_steller"):
    for _key in ("steller_dir", "steller_te_fasta"):
        if not os.path.exists(config.get(_key, "")):
            raise WorkflowError(f"config '{_key}' is missing or not downloaded.\n"
                                "Run: bash data/get_extra_annotation_data.sh")


rule steller_te:
    input:
        bam = te_input_bam,
        reference = config["reference"],
        te_fasta = config.get("steller_te_fasta", [])
    output:
        vcf = "{output_dir}/variants/steller_{sample}/{sample}_steller_repeats.vcf"
    params:
        script = os.path.join(config.get("steller_dir", ""), "steller", "steller.py"),
        prefix = "{output_dir}/variants/steller_{sample}/{sample}_steller",
        sr = config.get("steller_min_support", 3),
        maxreads = config.get("steller_max_reads", 100),
        mq = config.get("steller_min_mapq", 20)
    conda:
        "../envs/steller.yaml"
    log:
        "{output_dir}/logs/steller_{sample}.log"
    resources:
        threads = lambda wildcards, attempt: 2,
        time_hrs = lambda wildcards, attempt: attempt * 1,
        mem_gb = lambda wildcards, attempt: 2 + (attempt * 2)
    message:
        "Detecting non-reference TE insertions with sTELLeR for {input.bam} ..."
    shell:
        """
        test -s {input.bam}.bai || samtools index -@ {resources.threads} {input.bam} >{log} 2>&1
        mkdir -p $(dirname {params.prefix})
        python {params.script} --ref {input.reference} --TE_fasta {input.te_fasta} --bam {input.bam} \
            --style pb --sr {params.sr} -mr {params.maxreads} -m {params.mq} \
            -o {params.prefix} >>{log} 2>&1
        """


rule steller_te_filter:
    """High-confidence TE set + annotated table with the reason for every pass/fail."""
    input:
        vcf = "{output_dir}/variants/steller_{sample}/{sample}_steller_repeats.vcf",
        bam = te_input_bam,
        sv = lambda wc: te_sv_vcf(wc) if config.get("te_use_sv_vcf") else [],
        exclude = config.get("te_exclude_bed") or []
    output:
        vcf = "{output_dir}/annotated_variants/steller_{sample}/{sample}_steller_te.high_confidence.vcf.gz",
        tbi = "{output_dir}/annotated_variants/steller_{sample}/{sample}_steller_te.high_confidence.vcf.gz.tbi",
        tsv = "{output_dir}/annotated_variants/steller_{sample}/{sample}_steller_te.annotated.tsv"
    params:
        script = os.path.join(workflow.basedir, "scripts", "filter_te_calls.py"),
        prefix = "{output_dir}/annotated_variants/steller_{sample}/{sample}_steller_te",
        sex = te_sample_sex,
        sv_opt = lambda wc, input: f"--sv-vcf {input.sv}" if input.sv else "",
        excl_opt = lambda wc, input: f"--exclude-bed {input.exclude}" if input.exclude else "",
        extra = config.get("te_filter_extra", "")
    conda:
        "../envs/steller.yaml"
    log:
        "{output_dir}/logs/steller_filter_{sample}.log"
    resources:
        threads = lambda wildcards, attempt: 1,
        time_hrs = lambda wildcards, attempt: attempt * 1,
        mem_gb = lambda wildcards, attempt: 2 * attempt
    message:
        "Filtering sTELLeR TE calls for {input.bam} ..."
    shell:
        """
        test -s {input.bam}.bai || samtools index {input.bam} >{log} 2>&1
        python {params.script} {input.vcf} {params.prefix} {input.bam} {params.sex} \
            {params.sv_opt} {params.excl_opt} {params.extra} 2>>{log}
        bgzip -f {params.prefix}.high_confidence.vcf >>{log} 2>&1
        tabix -f -p vcf {output.vcf} >>{log} 2>&1
        """