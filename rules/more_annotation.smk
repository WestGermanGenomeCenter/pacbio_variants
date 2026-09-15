# this is ai-slop, needs to be made useful first:

#options so far not realized: 
# echtvar: no conda install, only gnomad data - that is already part of vep
# alphagenome: no conda, gpu needed, no real usability from real users yet
# vep plugins: cadd, alphamissense, spliceai (here snvs only and not indels due to licence): nice idea, not usable on the hpc due to only-online plugin install. data is available and usable, though 
# annovar: older, also only registration-only access to main executable, so no conda aswell 
# svdb: merging svs but not snps? not for now 
#
# ---------------------------------------------------------------------------
# STRANGER — tandem repeat expansion annotation for TRGT output
# Annotates TRGT VCF with pathogenic repeat thresholds from the STRchive catalog.
# You already run TRGT — this closes the clinical interpretation gap on its output.
# Catalog bundled with stranger; runs fully offline.
# ---------------------------------------------------------------------------
rule stranger_trgt:
    input:
        trgt_vcf="{output_dir}/variants/trgt_{sample}/{sample}.vcf.gz"
    output:
        annotated_vcf="{output_dir}/annotated_variants/stranger_{sample}/{sample}_trgt_stranger_annotated.vcf.gz",
        annotated_tbi="{output_dir}/annotated_variants/stranger_{sample}/{sample}_trgt_stranger_annotated.vcf.gz.tbi"
    conda:
        "../envs/stranger.yaml"
    log:
        "{output_dir}/logs/stranger_{sample}.log"
    resources:
        threads=lambda wildcards, attempt: 2,
        time_hrs=lambda wildcards, attempt: attempt * 1,
        mem_gb=lambda wildcards, attempt: 4 + (attempt * 2)
    params:
        # Optional: path to a custom STRchive-based catalog in stranger format.
        # Leave empty to use the bundled catalog (covers ~60 known disease loci).
        # database is apparently included in tool
    message:
        "Annotating TRGT repeat expansions with stranger for sample {input}..."
    shell:
        """

        stranger --trgt {input.trgt_vcf}  2>{log} | bgzip -c >{output.annotated_vcf}
        tabix -p vcf {output.annotated_vcf} >>{log} 2>&1
        """


# more to add:
#  https://pwwang.github.io/vcfstats/cliargs/