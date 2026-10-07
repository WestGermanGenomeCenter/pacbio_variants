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
        catalog=config["trgt_catalog_file"],

        # Optional: path to a custom STRchive-based catalog in stranger format.
        # Leave empty to use the bundled catalog (covers ~60 known disease loci).
        # database is apparently included in tool
    message:
        "Annotating TRGT repeat expansions with stranger for sample {input}..."
    shell:
        """

        stranger --trgt {input.trgt_vcf} -f {params.catalog} 2>{log} | bgzip -c >{output.annotated_vcf}
        tabix -p vcf {output.annotated_vcf} >>{log} 2>&1
        """


# more to add:
#  https://pwwang.github.io/vcfstats/cliargs/
rule colorsdb_sv:
    input:
        svs_phased="{output_dir}/variants/longphase_{sample}/{sample}_phased_SV.vcf",
        phased_cnv_and_svs="{output_dir}/variants/sawfish_phased_{sample}/{sample}_genotyped.sv.vcf.gz",
        db=config.get("colorsdb_sv_vcf", [])

    output:
        vcf1="{output_dir}/annotated_variants/colorsdb_longphase_{sample}/{sample}_longphase_colorsdb.vcf.gz",
        tbi1="{output_dir}/annotated_variants/colorsdb_longphase_{sample}/{sample}_longphase_colorsdb.vcf.gz.tbi",
        vcf2="{output_dir}/annotated_variants/colorsdb_sawfish_phased_{sample}/{sample}_sawfish_phased_colorsdb.vcf.gz",
        tbi2="{output_dir}/annotated_variants/colorsdb_sawfish_phased_{sample}/{sample}_sawfish_phased_colorsdb.vcf.gz.tbi"

    params:
        overlap=config.get("colorsdb_sv_overlap", 0.6),
        bnd=config.get("colorsdb_sv_bnd_distance", 10000)

    conda:
        "../envs/svdb.yaml"

    log:
        config["output_dir"] + "/logs/colorsdb_sv_{sample}.log"

    resources:
        threads=lambda wildcards, attempt: 1,
        time_hrs=lambda wildcards, attempt: attempt * 1,
        mem_gb=lambda wildcards, attempt: 4 * attempt

    message:
        "Adding CoLoRSdb frequencies to SVs of {wildcards.sample}..."

    shell:
        """
        svdb --query \
            --query_vcf {input.svs_phased} \
            --db {input.db} \
            --in_occ AC \
            --in_frq AF \
            --out_occ colorsdb_ac \
            --out_frq colorsdb_af \
            --overlap {params.overlap} \
            --bnd_distance {params.bnd} \
            > {output.vcf1} \
            2> {log}

        tabix -p vcf {output.vcf1} \
            >> {log} 2>&1

        svdb --query \
            --query_vcf {input.phased_cnv_and_svs} \
            --db {input.db} \
            --in_occ AC \
            --in_frq AF \
            --out_occ colorsdb_ac \
            --out_frq colorsdb_af \
            --overlap {params.overlap} \
            --bnd_distance {params.bnd} \
            > {output.vcf2} \
            2>> {log}

        tabix -p vcf {output.vcf2} \
            >> {log} 2>&1
        """

































