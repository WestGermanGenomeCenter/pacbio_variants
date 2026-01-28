# Guide: Converting PacBio Workflow to Support Trios and Multi-Sample Analyses

## Overview
This guide shows how to modify the WestGermanGenomeCenter PacBio variants workflow to support:
- **Trios** (mother-father-child)
- **Quadruplets** (mother-father-child1-child2)
- **Mother-child** (affected mother - unaffected child or vice versa)
- Single samples (existing functionality preserved)

## 1. Modified Samplesheet Structure

### New `samplesheet.csv` Format

```csv
sample_name,bam_file,family_id,relation,affected_status
mother_001,mother_001.bam,family_001,mother,unaffected
father_001,father_001.bam,family_001,father,unaffected
child_001,child_001.bam,family_001,child,affected
child_002,child_002.bam,family_001,child,affected
mother_002,mother_002.bam,family_002,mother,affected
child_003,child_003.bam,family_002,child,unaffected
singleton_001,singleton_001.bam,singleton_001,proband,affected
```

**Column Definitions:**
- `sample_name`: Unique identifier for each sample
- `bam_file`: Input BAM filename (must exist in root folder)
- `family_id`: Groups related samples (use sample_name for singletons)
- `relation`: One of {mother, father, child, proband}
- `affected_status`: {affected, unaffected}

### Alternative PED Format (Optional)

For compatibility with standard trio tools, also support PED format:

```csv
family_id,sample_name,paternal_id,maternal_id,sex,phenotype,bam_file
family_001,father_001,0,0,1,1,father_001.bam
family_001,mother_001,0,0,2,1,mother_001.bam
family_001,child_001,father_001,mother_001,1,2,child_001.bam
family_001,child_002,father_001,mother_001,2,2,child_002.bam
```

## 2. Modified Config.yaml

Add new sections to `config.yaml`:

```yaml
# Existing config options...

# New: Family analysis options
family_analysis:
  enable_trio_calling: True
  enable_denovo_detection: True
  enable_inheritance_filtering: True
  enable_compound_het_detection: True

# Tool-specific trio options
use_deeptrio: True  # Use DeepTrio instead of DeepVariant for trios
use_glnexus: True   # Joint calling for families
use_whatshap_pedigree: True  # Pedigree-aware phasing

# De novo variant calling
denovo_calling:
  min_child_depth: 10
  min_parent_depth: 10
  max_parent_alt_freq: 0.05
  min_child_alt_freq: 0.25

# Inheritance patterns to check
inheritance_modes:
  - de_novo
  - autosomal_recessive
  - autosomal_dominant
  - x_linked
  - compound_heterozygous
```

## 3. New Snakemake Rules

### 3.1 Create `rules/common_family.smk`

```python
"""
Common functions for family-based analysis
"""
import pandas as pd

def load_samplesheet(samplesheet_path):
    """Load and validate samplesheet"""
    df = pd.read_csv(samplesheet_path)
    required_cols = ['sample_name', 'bam_file', 'family_id', 'relation', 'affected_status']
    
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")
    
    return df

def get_families(samplesheet_df):
    """Get list of all family IDs"""
    return samplesheet_df['family_id'].unique().tolist()

def get_family_members(samplesheet_df, family_id):
    """Get all samples in a family"""
    return samplesheet_df[samplesheet_df['family_id'] == family_id]['sample_name'].tolist()

def is_trio(samplesheet_df, family_id):
    """Check if family is a trio (has mother, father, child)"""
    family_df = samplesheet_df[samplesheet_df['family_id'] == family_id]
    relations = set(family_df['relation'].tolist())
    return {'mother', 'father', 'child'}.issubset(relations)

def is_quad(samplesheet_df, family_id):
    """Check if family is a quadruplet"""
    family_df = samplesheet_df[samplesheet_df['family_id'] == family_id]
    children = family_df[family_df['relation'] == 'child']
    return is_trio(samplesheet_df, family_id) and len(children) >= 2

def is_duo(samplesheet_df, family_id):
    """Check if family is a mother-child or father-child duo"""
    family_df = samplesheet_df[samplesheet_df['family_id'] == family_id]
    relations = set(family_df['relation'].tolist())
    return (len(family_df) == 2 and 
            'child' in relations and 
            ('mother' in relations or 'father' in relations))

def get_parent_samples(samplesheet_df, family_id):
    """Get parent sample names for a family"""
    family_df = samplesheet_df[samplesheet_df['family_id'] == family_id]
    parents = family_df[family_df['relation'].isin(['mother', 'father'])]
    return parents['sample_name'].tolist()

def get_child_samples(samplesheet_df, family_id):
    """Get child sample names for a family"""
    family_df = samplesheet_df[samplesheet_df['family_id'] == family_id]
    children = family_df[family_df['relation'] == 'child']
    return children['sample_name'].tolist()

def get_sample_role(samplesheet_df, sample_name):
    """Get role of a sample (mother/father/child/proband)"""
    row = samplesheet_df[samplesheet_df['sample_name'] == sample_name]
    return row['relation'].values[0] if len(row) > 0 else 'unknown'

def create_ped_file(samplesheet_df, family_id, output_path):
    """Create PED file for a family"""
    family_df = samplesheet_df[samplesheet_df['family_id'] == family_id]
    
    with open(output_path, 'w') as f:
        for _, row in family_df.iterrows():
            # PED format: FamilyID IndividualID PaternalID MaternalID Sex Phenotype
            father_id = '0'
            mother_id = '0'
            
            if row['relation'] == 'child':
                parents = family_df[family_df['relation'].isin(['mother', 'father'])]
                for _, parent in parents.iterrows():
                    if parent['relation'] == 'father':
                        father_id = parent['sample_name']
                    elif parent['relation'] == 'mother':
                        mother_id = parent['sample_name']
            
            sex = '0'  # Unknown - would need to be added to samplesheet
            phenotype = '2' if row['affected_status'] == 'affected' else '1'
            
            f.write(f"{family_id}\t{row['sample_name']}\t{father_id}\t{mother_id}\t{sex}\t{phenotype}\n")
```

### 3.2 Create `rules/deeptrio.smk`

```python
"""
DeepTrio variant calling for trios
"""

rule deeptrio_make_examples:
    input:
        bam_child=lambda wildcards: get_bam_for_sample(wildcards.child),
        bam_parent1=lambda wildcards: get_bam_for_sample(get_parent_samples(SAMPLESHEET, wildcards.family_id)[0]),
        bam_parent2=lambda wildcards: get_bam_for_sample(get_parent_samples(SAMPLESHEET, wildcards.family_id)[1]),
        bai_child=lambda wildcards: get_bam_for_sample(wildcards.child) + ".bai",
        bai_parent1=lambda wildcards: get_bam_for_sample(get_parent_samples(SAMPLESHEET, wildcards.family_id)[0]) + ".bai",
        bai_parent2=lambda wildcards: get_bam_for_sample(get_parent_samples(SAMPLESHEET, wildcards.family_id)[1]) + ".bai",
        ref=config["reference"],
        fai=config["reference"] + ".fai"
    output:
        examples_child=temp(directory("deeptrio_temp/{family_id}/{child}/make_examples")),
        examples_parent1=temp(directory("deeptrio_temp/{family_id}/{child}/make_examples_parent1")),
        examples_parent2=temp(directory("deeptrio_temp/{family_id}/{child}/make_examples_parent2"))
    params:
        sample_child="{child}",
        sample_parent1=lambda wildcards: get_parent_samples(SAMPLESHEET, wildcards.family_id)[0],
        sample_parent2=lambda wildcards: get_parent_samples(SAMPLESHEET, wildcards.family_id)[1],
        extra=config.get("deeptrio_extra", "")
    threads: 16
    resources:
        mem_mb=64000,
        runtime=1440
    log:
        "logs/deeptrio/{family_id}/{child}_make_examples.log"
    conda:
        "../envs/deeptrio.yaml"
    shell:
        """
        deeptrio make_examples \
            --mode calling \
            --ref {input.ref} \
            --reads_parent1 {input.bam_parent1} \
            --reads_parent2 {input.bam_parent2} \
            --reads {input.bam_child} \
            --examples_parent1 {output.examples_parent1} \
            --examples_parent2 {output.examples_parent2} \
            --examples {output.examples_child} \
            --sample_name_parent1 {params.sample_parent1} \
            --sample_name_parent2 {params.sample_parent2} \
            --sample_name {params.sample_child} \
            --num_shards {threads} \
            {params.extra} \
            2>&1 | tee {log}
        """

rule deeptrio_call_variants:
    input:
        examples_child="deeptrio_temp/{family_id}/{child}/make_examples",
        examples_parent1="deeptrio_temp/{family_id}/{child}/make_examples_parent1",
        examples_parent2="deeptrio_temp/{family_id}/{child}/make_examples_parent2"
    output:
        tfrecords_child=temp(directory("deeptrio_temp/{family_id}/{child}/call_variants")),
        tfrecords_parent1=temp(directory("deeptrio_temp/{family_id}/{child}/call_variants_parent1")),
        tfrecords_parent2=temp(directory("deeptrio_temp/{family_id}/{child}/call_variants_parent2"))
    params:
        model_parent="gs://deepvariant/models/DeepTrio/1.6.0/DeepTrio-inception_v3-1.6.0+data-wgs_parent",
        model_child="gs://deepvariant/models/DeepTrio/1.6.0/DeepTrio-inception_v3-1.6.0+data-wgs_child"
    threads: 4
    resources:
        mem_mb=16000,
        runtime=720,
        gpu=1
    log:
        "logs/deeptrio/{family_id}/{child}_call_variants.log"
    conda:
        "../envs/deeptrio.yaml"
    shell:
        """
        # Call variants for parents
        deeptrio call_variants \
            --examples {input.examples_parent1} \
            --checkpoint {params.model_parent} \
            --outfile {output.tfrecords_parent1} &
        
        deeptrio call_variants \
            --examples {input.examples_parent2} \
            --checkpoint {params.model_parent} \
            --outfile {output.tfrecords_parent2} &
        
        # Call variants for child
        deeptrio call_variants \
            --examples {input.examples_child} \
            --checkpoint {params.model_child} \
            --outfile {output.tfrecords_child}
        
        wait
        """

rule deeptrio_postprocess:
    input:
        tfrecords_child="deeptrio_temp/{family_id}/{child}/call_variants",
        tfrecords_parent1="deeptrio_temp/{family_id}/{child}/call_variants_parent1",
        tfrecords_parent2="deeptrio_temp/{family_id}/{child}/call_variants_parent2",
        ref=config["reference"]
    output:
        vcf_child="variants/deeptrio_{family_id}/{child}.vcf.gz",
        vcf_parent1="variants/deeptrio_{family_id}/{child}_parent1.vcf.gz",
        vcf_parent2="variants/deeptrio_{family_id}/{child}_parent2.vcf.gz",
        gvcf_child="variants/deeptrio_{family_id}/{child}.g.vcf.gz",
        gvcf_parent1="variants/deeptrio_{family_id}/{child}_parent1.g.vcf.gz",
        gvcf_parent2="variants/deeptrio_{family_id}/{child}_parent2.g.vcf.gz"
    params:
        sample_child="{child}",
        sample_parent1=lambda wildcards: get_parent_samples(SAMPLESHEET, wildcards.family_id)[0],
        sample_parent2=lambda wildcards: get_parent_samples(SAMPLESHEET, wildcards.family_id)[1]
    threads: 1
    resources:
        mem_mb=8000,
        runtime=120
    log:
        "logs/deeptrio/{family_id}/{child}_postprocess.log"
    conda:
        "../envs/deeptrio.yaml"
    shell:
        """
        deeptrio postprocess_variants \
            --ref {input.ref} \
            --infile {input.tfrecords_child} \
            --outfile {output.vcf_child} \
            --sample_name {params.sample_child} \
            --gvcf_outfile {output.gvcf_child} &
        
        deeptrio postprocess_variants \
            --ref {input.ref} \
            --infile {input.tfrecords_parent1} \
            --outfile {output.vcf_parent1} \
            --sample_name {params.sample_parent1} \
            --gvcf_outfile {output.gvcf_parent1} &
        
        deeptrio postprocess_variants \
            --ref {input.ref} \
            --infile {input.tfrecords_parent2} \
            --outfile {output.vcf_parent2} \
            --sample_name {params.sample_parent2} \
            --gvcf_outfile {output.gvcf_parent2}
        
        wait
        """
```

### 3.3 Create `rules/glnexus.smk`

```python
"""
GLnexus joint calling for families
"""

rule glnexus_family:
    input:
        gvcfs=lambda wildcards: expand(
            "variants/deeptrio_{family_id}/{sample}.g.vcf.gz",
            family_id=wildcards.family_id,
            sample=get_family_members(SAMPLESHEET, wildcards.family_id)
        )
    output:
        bcf="variants/glnexus_{family_id}/joint_calls.bcf",
        vcf="variants/glnexus_{family_id}/joint_calls.vcf.gz"
    params:
        config_preset="DeepVariant_unfiltered",
        extra=""
    threads: 16
    resources:
        mem_mb=128000,
        runtime=480
    log:
        "logs/glnexus/{family_id}_joint_calling.log"
    conda:
        "../envs/glnexus.yaml"
    shell:
        """
        glnexus_cli \
            --config {params.config_preset} \
            --threads {threads} \
            --mem-gbytes {resources.mem_mb}/1024 \
            {input.gvcfs} \
            > {output.bcf} \
            2> {log}
        
        bcftools view {output.bcf} | \
            bcftools sort -Oz -o {output.vcf}
        
        bcftools index -t {output.vcf}
        """
```

### 3.4 Create `rules/denovo_detection.smk`

```python
"""
De novo variant detection in trios
"""

rule denovo_filter:
    input:
        vcf="variants/glnexus_{family_id}/joint_calls.vcf.gz",
        ped="pedigrees/{family_id}.ped"
    output:
        vcf="variants/denovo_{family_id}/{child}_denovo.vcf.gz",
        tsv="variants/denovo_{family_id}/{child}_denovo.tsv"
    params:
        min_child_dp=config.get("denovo_calling", {}).get("min_child_depth", 10),
        min_parent_dp=config.get("denovo_calling", {}).get("min_parent_depth", 10),
        max_parent_af=config.get("denovo_calling", {}).get("max_parent_alt_freq", 0.05),
        min_child_af=config.get("denovo_calling", {}).get("min_child_alt_freq", 0.25),
        child_sample="{child}"
    threads: 2
    resources:
        mem_mb=8000,
        runtime=60
    log:
        "logs/denovo/{family_id}/{child}_filter.log"
    conda:
        "../envs/denovo.yaml"
    script:
        "../scripts/filter_denovo.py"
```

### 3.5 Create `rules/whatshap_pedigree.smk`

```python
"""
Pedigree-aware phasing with WhatsHap
"""

rule whatshap_phase_pedigree:
    input:
        vcf="variants/glnexus_{family_id}/joint_calls.vcf.gz",
        bams=lambda wildcards: expand(
            "{sample}.bam",
            sample=get_family_members(SAMPLESHEET, wildcards.family_id)
        ),
        ped="pedigrees/{family_id}.ped",
        ref=config["reference"]
    output:
        vcf="variants/whatshap_pedigree_{family_id}/phased.vcf.gz",
        stats="variants/whatshap_pedigree_{family_id}/phasing_stats.txt"
    params:
        extra="--indels --distrust-genotypes"
    threads: 8
    resources:
        mem_mb=32000,
        runtime=480
    log:
        "logs/whatshap_pedigree/{family_id}_phase.log"
    conda:
        "../envs/whatshap.yaml"
    shell:
        """
        whatshap phase \
            --reference {input.ref} \
            --ped {input.ped} \
            --output {output.vcf} \
            --output-read-list {output.stats} \
            {params.extra} \
            {input.vcf} \
            {input.bams} \
            2>&1 | tee {log}
        
        tabix -p vcf {output.vcf}
        """

rule whatshap_stats_pedigree:
    input:
        vcf="variants/whatshap_pedigree_{family_id}/phased.vcf.gz"
    output:
        tsv="variants/whatshap_pedigree_{family_id}/stats.tsv",
        blocklist="variants/whatshap_pedigree_{family_id}/blocks.txt"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30
    log:
        "logs/whatshap_pedigree/{family_id}_stats.log"
    conda:
        "../envs/whatshap.yaml"
    shell:
        """
        whatshap stats \
            --tsv {output.tsv} \
            --block-list {output.blocklist} \
            {input.vcf} \
            2>&1 | tee {log}
        """
```

### 3.6 Create `rules/inheritance_filtering.smk`

```python
"""
Inheritance pattern filtering
"""

rule filter_recessive:
    input:
        vcf="variants/whatshap_pedigree_{family_id}/phased.vcf.gz",
        ped="pedigrees/{family_id}.ped"
    output:
        vcf="variants/inheritance_{family_id}/{child}_autosomal_recessive.vcf.gz"
    params:
        affected_sample="{child}"
    threads: 2
    resources:
        mem_mb=8000,
        runtime=60
    log:
        "logs/inheritance/{family_id}/{child}_recessive.log"
    conda:
        "../envs/inheritance.yaml"
    script:
        "../scripts/filter_recessive.py"

rule filter_dominant:
    input:
        vcf="variants/whatshap_pedigree_{family_id}/phased.vcf.gz",
        ped="pedigrees/{family_id}.ped"
    output:
        vcf="variants/inheritance_{family_id}/{child}_autosomal_dominant.vcf.gz"
    params:
        affected_sample="{child}"
    threads: 2
    resources:
        mem_mb=8000,
        runtime=60
    log:
        "logs/inheritance/{family_id}/{child}_dominant.log"
    conda:
        "../envs/inheritance.yaml"
    script:
        "../scripts/filter_dominant.py"

rule find_compound_hets:
    input:
        vcf="variants/whatshap_pedigree_{family_id}/phased.vcf.gz",
        ped="pedigrees/{family_id}.ped"
    output:
        vcf="variants/inheritance_{family_id}/{child}_compound_het.vcf.gz",
        tsv="variants/inheritance_{family_id}/{child}_compound_het_pairs.tsv"
    params:
        affected_sample="{child}",
        gene_annotations=config.get("gene_annotations", "")
    threads: 4
    resources:
        mem_mb=16000,
        runtime=120
    log:
        "logs/inheritance/{family_id}/{child}_compound_het.log"
    conda:
        "../envs/inheritance.yaml"
    script:
        "../scripts/find_compound_hets.py"
```

## 4. Supporting Python Scripts

### 4.1 `scripts/filter_denovo.py`

```python
#!/usr/bin/env python3
"""
Filter for de novo variants in trios
"""
import sys
import cyvcf2
from cyvcf2 import VCF, Writer

def get_parent_indices(vcf, ped_file):
    """Parse PED file to find parent sample indices"""
    # Implementation to parse PED and map to VCF samples
    pass

def is_denovo(variant, child_idx, parent_indices, params):
    """Check if variant is de novo"""
    child_gt = variant.genotypes[child_idx]
    
    # Child must be het or hom_alt
    if child_gt[0] == 0 and child_gt[1] == 0:
        return False
    
    # Check child DP and AF
    child_dp = variant.format('DP')[child_idx][0]
    if child_dp < params['min_child_dp']:
        return False
    
    child_ad = variant.format('AD')[child_idx]
    if len(child_ad) < 2:
        return False
    child_af = child_ad[1] / sum(child_ad) if sum(child_ad) > 0 else 0
    if child_af < params['min_child_af']:
        return False
    
    # Check both parents
    for parent_idx in parent_indices:
        parent_gt = variant.genotypes[parent_idx]
        parent_dp = variant.format('DP')[parent_idx][0]
        
        if parent_dp < params['min_parent_dp']:
            continue
        
        # Parent should be ref/ref or very low AF
        parent_ad = variant.format('AD')[parent_idx]
        if len(parent_ad) < 2:
            continue
        parent_af = parent_ad[1] / sum(parent_ad) if sum(parent_ad) > 0 else 0
        
        if parent_af > params['max_parent_af']:
            return False
    
    return True

def main(snakemake):
    vcf = VCF(snakemake.input.vcf)
    ped_file = snakemake.input.ped
    child_sample = snakemake.params.child_sample
    
    # Get indices
    child_idx = vcf.samples.index(child_sample)
    parent_indices = get_parent_indices(vcf, ped_file)
    
    params = {
        'min_child_dp': snakemake.params.min_child_dp,
        'min_parent_dp': snakemake.params.min_parent_dp,
        'max_parent_af': snakemake.params.max_parent_af,
        'min_child_af': snakemake.params.min_child_af
    }
    
    # Write output
    w = Writer(snakemake.output.vcf, vcf)
    tsv_out = open(snakemake.output.tsv, 'w')
    tsv_out.write("CHROM\tPOS\tREF\tALT\tGENE\tCHILD_GT\tCHILD_DP\tCHILD_AF\n")
    
    for variant in vcf:
        if is_denovo(variant, child_idx, parent_indices, params):
            w.write_record(variant)
            
            # Write to TSV
            child_ad = variant.format('AD')[child_idx]
            child_af = child_ad[1] / sum(child_ad) if sum(child_ad) > 0 else 0
            tsv_out.write(f"{variant.CHROM}\t{variant.POS}\t{variant.REF}\t{variant.ALT[0]}\t")
            tsv_out.write(f".\t{variant.genotypes[child_idx][:2]}\t")
            tsv_out.write(f"{variant.format('DP')[child_idx][0]}\t{child_af:.3f}\n")
    
    w.close()
    tsv_out.close()
    vcf.close()

if __name__ == '__main__':
    main(snakemake)
```

## 5. Modified Main Snakefile

Update the main `Snakefile` to include family logic:

```python
import pandas as pd
from pathlib import Path

# Load configuration
configfile: "config.yaml"

# Load samplesheet
SAMPLESHEET = pd.read_csv(config["samplesheet"])

# Include common family functions
include: "rules/common_family.smk"

# Get families and samples
FAMILIES = get_families(SAMPLESHEET)
SAMPLES = SAMPLESHEET['sample_name'].tolist()

# Categorize families
TRIOS = [f for f in FAMILIES if is_trio(SAMPLESHEET, f)]
QUADS = [f for f in FAMILIES if is_quad(SAMPLESHEET, f)]
DUOS = [f for f in FAMILIES if is_duo(SAMPLESHEET, f)]
SINGLETONS = [f for f in FAMILIES if f not in TRIOS + DUOS]

# Include existing rules
include: "rules/alignment.smk"
include: "rules/deepvariant.smk"  # For singletons
include: "rules/sniffles.smk"
include: "rules/trgt.smk"
# ... other existing tools

# Include new family-based rules
if config.get("use_deeptrio", False):
    include: "rules/deeptrio.smk"

if config.get("use_glnexus", False):
    include: "rules/glnexus.smk"

if config.get("family_analysis", {}).get("enable_denovo_detection", False):
    include: "rules/denovo_detection.smk"

if config.get("use_whatshap_pedigree", False):
    include: "rules/whatshap_pedigree.smk"

if config.get("family_analysis", {}).get("enable_inheritance_filtering", False):
    include: "rules/inheritance_filtering.smk"

# Define target outputs
rule all:
    input:
        # Singleton outputs
        expand("variants/deepvariant_{sample}/{sample}.vcf.gz", 
               sample=[s for s in SAMPLES if s in SINGLETONS]),
        
        # Trio outputs (if enabled)
        expand("variants/deeptrio_{family}/{child}.vcf.gz",
               family=TRIOS,
               child=lambda wildcards: get_child_samples(SAMPLESHEET, wildcards.family))
               if config.get("use_deeptrio", False) else [],
        
        # Joint calling outputs
        expand("variants/glnexus_{family}/joint_calls.vcf.gz",
               family=TRIOS + QUADS + DUOS)
               if config.get("use_glnexus", False) else [],
        
        # De novo outputs
        expand("variants/denovo_{family}/{child}_denovo.vcf.gz",
               family=TRIOS + QUADS,
               child=lambda wildcards: get_child_samples(SAMPLESHEET, wildcards.family))
               if config.get("family_analysis", {}).get("enable_denovo_detection", False) else [],
        
        # Existing outputs for all samples
        expand("variants/sniffles_{sample}/{sample}.vcf.gz", sample=SAMPLES),
        # ... other tools
```

## 6. Conda Environment Files

### `envs/deeptrio.yaml`

```yaml
name: deeptrio
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - python=3.9
  - google-deepvariant=1.6.0
  - bcftools=1.18
  - tabix=1.18
```

### `envs/glnexus.yaml`

```yaml
name: glnexus
channels:
  - bioconda
  - conda-forge
dependencies:
  - glnexus=1.4.1
  - bcftools=1.18
  - tabix=1.18
```

## 7. Tool Compatibility Matrix

| Tool | Trio Support | Implementation |
|------|-------------|----------------|
| **DeepVariant** | ❌ (use DeepTrio) | Single samples only |
| **DeepTrio** | ✅ Native | Joint calling with pedigree info |
| **GLnexus** | ✅ Joint calling | Merges gVCFs from all family members |
| **Sniffles2** | ⚠️ Partial | Run per-sample, then merge with SURVIVOR |
| **TRGT** | ⚠️ Partial | Run per-sample, manual inheritance checks |
| **HiFiCNV** | ❌ | Single sample CNV calling |
| **Paraphase** | ❌ | Single sample gene calling |
| **WhatsHap** | ✅ Pedigree mode | Uses pedigree info for better phasing |
| **LongPhase** | ⚠️ Limited | Can use multiple BAMs but no pedigree |
| **NanoCaller** | ❌ (use DeepTrio) | Single samples only |
| **BCFtools** | ⚠️ Manual | Can call on all samples, manual filtering |

## 8. Additional Modifications Needed

### 8.1 Update Sniffles for Family Merging

```python
# Add to rules/sniffles.smk

rule sniffles_merge_family:
    input:
        vcfs=lambda wildcards: expand(
            "variants/sniffles_{sample}/{sample}.vcf.gz",
            sample=get_family_members(SAMPLESHEET, wildcards.family_id)
        )
    output:
        vcf="variants/sniffles_merged_{family_id}/merged.vcf.gz"
    params:
        min_support=2
    threads: 4
    resources:
        mem_mb=16000,
        runtime=120
    log:
        "logs/sniffles/{family_id}_merge.log"
    conda:
        "../envs/sniffles.yaml"
    shell:
        """
        sniffles \
            --input {input.vcfs} \
            --vcf {output.vcf} \
            --threads {threads} \
            --min-support {params.min_support} \
            2>&1 | tee {log}
        
        bcftools index -t {output.vcf}
        """
```

### 8.2 Create PED Files

```python
# Add to main Snakefile or rules/common_family.smk

rule create_ped_files:
    input:
        samplesheet=config["samplesheet"]
    output:
        ped="pedigrees/{family_id}.ped"
    params:
        family_id="{family_id}"
    run:
        create_ped_file(SAMPLESHEET, params.family_id, output.ped)
```

### 8.3 Modified Target Rule with All Outputs

```python
# Complete rule all with all possible outputs

def get_all_outputs(config, samplesheet):
    """Generate all output files based on config and samplesheet"""
    outputs = []
    
    families = get_families(samplesheet)
    samples = samplesheet['sample_name'].tolist()
    trios = [f for f in families if is_trio(samplesheet, f)]
    quads = [f for f in families if is_quad(samplesheet, f)]
    duos = [f for f in families if is_duo(samplesheet, f)]
    singletons = [s for f in families for s in get_family_members(samplesheet, f) 
                  if f not in trios + duos]
    
    # Single sample variant calling
    for sample in singletons:
        outputs.append(f"variants/deepvariant_{sample}/{sample}.vcf.gz")
    
    # Trio-based calling
    if config.get("use_deeptrio", False):
        for family in trios + quads:
            for child in get_child_samples(samplesheet, family):
                outputs.append(f"variants/deeptrio_{family}/{child}.vcf.gz")
    
    # Joint calling
    if config.get("use_glnexus", False):
        for family in trios + quads + duos:
            outputs.append(f"variants/glnexus_{family}/joint_calls.vcf.gz")
    
    # De novo detection
    if config.get("family_analysis", {}).get("enable_denovo_detection", False):
        for family in trios + quads:
            for child in get_child_samples(samplesheet, family):
                outputs.append(f"variants/denovo_{family}/{child}_denovo.vcf.gz")
    
    # Pedigree phasing
    if config.get("use_whatshap_pedigree", False):
        for family in trios + quads + duos:
            outputs.append(f"variants/whatshap_pedigree_{family}/phased.vcf.gz")
    
    # Inheritance filtering
    if config.get("family_analysis", {}).get("enable_inheritance_filtering", False):
        for family in trios + quads:
            for child in get_child_samples(samplesheet, family):
                outputs.append(f"variants/inheritance_{family}/{child}_autosomal_recessive.vcf.gz")
                outputs.append(f"variants/inheritance_{family}/{child}_autosomal_dominant.vcf.gz")
                outputs.append(f"variants/inheritance_{family}/{child}_compound_het.vcf.gz")
    
    # SV calling for all samples
    for sample in samples:
        outputs.append(f"variants/sniffles_{sample}/{sample}.vcf.gz")
    
    # SV merging for families
    for family in families:
        if len(get_family_members(samplesheet, family)) > 1:
            outputs.append(f"variants/sniffles_merged_{family}/merged.vcf.gz")
    
    # TR calling for all samples
    for sample in samples:
        outputs.append(f"variants/trgt_{sample}/{sample}.vcf.gz")
    
    # CNV calling for all samples
    for sample in samples:
        outputs.append(f"variants/hificnv_{sample}/{sample}.vcf.gz")
    
    # PED files for families
    for family in trios + quads + duos:
        outputs.append(f"pedigrees/{family}.ped")
    
    return outputs

rule all:
    input:
        get_all_outputs(config, SAMPLESHEET)
```

## 9. Usage Examples

### Example 1: Trio Analysis

**Samplesheet:**
```csv
sample_name,bam_file,family_id,relation,affected_status
father_001,father.bam,trio_001,father,unaffected
mother_001,mother.bam,trio_001,mother,unaffected
child_001,child.bam,trio_001,child,affected
```

**Config:**
```yaml
use_deeptrio: True
use_glnexus: True
use_whatshap_pedigree: True
family_analysis:
  enable_denovo_detection: True
  enable_inheritance_filtering: True
```

**Outputs:**
- De novo variants: `variants/denovo_trio_001/child_001_denovo.vcf.gz`
- Recessive candidates: `variants/inheritance_trio_001/child_001_autosomal_recessive.vcf.gz`
- Phased variants: `variants/whatshap_pedigree_trio_001/phased.vcf.gz`

### Example 2: Mother-Child Duo (Affected Mother)

**Samplesheet:**
```csv
sample_name,bam_file,family_id,relation,affected_status
mother_002,mother.bam,duo_001,mother,affected
child_002,child.bam,duo_001,child,unaffected
```

**Analysis:** Look for variants present in mother but absent in child (protective variants in child or pathogenic in mother)

### Example 3: Quadruplet

**Samplesheet:**
```csv
sample_name,bam_file,family_id,relation,affected_status
father_003,father.bam,quad_001,father,unaffected
mother_003,mother.bam,quad_001,mother,unaffected
child_003a,child1.bam,quad_001,child,affected
child_003b,child2.bam,quad_001,child,affected
```

**Analysis:** De novo variants shared by both affected children

## 10. Testing Strategy

1. **Test with singletons first** - ensure backward compatibility
2. **Add one trio** - validate trio calling pipeline
3. **Add duo** - verify duo-specific filtering
4. **Add quadruplet** - test multi-child logic
5. **Mix all types** - validate complete workflow

## 11. Summary of Key Changes

### Files to Create:
1. `rules/common_family.smk` - Family utility functions
2. `rules/deeptrio.smk` - DeepTrio variant calling
3. `rules/glnexus.smk` - Joint calling
4. `rules/denovo_detection.smk` - De novo filtering
5. `rules/whatshap_pedigree.smk` - Pedigree phasing
6. `rules/inheritance_filtering.smk` - Inheritance pattern filters
7. `scripts/filter_denovo.py` - De novo filtering script
8. `scripts/filter_recessive.py` - Recessive filtering
9. `scripts/filter_dominant.py` - Dominant filtering
10. `scripts/find_compound_hets.py` - Compound het detection
11. `envs/deeptrio.yaml` - DeepTrio environment
12. `envs/glnexus.yaml` - GLnexus environment

### Files to Modify:
1. `Snakefile` - Add family logic, new includes, modified rule all
2. `config.yaml` - Add family analysis options
3. `samplesheet.csv` - Add family_id, relation, affected_status columns
4. `rules/sniffles.smk` - Add family merging rule

### Backward Compatibility:
- All existing single-sample workflows continue to work
- Family-based analysis only triggers when family_id groups samples
- Tools without family support continue running per-sample

This design maintains flexibility while adding powerful family-based analysis capabilities!