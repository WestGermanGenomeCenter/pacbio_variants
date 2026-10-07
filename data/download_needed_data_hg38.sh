```bash
#!/usr/bin/env bash
# =============================================================================
# Download / build reference data for the pb_variants pipeline
#
# DEFAULT:
#   All data are placed in the CURRENT DIRECTORY:
#
#       ./gnomad.v4.1.sv.sites.vcf.gz
#       ./small_variants/00-All.vcf
#       ./vep_hg38/
#       ./colorsdb/
#       ...
#
# EXAMPLE:
#
#   bash download_pb_variants_data.sh
#
# To use the actual pb_variants data directory:
#
#   DATA_DIR=/gpfs/project/projects/bmfz_gtl/software/pb_variants/data \
#       bash download_pb_variants_data.sh
#
# The destination is configurable:
#
#   DATA_DIR=/some/other/path bash download_pb_variants_data.sh
#
# IMPORTANT:
#   - All automatically downloaded resources are intended for hg38/GRCh38.
#   - Do NOT replace them with hg19/GRCh37 resources.
#   - Existing files are not downloaded again.
#   - Downloads use temporary files and are only moved into place after
#     successful completion.
#   - Large resources can require substantial disk space.
#
# Resources NOT downloaded automatically:
#
#   barcodes_lima_fasta
#       Lab/run-specific PacBio Lima barcode FASTA.
#
#   kraken2_db_folder
#       Existing cluster-wide Kraken2 database.
#
#   deepvariant_checkpoint_dir
#       Not currently required by the pipeline.
#
#   repeat_catalog_v1...
#       Historical TRGT v1 catalog. The exact upstream location needs to be
#       established before substituting another catalog.
#
# =============================================================================

set -euo pipefail


###############################################################################
# CONFIGURATION
###############################################################################

# Default destination is the current directory.
#
# Override with:
#
#   DATA_DIR=/gpfs/.../data bash download_pb_variants_data.sh
#
DATA_DIR="${DATA_DIR:-.}"

# AnnotSV version expected by the existing pipeline path.
#
# This is intentionally separate from the conda environment.
ANNOTSV_VERSION="${ANNOTSV_VERSION:-3.5.3}"

# VEP release.
#
# IMPORTANT:
# Keep this pinned for reproducibility. Change deliberately if the pipeline
# has been validated against another Ensembl/VEP release.
VEP_RELEASE="${VEP_RELEASE:-112}"


###############################################################################
# NORMALISE DATA_DIR
###############################################################################

# Convert a relative DATA_DIR to an absolute path.
#
# This is useful because the script changes directories during installation.
DATA_DIR="$(cd "$DATA_DIR" 2>/dev/null && pwd || {
    mkdir -p "$DATA_DIR"
    cd "$DATA_DIR"
    pwd
})"

echo
echo "======================================================================"
echo "pb_variants reference-data setup"
echo "======================================================================"
echo
echo "DATA_DIR:"
echo "  $DATA_DIR"
echo


###############################################################################
# REQUIREMENTS
###############################################################################

require_command()
{
    local cmd="$1"

    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo
        echo "ERROR: required command not found: $cmd"
        echo
        exit 1
    fi
}


# wget/curl are interchangeable for most downloads.
if ! command -v wget >/dev/null 2>&1 && \
   ! command -v curl >/dev/null 2>&1; then
    echo "ERROR: either wget or curl is required."
    exit 1
fi


###############################################################################
# DOWNLOAD HELPER
###############################################################################
#
# Usage:
#
#   fetch URL OUTPUT
#
# The download first goes to:
#
#   OUTPUT.tmp
#
# and is moved into place only after a successful download.
#
# This prevents a partially downloaded file from being mistaken for a
# completed reference resource on the next run.
#
###############################################################################

fetch()
{
    local url="$1"
    local out="$2"
    local tmp="${out}.tmp"

    if [[ -s "$out" ]]; then
        echo "[SKIP] $out already exists"
        return 0
    fi

    mkdir -p "$(dirname "$out")"

    echo
    echo "[GET ] $url"
    echo "[TO  ] $out"

    rm -f "$tmp"

    if command -v wget >/dev/null 2>&1; then
        wget \
            --continue \
            --tries=3 \
            -O "$tmp" \
            "$url"
    else
        curl \
            -L \
            --fail \
            --retry 3 \
            -o "$tmp" \
            "$url"
    fi

    if [[ ! -s "$tmp" ]]; then
        echo "ERROR: download produced an empty file:"
        echo "       $tmp"
        rm -f "$tmp"
        exit 1
    fi

    mv "$tmp" "$out"

    echo "[ OK ] $out"
}


###############################################################################
# 1. gnomAD SV v4.1
###############################################################################
#
# Pipeline config:
#
#   gnomad_svdb:
#     .../gnomad.v4.1.sv.sites.vcf.gz
#
# gnomAD provides a bgzipped VCF. The pipeline queries it by genomic region,
# so a tabix index is required.
#
###############################################################################

echo
echo "======================================================================"
echo "1/10  gnomAD SV v4.1"
echo "======================================================================"

GNOMAD_SV="$DATA_DIR/gnomad.v4.1.sv.sites.vcf.gz"

fetch \
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz" \
    "$GNOMAD_SV"


# Prefer the published index if available.
#
# If it cannot be downloaded, create the index locally with tabix.
if [[ ! -s "${GNOMAD_SV}.tbi" ]]; then

    echo "[IDX ] Downloading gnomAD tabix index"

    if ! fetch \
        "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz.tbi" \
        "${GNOMAD_SV}.tbi"
    then
        echo "[IDX ] Published index unavailable."
        echo "[IDX ] Creating index with tabix."

        require_command tabix
        tabix -p vcf "$GNOMAD_SV"
    fi
fi


###############################################################################
# 2. AnnotSV 3.5.3
###############################################################################
#
# Pipeline config:
#
#   annotsv_data_dir:
#     .../annotsv_data/annotsv_353/AnnotSV
#
# AnnotSV is deliberately installed independently from conda.
#
# The official AnnotSV installation separates:
#
#   make install
#
# from:
#
#   make install-human-annotation
#
# The latter downloads the human annotation database.
#
# We install into:
#
#   DATA_DIR/annotsv_data/annotsv_353/AnnotSV
#
###############################################################################

echo
echo "======================================================================"
echo "2/10  AnnotSV ${ANNOTSV_VERSION}"
echo "======================================================================"

ANNOTSV_ROOT="$DATA_DIR/annotsv_data/annotsv_353"
ANNOTSV_DIR="$ANNOTSV_ROOT/AnnotSV"

if [[ -x "$ANNOTSV_DIR/bin/AnnotSV" ]]; then

    echo "[SKIP] AnnotSV already installed:"
    echo "       $ANNOTSV_DIR"

else

    require_command git
    require_command make

    mkdir -p "$ANNOTSV_ROOT"

    ANNOTSV_TMP="$ANNOTSV_ROOT/AnnotSV-source.tmp"

    rm -rf "$ANNOTSV_TMP"

    echo "[GET ] AnnotSV ${ANNOTSV_VERSION}"

    git clone \
        --branch "${ANNOTSV_VERSION}" \
        --depth 1 \
        "https://github.com/lgmgeo/AnnotSV.git" \
        "$ANNOTSV_TMP"

    cd "$ANNOTSV_TMP"

    echo
    echo "[INST] Installing AnnotSV"
    echo

    make \
        PREFIX="$ANNOTSV_DIR" \
        install

    echo
    echo "[INST] Installing AnnotSV human annotations"
    echo

    make \
        PREFIX="$ANNOTSV_DIR" \
        install-human-annotation

    cd "$DATA_DIR"

    rm -rf "$ANNOTSV_TMP"

    echo
    echo "[ OK ] AnnotSV installed:"
    echo "       $ANNOTSV_DIR"
fi


###############################################################################
# 3. TRGT repeat catalog
###############################################################################
#
# Pipeline config expects EXACTLY:
#
#   repeat_catalog_v1.hg38.1_to_1000bp_motifs.TRGT.bed
#
# IMPORTANT:
#
# The current PacBio TRGT repository now contains newer GRCh38 repeat
# catalogs. The historical "v1 / 1-to-1000bp motifs" file should therefore
# NOT be silently replaced by a current catalog.
#
# We check whether it already exists.
#
###############################################################################

echo
echo "======================================================================"
echo "3/10  TRGT repeat catalog"
echo "======================================================================"

TRGT_REPEAT="$DATA_DIR/repeat_catalog_v1.hg38.1_to_1000bp_motifs.TRGT.bed"

if [[ -s "$TRGT_REPEAT" ]]; then

    echo "[ OK ] Existing TRGT catalog:"
    echo "       $TRGT_REPEAT"

else

    echo
    echo "[WARN] Historical TRGT v1 catalog is missing:"
    echo
    echo "       $TRGT_REPEAT"
    echo
    echo "The pipeline explicitly requests this historical catalog."
    echo "The current PacBio TRGT repository provides newer GRCh38 BED"
    echo "catalogs, so this script will NOT silently substitute one."
    echo
    echo "Obtain the exact v1 catalog and place it at:"
    echo
    echo "       $TRGT_REPEAT"
    echo

fi


###############################################################################
# 4. dbSNP 00-All.vcf
###############################################################################
#
# Pipeline config expects:
#
#   small_variants/00-All.vcf
#
# NCBI distributes the file compressed as:
#
#   00-All.vcf.gz
#
# We keep the compressed original and additionally create the uncompressed
# VCF expected by the pipeline.
#
###############################################################################

echo
echo "======================================================================"
echo "4/10  dbSNP 00-All.vcf"
echo "======================================================================"

DBSNP_DIR="$DATA_DIR/small_variants"

DBSNP_GZ="$DBSNP_DIR/00-All.vcf.gz"
DBSNP_VCF="$DBSNP_DIR/00-All.vcf"

mkdir -p "$DBSNP_DIR"

if [[ -s "$DBSNP_VCF" ]]; then

    echo "[SKIP] $DBSNP_VCF already exists"

else

    fetch \
        "https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz" \
        "$DBSNP_GZ"

    echo "[GUNZ] Creating:"
    echo "       $DBSNP_VCF"

    gzip -dc "$DBSNP_GZ" > "${DBSNP_VCF}.tmp"

    mv \
        "${DBSNP_VCF}.tmp" \
        "$DBSNP_VCF"

fi


###############################################################################
# 5. 1000 Genomes GRCh38 SV annotation
###############################################################################
#
# Pipeline config expects:
#
#   ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf
#
# Download compressed, then unpack.
#
###############################################################################

echo
echo "======================================================================"
echo "5/10  1000 Genomes GRCh38 SV VCF"
echo "======================================================================"

SV_ANNOT="$DATA_DIR/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf"

if [[ -s "$SV_ANNOT" ]]; then

    echo "[SKIP] $SV_ANNOT already exists"

else

    fetch \
        "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/supporting/GRCh38_positions/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf.gz" \
        "${SV_ANNOT}.gz"

    echo "[GUNZ] Creating:"
    echo "       $SV_ANNOT"

    gzip -dc "${SV_ANNOT}.gz" > "${SV_ANNOT}.tmp"

    mv \
        "${SV_ANNOT}.tmp" \
        "$SV_ANNOT"

fi


###############################################################################
# 6. Ensembl VEP GRCh38 cache
###############################################################################
#
# Pipeline config:
#
#   vep_cache_dir:
#     .../vep_hg38
#
# The VEP release is pinned through VEP_RELEASE.
#
# Default:
#
#   VEP_RELEASE=112
#
# To use another release:
#
#   VEP_RELEASE=116 DATA_DIR=. bash download_pb_variants_data.sh
#
# The cache itself is release-specific, so changing the release changes the
# biological annotation resource and should therefore be deliberate.
#
###############################################################################

echo
echo "======================================================================"
echo "6/10  Ensembl VEP GRCh38 cache"
echo "======================================================================"

VEP_CACHE="$DATA_DIR/vep_hg38"

mkdir -p "$VEP_CACHE"

# Look for an existing VEP cache rather than merely checking that the
# directory exists.
if find "$VEP_CACHE" -type f -print -quit 2>/dev/null | grep -q .; then

    echo "[SKIP] VEP cache directory is not empty:"
    echo "       $VEP_CACHE"

else

    VEP_TAR="$DATA_DIR/homo_sapiens_merged_vep_${VEP_RELEASE}_GRCh38.tar.gz"

    fetch \
        "https://ftp.ensembl.org/pub/release-${VEP_RELEASE}/variation/vep/homo_sapiens_merged_vep_${VEP_RELEASE}_GRCh38.tar.gz" \
        "$VEP_TAR"

    echo "[UNPK] VEP cache"

    tar \
        -xzf "$VEP_TAR" \
        -C "$VEP_CACHE"

    rm -f "$VEP_TAR"

fi


###############################################################################
# 7. SVTopo / HiFiCNV excluded regions
###############################################################################
#
# Pipeline config:
#
#   svtopo_exclude_regions:
#     .../cnv.excluded_regions.hg38.bed.gz
#
# The resource is a GRCh38 excluded-region BED.
#
###############################################################################

echo
echo "======================================================================"
echo "7/10  GRCh38 CNV excluded regions"
echo "======================================================================"

SVTOPO="$DATA_DIR/cnv.excluded_regions.hg38.bed.gz"

fetch \
    "https://raw.githubusercontent.com/PacificBiosciences/reference_genomes/main/hg38/cnv_excluded_regions/cnv.excluded_regions.bed.gz" \
    "$SVTOPO" || {

    echo
    echo "[WARN] Current PacBio reference-genomes URL did not resolve."
    echo
    echo "The pipeline requires:"
    echo "  $SVTOPO"
    echo
    echo "Check the current PacBio reference-genomes repository before"
    echo "substituting another excluded-region resource."
    echo

}


###############################################################################
# 8. Stranger GRCh38 repeat catalog
###############################################################################
#
# Pipeline config expects the LOCAL filename:
#
#   stranger_repeat_catalog_grch38.json
#
# The Stranger repository calls the upstream file:
#
#   variant_catalog_grch38.json
#
# We therefore download the upstream file but save it using the filename
# expected by this pipeline.
#
###############################################################################

echo
echo "======================================================================"
echo "8/10  Stranger GRCh38 catalog"
echo "======================================================================"

STRANGER="$DATA_DIR/stranger_repeat_catalog_grch38.json"

fetch \
    "https://raw.githubusercontent.com/moonso/stranger/main/stranger/resources/variant_catalog_grch38.json" \
    "$STRANGER"


###############################################################################
# 9. CoLoRSdb GRCh38 v1.2.0
###############################################################################
#
# Pipeline config:
#
#   colorsdb_sv_vcf:
#     .../colorsdb/CoLoRSdb.GRCh38.v1.2.0.pbsv.jasmine.vcf.gz
#
# Keep this as the published bgzipped VCF.
#
# A versioned Zenodo release is preferred over an unversioned external
# download server for reproducibility.
#
###############################################################################

echo
echo "======================================================================"
echo "9/10  CoLoRSdb GRCh38 v1.2.0"
echo "======================================================================"

COLORS_DIR="$DATA_DIR/colorsdb"

COLORS_VCF="$COLORS_DIR/CoLoRSdb.GRCh38.v1.2.0.pbsv.jasmine.vcf.gz"

mkdir -p "$COLORS_DIR"

fetch \
    "https://zenodo.org/records/14814308/files/CoLoRSdb.GRCh38.v1.2.0.pbsv.jasmine.vcf.gz?download=1" \
    "$COLORS_VCF"

# Download the published index as well.
if [[ ! -s "${COLORS_VCF}.tbi" ]]; then

    fetch \
        "https://zenodo.org/records/14814308/files/CoLoRSdb.GRCh38.v1.2.0.pbsv.jasmine.vcf.gz.tbi?download=1" \
        "${COLORS_VCF}.tbi"

fi


###############################################################################
# 10. DeepVariant SIF
###############################################################################
#
# Pipeline config:
#
#   sif_image_deepvariant:
#     .../deepvariant_1.9.0.sif
#
# This is OPTIONAL.
#
# It is only required when DeepVariant is run locally rather than through
# the HPC-provided installation/module.
#
# We do not require Singularity/Apptainer just to run this script.
#
###############################################################################

echo
echo "======================================================================"
echo "10/10  DeepVariant 1.9.0 SIF (optional)"
echo "======================================================================"

DEEPVARIANT_SIF="$DATA_DIR/deepvariant_1.9.0.sif"

if [[ -s "$DEEPVARIANT_SIF" ]]; then

    echo "[SKIP] $DEEPVARIANT_SIF already exists"

else

    if command -v singularity >/dev/null 2>&1; then

        echo "[GET ] DeepVariant 1.9.0 using Singularity"

        singularity pull \
            --name "$DEEPVARIANT_SIF" \
            "docker://google/deepvariant:1.9.0"

    elif command -v apptainer >/dev/null 2>&1; then

        echo "[GET ] DeepVariant 1.9.0 using Apptainer"

        apptainer pull \
            --name "$DEEPVARIANT_SIF" \
            "docker://google/deepvariant:1.9.0"

    else

        echo
        echo "[INFO] Singularity/Apptainer not found."
        echo
        echo "DeepVariant SIF was NOT downloaded."
        echo
        echo "If required later, run on a machine with Singularity:"
        echo
        echo "  singularity pull --name \"$DEEPVARIANT_SIF\" \\"
        echo "      docker://google/deepvariant:1.9.0"
        echo

    fi
fi


###############################################################################
# PIPELINE-SPECIFIC / CLUSTER-SPECIFIC RESOURCES
###############################################################################

echo
echo "======================================================================"
echo "Resources intentionally NOT downloaded"
echo "======================================================================"


###############################################################################
# PacBio Lima barcode FASTA
###############################################################################

BARCODE_FASTA="$DATA_DIR/../input/pacbio_barcodes.fasta"

echo
echo "1. PacBio Lima barcode FASTA"
echo
echo "Expected:"
echo "  $BARCODE_FASTA"
echo
echo "This is sequencing-run / lab-specific and must not be replaced by"
echo "an arbitrary public barcode file."


###############################################################################
# Kraken2
###############################################################################

KRAKEN_DB="/gpfs/project/databases/Kraken2-2022-09-28/kraken_db_plus"

echo
echo "2. Kraken2 database"
echo
echo "Configured database:"
echo "  $KRAKEN_DB"
echo
echo "This is a cluster-wide pre-built database."
echo "The script does not download or modify it."


###############################################################################
# DeepVariant checkpoint
###############################################################################

echo
echo "3. DeepVariant checkpoint"
echo
echo "Configured path:"
echo "  $DATA_DIR/deepvariant_pacbio_190_model"
echo
echo "Not currently needed."
echo "No download is performed."


###############################################################################
# Historical TRGT catalog reminder
###############################################################################

echo
echo "4. Historical TRGT catalog"
echo
echo "Required if the pipeline uses:"
echo "  $TRGT_REPEAT"
echo
echo "No newer catalog is substituted automatically."


###############################################################################
# FINAL CHECK
###############################################################################

echo
echo "======================================================================"
echo "FINAL CHECK"
echo "======================================================================"

check_file()
{
    local f="$1"

    if [[ -s "$f" ]]; then
        echo "[ OK ] $f"
    else
        echo "[MISS] $f"
    fi
}

check_dir()
{
    local d="$1"

    if [[ -d "$d" ]]; then
        echo "[ OK ] $d/"
    else
        echo "[MISS] $d/"
    fi
}


echo
echo "--- Core files ---"

check_file "$DATA_DIR/gnomad.v4.1.sv.sites.vcf.gz"
check_file "$DATA_DIR/gnomad.v4.1.sv.sites.vcf.gz.tbi"

check_dir "$DATA_DIR/annotsv_data/annotsv_353/AnnotSV"

check_file "$DATA_DIR/repeat_catalog_v1.hg38.1_to_1000bp_motifs.TRGT.bed"

check_file "$DATA_DIR/small_variants/00-All.vcf"

check_file \
    "$DATA_DIR/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf"

check_dir "$DATA_DIR/vep_hg38"

check_file "$DATA_DIR/cnv.excluded_regions.hg38.bed.gz"

check_file "$DATA_DIR/stranger_repeat_catalog_grch38.json"

check_file \
    "$DATA_DIR/colorsdb/CoLoRSdb.GRCh38.v1.2.0.pbsv.jasmine.vcf.gz"

check_file \
    "$DATA_DIR/colorsdb/CoLoRSdb.GRCh38.v1.2.0.pbsv.jasmine.vcf.gz.tbi"

check_file "$DATA_DIR/deepvariant_1.9.0.sif"


echo
echo "--- Cluster/lab-specific resources ---"

if [[ -s "$BARCODE_FASTA" ]]; then
    echo "[ OK ] $BARCODE_FASTA"
else
    echo "[----] $BARCODE_FASTA"
fi

if [[ -d "$KRAKEN_DB" ]]; then
    echo "[ OK ] $KRAKEN_DB"
else
    echo "[----] $KRAKEN_DB"
fi


###############################################################################
# SUMMARY
###############################################################################

echo
echo "======================================================================"
echo "SETUP FINISHED"
echo "======================================================================"
echo
echo "Reference-data directory:"
echo "  $DATA_DIR"
echo
echo "Run again safely with:"
echo
echo "  DATA_DIR=\"$DATA_DIR\" bash $0"
echo
echo "Existing resources will be skipped."
echo
echo "IMPORTANT:"
echo "  Review any [MISS] entries above before running the pipeline."
echo
```
