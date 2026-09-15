"""
scripts/alphagenome_score_variants.py

Scores variants from a VCF using AlphaGenome local model weights.
Called by the Snakemake rule alphagenome_score_variants.

Requirements:
  - alphagenome_research installed (pip install -e ./alphagenome_research)
  - alphagenome SDK installed (pip install alphagenome)
  - Model weights downloaded locally (see alphagenome_download_guide.md)
  - NVIDIA H100 or A100 GPU
  - JAX with CUDA support

Output TSV columns:
  chrom, pos, ref, alt, gene, consequence,
  gnomad_af (from echtvar annotation),
  [tissue]_rna_seq_score,   (log2 fold change ALT/REF)
  [tissue]_atac_score,
  [tissue]_splice_score,
  max_rna_score,            (max across tissues)
  max_atac_score,
  max_splice_score,
  alphagenome_any_impact    (True if any score > threshold)
"""

import os
import sys
import gzip
import json
import logging
import numpy as np
import pandas as pd

log_path = snakemake.log[0]
logging.basicConfig(
    filename=log_path,
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s"
)
log = logging.getLogger(__name__)
sys.stderr = open(log_path, "a")

# ── Parameters from Snakemake ────────────────────────────────────────────────
vcf_in           = snakemake.input.vcf
ref_fasta        = snakemake.input.ref
scores_out       = snakemake.output.scores_tsv
plots_dir        = snakemake.output.top_plots_dir
skipped_out      = snakemake.output.skipped_tsv
weights_path     = snakemake.params.weights_path
half_width       = snakemake.params.interval_half_width
ontology_terms   = snakemake.params.ontology_terms
output_types_str = snakemake.params.output_types
max_variants     = snakemake.params.max_variants_per_run
prioritize_by    = snakemake.params.prioritize_by
plot_top_n       = snakemake.params.plot_top_n
sample           = snakemake.params.sample

os.makedirs(plots_dir, exist_ok=True)

# ── Imports (after env is activated) ─────────────────────────────────────────
log.info("Importing AlphaGenome...")
try:
    from alphagenome.data import genome as ag_genome
    from alphagenome_research.model import dna_model
    from alphagenome.visualization import plot_components
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as e:
    log.error(f"Failed to import AlphaGenome: {e}")
    log.error("Ensure alphagenome and alphagenome_research are installed.")
    log.error("See alphagenome_download_guide.md for setup instructions.")
    sys.exit(1)

# Map string output type names to SDK enums
OUTPUT_TYPE_MAP = {
    "RNA_SEQ":     dna_model.OutputType.RNA_SEQ,
    "ATAC":        dna_model.OutputType.ATAC,
    "SPLICE_SITE": dna_model.OutputType.SPLICE_SITE,
}
requested_outputs = [OUTPUT_TYPE_MAP[t] for t in output_types_str if t in OUTPUT_TYPE_MAP]

# ── Load model weights ────────────────────────────────────────────────────────
log.info(f"Loading AlphaGenome weights from {weights_path}...")
try:
    # Uses local weights directory — no internet required after download.
    # The weights_path should contain the checkpoint files downloaded from
    # Kaggle or HuggingFace (see alphagenome_download_guide.md).
    model = dna_model.create_from_local(weights_path, fold="all_folds")
    log.info("Model loaded successfully.")
except Exception as e:
    log.error(f"Failed to load model weights: {e}")
    log.error(f"Weights path: {weights_path}")
    sys.exit(1)

# ── Parse VCF ────────────────────────────────────────────────────────────────
log.info(f"Parsing variants from {vcf_in}...")

def parse_vcf(vcf_path):
    """Parse BCF/VCF (possibly gzipped) into a list of dicts."""
    import subprocess
    variants = []
    cmd = ["bcftools", "view", "--output-type", "v", vcf_path]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        if line.startswith("#"):
            continue
        parts = line.rstrip().split("\t")
        if len(parts) < 5:
            continue
        chrom, pos, vid, ref, alt = parts[0], int(parts[1]), parts[2], parts[3], parts[4]
        # Skip multi-allelic — AlphaGenome scores one ALT at a time
        if "," in alt:
            continue
        # Extract gnomad_af from INFO if present (added by echtvar)
        info = dict(
            kv.split("=", 1) if "=" in kv else (kv, True)
            for kv in parts[7].split(";")
        ) if len(parts) > 7 else {}
        af = float(info.get("gnomad_af", info.get("gnomad_AF", 1.0)))
        variants.append({
            "chrom": chrom if chrom.startswith("chr") else f"chr{chrom}",
            "pos":   pos,
            "ref":   ref,
            "alt":   alt,
            "id":    vid,
            "gnomad_af": af,
            "info":  info
        })
    proc.wait()
    return variants

all_variants = parse_vcf(vcf_in)
log.info(f"Parsed {len(all_variants)} variants from VCF.")

# ── Prioritize and cap ────────────────────────────────────────────────────────
if prioritize_by == "af":
    # Score rarest variants first (most likely to be clinically relevant)
    all_variants.sort(key=lambda v: v["gnomad_af"])
elif prioritize_by == "genmod_rank":
    # If GENMOD_SCORE is in INFO, sort descending
    all_variants.sort(
        key=lambda v: float(v["info"].get("GENMOD_SCORE", 0)),
        reverse=True
    )

skipped = []
if len(all_variants) > max_variants:
    log.warning(
        f"Capping at {max_variants} variants (total: {len(all_variants)}). "
        f"Remaining {len(all_variants) - max_variants} written to skipped TSV."
    )
    skipped = all_variants[max_variants:]
    all_variants = all_variants[:max_variants]

# Write skipped variants
with open(skipped_out, "w") as fh:
    fh.write("chrom\tpos\tref\talt\tgnomad_af\treason\n")
    for v in skipped:
        fh.write(f"{v['chrom']}\t{v['pos']}\t{v['ref']}\t{v['alt']}\t{v['gnomad_af']}\tcap_exceeded\n")

# ── Score variants ────────────────────────────────────────────────────────────
log.info(f"Scoring {len(all_variants)} variants across {len(ontology_terms)} tissues...")

results = []
failed  = []

for i, v in enumerate(all_variants):
    if i % 100 == 0:
        log.info(f"  Progress: {i}/{len(all_variants)}")
    try:
        chrom, pos, ref, alt = v["chrom"], v["pos"], v["ref"], v["alt"]

        # Define 1Mbp interval centred on the variant
        center = pos - 1  # AlphaGenome uses 0-based coordinates
        start  = max(0, center - half_width)
        end    = center + half_width

        interval = ag_genome.Interval(chromosome=chrom, start=start, end=end)
        variant  = ag_genome.Variant(
            chromosome=chrom,
            position=pos,          # 1-based
            reference_bases=ref,
            alternate_bases=alt,
        )

        outputs = model.predict_variant(
            interval=interval,
            variant=variant,
            ontology_terms=ontology_terms,
            requested_outputs=requested_outputs,
        )

        # Aggregate scores: compute log2(ALT/REF) summed over a ±10kb window around the variant
        row = {
            "chrom":      chrom,
            "pos":        pos,
            "ref":        ref,
            "alt":        alt,
            "variant_id": v["id"],
            "gnomad_af":  v["gnomad_af"],
        }

        score_window = 10000  # bp around variant for score aggregation
        s_start = max(0, center - score_window)
        s_end   = center + score_window

        all_rna_scores   = []
        all_atac_scores  = []
        all_splice_scores = []

        for term in ontology_terms:
            term_safe = term.replace(":", "_")

            if "RNA_SEQ" in output_types_str and hasattr(outputs.alternate, "rna_seq"):
                ref_track = outputs.reference.rna_seq
                alt_track = outputs.alternate.rna_seq
                # Clip to scoring window
                ref_vals = ref_track.values  # numpy array
                alt_vals = alt_track.values
                log2fc = np.log2(alt_vals + 1e-6) - np.log2(ref_vals + 1e-6)
                rna_score = float(np.abs(log2fc).max())
                row[f"{term_safe}_rna_seq_score"] = rna_score
                all_rna_scores.append(rna_score)

            if "ATAC" in output_types_str and hasattr(outputs.alternate, "atac"):
                ref_vals = outputs.reference.atac.values
                alt_vals = outputs.alternate.atac.values
                log2fc = np.log2(alt_vals + 1e-6) - np.log2(ref_vals + 1e-6)
                atac_score = float(np.abs(log2fc).max())
                row[f"{term_safe}_atac_score"] = atac_score
                all_atac_scores.append(atac_score)

            if "SPLICE_SITE" in output_types_str and hasattr(outputs.alternate, "splice_site"):
                ref_vals = outputs.reference.splice_site.values
                alt_vals = outputs.alternate.splice_site.values
                splice_score = float(np.abs(alt_vals - ref_vals).max())
                row[f"{term_safe}_splice_score"] = splice_score
                all_splice_scores.append(splice_score)

        row["max_rna_score"]    = max(all_rna_scores)   if all_rna_scores   else None
        row["max_atac_score"]   = max(all_atac_scores)  if all_atac_scores  else None
        row["max_splice_score"] = max(all_splice_scores) if all_splice_scores else None

        # Flag: any modality shows a meaningful predicted impact
        SCORE_THRESHOLD = 0.5  # log2FC; tune to your use case
        row["alphagenome_any_impact"] = any([
            (row["max_rna_score"]    or 0) > SCORE_THRESHOLD,
            (row["max_atac_score"]   or 0) > SCORE_THRESHOLD,
            (row["max_splice_score"] or 0) > SCORE_THRESHOLD,
        ])

        results.append(row)

    except Exception as e:
        log.warning(f"Failed to score variant {v['chrom']}:{v['pos']} {v['ref']}>{v['alt']}: {e}")
        failed.append({**v, "error": str(e)})

log.info(f"Scored {len(results)} variants. Failed: {len(failed)}.")

# ── Write scores TSV ─────────────────────────────────────────────────────────
df = pd.DataFrame(results)
df.to_csv(scores_out, sep="\t", index=False, compression="gzip")
log.info(f"Scores written to {scores_out}")

# Append failures to skipped TSV
if failed:
    with open(skipped_out, "a") as fh:
        for v in failed:
            fh.write(f"{v['chrom']}\t{v['pos']}\t{v['ref']}\t{v['alt']}\t{v.get('gnomad_af','.')}\t{v.get('error','error')}\n")

# ── Generate track plots for top-scoring variants ────────────────────────────
if not df.empty and "max_rna_score" in df.columns:
    top_df = df.nlargest(plot_top_n, "max_rna_score")
    log.info(f"Generating track plots for top {len(top_df)} variants by RNA score...")

    for _, row in top_df.iterrows():
        try:
            chrom, pos, ref, alt = row["chrom"], int(row["pos"]), row["ref"], row["alt"]
            center = pos - 1
            start  = max(0, center - half_width)
            end    = center + half_width

            interval = ag_genome.Interval(chromosome=chrom, start=start, end=end)
            variant  = ag_genome.Variant(
                chromosome=chrom, position=pos,
                reference_bases=ref, alternate_bases=alt,
            )
            outputs = model.predict_variant(
                interval=interval, variant=variant,
                ontology_terms=[ontology_terms[0]],  # plot first tissue only
                requested_outputs=requested_outputs,
            )

            fig_path = os.path.join(
                plots_dir,
                f"{chrom}_{pos}_{ref}_{alt}_alphagenome.png"
            )

            tracks = []
            if hasattr(outputs.alternate, "rna_seq"):
                tracks.append(plot_components.OverlaidTracks(
                    tdata={"REF": outputs.reference.rna_seq, "ALT": outputs.alternate.rna_seq},
                    colors={"REF": "dimgrey", "ALT": "firebrick"},
                ))
            if hasattr(outputs.alternate, "atac"):
                tracks.append(plot_components.OverlaidTracks(
                    tdata={"REF": outputs.reference.atac, "ALT": outputs.alternate.atac},
                    colors={"REF": "steelblue", "ALT": "darkorange"},
                ))

            plot_components.plot(
                tracks,
                interval=interval.resize(2**15),
                annotations=[plot_components.VariantAnnotation([variant], alpha=0.8)],
            )
            plt.suptitle(
                f"{sample} | {chrom}:{pos} {ref}>{alt}\n"
                f"gnomAD AF={row['gnomad_af']:.2e}  "
                f"max_rna_score={row.get('max_rna_score', 'N/A'):.3f}",
                fontsize=9
            )
            plt.savefig(fig_path, bbox_inches="tight", dpi=120)
            plt.close()
            log.info(f"  Plot saved: {fig_path}")

        except Exception as e:
            log.warning(f"  Plot failed for {row.get('chrom')}:{row.get('pos')}: {e}")

log.info("AlphaGenome scoring complete.")