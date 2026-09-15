# Getting started with your `pacbio_variants` results

You ran the [pacbio_variants](https://github.com/WestGermanGenomeCenter/pacbio_variants) pipeline on **5 human samples (hg38)**. This guide walks you through what to look at first, in a sensible order, so you're not staring at 20,000 output files wondering where to start.

Your 5 samples (barcodes): `bc2068`, `bc2069`, `bc2070`, `bc2071`, `bc2072`.

> Tip: throughout the output folder, files are simply named `all_smrtcells_1519_<barcode>_...`. Just swap in the barcode of the sample you care about.

---

## 0. Start here: did the run actually go well?

Open, in your browser:

```
multiqc_report.html
```

This is the single best "is everything ok" file. [MultiQC](https://multiqc.info/) collects QC metrics from mosdepth (coverage) and, if enabled, Kraken2 (contamination) across all 5 samples into one report with plots you can compare side by side.

What to look for:
- **Coverage** (from mosdepth) — is it roughly what you expect for HiFi WGS (commonly ~20–30x)? Are all 5 samples similar, or is one an outlier?
- **Contamination** (if Kraken2 was enabled) — do the samples look "clean" (mostly human), or is there a large fraction of non-human reads?

Per-sample raw numbers, if you want them, live in:
```
mosdepth/all_smrtcells_1519_<barcode>.mosdepth.summary.txt
```

Also worth a quick look — the DAG/rulegraph that shows exactly what ran for this execution:
```
pb_variants_rulegraph.2026_09_14_10_11_AM.pdf
```
...and your exact settings for this run are archived here, so you always know what you did:
```
config_2026_09_14_10_11_AM.yaml
samplesheet_2026_09_14_10_11_AM.csv
```

**Rule of thumb:** if `multiqc_report.html` exists and coverage looks sane, the pipeline finished successfully — the [README's own troubleshooting section](https://github.com/WestGermanGenomeCenter/pacbio_variants#supervising-the-run) uses exactly this file as the "pipeline is done" signal.

### Two things in the MultiQC report that can look confusing at first

**1. Why you see the same variants twice — once "whatshap", once "longphase":**
Both DeepVariant SNPs/indels and Sniffles SVs get phased by *two different tools* — [WhatsHap](https://whatshap.readthedocs.io/) and [LongPhase](https://github.com/twolinin/longphase) — run independently on the same underlying calls. So it's normal, not a bug, that MultiQC (and the `annotated_variants/` folder) shows what looks like duplicate entries per sample: `vep_whatshap_...` and `vep_longphase_...` are usually annotating the *same input variants*, just carrying two independent phasing results. Pick whichever phasing tool's output you trust more (or compare them), rather than treating both as separate variant sets.

**2. What the "overlap"/comparison numbers (truvari, rtg-tools) actually mean:**
The pipeline also cross-checks its callers against each other:
- **SVs:** [truvari](https://github.com/ACEnglish/truvari) compares sawfish calls against Sniffles, with **Sniffles treated as the "truth" set**.
- **SNPs/indels:** [rtg-tools vcfeval](https://github.com/RealTimeGenomics/rtg-tools) compares NanoCaller calls against DeepVariant (or bcftools), with **DeepVariant treated as the "truth" set**.

Important: neither DeepVariant nor Sniffles is actual ground truth (e.g. from an orthogonal validated dataset) — they're just the more established/mature caller in each pair, used as a convenient reference point so the second caller's calls can be scored as overlapping/matching or not. Treat these numbers as "how well do the two callers agree with each other", not "how accurate are the calls".

---

## 1. Look at the small variants (SNPs/indels) — DeepVariant

DeepVariant is the primary small-variant caller here. For each sample:

```
variants/deepvariant_all_smrtcells_1519_<barcode>/all_smrtcells_1519_<barcode>_variants.vcf.gz
```

Before opening the VCF itself, open DeepVariant's own built-in HTML report — no extra tool needed, just double-click/open in a browser:

```
variants/deepvariant_all_smrtcells_1519_<barcode>/all_smrtcells_1519_<barcode>_variants.visual_report.html
```

This gives you a quick visual summary (quality distributions, Ti/Tv ratio, indel size spectrum, etc.) per sample — a good gut-check before digging into individual variants.

The calls then get **phased** (i.e., assigned to a parental chromosome copy) in two stages:
- `variants/longphase_all_smrtcells_1519_<barcode>/all_smrtcells_1519_<barcode>_phased.vcf.gz` — after [LongPhase](https://github.com/twolinin/longphase)
- `variants/whatshap_all_smrtcells_1519_<barcode>/all_smrtcells_1519_<barcode>_phased.vcf.gz` — after [WhatsHap](https://whatshap.readthedocs.io/)

These same directories also contain the **haplotagged BAM files** (reads tagged with which haplotype they belong to) — useful later for IGV:
```
variants/whatshap_all_smrtcells_1519_<barcode>/all_smrtcells_1519_<barcode>_haplotaged.bam
```

There's a second, independent small-variant caller in the pipeline too — **NanoCaller** (`variants/nanocaller_all_smrtcells_1519_<barcode>/`) — useful as a cross-check against DeepVariant if you want extra confidence in a specific call.

---

## 2. Look at the structural variants (SVs) — Sniffles

For structural variants (deletions, insertions, duplications, inversions, translocations), start with [Sniffles2](https://github.com/fritzsedlazeck/Sniffles):

```
variants/sniffles_all_smrtcells_1519_<barcode>/all_smrtcells_1519_<barcode>_svs.vcf
```

Phased by LongPhase into:
```
variants/longphase_all_smrtcells_1519_<barcode>/all_smrtcells_1519_<barcode>_phased_SV.vcf.gz
```

There's a second, independent SV/CNV caller you can cross-check against: **sawfish**, which also reports its own phasing:
```
variants/sawfish_phased_all_smrtcells_1519_<barcode>/
```
It also gives you a per-sample copy-number bedgraph track (handy in IGV):
```
variants/sawfish_phased_all_smrtcells_1519_<barcode>/samples/<sample_id>/copynum.bedgraph
```

Other specialised variant types worth knowing about, each in its own `variants/<tool>_.../` folder:
- **HiFiCNV** — copy-number variants, with a bigWig depth track (`*.depth.bw`) great for IGV
- **TRGT** — tandem repeat genotyping
- **Paraphase** — variants in segmentally-duplicated/paralogous genes
- **Mitorsaw** — mitochondrial variants (hg38 only), which also ships a ready-made IGV alignment (`mito_igv_custom/custom_alignments.bam`)

---

## 3. Check whether annotation actually helped

If you enabled annotation in `config.yaml` (it looks like you did, based on your output), every variant type gets run through an annotator, so you get gene names, predicted consequences (missense/nonsense/etc.), population frequencies, and clinical significance where available, instead of bare coordinates.

**SNPs/indels:**
- [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) (Ensembl Variant Effect Predictor): `annotated_variants/vep_whatshap_all_smrtcells_1519_<barcode>/all_smrtcells_1519_<barcode>_vep_whatshap_annotated.vcf`
  - VEP also writes a handy standalone HTML summary — open this first, it's much smaller and faster than the full VCF: `..._vep_whatshap_annotated.vcf_summary.html`
- [SnpSift](https://pcingola.github.io/SnpSift/) adds further annotation on top: `annotated_variants/snpsift_whatshap_all_smrtcells_1519_<barcode>/all_smrtcells_1519_<barcode>_snpsift_whatshap_annotated.vcf`

**SVs/CNVs:**
- [AnnotSV](https://lbgi.fr/AnnotSV/) for sawfish SVs, output as a **TSV table** (opens fine in Excel/pandas, easier to scan than a VCF): `annotated_variants/annotsv_sawfish_all_smrtcells_1519_<barcode>/all_smrtcells_1519_<barcode>_genotyped.sv.annotated.tsv`
- [sansa](https://github.com/isbrandtb/sansa) for both sawfish and sniffles SVs, output as annotated BCF + a friendlier `.tsv.gz`: `annotated_variants/sansa_svs_sniffles_all_smrtcells_1519_<barcode>/all_smrtcells_1519_<barcode>_sniffles_longphase_annotated.tsv.gz`

**Did annotation help?** A quick sanity check: open one of the `_summary.html` VEP files, or the AnnotSV `.tsv`, and see whether known genes/consequences show up for a couple of variants you already spotted in the raw VCF. If the raw VCF just gave you `chr1:12345 A>G` and the annotated file now tells you it's a missense variant in a named gene — annotation is doing its job.

---

## 4. Open things up visually — IGV

For visual inspection, [IGV (Integrative Genomics Viewer)](https://igv.org/) — either the [desktop app](https://igv.org/doc/desktop/) or [igv.js](https://igv.org/doc/igvjs/) in-browser — is the standard tool. Load the **hg38** reference genome, then add tracks:

| What to look at | File(s) to load |
|---|---|
| Raw alignments | `bams/all_smrtcells_1519_<barcode>_aligned.bam` (+ its `.bai`) |
| Haplotype-tagged alignments (color reads by haplotype) | `variants/whatshap_all_smrtcells_1519_<barcode>/all_smrtcells_1519_<barcode>_haplotaged.bam` |
| SNPs/indels | `variants/deepvariant_all_smrtcells_1519_<barcode>/all_smrtcells_1519_<barcode>_variants.vcf.gz` |
| Structural variants | `variants/sniffles_all_smrtcells_1519_<barcode>/all_smrtcells_1519_<barcode>_svs.vcf` |
| Copy number depth track | `variants/hificnv_all_smrtcells_1519_<barcode>/all_smrtcells_1519_<barcode>.<sample_id>.depth.bw` |
| Copy number segments | `variants/sawfish_phased_all_smrtcells_1519_<barcode>/samples/<sample_id>/copynum.bedgraph` |
| Methylation (5mC) tracks | `variants/cpg_tools_all_smrtcells_1519_<barcode>/all_smrtcells_1519_<barcode>.hap1.bed.gz` / `.hap2.bed.gz` — load together with the haplotagged BAM above |
| Mitochondrial variants + reads | `variants/mitorsaw_all_smrtcells_1519_<barcode>/mito_igv_custom/custom_alignments.bam` + `all_smrtcells_1519_<barcode>_mitochondiral_variants.vcf.gz` |

A good first exercise: pick one interesting variant from the DeepVariant or Sniffles VCF, navigate IGV to that coordinate, and visually confirm the read support looks convincing (enough reads, consistent across both haplotypes if heterozygous, no obvious mapping artifacts).

For **Paraphase** results specifically, the pipeline can generate its own visualization via [paraviewer](https://github.com/PacificBiosciences/paraphase) if enabled in the config — check `variants/paraphase_all_smrtcells_1519_<barcode>/` for those plots before resorting to plain IGV, since paralog regions are easier to misread in a generic viewer.

---

## Suggested order of operations (TL;DR)

1. **`multiqc_report.html`** — confirm coverage/QC looks good for all 5 samples.
2. **DeepVariant `visual_report.html`** — quick per-sample SNP/indel QC.
3. **Sniffles `_svs.vcf`** — skim the structural variants.
4. **VEP `_summary.html`** and **AnnotSV `.tsv`** — check annotation is adding real gene/consequence information.
5. **IGV** — load reference (hg38) + haplotagged BAM + your VCFs of interest, and visually inspect a handful of variants that look interesting or unexpected.

## Helpful links

- [MultiQC docs](https://multiqc.info/docs/)
- [DeepVariant](https://github.com/google/deepvariant)
- [Sniffles2](https://github.com/fritzsedlazeck/Sniffles)
- [sawfish](https://github.com/PacificBiosciences/sawfish)
- [WhatsHap](https://whatshap.readthedocs.io/)
- [LongPhase](https://github.com/twolinin/longphase)
- [VEP (Ensembl)](https://www.ensembl.org/info/docs/tools/vep/index.html)
- [SnpSift](https://pcingola.github.io/SnpSift/)
- [AnnotSV](https://lbgi.fr/AnnotSV/)
- [IGV desktop](https://igv.org/doc/desktop/) / [igv.js](https://igv.org/doc/igvjs/)
- [pacbio_variants pipeline README](https://github.com/WestGermanGenomeCenter/pacbio_variants)

---

*Note: several output files (BAMs, VEP-annotated VCFs, cpg_tools bed files) are tens of GB — if you're viewing them remotely, use `tabix`/`samtools view region` to pull small regions instead of downloading whole files, or point IGV directly at files on the server if it supports remote/streamed access.*