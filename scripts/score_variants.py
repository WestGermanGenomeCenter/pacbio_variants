#!/usr/bin/env python3
"""
prioritize_snvs.py (v2.1) - shortlist of plausible SNVs/indels from VEP-annotated
DeepVariant and/or NanoCaller VCFs. Stdlib only.

Per variant:  score = impact * rarity * confidence * bonuses
  impact      VEP IMPACT on the PRIMARY transcript, chosen VEP-style (MANE Select / canonical
              protein-coding first; a non-canonical coding transcript only if more severe;
              NMD / retained-intron / non-coding only as a last resort). A more severe
              consequence on another protein-coding transcript is reported, not scored.
              ClinVar P/LP sets a floor of 0.9 (HIGH = 1.0); ClinVar benign is down-weighted.
  rarity      from MAX_AF (fallback gnomADg/gnomADe/1000G AF); no frequency = rare, but
              capped at 0.8 if the variant has an rsID; P/LP variants that are common are
              kept with a low rarity floor and flagged
  confidence  mean over callers that PASS: QUAL (saturating, per-caller scale), VAF consistency
              with the genotype (two-sided binomial p-value, no penalty above p = 0.05),
              depth, strand balance (if ADF/ADR exist); x0.75 if only one of two callers
  bonus       x1.25 if a rare high/moderate variant has a partner in the same gene
              on the other haplotype (phased, same phase set, "in trans")

Also: only chr1-22 and chrX are kept (alt/unplaced contigs, chrY and chrM skipped unless added with
--include-contigs); non-coding and pseudogene transcripts and coding genes without a VEP GENE_PHENO
flag are down-weighted (gene_rel); --mask-bed regions (segdups, low mappability) halve confidence.

RAM: records that cannot reach the shortlist (impact below --min-impact, AF above
--max-af, unless ClinVar P/LP) are skipped while streaming. Relax these two options
only if you really need --output-all to contain everything.

Usage:
  prioritize_snvs_v2.py --deepvariant dv.vcf --nanocaller nc.vcf \
      --reference reference.fasta --output shortlist.tsv --output-all all.tsv
"""
import argparse
import bisect
import csv
import gzip
import math
import os
import re
import resource
import sys
import time
from collections import Counter, defaultdict
from itertools import combinations

# ---- tunable constants (placeholders, to be calibrated) -------------------------
IMPACT_RANK = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}
IMPACT_SCORE = {"HIGH": 1.0, "MODERATE": 0.6, "LOW": 0.1, "MODIFIER": 0.0}
CLINVAR_PLP_SCORE = 0.9      # floor for ClinVar P/LP (HIGH = 1.0)
PLP_COMMON_RARITY_FLOOR = 0.15  # P/LP variants above --max-af are kept with this rarity
RSID_NO_AF_RARITY = 0.8       # rsID present but no population frequency: not "maximally rare"
ISOFORM_ONLY_REL = 0.3        # effect only on a retained-intron / NMD / similar isoform
ISOFORM_ONLY_BIOTYPES = {"retained_intron", "nonsense_mediated_decay", "processed_transcript",
                         "non_stop_decay", "protein_coding_CDS_not_defined", "TEC"}
VAF_P_FULL = 0.05            # binomial p at/above which the VAF term is 1.0
HOM_EXPECTED_VAF = 0.97
QUAL_FULL = {"dv": 30.0, "nc": 30.0}   # QUAL at/above which the QUAL term saturates (to be calibrated, e.g. GIAB)
CLINVAR_BENIGN_FACTOR = 0.1
COMPHET_TRANS_BONUS = 1.25
COMPHET_MAX_AF = 0.01          # partners must be rarer than this
COMPHET_MIN_IMPACT_SCORE = 0.6 # ... and at least MODERATE (or ClinVar P/LP)
COMPHET_MAX_GENE_N = 30        # skip pairing in genes with more candidates
FREQ_FALLBACK = ("gnomADg_AF", "gnomADe_AF", "AF")
CONF_WEIGHTS = {"qual": 0.35, "vaf": 0.25, "depth": 0.25, "strand": 0.15}
SINGLE_CALLER_FACTOR = 0.75
RARITY_COMMON = 0.05   # AF at/above this -> rarity 0
RARITY_FLOOR = 1e-5    # AF at/below this -> rarity 1
NONCODING_REL = 0.15   # gene relevance of lncRNA/other non-coding transcripts without a phenotype link
PSEUDOGENE_REL = 0.10
CODING_NO_PHENO_REL = 0.6   # protein-coding gene with no GENE_PHENO flag from VEP
MASK_FACTOR = 0.5      # confidence x this inside a --mask-bed region (segdups, low mappability, ...)
PROBLEM_GENE_FACTOR = 0.5   # confidence x this in repeat-rich / paralogous gene families
CLUSTER_FACTOR = 0.7        # confidence x this for a variant inside a dense cluster of candidates
# calls in these families are dominated by mapping artefacts (D4Z4, kringle IV, mucin VNTRs, ...)
ARTIFACT_GENE_RE = re.compile(r"^(DUX4.*|LPA|ZAN|MUC\d+[A-Z]?|NBPF\d.*|KIR\d.*|HLA-.*|OR\d+[A-Z]+\d+.*|"
                              r"PRAMEF\d.*|GOLGA6L\d.*|TBC1D3.*|USP17L.*|FRG2.*)$")
# real disease genes with close paralogs/pseudogenes: not penalised, but confirm (e.g. paraphase)
PARALOG_GENES = {"SMN1", "SMN2", "CYP21A2", "CYP2D6", "GBA", "GBA1", "STRC", "PMS2", "NCF1", "HBA1",
                 "HBA2", "IKBKG", "PKD1", "SBDS", "CHEK2", "OTOA"}
PRIMARY_RE = re.compile(r"^(chr)?([1-9]|1[0-9]|2[0-2]|X)$")
KEEP_CSQ = ("SYMBOL", "Gene", "Feature", "BIOTYPE", "Consequence", "IMPACT", "HGVSc",
            "HGVSp", "CLIN_SIG", "SIFT", "PolyPhen", "MANE_SELECT", "Existing_variation",
            "GENE_PHENO")


# ---------------------------------------------------------------- helpers ----------
def to_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def int_list(s):
    if not s:
        return None
    try:
        return [int(x) for x in s.split(",")]
    except ValueError:
        return None


def fmt(x, nd=3):
    if x is None:
        return ""
    return f"{x:.{nd}g}" if isinstance(x, float) else str(x)


def fmt_time(s):
    s = int(s)
    return f"{s // 3600}:{s % 3600 // 60:02d}:{s % 60:02d}" if s >= 3600 else f"{s // 60}:{s % 60:02d}"


def peak_ram_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024  # Linux: KB


class Progress:
    """Progress by bytes of the (compressed) input file; \r on a terminal, lines in logs."""

    def __init__(self, label, total_bytes, interval):
        self.label, self.total = label, max(total_bytes, 1)
        self.tty = sys.stderr.isatty()
        self.interval = interval if interval is not None else (1.0 if self.tty else 30.0)
        self.t0 = self.last = time.monotonic()

    def tick(self, lines, pos, force=False):
        now = time.monotonic()
        if not force and now - self.last < self.interval:
            return
        self.last = now
        pct = 100.0 * pos / self.total
        el = now - self.t0
        eta = el * (100 - pct) / pct if pct > 0 else 0
        msg = (f"[{self.label}] {pct:5.1f}% | {lines:,} lines | "
               f"{lines / max(el, 1e-9) / 1000:.0f}k lines/s | ETA {fmt_time(eta)} | "
               f"peak RAM {peak_ram_gb():.1f} GB")
        sys.stderr.write(("\r" + msg + "   ") if self.tty else (msg + "\n"))
        sys.stderr.flush()

    def done(self, lines, pos):
        self.tick(lines, pos, force=True)
        if self.tty:
            sys.stderr.write("\n")


class Fasta:
    """Random access to an indexed FASTA (needs reference.fasta.fai, see samtools faidx)."""

    def __init__(self, path):
        if not os.path.exists(path + ".fai"):
            sys.exit(f"ERROR: {path}.fai not found - run: samtools faidx {path}")
        self.fh = open(path, "rb")
        self.idx = {}
        with open(path + ".fai") as f:
            for line in f:
                name, length, off, lb, lw = line.split("\t")[:5]
                self.idx[name] = (int(length), int(off), int(lb), int(lw))

    def fetch(self, chrom, start, end):
        """0-based half-open."""
        length, off, lb, lw = self.idx[chrom]
        end = min(end, length)
        if start < 0 or start >= end:
            return ""
        first = off + (start // lb) * lw + start % lb
        last = off + ((end - 1) // lb) * lw + (end - 1) % lb
        self.fh.seek(first)
        return self.fh.read(last - first + 1).decode().replace("\n", "").replace("\r", "").upper()


def normalize(chrom, pos, ref, alt, fasta):
    """Trim and left-align (Tan et al. 2015). Without a reference: trim only."""
    if len(ref) == 1 and len(alt) == 1:
        return pos, ref, alt
    p, r, al = pos, ref.upper(), alt.upper()
    if fasta is not None and chrom in fasta.idx:
        changed = True
        while changed:
            changed = False
            if r and al and r[-1] == al[-1]:
                r, al, changed = r[:-1], al[:-1], True
            if (not r or not al) and p > 1:
                base = fasta.fetch(chrom, p - 2, p - 1)
                if not base:
                    break
                r, al, p, changed = base + r, base + al, p - 1, True
        if not r or not al:
            return pos, ref, alt
    else:
        while len(r) > 1 and len(al) > 1 and r[-1] == al[-1]:
            r, al = r[:-1], al[:-1]
    while len(r) > 1 and len(al) > 1 and r[0] == al[0]:
        r, al, p = r[1:], al[1:], p + 1
    return p, r, al


def vep_allele(ref, alt):
    """Allele string as VEP writes it in CSQ (first shared base trimmed for indels)."""
    if (len(ref) > 1 or len(alt) > 1) and ref[0] == alt[0]:
        return alt[1:] or "-"
    return alt


def clinvar_class(clin_sig):
    toks = {t.lower() for t in re.split(r"[&/,]", clin_sig or "") if t}
    if not toks:
        return ""
    if any(t.startswith("conflicting") for t in toks):
        return "conflicting"
    plp = bool(toks & {"pathogenic", "likely_pathogenic"})
    ben = bool(toks & {"benign", "likely_benign"})
    if plp and ben:
        return "conflicting"
    if plp:
        return "P/LP"
    if ben:
        return "benign"
    if "uncertain_significance" in toks:
        return "VUS"
    return ""


def pick_csq(csq, ix, alleles):
    """Transcript for this allele, VEP-style: MANE Select / canonical protein-coding first; a
    non-canonical protein-coding transcript only if it is more severe; NMD, retained-intron and
    non-coding transcripts only if nothing else overlaps the allele."""
    ia, iimp, imane = ix.get("Allele", -1), ix.get("IMPACT", -1), ix.get("MANE_SELECT", -1)
    ican, ibio = ix.get("CANONICAL", -1), ix.get("BIOTYPE", -1)

    def g(vals, i):
        return vals[i] if 0 <= i < len(vals) else ""

    ents = []
    for entry in csq.split(","):
        vals = entry.split("|")
        if g(vals, ia) not in alleles:
            continue
        coding = g(vals, ibio) == "protein_coding"
        t1 = bool(g(vals, imane)) or (g(vals, ican) == "YES" and coding)
        ents.append((IMPACT_RANK.get(g(vals, iimp), -1), t1, coding, vals))
    if not ents:
        return None
    
    t1s = [e for e in ents if e[1]]
    cods = [e for e in ents if e[2]]
    non_cods = [e for e in ents if not e[2]]
    
    best_t1 = max(t1s, key=lambda e: (e[0], bool(g(e[3], imane)))) if t1s else None
    if best_t1 and best_t1[0] >= IMPACT_RANK["MODERATE"]:
        pick = best_t1
    elif cods:
        best_c = max(cods, key=lambda e: (e[0], e[1]))
        pick = best_c if (best_t1 is None or best_c[0] > best_t1[0]) else best_t1
    else:
        pick = max(ents, key=lambda e: (e[0], e[1]))
    
    # Track worst entries separately to avoid silently dropping severe isoform-only hits
    worst_coding = max(cods, key=lambda e: e[0]) if cods else None
    worst_isoform = max(non_cods, key=lambda e: e[0]) if non_cods else None
    
    best = pick[3]
    out = {k: g(best, ix.get(k, -1)) for k in KEEP_CSQ}
    out["other_tx"] = ""
    out["non_canon"] = not pick[1]
    
    if worst_coding and worst_coding[0] > pick[0] and worst_coding[0] >= IMPACT_RANK["MODERATE"]:
        w = worst_coding[3]
        out["other_tx"] = (f"also {g(w, iimp)} {g(w, ix.get('Consequence', -1))} on "
                           f"{g(w, ix.get('Feature', -1))} ({g(w, ibio)})")
    
    if worst_isoform and worst_isoform[0] > pick[0] and worst_isoform[0] >= IMPACT_RANK["MODERATE"]:
        w = worst_isoform[3]
        note = f"also {g(w, iimp)} {g(w, ix.get('Consequence', -1))} on {g(w, ix.get('Feature', -1))} ({g(w, ibio)}): isoform-only effect, down-weighted"
        if out["other_tx"]:
            out["other_tx"] += "; " + note
        else:
            out["other_tx"] = note

    freq = to_float(g(best, ix.get("MAX_AF", -1)))
    if freq is None:
        fv = [to_float(g(best, ix.get(k, -1))) for k in FREQ_FALLBACK]
        fv = [v for v in fv if v is not None]
        freq = max(fv) if fv else None
    out["max_af"] = freq
    out["cv"] = clinvar_class(out["CLIN_SIG"])
    return out


class Call:
    __slots__ = ("qual", "filt", "status", "gt", "hom", "dp", "vaf", "alt_f", "alt_r",
                 "ps", "hap", "qual_pct", "alt_n", "ref_n", "tot_n", "vaf_p")


def better(new, old):
    if (new.status == "PASS") != (old.status == "PASS"):
        return new.status == "PASS"
    return (new.qual or 0) > (old.qual or 0)


# ---------------------------------------------------------------- loading ----------
def load_caller(path, label, a, fasta, stats, annot):
    """Stream one annotated VCF. Keeps compact Call objects (PASS and failed) for sites
    that are still candidates after the impact/AF pre-filter; annotation is stored once
    per normalised variant in `annot`."""
    calls, ix = {}, {}
    raw = open(path, "rb")
    fh = gzip.GzipFile(fileobj=raw) if path.endswith(".gz") else raw
    prog = Progress(label, os.path.getsize(path), a.progress_interval)
    min_rank = IMPACT_RANK[a.min_impact]
    tokens = [f"|{k}|" for k, v in IMPACT_RANK.items() if v >= min_rank] + ["athogenic"]
    n = 0
    for bline in fh:
        n += 1
        if not n & 1023:
            prog.tick(n, raw.tell())
        if bline.startswith(b"#"):
            if bline.startswith(b"##INFO=<ID=CSQ,"):
                m = re.search(r"Format:\s*([^\">]+)", bline.decode("utf-8", "replace"))
                if m:
                    ix = {name: i for i, name in enumerate(m.group(1).strip().split("|"))}
            continue
        f = bline.decode("utf-8", "replace").rstrip("\n").split("\t")
        if len(f) < 10:
            continue
        chrom, pos, _id, ref, alts, qual, filt, info = f[:8]
        stats[(label, "records read")] += 1
        if not contig_ok(chrom, a.include_contigs):
            stats[(label, "skipped: contig not chr1-22/X (alt, unplaced, chrY, chrM)")] += 1
            continue
        csq = next((x[4:] for x in info.split(";") if x.startswith("CSQ=")), None)
        if csq is None or not ix:
            stats[(label, "skipped: no CSQ annotation")] += 1
            continue
        if min_rank > 0 and not any(t in csq for t in tokens):
            stats[(label, "skipped by pre-filter (impact below --min-impact)")] += 1
            continue
        d = dict(zip(f[8].split(":"), f[9].split(":")))
        gt = d.get("GT", "./.")
        gt_alleles = re.split(r"[/|]", gt)
        pos = int(pos)
        for i, alt in enumerate(alts.split(",")):
            if alt.startswith("<") or alt in ("*", "."):
                continue
            c = pick_csq(csq, ix, {vep_allele(ref, alt), alt})
            if c is None:
                stats[(label, "skipped: no CSQ entry for allele")] += 1
                continue
            if (IMPACT_RANK.get(c["IMPACT"], -1) < min_rank and c["cv"] != "P/LP") or \
               (c["max_af"] is not None and c["max_af"] >= a.max_af and c["cv"] != "P/LP"):
                stats[(label, "skipped by pre-filter (impact/AF)")] += 1
                continue
            p2, r2, a2 = normalize(chrom, pos, ref, alt, fasta)
            if (p2, r2, a2) != (pos, ref, alt):
                stats[(label, "normalised (pos/alleles changed)")] += 1
            key = (chrom, p2, r2, a2)

            call = Call()
            call.qual, call.filt, call.gt = to_float(qual), filt, gt
            call.hom = all(g_ == str(i + 1) for g_ in gt_alleles)   # this alt on every copy (False for ./. and 0/0)
            if filt not in ("PASS", "."):
                call.status = filt
            elif str(i + 1) not in gt_alleles:
                call.status = "no_alt_GT"
            else:
                call.status = "PASS"
            ad = int_list(d.get("AD"))
            adf, adr = int_list(d.get("ADF")), int_list(d.get("ADR"))
            if ad and len(ad) > i + 1 and sum(ad) > 0:
                call.vaf = ad[i + 1] / sum(ad)
            else:
                v = (d.get("VAF") or d.get("VF") or "").split(",")
                call.vaf = to_float(v[i]) if len(v) > i else None
            call.dp = int(d["DP"]) if d.get("DP", ".").isdigit() else None
            call.alt_f = adf[i + 1] if adf and len(adf) > i + 1 else None
            call.alt_r = adr[i + 1] if adr and len(adr) > i + 1 else None
            alt_n = ad[i + 1] if ad and len(ad) > i + 1 else None
            if alt_n is None and call.alt_f is not None and call.alt_r is not None:
                alt_n = call.alt_f + call.alt_r
            call.alt_n = alt_n
            call.ref_n = ad[0] if ad else None
            call.tot_n = sum(ad) if ad and sum(ad) > 0 else call.dp
            call.vaf_p = None
            if call.status == "PASS":
                if call.dp is not None and call.dp < a.min_dp:
                    call.status = "low_dp"
                elif alt_n is not None and alt_n < a.min_alt_reads:
                    call.status = "low_alt_reads"
                elif not call.hom and call.vaf is not None and call.vaf < a.min_vaf:
                    call.status = "low_vaf"
            ps = d.get("PS", ".")
            call.ps = None if ps in (".", "") else ps
            call.hap = None
            if "|" in gt and len(gt_alleles) == 2 and call.ps:
                if gt_alleles == ["0", str(i + 1)]:
                    call.hap = 1
                elif gt_alleles == [str(i + 1), "0"]:
                    call.hap = 0
            call.qual_pct = None

            old = calls.get(key)
            if old is None or better(call, old):
                calls[key] = call
            if key not in annot:
                annot[key] = c
            stats[(label, "kept PASS" if call.status == "PASS"
                   else "kept failed (carried along)")] += 1
    prog.done(n, raw.tell())
    raw.close()
    if not ix:
        print(f"WARNING: no CSQ header in {path}; nothing could be annotated", file=sys.stderr)
    return calls


def add_qual_terms(calls, label, full):
    """QUAL -> 0..1, saturating at `full` (per caller). Not a percentile: a weak candidate set
    should not manufacture zeros, and callers' QUAL scales differ."""
    for c in calls.values():
        c.qual_pct = min(c.qual / full, 1.0) if (c.status == "PASS" and c.qual is not None) else None


# ---------------------------------------------------------------- scoring ----------
def binom_two_sided_p(k, n, p):
    """Exact two-sided binomial p-value (sum of outcomes no more likely than k)."""
    if n <= 0:
        return 1.0
    p = min(max(p, 1e-6), 1 - 1e-6)
    lg, lp, lq = math.lgamma, math.log(p), math.log(1 - p)

    def logpmf(x):
        return lg(n + 1) - lg(x + 1) - lg(n - x + 1) + x * lp + (n - x) * lq

    ref = logpmf(k) + 1e-9
    return min(1.0, sum(math.exp(l) for l in (logpmf(x) for x in range(n + 1)) if l <= ref))


def caller_confidence(c, full_dp, vaf_model="binomial"):
    terms = {}
    if c.qual_pct is not None:
        terms["qual"] = c.qual_pct
    if vaf_model == "binomial" and c.alt_n is not None:
        # For multi-allelic sites, tot_n includes other alts. 
        # Use ref + this alt as n if ref_n is available to avoid spuriously penalizing.
        n_binom = (c.ref_n + c.alt_n) if c.ref_n is not None else c.tot_n
        if n_binom > 0:
            c.vaf_p = binom_two_sided_p(c.alt_n, n_binom, HOM_EXPECTED_VAF if c.hom else 0.5)
            vaf_conf = min(1.0, c.vaf_p / VAF_P_FULL)
            # Depth-aware caution: if n is very low, cap confidence to avoid false certainty
            if n_binom < 10:
                vaf_conf = min(vaf_conf, 0.7)
            terms["vaf"] = vaf_conf
    elif vaf_model == "binomial" and c.vaf is not None:
        terms["vaf"] = max(0.0, 1 - abs(c.vaf - (1.0 if c.hom else 0.5)) / 0.5)
    if c.dp is not None:
        terms["depth"] = min(c.dp / full_dp, 1.0)
    if c.alt_f is not None and c.alt_r is not None and c.alt_f + c.alt_r > 0:
        terms["strand"] = min(1.0, 2 * min(c.alt_f, c.alt_r) / (c.alt_f + c.alt_r))
    if not terms:
        return None
    w = sum(CONF_WEIGHTS[k] for k in terms)
    return sum(CONF_WEIGHTS[k] * v for k, v in terms.items()) / w


def contig_ok(chrom, extra):
    """Primary autosomes + X only; alt/random/unplaced contigs, chrY and chrM are skipped unless added."""
    return PRIMARY_RE.match(chrom) is not None or chrom in extra


def gene_relevance(csq):
    """(factor, note): non-coding and pseudogene transcripts rarely explain a phenotype."""
    bio, pheno = csq.get("BIOTYPE", ""), csq.get("GENE_PHENO") == "1"
    if not bio:
        return 1.0, ""
    if bio == "protein_coding":
        return (1.0, "") if pheno else (CODING_NO_PHENO_REL, "no known gene-phenotype link")
    if bio in ISOFORM_ONLY_BIOTYPES:      # gene-level phenotype flags do not rescue these
        return ISOFORM_ONLY_REL, f"effect only on a {bio} isoform"
    if "pseudogene" in bio:
        return PSEUDOGENE_REL, f"{bio} (not a functional gene)"
    return (1.0, "") if pheno else (NONCODING_REL, f"non-coding transcript ({bio})")


def load_masks(specs):
    """--mask-bed NAME=PATH (BED, optionally gzipped) -> [(name, {chrom: (starts, ends)})]"""
    masks = []
    for spec in specs:
        name, _, path = spec.partition("=") if "=" in spec else (os.path.basename(spec), "", spec)
        by = defaultdict(list)
        with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
            for line in fh:
                f = line.split("\t")
                if len(f) >= 3 and not line.startswith(("#", "track", "browser")):
                    by[f[0]].append((int(f[1]), int(f[2])))
        idx = {}
        for chrom, iv in by.items():
            iv.sort()
            merged = [list(iv[0])]
            for st, en in iv[1:]:
                if st <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], en)
                else:
                    merged.append([st, en])
            idx[chrom] = ([m[0] for m in merged], [m[1] for m in merged])
        masks.append((name, idx))
    return masks


def mask_hits(masks, chrom, pos):
    hits = []
    for name, idx in masks:
        if chrom in idx:
            starts, ends = idx[chrom]
            i = bisect.bisect_right(starts, pos - 1) - 1
            if i >= 0 and pos - 1 < ends[i]:
                hits.append(name)
    return hits


def rarity_score(freq):
    f = RARITY_FLOOR if not freq else max(freq, RARITY_FLOOR)
    num = math.log10(RARITY_COMMON) - math.log10(f)
    den = math.log10(RARITY_COMMON) - math.log10(RARITY_FLOOR)
    return min(1.0, max(0.0, num / den))


def pair_status(A, B, order):
    """cis / trans if both variants are phased in the same caller and phase set."""
    for lab in order:
        ca, cb = A["passed"].get(lab), B["passed"].get(lab)
        if ca and cb and ca.ps and ca.ps == cb.ps and ca.hap is not None and cb.hap is not None:
            return "trans" if ca.hap != cb.hap else "cis"
    return "unphased"


def flag_clusters(items, window, min_n):
    """Dense groups of candidates (clustered indels/SNVs) are typical of repeat/mapping artefacts."""
    if min_n <= 0:
        return 0
    by_chrom = defaultdict(list)
    for it in items:
        by_chrom[it["row"]["chrom"]].append((it["row"]["pos"], it))
    n = 0
    for lst in by_chrom.values():
        lst.sort(key=lambda x: x[0])
        pos = [p for p, _ in lst]
        for p, it in lst:
            near = bisect.bisect_right(pos, p + window) - bisect.bisect_left(pos, p - window)
            if near >= min_n:
                it["clustered"] = True
                it["score"] *= CLUSTER_FACTOR
                it["why"].append(f"{near - 1} other candidate variants within {window} bp "
                                 "(repeat/mapping artefact?)")
                n += 1
    return n


def compound_het(items, order):
    by_gene = defaultdict(list)
    for it in items:
        if it["imp"] >= COMPHET_MIN_IMPACT_SCORE and it["rel"] >= CODING_NO_PHENO_REL \
                and (it["af"] is None or it["af"] < COMPHET_MAX_AF) \
                and it["gene"] and any(not c.hom for c in it["passed"].values()):
            by_gene[it["gene"]].append(it)
    n_trans = 0
    for gene, lst in by_gene.items():
        for it in lst:
            it["row"]["gene_n_rare"] = len(lst)
        if len(lst) < 2 or len(lst) > COMPHET_MAX_GENE_N:
            continue
        stat = defaultdict(set)
        for A, B in combinations(lst, 2):
            s = pair_status(A, B, order)
            stat[id(A)].add((s, B["row"]["locus"]))
            stat[id(B)].add((s, A["row"]["locus"]))
        for it in lst:
            kinds = {s for s, _ in stat[id(it)]}
            label = "trans" if "trans" in kinds else "unphased" if "unphased" in kinds else "cis"
            it["row"]["compound_het"] = label
            if label == "trans":
                partners = sorted(loc for s, loc in stat[id(it)] if s == "trans")[:3]
                it["why"].append("compound-het in trans with " + ",".join(partners))
                if not it.get("clustered"):
                    it["score"] *= COMPHET_TRANS_BONUS
                n_trans += 1
            elif label == "unphased":
                it["why"].append("2+ rare variants in gene (phase unknown)")
    return n_trans


# ---------------------------------------------------------------- main -------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--deepvariant", help="VEP-annotated DeepVariant VCF(.gz)")
    ap.add_argument("--nanocaller", help="VEP-annotated NanoCaller VCF(.gz)")
    ap.add_argument("--reference", help="indexed reference FASTA (.fai) for left-aligning indels")
    ap.add_argument("--output", required=True, help="shortlist TSV")
    ap.add_argument("--output-all", help="optional TSV with every scored variant")
    ap.add_argument("--top-n", type=int, default=100)
    ap.add_argument("--full-credit-dp", type=int, default=20,
                    help="depth at which the depth term reaches 1.0 (default 20)")
    ap.add_argument("--min-impact", default="LOW", choices=list(IMPACT_RANK),
                    help="skip records whose most severe consequence is below this "
                         "(ClinVar P/LP always kept). Default LOW")
    ap.add_argument("--max-af", type=float, default=RARITY_COMMON,
                    help="skip records with MAX_AF at/above this (default 0.05; rarity is 0 above it anyway)")
    ap.add_argument("--min-dp", type=int, default=8, help="calls below this depth do not count as PASS (8; 0 = off)")
    ap.add_argument("--min-alt-reads", type=int, default=4,
                    help="calls with fewer alt reads do not count as PASS (4; 0 = off)")
    ap.add_argument("--min-vaf", type=float, default=0.2,
                    help="het calls below this VAF do not count as PASS (0.2; 0 = off)")
    ap.add_argument("--vaf-model", choices=["binomial", "off"], default="binomial",
                    help="VAF-vs-genotype term: exact binomial consistency (default) or off "
                         "(use off for tumour/mosaic samples)")
    ap.add_argument("--qual-full-dv", type=float, default=QUAL_FULL["dv"], help="DeepVariant QUAL that saturates the QUAL term (30)")
    ap.add_argument("--qual-full-nc", type=float, default=QUAL_FULL["nc"], help="NanoCaller QUAL that saturates the QUAL term (30)")
    ap.add_argument("--cluster-window", type=int, default=2000,
                    help="window (bp) for the dense-cluster check (2000)")
    ap.add_argument("--cluster-min", type=int, default=3,
                    help="candidates within the window (incl. itself) that make a cluster (3; 0 = off)")
    ap.add_argument("--include-contigs", default="",
                    help="comma list of extra contigs to keep, e.g. chrY,chrM (default: chr1-22 and chrX only)")
    ap.add_argument("--mask-bed", action="append", default=[], metavar="NAME=PATH",
                    help="BED of problematic regions (segdups, low mappability...); variants inside get "
                         "confidence x0.5 (repeatable)")
    ap.add_argument("--progress-interval", type=float, default=None,
                    help="seconds between progress lines (default 1 on a terminal, 30 in logs)")
    a = ap.parse_args()
    a.include_contigs = {c for c in a.include_contigs.split(",") if c}
    a.qual_full = {"dv": a.qual_full_dv, "nc": a.qual_full_nc}
    masks = load_masks(a.mask_bed)

    inputs = [(lab, p) for lab, p in (("dv", a.deepvariant), ("nc", a.nanocaller)) if p]
    if not inputs:
        ap.error("give at least one of --deepvariant / --nanocaller")
    fasta = Fasta(a.reference) if a.reference else None
    if fasta is None:
        print("NOTE: no --reference given; indels are only trimmed, not left-aligned, "
              "so the same indel written differently by two callers will not merge.",
              file=sys.stderr)

    stats, callers, annot = Counter(), {}, {}
    for lab, path in inputs:
        print(f"[{lab}] reading {path}", file=sys.stderr)
        callers[lab] = load_caller(path, lab, a, fasta, stats, annot)
        add_qual_terms(callers[lab], lab, a.qual_full.get(lab, 30.0))
    order = list(callers)

    keys = set().union(*callers.values())
    print(f"[merge] scoring {len(keys):,} candidate sites", file=sys.stderr)
    items = []
    for key in keys:
        chrom, pos, ref, alt = key
        present = {lab: c[key] for lab, c in callers.items() if key in c}
        passed = {lab: c for lab, c in present.items() if c.status == "PASS"}
        if not passed:
            stats[("merge", "dropped: only failed calls, no PASS")] += 1
            continue
        csq = annot.get(key, {})
        impact, cv, af = csq.get("IMPACT", ""), csq.get("cv", ""), csq.get("max_af")
        imp = IMPACT_SCORE.get(impact, 0.0)
        if cv == "P/LP":
            imp = max(imp, CLINVAR_PLP_SCORE)
        if cv == "benign":
            imp *= CLINVAR_BENIGN_FACTOR
        rar = rarity_score(af)
        rar_note = ""
        if cv == "P/LP" and af is not None and af >= a.max_af:
            rar = max(rar, PLP_COMMON_RARITY_FLOOR)
            rar_note = f"ClinVar P/LP but common (AF {af:.2g}): check review status/penetrance"
        elif not af and re.search(r"(^|&)rs\d+", csq.get("Existing_variation", "")):
            rar = min(rar, RSID_NO_AF_RARITY)
            rar_note = "has an rsID but no population frequency"
        rel, rel_note = gene_relevance(csq)
        confs = {lab: caller_confidence(c, a.full_credit_dp, a.vaf_model) for lab, c in passed.items()}
        vals = [v for v in confs.values() if v is not None]
        conf = sum(vals) / len(vals) if vals else 0.0
        single = len(callers) > 1 and len(passed) == 1
        if single:
            conf *= SINGLE_CALLER_FACTOR
        sym = csq.get("SYMBOL", "")
        prob_note = ""
        if ARTIFACT_GENE_RE.match(sym):
            conf *= PROBLEM_GENE_FACTOR
            prob_note = f"{sym}: repeat-rich/paralogous gene family, calls are often mapping artefacts"
        elif sym in PARALOG_GENES:
            prob_note = f"{sym} has close paralogs/pseudogenes: confirm (e.g. paraphase)"
        mhits = mask_hits(masks, chrom, pos) if masks else []
        if mhits:
            conf *= MASK_FACTOR

        gts = {frozenset(re.split(r"[/|]", c.gt)) for c in passed.values()}
        gt_match = "" if len(passed) < 2 else ("yes" if len(gts) == 1 else "NO")

        why = []
        if impact:
            why.append(f"{impact} {csq.get('Consequence', '')}")
        if cv:
            why.append(f"ClinVar {cv}")
        if rel_note:
            why.append(rel_note)
        if csq.get("non_canon") and csq.get("IMPACT") in ("HIGH", "MODERATE"):
            why.append("effect only on a non-canonical transcript")
        if csq.get("other_tx"):
            why.append("more severe on another transcript: " + csq["other_tx"])
        if prob_note:
            why.append(prob_note)
        why.append(rar_note or ("absent from gnomAD/1000G" if af is None else f"max AF {af:.2g}"))
        if "COSV" in csq.get("Existing_variation", ""):
            why.append("COSMIC-listed (somatic hotspot? check sample type/mosaicism)")
        if mhits:
            why.append("in problematic region: " + ",".join(mhits))
        why.append("called by " + "+".join(l.upper() for l in passed) + (" (only one)" if single else ""))
        if gt_match == "NO":
            why.append("genotype disagrees between callers")
        if any(c.hom for c in passed.values()) and imp >= COMPHET_MIN_IMPACT_SCORE and rel >= CODING_NO_PHENO_REL \
                and (af is None or af < COMPHET_MAX_AF):
            why.append("homozygous (recessive candidate)")
        for lab, c in present.items():
            L = lab.upper()
            if c.status != "PASS":
                why.append(f"{L} {c.status} (GT {c.gt}, DP {c.dp}, VAF {fmt(c.vaf, 2)})")
                continue
            if c.dp is not None and c.dp < a.full_credit_dp / 2:
                why.append(f"low depth {L} (DP {c.dp})")
            if c.alt_f is not None and c.alt_r is not None:
                t = c.alt_f + c.alt_r
                if t >= 4 and min(c.alt_f, c.alt_r) / t < 0.1:
                    why.append(f"strand-biased alt reads {L}")
            if c.vaf_p is not None and c.vaf_p < VAF_P_FULL:
                obs_vaf = c.vaf if c.vaf is not None else (c.alt_n / c.tot_n if c.tot_n else None)
                vaf_str = f"{obs_vaf:.2f}" if obs_vaf is not None else "NA"
                why.append(f"VAF {vaf_str} atypical for genotype {c.gt} on {L} (binomial p={c.vaf_p:.2g})")

        row = {"locus": f"{chrom}:{pos}", "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
               "symbol": csq.get("SYMBOL") or csq.get("Gene", ""),
               "transcript": csq.get("Feature", ""), "biotype": csq.get("BIOTYPE", ""),
               "consequence": csq.get("Consequence", ""), "impact": impact,
               "hgvsc": csq.get("HGVSc", ""), "hgvsp": csq.get("HGVSp", ""),
               "clinvar": csq.get("CLIN_SIG", ""), "sift": csq.get("SIFT", ""),
               "polyphen": csq.get("PolyPhen", ""), "gene_pheno": csq.get("GENE_PHENO", ""),
               "existing_variation": csq.get("Existing_variation", ""),
               "max_af": fmt(af), "called_by": "+".join(l.upper() for l in passed),
               "gt_match": gt_match, "gene_n_rare": "", "compound_het": "",
               "gene_rel": fmt(rel, 2), "masks": ",".join(mhits),
               "cosmic": "yes" if "COSV" in csq.get("Existing_variation", "") else "",
               "impact_score": fmt(imp, 2), "rarity_score": fmt(rar, 2), "conf_score": fmt(conf, 2)}
        for lab in callers:
            c = present.get(lab)
            row[f"{lab}_status"] = c.status if c else "not_seen"
            row[f"{lab}_qual"] = fmt(c.qual) if c else ""
            row[f"{lab}_gt"] = c.gt if c else ""
            row[f"{lab}_dp"] = fmt(c.dp) if c else ""
            row[f"{lab}_vaf"] = fmt(c.vaf, 2) if c else ""
            row[f"{lab}_alt_fwd_rev"] = f"{c.alt_f}/{c.alt_r}" if c and c.alt_f is not None else ""
            row[f"{lab}_conf"] = fmt(confs.get(lab), 2)
        ps = next((c.ps for c in passed.values() if c.ps), "")
        row["phase_set"] = ps or ""
        items.append({"row": row, "passed": passed, "imp": imp, "af": af, "why": why, "rel": rel,
                      "gene": csq.get("Gene") or csq.get("SYMBOL", ""),
                      "score": imp * rel * rar * conf})

    n_clustered = flag_clusters(items, a.cluster_window, a.cluster_min)
    n_trans = compound_het(items, order)
    for it in items:
        it["row"]["reasons"] = "; ".join(it["why"])
        it["row"]["_score"] = it["score"]
        it["row"]["score"] = fmt(it["score"], 3)
    rows = [it["row"] for it in items]
    rows.sort(key=lambda r: (-r["_score"], r["called_by"] != "DV+NC", r["chrom"], r["pos"]))
    for i, r in enumerate(rows, 1):
        r["rank"] = i

    cols = ["rank", "locus", "chrom", "pos", "ref", "alt", "symbol", "transcript", "biotype",
            "consequence", "impact", "hgvsc", "hgvsp", "clinvar", "sift", "polyphen",
            "gene_pheno", "existing_variation", "max_af", "called_by", "gt_match"]
    for lab in callers:
        cols += [f"{lab}_{x}" for x in ("status", "qual", "gt", "dp", "vaf", "alt_fwd_rev", "conf")]
    cols += ["phase_set", "gene_n_rare", "compound_het", "gene_rel", "masks", "cosmic", "impact_score", "rarity_score",
             "conf_score", "score", "reasons"]

    def write(path, subset):
        with open(path, "w", newline="") as out:
            w = csv.DictWriter(out, fieldnames=cols, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(subset)

    shortlist = [r for r in rows if r["_score"] > 0][: a.top_n]
    write(a.output, shortlist)
    if a.output_all:
        write(a.output_all, rows)

    print("\n=== summary ===", file=sys.stderr)
    for (lab, what), n in sorted(stats.items()):
        print(f"[{lab}] {what}: {n:,}", file=sys.stderr)
    print(f"candidate sites scored: {len(rows):,}; compound-het (trans) variants: {n_trans}; "
          f"down-weighted as dense clusters: {n_clustered}",
          file=sys.stderr)
    print("shortlist by impact:", dict(Counter(r["impact"] for r in shortlist)), file=sys.stderr)
    print("shortlist by caller:", dict(Counter(r["called_by"] for r in shortlist)), file=sys.stderr)
    print(f"wrote {len(shortlist)} rows -> {a.output}; peak RAM {peak_ram_gb():.1f} GB", file=sys.stderr)


if __name__ == "__main__":
    main()