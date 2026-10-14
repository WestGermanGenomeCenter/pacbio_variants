#!/usr/bin/env python3
"""
prioritize_snvs.py (v2) - shortlist of plausible SNVs/indels from VEP-annotated
DeepVariant and/or NanoCaller VCFs. Stdlib only.

Per variant:  score = impact * rarity * confidence * bonuses
  impact      VEP IMPACT of the most severe transcript (ClinVar P/LP overrides,
              ClinVar benign is down-weighted)
  rarity      from MAX_AF (fallback gnomADg/gnomADe/1000G AF); absent = rare
  confidence  mean over callers that PASS: QUAL percentile, VAF vs genotype, depth,
              strand balance (if ADF/ADR exist); x0.75 if only one of two callers
  bonus       x1.25 if a rare high/moderate variant has a partner in the same gene
              on the other haplotype (phased, same phase set, "in trans")

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
CLINVAR_PLP_SCORE = 1.2
CLINVAR_BENIGN_FACTOR = 0.1
COMPHET_TRANS_BONUS = 1.25
COMPHET_MAX_AF = 0.01          # partners must be rarer than this
COMPHET_MIN_IMPACT_SCORE = 0.6 # ... and at least MODERATE (or ClinVar P/LP)
COMPHET_MAX_GENE_N = 300       # skip pairing in genes with more candidates
FREQ_FALLBACK = ("gnomADg_AF", "gnomADe_AF", "AF")
CONF_WEIGHTS = {"qual": 0.35, "vaf": 0.25, "depth": 0.25, "strand": 0.15}
SINGLE_CALLER_FACTOR = 0.75
RARITY_COMMON = 0.05   # AF at/above this -> rarity 0
RARITY_FLOOR = 1e-5    # AF at/below this -> rarity 1
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
    """Progress by bytes of the (compressed) input file; \\r on a terminal, lines in logs."""

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
    """Most severe transcript for this allele (tie-break: MANE Select, canonical, coding)."""
    ia, iimp, imane = ix.get("Allele", -1), ix.get("IMPACT", -1), ix.get("MANE_SELECT", -1)
    ican, ibio = ix.get("CANONICAL", -1), ix.get("BIOTYPE", -1)

    def g(vals, i):
        return vals[i] if 0 <= i < len(vals) else ""

    best, best_key = None, None
    for entry in csq.split(","):
        vals = entry.split("|")
        if g(vals, ia) not in alleles:
            continue
        key = (IMPACT_RANK.get(g(vals, iimp), -1), bool(g(vals, imane)),
               g(vals, ican) == "YES", g(vals, ibio) == "protein_coding")
        if best_key is None or key > best_key:
            best, best_key = vals, key
    if best is None:
        return None
    out = {k: g(best, ix.get(k, -1)) for k in KEEP_CSQ}
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
                 "ps", "hap", "qual_pct")


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
               (c["max_af"] is not None and c["max_af"] >= a.max_af):
                stats[(label, "skipped by pre-filter (impact/AF)")] += 1
                continue
            p2, r2, a2 = normalize(chrom, pos, ref, alt, fasta)
            if (p2, r2, a2) != (pos, ref, alt):
                stats[(label, "normalised (pos/alleles changed)")] += 1
            key = (chrom, p2, r2, a2)

            call = Call()
            call.qual, call.filt, call.gt = to_float(qual), filt, gt
            call.hom = len(set(gt_alleles)) == 1
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


def add_qual_percentiles(calls):
    passed = [c for c in calls.values() if c.status == "PASS"]
    quals = sorted(c.qual for c in passed if c.qual is not None)
    for c in passed:
        c.qual_pct = (bisect.bisect_left(quals, c.qual) / (len(quals) - 1)
                      if c.qual is not None and len(quals) > 1 else None)


# ---------------------------------------------------------------- scoring ----------
def caller_confidence(c, full_dp):
    terms = {}
    if c.qual_pct is not None:
        terms["qual"] = c.qual_pct
    if c.vaf is not None:
        terms["vaf"] = max(0.0, 1 - abs(c.vaf - (1.0 if c.hom else 0.5)) / 0.5)
    if c.dp is not None:
        terms["depth"] = min(c.dp / full_dp, 1.0)
    if c.alt_f is not None and c.alt_r is not None and c.alt_f + c.alt_r > 0:
        terms["strand"] = min(1.0, 2 * min(c.alt_f, c.alt_r) / (c.alt_f + c.alt_r))
    if not terms:
        return None
    w = sum(CONF_WEIGHTS[k] for k in terms)
    return sum(CONF_WEIGHTS[k] * v for k, v in terms.items()) / w


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


def compound_het(items, order):
    by_gene = defaultdict(list)
    for it in items:
        if it["imp"] >= COMPHET_MIN_IMPACT_SCORE and (it["af"] is None or it["af"] < COMPHET_MAX_AF) \
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
    ap.add_argument("--progress-interval", type=float, default=None,
                    help="seconds between progress lines (default 1 on a terminal, 30 in logs)")
    a = ap.parse_args()

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
        add_qual_percentiles(callers[lab])
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
        imp = CLINVAR_PLP_SCORE if cv == "P/LP" else IMPACT_SCORE.get(impact, 0.0)
        if cv == "benign":
            imp *= CLINVAR_BENIGN_FACTOR
        rar = rarity_score(af)
        confs = {lab: caller_confidence(c, a.full_credit_dp) for lab, c in passed.items()}
        vals = [v for v in confs.values() if v is not None]
        conf = sum(vals) / len(vals) if vals else 0.0
        single = len(callers) > 1 and len(passed) == 1
        if single:
            conf *= SINGLE_CALLER_FACTOR

        gts = {frozenset(re.split(r"[/|]", c.gt)) for c in passed.values()}
        gt_match = "" if len(passed) < 2 else ("yes" if len(gts) == 1 else "NO")

        why = []
        if impact:
            why.append(f"{impact} {csq.get('Consequence', '')}")
        if cv:
            why.append(f"ClinVar {cv}")
        why.append("absent from gnomAD/1000G" if af is None else f"max AF {af:.2g}")
        why.append("called by " + "+".join(l.upper() for l in passed) + (" (only one)" if single else ""))
        if gt_match == "NO":
            why.append("genotype disagrees between callers")
        if any(c.hom for c in passed.values()) and imp >= COMPHET_MIN_IMPACT_SCORE \
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
            if c.vaf is not None and abs(c.vaf - (1.0 if c.hom else 0.5)) > 0.25:
                why.append(f"VAF off expectation {L} ({c.vaf:.2f})")

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
        items.append({"row": row, "passed": passed, "imp": imp, "af": af, "why": why,
                      "gene": csq.get("Gene") or csq.get("SYMBOL", ""),
                      "score": imp * rar * conf})

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
    cols += ["phase_set", "gene_n_rare", "compound_het", "impact_score", "rarity_score",
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
    print(f"candidate sites scored: {len(rows):,}; compound-het (trans) variants: {n_trans}",
          file=sys.stderr)
    print("shortlist by impact:", dict(Counter(r["impact"] for r in shortlist)), file=sys.stderr)
    print("shortlist by caller:", dict(Counter(r["called_by"] for r in shortlist)), file=sys.stderr)
    print(f"wrote {len(shortlist)} rows -> {a.output}; peak RAM {peak_ram_gb():.1f} GB", file=sys.stderr)


if __name__ == "__main__":
    main()