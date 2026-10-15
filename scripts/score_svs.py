#!/usr/bin/env python3
"""
prioritize_svs.py - shortlist of plausible SVs from AnnotSV TSVs of Sniffles and/or
Sawfish calls (same structure as prioritize_snvs_v2.py). Stdlib only.

Per SV:  score = impact * rarity * confidence
  impact      max of: gene effect x gene relevance (CDS/exon overlap by SV type, weighted
              by OMIM/GenCC/HI/pLI/LOEUF), known pathogenic overlap (AnnotSV P_*_source),
              AnnotSV ranking score; AnnotSV class 1-2 (score <= -0.9) is down-weighted
  rarity      max AF of AnnotSV benign-SV overlaps (B_*_AFmax, matching SV type) and any
              SVDB/CoLoRSdb-style INFO keys; absent = rare
  confidence  mean over callers that PASS: QUAL percentile, alt-read support, VAF vs
              genotype, depth, PRECISE flag; x0.8 in segdups/gaps, x0.6 in blacklist;
              x0.75 if only one of two supplied callers has it

Sniffles and Sawfish SVs are matched by type, reciprocal overlap (DEL/DUP/INV) or
breakpoint distance + size ratio (INS/BND). Failed (non-PASS) records are carried along.

Usage:
  prioritize_svs.py --sniffles sn.annotated.tsv --sawfish sf.annotated.tsv \
      --output sv_shortlist.tsv --output-all sv_all_scored.tsv --top-n 100
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

# ---- tunable constants (placeholders, to be calibrated) -------------------------
CONF_WEIGHTS = {"qual": 0.20, "support": 0.20, "vaf": 0.15, "depth": 0.10, "precise": 0.10,
                "strand": 0.15, "covchg": 0.10, "cn": 0.15}
SINGLE_CALLER_FACTOR = 0.75
SEGDUP_FACTOR, GAP_FACTOR, BLACKLIST_FACTOR = 0.8, 0.8, 0.6
KNOWN_PATH_STRONG, KNOWN_PATH_PARTIAL = 0.8, 0.5
ANNOTSV_BENIGN_FACTOR = 0.1
DEPTH_ANOMALY_FACTOR = 0.5   # flank depth far above the sample median: collapsed repeat / segdup?
FAILED_DEPTH_FACTOR = 0.7    # the other caller filtered this locus for excessive depth
LARGE_SV_FACTOR = 0.5        # arm-scale events: aneuploidy / mosaicism / segmentation artefact?
RECESSIVE_HET_FACTOR = 0.5   # gene relevance x this for a het call in a recessive-only gene
RARITY_COMMON = 0.05
RARITY_FLOOR = 1e-4
GENCC_SCORE = {"definitive": 1.0, "strong": 0.9, "moderate": 0.7, "supportive": 0.5,
               "limited": 0.3, "disputed": 0.1}
KIND_P = {"DEL": ["loss"], "DUP": ["gain"], "INS": ["ins"], "CNV": ["loss", "gain"]}
KIND_B = {"DEL": ["loss"], "DUP": ["gain"], "INS": ["ins"], "INV": ["inv"], "CNV": ["loss", "gain"]}

FULL_COLS = ["CytoBand", "Gene_name", "Gene_count", "RE_gene", "P_gain_source", "P_gain_phen",
             "P_loss_source", "P_loss_phen", "P_ins_source", "P_ins_phen", "po_P_gain_source",
             "po_P_gain_percent", "po_P_loss_source", "po_P_loss_percent", "B_gain_AFmax", "B_loss_AFmax", "B_ins_AFmax", "B_inv_AFmax",
             "SegDup_left", "SegDup_right", "Gap_left", "Gap_right", "ENCODE_blacklist_left",
             "ENCODE_blacklist_right", "AnnotSV_ranking_score", "AnnotSV_ranking_criteria",
             "ACMG_class"]
SPLIT_COLS = ["Gene_name", "Tx", "Overlapped_CDS_length", "Overlapped_CDS_percent", "Frameshift",
              "Exon_count", "Location", "Location2", "Dist_nearest_SS", "HI", "TS", "GnomAD_pLI",
              "LOEUF_bin", "GenCC_classification", "GenCC_moi", "OMIM_inheritance", "OMIM_morbid",
              "OMIM_phenotype"]
# AnnotSV 3.5 column order; used to repair header names that lost their tab (e.g. a
# header token like "B_loss_sourceB_loss_coord") and to sanity-check the column count.
KNOWN_COLS = set("""AnnotSV_ID SV_chrom SV_start SV_end SV_length SV_type Samples_ID ID REF ALT QUAL
FILTER INFO FORMAT Annotation_mode CytoBand Gene_name Closest_left Closest_right Gene_count Tx
Tx_version Tx_start Tx_end Overlapped_tx_length Overlapped_CDS_length Overlapped_CDS_percent
Frameshift Exon_count Location Location2 Dist_nearest_SS Nearest_SS_type Intersect_start
Intersect_end RE_gene P_gain_phen P_gain_hpo P_gain_source P_gain_coord P_loss_phen P_loss_hpo
P_loss_source P_loss_coord P_ins_phen P_ins_hpo P_ins_source P_ins_coord po_P_gain_phen
po_P_gain_hpo po_P_gain_source po_P_gain_coord po_P_gain_percent po_P_loss_phen po_P_loss_hpo
po_P_loss_source po_P_loss_coord po_P_loss_percent P_snvindel_nb P_snvindel_phen B_gain_source
B_gain_coord B_gain_AFmax B_loss_source B_loss_coord B_loss_AFmax B_ins_source B_ins_coord
B_ins_AFmax B_inv_source B_inv_coord B_inv_AFmax po_B_gain_allG_source po_B_gain_allG_coord
po_B_gain_someG_source po_B_gain_someG_coord po_B_loss_allG_source po_B_loss_allG_coord
po_B_loss_someG_source po_B_loss_someG_coord GC_content_left GC_content_right Repeat_coord_left
Repeat_type_left Repeat_coord_right Repeat_type_right Gap_left Gap_right SegDup_left SegDup_right
ENCODE_blacklist_left ENCODE_blacklist_characteristics_left ENCODE_blacklist_right
ENCODE_blacklist_characteristics_right ACMG HI TS DDD_HI_percent ExAC_delZ ExAC_dupZ ExAC_cnvZ
ExAC_synZ ExAC_misZ GenCC_disease GenCC_moi GenCC_classification GenCC_pmid NCBI_gene_ID
HGNC_gene_ID OMIM_ID OMIM_phenotype OMIM_inheritance OMIM_morbid OMIM_morbid_candidate LOEUF_bin
GnomAD_pLI ExAC_pLI AnnotSV_ranking_score AnnotSV_ranking_criteria ACMG_class""".split())
MATE_RE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")


# ---------------------------------------------------------------- helpers ----------
def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def to_int(x):
    v = num(x)
    return int(v) if v is not None and v == v else None


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


def short(s, n):
    s = s or ""
    return s if len(s) <= n else s[: n - 1] + "…"


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
               f"{lines / max(el, 1e-9) / 1000:.1f}k lines/s | ETA {fmt_time(eta)} | "
               f"peak RAM {peak_ram_gb():.1f} GB")
        sys.stderr.write(("\r" + msg + "   ") if self.tty else (msg + "\n"))
        sys.stderr.flush()

    def done(self, lines, pos):
        self.tick(lines, pos, force=True)
        if self.tty:
            sys.stderr.write("\n")


def parse_info(s):
    d = {}
    for kv in s.split(";"):
        if kv:
            k, sep, v = kv.partition("=")
            d[k] = v if sep else True
    return d


def norm_type(t, alt):
    t = (t or "").upper()
    if not t:
        m = re.match(r"<([^>]+)>", alt or "")
        t = m.group(1).upper() if m else ""
    t = t.split(":")[0]
    if t == "CN0":
        return "DEL"
    if re.fullmatch(r"CN\d+", t):
        return "DUP" if int(t[2:]) > 2 else "DEL"
    return "BND" if t in ("TRA", "BND") else t


def repair_header(cols):
    """Split header tokens that are two known AnnotSV names glued together."""
    out, fixed = [], []
    for c in cols:
        if c in KNOWN_COLS or len(c) < 8:
            out.append(c)
            continue
        for i in range(1, len(c)):
            if c[:i] in KNOWN_COLS and c[i:] in KNOWN_COLS:
                out += [c[:i], c[i:]]
                fixed.append(c)
                break
        else:
            out.append(c)
    return out, fixed


def af_values(s):
    """AnnotSV AFmax cell -> list of fractions (AFmax is a fraction, e.g. 0.4511; values > 1 are
    treated as percent just in case)."""
    out = []
    for tok in re.split(r"[;,|]", s or ""):
        v = num(tok)
        if v is not None:
            out.append(v / 100.0 if v > 1 else v)
    return out


# ---------------------------------------------------------------- gene logic -------
def real_sources(s):
    """Pathogenic-SV sources without the gene-level 'morbid:GENE' entries (handled by the gene logic)."""
    return [t for t in (s or "").split(";") if t and not t.lower().startswith("morbid:")]


def pathogenic_cover(sources, percents):
    """Max % of a (non-morbid) known pathogenic SV covered by the query SV."""
    toks, pcts = (sources or "").split(";"), [num(t) for t in (percents or "").split(";")]
    if len(toks) == len(pcts):
        vals = [p for t, p in zip(toks, pcts) if p is not None and not t.lower().startswith("morbid:")]
    else:
        vals = [p for p in pcts if p is not None] if real_sources(sources) else []
    return max(vals) if vals else None


def max_number(s):
    vals = [v for v in (num(t) for t in re.split(r"[;,|]", s or "")) if v is not None]
    return max(vals) if vals else None


def gene_effect(svtype, s):
    """How likely the SV disrupts this gene/transcript (0..1) from AnnotSV split-line fields."""
    cds = to_int(s.get("Overlapped_CDS_length")) or 0
    pct = num(s.get("Overlapped_CDS_percent")) or 0.0
    loc, loc2 = (s.get("Location") or "").lower(), s.get("Location2") or ""
    # Location2 "CDS" alone only says the intron lies within the CDS span; point SVs
    # (INS/BND) need an exon in Location to hit coding sequence
    in_cds = cds > 0 or ("exon" in loc and "CDS" in loc2)
    in_exon = in_cds or "exon" in loc
    fs = (s.get("Frameshift") or "").lower()
    dss = to_int(s.get("Dist_nearest_SS"))
    near_ss = dss is not None and dss <= 20
    whole_tx = "txstart-txend" in loc
    if svtype in ("DEL", "CNV"):
        if in_cds:
            if whole_tx and pct >= 99.9:
                return 1.0, "whole-gene deletion"
            return (1.0 if fs == "yes" or pct >= 50 else 0.8), "CDS deletion"
        if in_exon or whole_tx:
            return 0.35, "UTR/non-coding exon or whole non-coding transcript deletion"
        if near_ss:
            return 0.3, "deletion near splice site"
        return (0.1, "intronic deletion") if loc else (0.0, "")
    if svtype == "DUP":
        if in_cds:
            return (0.3, "whole-gene duplication") if pct >= 99.9 else (0.7, "partial CDS duplication")
        if in_exon:
            return 0.3, "UTR/non-coding exon duplication"
        return (0.05, "intronic duplication") if loc else (0.0, "")
    if svtype == "INS":
        if in_cds:
            return 0.9, "insertion in CDS"
        if in_exon:
            return 0.35, "insertion in UTR/non-coding exon"
        return (0.3, "insertion near splice site") if near_ss else ((0.1, "intronic insertion") if loc else (0.0, ""))
    if svtype in ("INV", "BND"):
        if in_cds:
            return 0.7, f"{svtype} breakpoint in CDS"
        if in_exon:
            return 0.3, f"{svtype} breakpoint in exon"
        return (0.15, f"{svtype} breakpoint in intron") if loc else (0.0, "")
    return 0.0, ""


def recessive_only(s):
    """True if GenCC/OMIM list only recessive inheritance for this gene (AR, XLR)."""
    for key in ("GenCC_moi", "OMIM_inheritance"):
        toks = {t.strip().upper() for t in re.split(r"[;,|/]", s.get(key) or "") if t.strip()}
        if toks:
            return bool(toks & {"AR", "XLR"}) and not (toks & {"AD", "XLD", "XL", "SD", "SMU", "DD", "IC"})
    return False


def gene_disease(s):
    """Gene relevance 0..1 from OMIM/GenCC/haploinsufficiency/constraint."""
    vals = []
    if (s.get("OMIM_morbid") or "").lower() == "yes":
        vals.append(0.8)
    for tok in re.split(r"[;,]", (s.get("GenCC_classification") or "").lower()):
        for k, v in GENCC_SCORE.items():
            if tok.strip().startswith(k):
                vals.append(v)
    hi = to_int(s.get("HI"))
    if hi == 3:
        vals.append(1.0)
    elif hi == 2:
        vals.append(0.6)
    pli, loeuf = num(s.get("GnomAD_pLI")), to_int(s.get("LOEUF_bin"))
    if pli is not None and pli >= 0.9:
        vals.append(0.5)
    if loeuf is not None and loeuf <= 1:
        vals.append(0.5)
    return max(vals) if vals else 0.0


# ---------------------------------------------------------------- loading ----------
class SVCall:
    __slots__ = ("label", "aid", "chrom", "start", "end", "svtype", "svlen", "mate", "id", "qual",
                 "filt", "status", "gt", "hom", "alt", "ref", "vaf", "dp", "precise", "strand", "ps", "hap",
                 "cn", "cn_term", "claim", "cov", "cov_flank", "cov_ratio", "covchg", "dgenes", "depth_x",
                 "qual_pct", "info", "full", "best", "best_val")


def build_call(f, G, sidx, label, aid, min_cov_len):
    c = SVCall()
    c.label, c.aid = label, aid
    c.chrom = G(f, "SV_chrom")
    c.start = to_int(G(f, "SV_start")) or 0
    c.end = to_int(G(f, "SV_end")) or c.start
    info = parse_info(G(f, "INFO"))
    ref, alt = G(f, "REF"), G(f, "ALT")
    c.svtype = norm_type(G(f, "SV_type") or str(info.get("SVTYPE", "")), alt)
    ln = to_int(G(f, "SV_length"))
    if ln is None:
        ln = to_int(str(info.get("SVLEN", "")).split(",")[0])
    if ln is None and c.svtype == "INS" and alt and not alt.startswith("<"):
        ln = len(alt) - len(ref)
    if ln is None and c.svtype not in ("INS", "BND"):
        ln = c.end - c.start
    c.svlen = abs(ln) if ln is not None and c.svtype != "BND" else None
    m = MATE_RE.search(alt or "")
    c.mate = (m.group(1), int(m.group(2))) if m else None
    c.id, c.qual, c.filt = G(f, "ID"), num(G(f, "QUAL")), G(f, "FILTER")

    d = dict(zip(G(f, "FORMAT").split(":"), f[sidx].split(":") if sidx < len(f) else []))
    c.gt = d.get("GT", "./.")
    gts = re.split(r"[/|]", c.gt)
    c.hom = len(set(gts)) == 1
    dv, dr = to_int(d.get("DV")), to_int(d.get("DR"))
    if dv is not None:
        c.alt, c.ref = dv, dr
    else:
        ad = int_list(d.get("AD"))
        if ad and len(ad) >= 2:
            c.alt, c.ref = sum(ad[1:]), ad[0]
        else:
            c.alt, c.ref = to_int(info.get("SUPPORT")), None
    if c.alt is not None and c.ref is not None and c.alt + c.ref > 0:
        c.vaf = c.alt / (c.alt + c.ref)
        c.dp = c.alt + c.ref
    else:
        c.vaf = num(d.get("AF") or d.get("VAF") or str(info.get("AF", "")).split(",")[0])
        c.dp = to_int(d.get("DP"))
    c.precise = 1.0 if "PRECISE" in info else (0.5 if "IMPRECISE" in info else None)
    # Sawfish CNV records: FORMAT GT:CN:CNQ, INFO SVCLAIM (D = depth, J = junction, DJ = both)
    cn = to_int(d.get("CN"))
    c.cn, c.cn_term = cn, None
    if cn is not None and c.svtype in ("DEL", "DUP", "CNV"):
        want = {"DEL": (0 if c.hom else 1), "DUP": (4 if c.hom else 3), "CNV": None}[c.svtype]
        right_dir = (cn < 2) if c.svtype == "DEL" else (cn > 2) if c.svtype == "DUP" else cn != 2
        c.cn_term = 1.0 if (want is not None and cn == want) else (0.7 if right_dir else 0.0)
    claim = info.get("SVCLAIM")
    c.claim = claim if isinstance(claim, str) else None
    # Sniffles2 COVERAGE: upstream, start, center, end, downstream
    cov = [num(x) for x in str(info.get("COVERAGE", "")).split(",")]
    c.cov, c.cov_flank, c.cov_ratio, c.covchg, c.depth_x = None, None, None, None, None
    if len(cov) == 5 and all(v is not None for v in cov):
        c.cov = ",".join(f"{v:g}" for v in cov)
        c.cov_flank = (cov[0] + cov[4]) / 2
        ref_flank = (min(cov[0], cov[4]) if c.svtype == "DEL" else
                     max(cov[0], cov[4]) if c.svtype == "DUP" else c.cov_flank)
        if ref_flank > 0:
            c.cov_ratio = cov[2] / ref_flank
            if c.svlen and c.svlen >= min_cov_len and c.svtype in ("DEL", "DUP"):
                r = c.cov_ratio
                if c.svtype == "DEL":
                    c.covchg = min(1.0, max(0.0, (1 - r) / (1.0 if c.hom else 0.5)))
                else:
                    c.covchg = min(1.0, max(0.0, (r - 1) / (1.0 if c.hom else 0.5)))
    c.dgenes = {}
    st = info.get("STRAND")  # Sniffles2: "+-" = supporting reads on both strands
    c.strand = (None if not isinstance(st, str) else 1.0 if ("+" in st and "-" in st)
                else 0.3 if st in ("+", "-") else None)
    ps = d.get("PS", ".")
    c.ps = None if ps in (".", "") else ps
    c.hap = None
    if "|" in c.gt and len(gts) == 2 and c.ps:
        c.hap = 1 if gts == ["0", "1"] else 0 if gts == ["1", "0"] else None
    if c.filt not in ("PASS", "."):
        c.status = c.filt
    elif all(g in ("0", ".") for g in gts):
        c.status = "no_alt_GT"
    else:
        c.status = "PASS"
    c.qual_pct, c.info, c.full, c.best, c.best_val = None, info, None, None, -1.0
    return c


def load_annotsv(path, label, a, stats):
    raw = open(path, "rb")
    fh = gzip.GzipFile(fileobj=raw) if path.endswith(".gz") else raw
    prog = Progress(label, os.path.getsize(path), a.progress_interval)
    ix, sidx, G, selfull, selsplit = None, None, None, [], []
    svs, n = {}, 0
    for bline in fh:
        n += 1
        if not n & 255:
            prog.tick(n, raw.tell())
        line = bline.decode("utf-8", "replace").rstrip("\n")
        if ix is None:
            cols, fixed = repair_header(line.split("\t"))
            if fixed:
                print(f"[{label}] NOTE: header tokens with a missing tab were split: "
                      f"{', '.join(fixed)}", file=sys.stderr)
            ix = {c: i for i, c in enumerate(cols)}
            if "Annotation_mode" not in ix or "FORMAT" not in ix or "AnnotSV_ID" not in ix:
                sys.exit(f"ERROR: {path} does not look like an AnnotSV TSV "
                         "(need AnnotSV_ID, FORMAT and Annotation_mode columns)")
            samples = cols[ix["FORMAT"] + 1: ix["Annotation_mode"]]
            if a.sample:
                if a.sample not in samples:
                    sys.exit(f"ERROR: sample '{a.sample}' not in {path}; found {samples[:5]}...")
                sidx = ix[a.sample]
            elif len(samples) == 1:
                sidx = ix[samples[0]]
            else:
                sys.exit(
                    f"ERROR: {path} has {len(samples):,} sample columns (first: {', '.join(samples[:3])}).\n"
                    "AnnotSV output for one sample's Sniffles/Sawfish VCF should have exactly one.\n"
                    "This looks like a cohort callset (e.g. 1000 Genomes) rather than your sample's calls -\n"
                    "check which VCF was handed to AnnotSV. If this is intended, pass --sample NAME.")
            missing = [c for c in FULL_COLS + SPLIT_COLS if c not in ix]
            if missing:
                print(f"[{label}] NOTE: columns not found (treated as empty): "
                      f"{', '.join(sorted(set(missing)))}", file=sys.stderr)
            G = lambda f, name: f[ix[name]] if name in ix and ix[name] < len(f) else ""
            selfull = [(c, ix[c]) for c in FULL_COLS if c in ix]
            selsplit = [(c, ix[c]) for c in SPLIT_COLS if c in ix]
            continue
        f = line.split("\t")
        if len(f) != len(ix):
            stats[(label, "lines whose column count differs from the header")] += 1
            if stats[(label, "lines whose column count differs from the header")] == 1:
                print(f"[{label}] WARNING: header has {len(ix)} columns but line {n} has {len(f)}; "
                      "column mapping may be shifted - check the TSV", file=sys.stderr)
        aid, mode = G(f, "AnnotSV_ID"), G(f, "Annotation_mode")
        if not aid:
            continue
        stats[(label, f"lines: {mode or '?'}")] += 1
        c = svs.get(aid)
        if c is None:
            c = svs[aid] = build_call(f, G, sidx, label, aid, a.min_cov_len)
        if mode == "full":
            c.full = {k: f[i] for k, i in selfull if i < len(f)}
        elif mode == "split":
            s = {k: f[i] for k, i in selsplit if i < len(f)}
            eff, desc = gene_effect(c.svtype, s)
            dis = gene_disease(s)
            rec = recessive_only(s) and not c.hom
            if rec:
                dis *= RECESSIVE_HET_FACTOR
            val = eff * (0.5 + 0.5 * dis)
            gene = s.get("Gene_name", "")
            if gene and eff >= 0.7 and dis >= 0.4:
                c.dgenes[gene] = max(c.dgenes.get(gene, (0, "")), (dis, "AR" if rec else ""))
            if val > c.best_val:
                c.best_val = val
                c.best = dict(s, effect=eff, disease=dis, desc=desc, rec_het=rec)
    prog.done(n, raw.tell())
    raw.close()

    calls = []
    for c in svs.values():
        if c.svlen is not None and c.svlen < a.min_svlen and c.svtype != "BND":
            stats[(label, f"skipped: shorter than --min-svlen {a.min_svlen}")] += 1
            continue
        calls.append(c)
        stats[(label, "kept PASS" if c.status == "PASS" else "kept failed (carried along)")] += 1
    # flank depth relative to the sample's typical depth (median over this caller's PASS calls)
    flanks = sorted(c.cov_flank for c in calls if c.status == "PASS" and c.cov_flank)
    ref_dp = a.expected_dp or (flanks[len(flanks) // 2] if flanks else None)
    if ref_dp:
        stats[(label, f"sample depth reference used for the depth-anomaly check: {ref_dp:g}x")] += 1
        for c in calls:
            if c.cov_flank:
                c.depth_x = c.cov_flank / ref_dp
    # QUAL percentile among PASS calls, only if QUAL is informative
    passed = [c for c in calls if c.status == "PASS" and c.qual is not None]
    quals = sorted(c.qual for c in passed)
    if len(set(quals)) >= 5:
        for c in passed:
            c.qual_pct = bisect.bisect_left(quals, c.qual) / (len(quals) - 1)
    else:
        print(f"[{label}] NOTE: QUAL has too few distinct values to be used", file=sys.stderr)
    return calls


# ---------------------------------------------------------------- matching ---------
def compatible(t1, t2):
    return t1 == t2 or {t1, t2} in ({"CNV", "DEL"}, {"CNV", "DUP"})


def compare(x, y, a):
    """Similarity (0..1) and metrics if x and y are plausibly the same event, else None."""
    if not compatible(x.svtype, y.svtype):
        return None
    ratio = (min(x.svlen, y.svlen) / max(x.svlen, y.svlen)
             if x.svlen and y.svlen else None)
    if x.svtype == "BND":
        d = abs(x.start - y.start)
        if d > a.max_dist or (x.mate and y.mate and (x.mate[0] != y.mate[0]
                                                      or abs(x.mate[1] - y.mate[1]) > a.max_dist)):
            return None
        return 1 - d / (a.max_dist + 1), None, ratio
    if x.svtype == "INS":
        d = abs(x.start - y.start)
        if d > a.max_dist or (ratio is not None and ratio < a.min_size_ratio):
            return None
        return 0.5 * (1 - d / (a.max_dist + 1)) + 0.5 * (ratio if ratio is not None else 0.5), None, ratio
    ov = min(x.end, y.end) - max(x.start, y.start)
    if ov <= 0:
        return None
    ro = min(ov / max(1, x.end - x.start), ov / max(1, y.end - y.start))
    return (ro, ro, ratio) if ro >= a.min_ro else None


def match_callers(A, B, a):
    """Greedy one-to-one matching; PASS-containing pairs win ties. Returns {i_in_A: (j_in_B, ro, ratio)}."""
    by_chrom = defaultdict(list)
    for j, c in enumerate(B):
        by_chrom[c.chrom].append(j)
    index = {}
    for chrom, js in by_chrom.items():
        js.sort(key=lambda j: B[j].start)
        index[chrom] = (js, [B[j].start for j in js], max(B[j].end - B[j].start for j in js))
    cand = []
    for i, x in enumerate(A):
        if x.chrom not in index:
            continue
        js, starts, maxspan = index[x.chrom]
        lo = bisect.bisect_left(starts, x.start - maxspan - a.max_dist)
        hi = bisect.bisect_right(starts, x.end + a.max_dist)
        for k in range(lo, hi):
            j = js[k]
            r = compare(x, B[j], a)
            if r:
                npass = (x.status == "PASS") + (B[j].status == "PASS")
                cand.append((npass, r[0], i, j, r[1], r[2]))
    cand.sort(reverse=True)
    pairs, usedA, usedB = {}, set(), set()
    for _, _, i, j, ro, ratio in cand:
        if i not in usedA and j not in usedB:
            pairs[i] = (j, ro, ratio)
            usedA.add(i)
            usedB.add(j)
    return pairs


# ---------------------------------------------------------------- scoring ----------
def caller_confidence(c, a):
    t = {}
    if c.qual_pct is not None:
        t["qual"] = c.qual_pct
    if c.alt is not None:
        t["support"] = min(c.alt / a.full_credit_support, 1.0)
    if c.vaf is not None:
        t["vaf"] = max(0.0, 1 - abs(c.vaf - (1.0 if c.hom else 0.5)) / 0.5)
    depth = c.cov_flank if c.cov_flank is not None else c.dp
    if depth is not None:
        t["depth"] = min(depth / a.full_credit_dp, 1.0)
    if c.covchg is not None:
        t["covchg"] = c.covchg
    if c.cn_term is not None:
        t["cn"] = c.cn_term
    if c.precise is not None:
        t["precise"] = c.precise
    if c.strand is not None:
        t["strand"] = c.strand
    if not t:
        return None
    w = sum(CONF_WEIGHTS[k] for k in t)
    return sum(CONF_WEIGHTS[k] * v for k, v in t.items()) / w


def rarity_score(freq):
    f = RARITY_FLOOR if not freq else max(freq, RARITY_FLOOR)
    num_ = math.log10(RARITY_COMMON) - math.log10(f)
    return min(1.0, max(0.0, num_ / (math.log10(RARITY_COMMON) - math.log10(RARITY_FLOOR))))


def db_frequency(c, db_keys, db_size):
    """Frequency/occurrence from SVDB/CoLoRSdb-style INFO keys of one call."""
    afs, occs, srcs = [], [], []
    for k in db_keys:
        v = c.info.get(k)
        if v is None or v is True:
            continue
        vals = [x for x in (num(t) for t in str(v).split(",")) if x is not None]
        if not vals:
            continue
        if re.search(r"(?i)occ$", k):
            occs.append((max(vals), k))
        else:
            afs.append((max(vals), k))
    if not afs and occs and db_size:
        afs = [(o / db_size, k) for o, k in occs]
    return (max(afs) if afs else (None, None)), (max(o for o, _ in occs) if occs else None)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sniffles", help="AnnotSV TSV of the Sniffles VCF")
    ap.add_argument("--sawfish", help="AnnotSV TSV of the Sawfish VCF")
    ap.add_argument("--sample", help="sample column name (only needed for multi-sample TSVs)")
    ap.add_argument("--output", required=True, help="shortlist TSV")
    ap.add_argument("--output-all", help="optional TSV with every scored SV")
    ap.add_argument("--top-n", type=int, default=100)
    ap.add_argument("--min-svlen", type=int, default=50, help="skip SVs shorter than this (default 50)")
    ap.add_argument("--min-ro", type=float, default=0.5, help="reciprocal overlap to match DEL/DUP/INV (0.5)")
    ap.add_argument("--max-dist", type=int, default=500, help="max breakpoint distance to match INS/BND (500)")
    ap.add_argument("--min-size-ratio", type=float, default=0.5, help="min size ratio to match INS (0.5)")
    ap.add_argument("--min-po-percent", type=float, default=50.0,
                    help="min %% of a known pathogenic SV that must be covered to count (50)")
    ap.add_argument("--min-cov-len", type=int, default=1000,
                    help="min SV length for the Sniffles COVERAGE change check (1000)")
    ap.add_argument("--single-exempt-len", type=int, default=100000,
                    help="SVs at least this long (or depth-only) are not penalised for a single caller (100000)")
    ap.add_argument("--expected-dp", type=float,
                    help="typical sample depth (default: median Sniffles flank depth)")
    ap.add_argument("--max-depth-x", type=float, default=2.5,
                    help="flank depth above this multiple of the sample depth counts as anomalous (2.5)")
    ap.add_argument("--large-sv-len", type=int, default=5000000,
                    help="SVs at least this long are flagged as arm-scale and down-weighted (5 Mb)")
    ap.add_argument("--max-depth-only", type=int, default=20,
                    help="max depth-only (SVCLAIM=D) CNV segments in the shortlist (20)")
    ap.add_argument("--cnv-warn-len", type=int, default=1000000,
                    help="warn if 3+ scored SVs are at least this long (1 Mb)")
    ap.add_argument("--full-credit-support", type=int, default=6, help="alt reads for full support credit (6)")
    ap.add_argument("--full-credit-dp", type=int, default=20, help="depth for full depth credit (20)")
    ap.add_argument("--db-key", action="append", default=[],
                    help="INFO key holding a population SV db frequency/occurrence (repeatable; "
                         "autodetected: keys ending in OCC/FRQ or containing 'colors')")
    ap.add_argument("--db-size", type=int, help="number of samples in the SV db (to turn *_OCC into a frequency)")
    ap.add_argument("--progress-interval", type=float, default=None)
    a = ap.parse_args()

    inputs = [(lab, p) for lab, p in (("sn", a.sniffles), ("sf", a.sawfish)) if p]
    if not inputs:
        ap.error("give at least one of --sniffles / --sawfish")

    stats, loaded = Counter(), {}
    for lab, path in inputs:
        print(f"[{lab}] reading {path}", file=sys.stderr)
        loaded[lab] = load_annotsv(path, lab, a, stats)
    order = list(loaded)

    # --- INFO keys that look like a population SV database
    keycount = Counter(k for cs in loaded.values() for c in cs[:2000] for k in c.info)
    db_keys = a.db_key or [k for k in keycount if re.search(r"(?i)(occ|frq)$", k) or "colors" in k.lower()]
    print(f"[info] INFO keys seen: {', '.join(sorted(keycount)) or '(none)'}", file=sys.stderr)
    print(f"[info] population-db keys used: {', '.join(db_keys) or '(none found; use --db-key)'}",
          file=sys.stderr)

    # --- match callers
    entries = []
    if len(order) == 2:
        A, B = loaded[order[0]], loaded[order[1]]
        pairs = match_callers(A, B, a)
        usedB = set()
        for i, x in enumerate(A):
            if i in pairs:
                j, ro, ratio = pairs[i]
                usedB.add(j)
                entries.append(({order[0]: x, order[1]: B[j]}, ro, ratio))
            else:
                entries.append(({order[0]: x}, None, None))
        entries += [({order[1]: y}, None, None) for j, y in enumerate(B) if j not in usedB]
        stats[("merge", "matched across callers")] = len(pairs)
    else:
        entries = [({order[0]: x}, None, None) for x in loaded[order[0]]]

    print(f"[merge] scoring {len(entries):,} SVs", file=sys.stderr)
    rows = []
    for present, ro, ratio in entries:
        passed = {lab: c for lab, c in present.items() if c.status == "PASS"}
        if not passed:
            stats[("merge", "dropped: only failed calls, no PASS")] += 1
            continue
        prim = next(iter(passed.values()))
        ann = next((c for c in [prim] + list(present.values()) if c.full is not None), prim)
        full = ann.full or {}
        best = ann.best
        t = prim.svtype

        # impact
        gene_term = best["effect"] * (0.5 + 0.5 * best["disease"]) if best else 0.0
        path, path_why = 0.0, ""
        for k in KIND_P.get(t, []):
            if real_sources(full.get(f"P_{k}_source")):
                path, path_why = KNOWN_PATH_STRONG, f"known pathogenic {k}: {short(full.get(f'P_{k}_phen'), 60)}"
            elif path < KNOWN_PATH_PARTIAL:
                cov = pathogenic_cover(full.get(f"po_P_{k}_source"), full.get(f"po_P_{k}_percent"))
                if cov is not None and cov >= a.min_po_percent:
                    path, path_why = KNOWN_PATH_PARTIAL, f"covers {cov:.0f}% of a known pathogenic {k}"
        sc = num(full.get("AnnotSV_ranking_score"))
        annot_term = max(0.0, min(1.0, sc)) if sc is not None else 0.0
        impact = max(gene_term, path, annot_term)
        if sc is not None and sc <= -0.9:
            impact *= ANNOTSV_BENIGN_FACTOR
        acmg = (full.get("ACMG_class") or "").replace("full=", "")

        # rarity
        cands = []
        for k in KIND_B.get(t, []):
            vals = af_values(full.get(f"B_{k}_AFmax"))
            if vals:
                cands.append((max(vals), f"AnnotSV B_{k}"))
        occ = None
        for c in present.values():
            (dbaf, dbk), o = db_frequency(c, db_keys, a.db_size)
            if dbaf is not None:
                cands.append((dbaf, dbk))
            if o is not None:
                occ = o if occ is None else max(occ, o)
        af, af_src = max(cands) if cands else (None, "")
        rar = rarity_score(af)

        # confidence
        confs = {lab: caller_confidence(c, a) for lab, c in passed.items()}
        vals = [v for v in confs.values() if v is not None]
        conf = sum(vals) / len(vals) if vals else 0.0
        ctx = 1.0
        if full.get("SegDup_left") or full.get("SegDup_right"):
            ctx *= SEGDUP_FACTOR
        if full.get("Gap_left") or full.get("Gap_right"):
            ctx *= GAP_FACTOR
        if full.get("ENCODE_blacklist_left") or full.get("ENCODE_blacklist_right"):
            ctx *= BLACKLIST_FACTOR
        # a lone call is only penalised if the other caller could plausibly have seen it:
        # depth-only (SVCLAIM=D) or Mb-scale CNVs are normally invisible to the other caller
        failed = {lab: c for lab, c in present.items() if c.status != "PASS"}
        exempt = (not failed) and (any(c.claim == "D" for c in passed.values())
                                   or (prim.svlen or 0) >= a.single_exempt_len)
        single = len(order) > 1 and len(passed) == 1 and not exempt
        conf *= ctx * (SINGLE_CALLER_FACTOR if single else 1.0)
        anomalies = []
        if any("maxscoringdepth" in c.status.lower() for c in failed.values()):
            conf *= FAILED_DEPTH_FACTOR
            anomalies.append("other caller filtered it for excessive depth (MaxScoringDepth)")
        dx = max((c.depth_x for c in passed.values() if c.depth_x), default=None)
        if dx is not None and dx > a.max_depth_x:
            conf *= DEPTH_ANOMALY_FACTOR
            anomalies.append(f"flank depth {dx:.1f}x the sample median: collapsed repeat/segdup artefact?")
        if (prim.svlen or 0) >= a.large_sv_len:
            conf *= LARGE_SV_FACTOR
            anomalies.append(f"arm-scale event ({prim.svlen / 1e6:.1f} Mb): check aneuploidy/mosaicism/segmentation")
        score = impact * rar * conf

        gts = {frozenset(re.split(r"[/|]", c.gt)) for c in passed.values()}
        gt_match = "" if len(passed) < 2 else ("yes" if len(gts) == 1 else "NO")
        depth_only = any(c.claim == "D" for c in passed.values())

        # reasons
        why = [f"{t} {prim.svlen:,} bp" if prim.svlen else t]
        if best and best["effect"] > 0:
            why.append(f"{best['desc']} in {best.get('Gene_name', '')}"
                       + (f" ({best.get('Location', '')})" if best.get("Location") else ""))
        if best and best["disease"] > 0:
            why.append(f"gene relevance {best['disease']:.1f}"
                       + (" (recessive-only gene, het call: needs a second hit)" if best.get("rec_het") else ""))
        dg = prim.dgenes
        if len(dg) > 1:
            why.append(f"{len(dg)} disease-associated genes hit")
        if path_why:
            why.append(path_why)
        if sc is not None:
            why.append(f"AnnotSV score {sc:.2f} (class {acmg or '?'})")
        why.append("no overlap with benign SV/db" if af is None else f"max AF {af:.2g} ({af_src})")
        why.append("called by " + "+".join(l.upper() for l in passed)
                   + (" (only one)" if single else " (only one; not penalised, depth-only/large)"
                      if len(order) > 1 and len(passed) == 1 else ""))
        for lab, c in present.items():
            if c.status != "PASS":
                why.append(f"{lab.upper()} {c.status} (GT {c.gt}, alt reads {c.alt})")
            elif c.alt is not None and c.alt < a.full_credit_support / 2:
                why.append(f"low support {lab.upper()} ({c.alt} alt reads)")
            if c.covchg is not None:
                why.append(f"coverage change {'supports' if c.covchg >= 0.6 else 'weak/absent'} "
                           f"{lab.upper()} (center/flank {c.cov_ratio:.2f})")
            if c.claim == "D":
                why.append(f"{lab.upper()}: depth-only evidence (SVCLAIM=D), breakpoints approximate")
            if c.cn is not None and c.cn_term is not None and c.cn_term < 1:
                why.append(f"{lab.upper()}: CN={c.cn} {'fits' if c.cn_term else 'contradicts'} {c.svtype} genotype {c.gt}")
        if gt_match == "NO":
            why.append("genotype disagrees between callers ("
                       + " vs ".join(f"{l.upper()} {c.gt}" for l, c in passed.items()) + ")")
        if ctx < 1:
            why.append("segdup/gap/blacklist context")
        why.extend(anomalies)

        genes = [g for g in (full.get("Gene_name") or "").split(";") if g]
        row = {
            "locus": f"{prim.chrom}:{prim.start}" + (f"-{prim.end}" if t not in ("INS", "BND") else ""),
            "chrom": prim.chrom, "start": prim.start, "end": prim.end, "svtype": t,
            "svlen": fmt(prim.svlen), "cytoband": full.get("CytoBand", ""),
            "genes": ";".join(genes[:6]) + (f";+{len(genes) - 6} more" if len(genes) > 6 else ""),
            "gene_count": full.get("Gene_count", ""),
            "best_gene": (best or {}).get("Gene_name", ""), "tx": (best or {}).get("Tx", ""),
            "location": ((best or {}).get("Location", "") + " " + (best or {}).get("Location2", "")).strip(),
            "cds_pct": (best or {}).get("Overlapped_CDS_percent", ""),
            "frameshift": (best or {}).get("Frameshift", ""),
            "gene_effect": fmt((best or {}).get("effect"), 2), "gene_relevance": fmt((best or {}).get("disease"), 2),
            "omim_phenotype": short((best or {}).get("OMIM_phenotype"), 80),
            "gencc": (best or {}).get("GenCC_classification", ""),
            "hi": (best or {}).get("HI", ""), "pli": (best or {}).get("GnomAD_pLI", ""),
            "loeuf_bin": (best or {}).get("LOEUF_bin", ""),
            "n_disease_genes": len(dg),
            "disease_genes": ";".join(f"{g}{'(AR)' if m[1] else ''}" for g, m in
                                      sorted(dg.items(), key=lambda kv: -kv[1][0])[:8]),
            "known_pathogenic": path_why, "annotsv_score": fmt(sc, 3), "annotsv_class": acmg,
            "annotsv_criteria": short(full.get("AnnotSV_ranking_criteria"), 200),
            "max_af": fmt(af), "af_source": af_src, "db_occ": fmt(occ),
            "segdup": "yes" if (full.get("SegDup_left") or full.get("SegDup_right")) else "",
            "called_by": "+".join(l.upper() for l in passed), "gt_match": gt_match,
            "evidence_class": "depth-only" if depth_only else "junction",
            "_annot": sc if sc is not None else 0.0, "_ndg": len(dg), "_svlen": prim.svlen or 0,
            "match_ro": fmt(ro, 2), "match_size_ratio": fmt(ratio, 2),
            "start_delta": (present[order[1]].start - present[order[0]].start) if len(present) == 2 else "",
            "end_delta": (present[order[1]].end - present[order[0]].end) if len(present) == 2 else "",
            "impact_score": fmt(impact, 2), "rarity_score": fmt(rar, 2), "conf_score": fmt(conf, 2),
            "reasons": "; ".join(why), "_score": score, "score": fmt(score, 3),
        }
        for lab in order:
            c = present.get(lab)
            row[f"{lab}_status"] = c.status if c else "not_seen"
            row[f"{lab}_id"] = c.id if c else ""
            row[f"{lab}_start"] = c.start if c else ""
            row[f"{lab}_end"] = c.end if c else ""
            row[f"{lab}_svlen"] = fmt(c.svlen) if c else ""
            row[f"{lab}_qual"] = fmt(c.qual) if c else ""
            row[f"{lab}_gt"] = c.gt if c else ""
            row[f"{lab}_alt_reads"] = fmt(c.alt) if c else ""
            row[f"{lab}_ref_reads"] = fmt(c.ref) if c else ""
            row[f"{lab}_vaf"] = fmt(c.vaf, 2) if c else ""
            row[f"{lab}_cn"] = fmt(c.cn) if c else ""
            row[f"{lab}_claim"] = (c.claim or "") if c else ""
            row[f"{lab}_cov"] = (c.cov or "") if c else ""
            row[f"{lab}_cov_ratio"] = fmt(c.cov_ratio, 2) if c else ""
            row[f"{lab}_conf"] = fmt(confs.get(lab), 2)
        row["phase_set"] = next((c.ps for c in passed.values() if c.ps), "")
        rows.append(row)

    rows.sort(key=lambda r: (-r["_score"], -r["_annot"], -r["_ndg"], "+" not in r["called_by"],
                             r["chrom"], r["start"]))
    for i, r in enumerate(rows, 1):
        r["rank"] = i

    cols = ["rank", "locus", "chrom", "start", "end", "svtype", "svlen", "cytoband", "genes", "gene_count",
            "best_gene", "tx", "location", "cds_pct", "frameshift", "gene_effect", "gene_relevance",
            "omim_phenotype", "gencc", "hi", "pli", "loeuf_bin", "n_disease_genes", "disease_genes",
            "known_pathogenic", "annotsv_score",
            "annotsv_class", "annotsv_criteria", "max_af", "af_source", "db_occ", "segdup", "called_by", "gt_match",
            "evidence_class", "match_ro", "match_size_ratio", "start_delta", "end_delta"]
    for lab in order:
        cols += [f"{lab}_{x}" for x in ("status", "id", "start", "end", "svlen", "qual", "gt",
                                         "alt_reads", "ref_reads", "vaf", "cn", "claim", "cov",
                                         "cov_ratio", "conf")]
    cols += ["phase_set", "impact_score", "rarity_score", "conf_score", "score", "reasons"]

    def write(path, subset):
        with open(path, "w", newline="") as out:
            w = csv.DictWriter(out, fieldnames=cols, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(subset)

    # depth-only CNV segments have no junction evidence and tend to flood the list: cap them
    shortlist, n_depth, n_depth_cut = [], 0, 0
    for r in rows:
        if r["_score"] <= 0 or len(shortlist) >= a.top_n:
            break
        if r["evidence_class"] == "depth-only":
            if n_depth >= a.max_depth_only:
                n_depth_cut += 1
                continue
            n_depth += 1
        shortlist.append(r)
    write(a.output, shortlist)
    if a.output_all:
        write(a.output_all, rows)

    print("\n=== summary ===", file=sys.stderr)
    for (lab, what), n in sorted(stats.items()):
        print(f"[{lab}] {what}: {n:,}", file=sys.stderr)
    print(f"SVs scored: {len(rows):,}", file=sys.stderr)
    if n_depth_cut:
        print(f"depth-only CNV segments held back by --max-depth-only {a.max_depth_only}: {n_depth_cut} "
              "(they are all in --output-all)", file=sys.stderr)
    big = sum(1 for r in rows if r["_score"] > 0 and r["_svlen"] >= a.cnv_warn_len)
    if big >= 3:
        print(f"WARNING: {big} scored SVs of >= {a.cnv_warn_len / 1e6:g} Mb. That is unusual for constitutional "
              "DNA; the sample may be a tumour/cell line or the CNV segmentation is noisy.", file=sys.stderr)
    print("shortlist by type:", dict(Counter(r["svtype"] for r in shortlist)), file=sys.stderr)
    print("shortlist by caller:", dict(Counter(r["called_by"] for r in shortlist)), file=sys.stderr)
    print(f"wrote {len(shortlist)} rows -> {a.output}; peak RAM {peak_ram_gb():.1f} GB", file=sys.stderr)


if __name__ == "__main__":
    main()