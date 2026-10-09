#!/usr/bin/env python3
"""
variant_overlap_report.py - per-sample caller comparison report (PDF + TSV tables).

Reads <variants-dir>/<tool>_<sample>/ and writes <outdir>/<sample>_variant_overlap.pdf with
  * overview tables (bcftools stats, SV summary)
  * SNV/indel comparison   deepvariant vs nanocaller  (Venn, types, spectrum, GT concordance, ...)
  * SV comparison          sawfish vs sniffles        (Venn, types, sizes, size agreement, ...)
  * one genome-wide track page per compared tool (calls along the genome, private vs shared,
    windowed counts, chromosome strip; chrM is drawn in the mitorsaw colour, scaled by its call count)
  * density heatmaps (chromosome x window) per tool and the log2 ratio between the two tools
  * TRGT page (repeat loci along the genome, STR_STATUS from stranger if available)
  * mitorsaw page (chrM variants, coverage, gene map)
  * paraphase page (copy numbers per region, regions with CN != 2)
Phased VCFs (whatshap, longphase) are not used: they are derived from the caller outputs.

SNVs/indels match on chrom:pos:ref:alt after splitting multiallelics (PASS only).
SVs match when type and chromosome agree and either both breakpoints are within --max-dist,
or the reciprocal overlap is >= --reciprocal-overlap (INS: position within --max-dist and size
ratio >= --reciprocal-overlap; BND: both breakpoints within --max-dist). Only chr1-22,X,Y,M are used.
Missing files/directories only produce warnings.

Conda dependencies:
  conda install -c conda-forge -c bioconda bcftools python pandas numpy matplotlib
"""
import argparse
import csv
import datetime as dt
import json
import math
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import textwrap
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

WARNINGS = []


def warn(msg):
    WARNINGS.append(msg)
    print(f"WARNING: {msg}", file=sys.stderr, flush=True)


def info(msg):
    print(msg, file=sys.stderr, flush=True)


try:
    import numpy as np
    import pandas as pd
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.lines import Line2D
    from matplotlib.patches import Circle, Rectangle
    from matplotlib.ticker import FuncFormatter
except ImportError as e:
    sys.exit(
        f"Missing python dependency ({e}).\n"
        "Install with: conda install -c conda-forge -c bioconda bcftools python pandas numpy matplotlib"
    )


def q(x):
    return shlex.quote(str(x))


# --------------------------------------------------------------------------------------
# Colours (edit freely). The mitorsaw colour is used for chrM in the chromosome strips.
# --------------------------------------------------------------------------------------
COLORS = {
    "deepvariant": "#E69F00", "nanocaller": "#0072B2",
    "sawfish": "#009E73", "sniffles": "#D55E00", "hificnv": "#C9B800",
    "shared": "#9B7FC7",
    "mitorsaw": "#2CA02C", "trgt": "#4C78C9", "paraphase": "#CC79A7", "other": "#666666",
}


def tool_color(name):
    return COLORS.get(name, COLORS["other"])


# --------------------------------------------------------------------------------------
# Where to look: (name, kind, directory, patterns); {s} = sample; first hit wins.
# A pattern ending in .vcf.gz falls back to the plain .vcf. Phased VCFs are NOT listed on purpose.
# --------------------------------------------------------------------------------------
REGISTRY = [
    ("deepvariant", "snv", "deepvariant_{s}", ["{s}_variants.vcf.gz"]),
    ("nanocaller", "snv", "nanocaller_{s}", ["{s}_nanocaller.vcf.gz"]),
    ("sawfish", "sv", "sawfish_phased_{s}", ["{s}_genotyped.sv.vcf.gz"]),
    ("sniffles", "sv", "sniffles_{s}", ["{s}_svs.vcf.gz"]),
    ("hificnv", "sv", "hificnv_{s}", ["*.vcf.gz"]),
    ("trgt", "special", "trgt_{s}", ["{s}.vcf.gz"]),
    ("mitorsaw", "special", "mitorsaw_{s}", ["*variants.vcf.gz"]),
    ("paraphase", "special", "paraphase_{s}", ["*.paraphase.json"]),
]
# stranger-annotated TRGT VCF, searched under --annotated-dir
STRANGER = ("stranger_{s}", ["{s}_trgt_stranger_annotated.vcf.gz"])

HEAT_BIN_SNV = 1_000_000
HEAT_BIN_SV = 5_000_000
SV_TYPES = ["DEL", "INS", "DUP", "INV", "BND"]
SV_MARKERS = {"DEL": "v", "INS": "^", "DUP": "s", "INV": "D", "BND": "x"}

PRIMARY = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
GRCH38 = dict(zip(PRIMARY, [
    248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636,
    138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345,
    83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415, 16569]))

# rCRS coordinates (GRCh38 chrM), 1-based; tRNAs omitted
MT_GENES = [("D-loop", 16024, 16569), ("D-loop", 1, 576), ("RNR1", 648, 1601), ("RNR2", 1671, 3229),
            ("ND1", 3307, 4262), ("ND2", 4470, 5511), ("CO1", 5904, 7445), ("CO2", 7586, 8269),
            ("ATP8", 8366, 8572), ("ATP6", 8527, 9207), ("CO3", 9207, 9990), ("ND3", 10059, 10404),
            ("ND4L", 10470, 10766), ("ND4", 10760, 12137), ("ND5", 12337, 14148),
            ("ND6", 14149, 14673), ("CYB", 14747, 15887)]


def canon(name):
    m = re.fullmatch(r"(?:chr)?([0-9]{1,2}|X|Y|M|MT)", str(name), re.I)
    if not m:
        return None
    g = m.group(1).upper()
    g = "M" if g == "MT" else g
    if g.isdigit() and not 1 <= int(g) <= 22:
        return None
    return "chr" + g


class Genome:
    """Concatenated genome axis. chrM is drawn with a minimum width so it stays visible."""

    def __init__(self, lengths_raw):
        canon_len = {}
        for raw, ln in lengths_raw.items():
            c = canon(raw)
            if c and c not in canon_len:
                canon_len[c] = int(ln)
        self.names = [c for c in PRIMARY if c in canon_len]
        self.len = np.array([canon_len[c] for c in self.names], dtype=np.int64)
        non_m = sum(canon_len[c] for c in self.names if c != "chrM")
        disp = [max(int(non_m * 0.012), canon_len[c]) if c == "chrM" else canon_len[c] for c in self.names]
        self.disp = np.array(disp, dtype=np.float64)
        self.off_arr = np.concatenate([[0.0], np.cumsum(self.disp)[:-1]])
        self.scale_arr = self.disp / self.len
        self.total = float(self.disp.sum())
        self.code = {c: i for i, c in enumerate(self.names)}
        self.raw2code = {}
        for c, i in self.code.items():
            for alias in (c, c[3:], "MT" if c == "chrM" else None, "chrMT" if c == "chrM" else None):
                if alias:
                    self.raw2code[alias] = i
        for raw in lengths_raw:
            c = canon(raw)
            if c in self.code:
                self.raw2code[raw] = self.code[c]

    def x(self, codes, pos):
        codes = np.asarray(codes, dtype=np.int64)
        return self.off_arr[codes] + np.asarray(pos, dtype=np.float64) * self.scale_arr[codes]

    def locate(self, xv):
        i = int(np.clip(np.searchsorted(self.off_arr, xv, side="right") - 1, 0, len(self.names) - 1))
        return self.names[i], (xv - self.off_arr[i]) / self.scale_arr[i]


def probe_header(path):
    r = subprocess.run(["bcftools", "view", "-h", str(path)], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stderr.strip()[-300:] or "cannot read header")
    fmt = set(re.findall(r"^##FORMAT=<ID=([^,>]+)", r.stdout, re.M))
    inf = set(re.findall(r"^##INFO=<ID=([^,>]+)", r.stdout, re.M))
    contigs = {m.group(1): int(m.group(2)) for m in re.finditer(r"^##contig=<ID=([^,>]+),length=(\d+)", r.stdout, re.M)}
    return fmt, inf, contigs


def load_lengths(args, vcfs):
    fai = args.fai or (f"{args.ref}.fai" if args.ref and Path(f"{args.ref}.fai").is_file() else None)
    if fai and Path(fai).is_file():
        lengths = {}
        for line in open(fai):
            f = line.split("\t")
            lengths[f[0]] = int(f[1])
        return lengths, f"fai {fai}"
    for v in vcfs:
        try:
            contigs = probe_header(v)[2]
        except Exception:
            continue
        if any(canon(c) for c in contigs):
            return contigs, f"VCF header of {Path(v).name}"
    warn("no contig lengths found (use --fai or --ref with .fai, or ##contig lines in the VCFs); "
         "using built-in GRCh38 lengths")
    return dict(GRCH38), "built-in GRCh38"


# --------------------------------------------------------------------------------------
# Discovery
# --------------------------------------------------------------------------------------
def find_file(sample, bases, dirpat, patterns):
    dname = dirpat.replace("{s}", sample)
    for base in bases:
        d = Path(base) / dname
        if not d.is_dir():
            continue
        for pat in patterns:
            for p in [pat] + ([pat[:-3]] if pat.endswith(".vcf.gz") else []):
                hits = sorted(d.glob(p.replace("{s}", sample)))
                if hits:
                    return hits[0]
    return None


# --------------------------------------------------------------------------------------
# bcftools streaming
# --------------------------------------------------------------------------------------
def stream_bcftools(cmd, columns, dtypes, chunksize=1_000_000):
    with tempfile.TemporaryFile() as err:
        proc = subprocess.Popen(["bash", "-o", "pipefail", "-c", cmd], stdout=subprocess.PIPE, stderr=err)
        try:
            try:
                reader = pd.read_csv(proc.stdout, sep="\t", header=None, names=columns, dtype=dtypes,
                                     chunksize=chunksize, keep_default_na=False, quoting=csv.QUOTE_NONE)
                for chunk in reader:
                    yield chunk
            except pd.errors.EmptyDataError:
                pass
        finally:
            proc.stdout.close()
            rc = proc.wait()
        if rc != 0:
            err.seek(0)
            raise RuntimeError(err.read().decode(errors="replace").strip()[-400:] or f"bcftools exit code {rc}")


# --------------------------------------------------------------------------------------
# Small variants (SNV / indel)
# --------------------------------------------------------------------------------------
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
SUBS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
SUBMAP = {}
for _r in "ACGT":
    for _a in "ACGT":
        if _r != _a:
            _rr, _aa = (_r, _a) if _r in "CT" else (COMP[_r], COMP[_a])
            SUBMAP[f"{_r}>{_a}"] = SUBS.index(f"{_rr}>{_aa}")
VTYPES = ["SNP", "insertion", "deletion", "MNP/other"]
GTCLS = ["het", "hom-alt", "other"]


def _vaf(col, mode):
    if mode in ("VAF", "AF"):
        return pd.to_numeric(col.str.split(",").str[0], errors="coerce").to_numpy("float32")
    if mode == "AD":
        parts = col.str.split(",", n=2, expand=True)
        if parts.shape[1] < 2:
            return np.full(len(col), np.nan, dtype="float32")
        r, a = pd.to_numeric(parts[0], errors="coerce"), pd.to_numeric(parts[1], errors="coerce")
        return (a / (r + a)).to_numpy("float32")
    return np.full(len(col), np.nan, dtype="float32")


def snv_table(path, args, genome):
    """One row per unique small variant: key, chromosome code, pos, type, substitution class, ..."""
    fmt_ids = probe_header(path)[0]
    mode = next((t for t in ("VAF", "AF") if t in fmt_ids), "AD" if "AD" in fmt_ids else None)
    vaf_expr = {"VAF": "[%VAF]", "AF": "[%AF]", "AD": "[%AD]"}.get(mode, ".")
    gt_expr = "[%GT]" if "GT" in fmt_ids else "."
    cmd = "bcftools norm -m -any "
    if args.ref:
        cmd += f"-f {q(args.ref)} -c w "
    cmd += f"-Ou {q(path)}"
    if not args.no_pass_filter:
        cmd += " | bcftools view -f .,PASS -Ou"
    cmd += f" | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t{gt_expr}\\t{vaf_expr}\\n'"
    cols = ["chrom", "pos", "ref", "alt", "qual", "gt", "vaf"]
    dtypes = {c: str for c in cols}
    dtypes["pos"] = "int64"
    parts = defaultdict(list)
    for ch in stream_bcftools(cmd, cols, dtypes):
        codes = ch["chrom"].map(genome.raw2code)
        keep = codes.notna().to_numpy() & ~ch["alt"].str.contains(r"^<|[\[\]]|^\*$", regex=True).to_numpy()
        ch, codes = ch[keep], codes[keep]
        if ch.empty:
            continue
        rl, al = ch["ref"].str.len().to_numpy(), ch["alt"].str.len().to_numpy()
        vtype = np.where((rl == 1) & (al == 1), 0, np.where(al > rl, 1, np.where(al < rl, 2, 3))).astype("int8")
        sub = (ch["ref"] + ">" + ch["alt"]).map(SUBMAP).fillna(-1).to_numpy()
        g = ch["gt"].str.replace("|", "/", regex=False)
        gtc = np.full(len(ch), 2, dtype="int8")
        gtc[g.isin(["0/1", "1/0"]).to_numpy()] = 0
        gtc[(g == "1/1").to_numpy()] = 1
        parts["key"].append(pd.util.hash_pandas_object(ch[["chrom", "pos", "ref", "alt"]], index=False).to_numpy("uint64"))
        parts["ch"].append(codes.to_numpy().astype("int8"))
        parts["pos"].append(ch["pos"].to_numpy().astype("int32"))
        parts["vtype"].append(vtype)
        parts["sub"].append(np.where(vtype == 0, sub, -1).astype("int8"))
        parts["ilen"].append(np.clip(al - rl, -50, 50).astype("int16"))
        parts["qual"].append(pd.to_numeric(ch["qual"], errors="coerce").to_numpy("float32"))
        parts["vaf"].append(_vaf(ch["vaf"], mode))
        parts["gt"].append(gtc)
    if not parts:
        return pd.DataFrame({k: [] for k in ("key", "ch", "pos", "vtype", "sub", "ilen", "qual", "vaf", "gt")})
    df = pd.DataFrame({k: np.concatenate(v) for k, v in parts.items()})
    return df.drop_duplicates("key").sort_values("key", ignore_index=True)


def compare_snv(A, B):
    _, ia, ib = np.intersect1d(A["key"].to_numpy(), B["key"].to_numpy(), assume_unique=True, return_indices=True)
    sa, sb = np.zeros(len(A), bool), np.zeros(len(B), bool)
    sa[ia], sb[ib] = True, True
    return dict(sa=sa, sb=sb, ia=ia, ib=ib)


def tstv(sub):
    ts, tv = np.isin(sub, (2, 4)).sum(), np.isin(sub, (0, 1, 3, 5)).sum()
    return ts / tv if tv else float("nan")


def bcftools_stats(path, args):
    flt = "" if args.no_pass_filter else "-f .,PASS "
    out = subprocess.run(f"bcftools stats -s - {flt}{q(path)}", shell=True, capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(out.stderr.strip()[-300:])
    res, seen = {}, False
    for line in out.stdout.splitlines():
        f = line.split("\t")
        if f[0] == "SN" and len(f) >= 4:
            res[f[2].strip().rstrip(":")] = int(f[3])
        elif f[0] == "TSTV" and len(f) >= 5:
            res["ts/tv"] = f[4]
        elif f[0] == "PSC" and not seen and len(f) >= 6:
            seen = True
            res["het"], res["homalt"] = int(f[5]), int(f[4])
    return res


# --------------------------------------------------------------------------------------
# Structural variants
# --------------------------------------------------------------------------------------
RE_SVTYPE = re.compile(r"(?:^|;)SVTYPE=([^;]+)")
RE_END = re.compile(r"(?:^|;)END=(\d+)")
RE_SVLEN = re.compile(r"(?:^|;)SVLEN=(-?\d+)")
RE_BND = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")


def parse_sv(chrom, pos, ref, alt, info_str):
    alt = alt.split(",")[0]
    m = RE_SVTYPE.search(info_str)
    svtype = m.group(1).split(":")[0].upper() if m else None
    bnd = RE_BND.search(alt)
    if svtype is None:
        if bnd:
            svtype = "BND"
        elif alt.startswith("<"):
            svtype = alt.strip("<>").split(":")[0].upper()
        else:
            svtype = "INS" if len(alt) > len(ref) else "DEL"
    if svtype == "BND":
        chr2, end = (bnd.group(1), int(bnd.group(2))) if bnd else (chrom, pos)
        return dict(chrom=chrom, chr2=chr2, type="BND", start=pos, end=end, size=0)
    size, end = 0, None
    ml, me = RE_SVLEN.search(info_str), RE_END.search(info_str)
    if ml:
        size = abs(int(ml.group(1)))
    if me:
        end = int(me.group(1))
    if svtype == "INS":
        if not size and not alt.startswith("<"):
            size = max(len(alt) - len(ref), 0)
        end = pos
    else:
        if not size and end is not None:
            size = max(end - pos, 0)
        if end is None:
            end = pos + size
        end = max(end, pos)
    return dict(chrom=chrom, chr2=chrom, type=svtype, start=pos, end=end, size=size)


def read_svs(path, args, genome):
    cmd = f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t%INFO\\n' {q(path)}"
    cols = ["chrom", "pos", "ref", "alt", "filt", "info"]
    rows, n_nonpass, n_small = [], 0, 0
    for ch in stream_bcftools(cmd, cols, {c: str for c in cols}):
        for chrom, pos, ref, alt, filt, inf in ch.itertuples(index=False):
            if chrom not in genome.raw2code:
                continue
            if not args.no_pass_filter and filt not in ("PASS", "."):
                n_nonpass += 1
                continue
            rec = parse_sv(chrom, int(pos), ref, alt, inf)
            if rec["type"] == "BND":
                if (rec["chrom"], rec["start"]) > (rec["chr2"], rec["end"]):
                    continue  # breakend pairs: keep one of the two records
            elif rec["size"] and rec["size"] < args.min_sv_size:
                n_small += 1
                continue
            rows.append(rec)
    df = pd.DataFrame(rows, columns=["chrom", "chr2", "type", "start", "end", "size"])
    df["ch"] = df["chrom"].map(genome.raw2code).astype("int8") if len(df) else pd.Series(dtype="int8")
    return df, n_nonpass, n_small


def sv_match(typ, sa, ea, za, sb, eb, zb, d, r):
    if typ == "BND":
        return abs(sa - sb) <= d and abs(ea - eb) <= d
    if typ == "INS":
        if abs(sa - sb) > d:
            return False
        return min(za, zb) / max(za, zb) >= r if za and zb else True
    if abs(sa - sb) <= d and abs(ea - eb) <= d:
        return True
    ov = min(ea, eb) - max(sa, sb)
    return ov > 0 and ov / max(ea - sa, eb - sb, 1) >= r


def sv_cluster(sv_by_tool, args):
    """All SV records of both tools with a cluster root and the membership mask of that cluster."""
    frames = []
    for i, n in enumerate(sv_by_tool):
        df = sv_by_tool[n].copy()
        df["tool"] = i
        frames.append(df)
    allsv = pd.concat(frames, ignore_index=True)
    parent = list(range(len(allsv)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    d, r = args.max_dist, args.reciprocal_overlap
    for (_, _, typ), g in allsv.groupby(["chrom", "chr2", "type"], sort=False):
        g = g.sort_values("start")
        idx, st, en = g.index.to_numpy(), g["start"].to_numpy(), g["end"].to_numpy()
        sz, tl = g["size"].to_numpy(), g["tool"].to_numpy()
        n = len(idx)
        for a in range(n):
            lim = (st[a] if typ == "BND" else en[a]) + d
            for b in range(a + 1, n):
                if st[b] > lim:
                    break
                if tl[a] != tl[b] and sv_match(typ, st[a], en[a], sz[a], st[b], en[b], sz[b], d, r):
                    ra, rb = find(idx[a]), find(idx[b])
                    if ra != rb:
                        parent[rb] = ra
    roots = np.fromiter((find(i) for i in range(len(allsv))), dtype=np.int64, count=len(allsv))
    bits = np.left_shift(1, allsv["tool"].to_numpy().astype(np.int64))
    mask = np.zeros(len(allsv), dtype=np.int64)
    np.bitwise_or.at(mask, roots, bits)
    allsv["root"], allsv["mask"] = roots, mask[roots]
    return allsv


def sv_clusters(allsv, n_tools=2):
    cl = allsv.groupby("root").agg(mask=("mask", "first"), type=("type", "first"))
    for i in range(n_tools):
        sub = allsv[allsv["tool"] == i].groupby("root").agg(size=("size", "median"), start=("start", "median"))
        cl[f"size_{i}"], cl[f"start_{i}"] = sub["size"], sub["start"]
    return cl


# --------------------------------------------------------------------------------------
# Special tools: TRGT, mitorsaw, paraphase
# --------------------------------------------------------------------------------------
def read_trgt(path, genome, args):
    fmt_ids = probe_header(path)[0]
    gt = "[%GT]" if "GT" in fmt_ids else "."
    al = "[%AL]" if "AL" in fmt_ids else "."
    mc = "[%MC]" if "MC" in fmt_ids else "."
    cmd = f"bcftools query -f '%CHROM\\t%POS\\t%INFO\\t{gt}\\t{al}\\t{mc}\\n' {q(path)}"
    cols = ["chrom", "pos", "info", "gt", "al", "mc"]
    dtypes = {c: str for c in cols}
    dtypes["pos"] = "int64"
    out = []
    for ch in stream_bcftools(cmd, cols, dtypes):
        codes = ch["chrom"].map(genome.raw2code)
        ch, codes = ch[codes.notna().to_numpy()], codes[codes.notna()]
        if ch.empty:
            continue
        inf = ch["info"]
        status = inf.str.extract(r"(?:^|;)STR_STATUS=([^;]+)")[0].fillna("")
        st = np.full(len(ch), -1, dtype="int8")
        st[status.str.contains("normal").to_numpy()] = 0
        st[status.str.contains("pre_mutation").to_numpy()] = 1
        st[status.str.contains("full_mutation").to_numpy()] = 2
        al_num = ch["al"].str.split(",", expand=True).apply(pd.to_numeric, errors="coerce")
        mc_num = ch["mc"].str.replace("_", ",", regex=False).str.split(",", expand=True).apply(pd.to_numeric, errors="coerce")
        g = ch["gt"].str.replace("|", "/", regex=False)
        gp = g.str.split("/", expand=True)
        het = (gp[0] != gp[1]).to_numpy() if gp.shape[1] > 1 else np.zeros(len(ch), bool)
        out.append(pd.DataFrame({
            "ch": codes.to_numpy().astype("int8"), "pos": ch["pos"].to_numpy(),
            "end": pd.to_numeric(inf.str.extract(r"(?:^|;)END=(\d+)")[0], errors="coerce").to_numpy(),
            "trid": inf.str.extract(r"(?:^|;)TRID=([^;]+)")[0].fillna("").to_numpy(),
            "disease": inf.str.extract(r"(?:^|;)Disease=([^;]+)")[0].fillna("").to_numpy(),
            "status": st, "al_max": al_num.max(axis=1).to_numpy(), "mc_max": mc_num.max(axis=1).to_numpy(),
            "called": ~g.isin([".", "./.", ".|.", ""]).to_numpy(), "het": het,
        }))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def read_mito(path, args):
    fmt_ids = probe_header(path)[0]
    af_tag = next((t for t in ("HF", "AF", "VAF", "HP") if t in fmt_ids), None)
    dp_tag = next((t for t in ("DP", "COV") if t in fmt_ids), None)
    gt = "[%GT]" if "GT" in fmt_ids else "."
    afx = f"[%{af_tag}]" if af_tag else ("[%AD]" if "AD" in fmt_ids else ".")
    dpx = f"[%{dp_tag}]" if dp_tag else "."
    cmd = f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t{gt}\\t{afx}\\t{dpx}\\n' {q(path)}"
    cols = ["chrom", "pos", "ref", "alt", "filt", "gt", "af", "dp"]
    parts = [c for c in stream_bcftools(cmd, cols, {k: str for k in cols})]
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame({c: [] for c in cols})
    if not args.no_pass_filter and len(df):
        df = df[df["filt"].isin(["PASS", "."])]
    mode = af_tag if af_tag else ("AD" if "AD" in fmt_ids else None)
    df = df.assign(
        pos=pd.to_numeric(df["pos"], errors="coerce"),
        af=_vaf(df["af"], "AF" if af_tag else mode) if len(df) else [],
        dp=pd.to_numeric(df["dp"].str.split(",").str[0], errors="coerce") if len(df) else [],
    )
    df.attrs["format_tags"] = ",".join(sorted(fmt_ids)) or "none"
    return df


def load_json(path):
    with open(path) as fh:
        return json.load(fh)


def longest_numeric_list(obj, min_len=1000):
    best = None
    stack = [obj]
    while stack:
        o = stack.pop()
        if isinstance(o, list):
            if len(o) >= min_len and all(isinstance(v, (int, float)) for v in o[:50]):
                if best is None or len(o) > len(best):
                    best = o
            else:
                stack.extend(o[:200])
        elif isinstance(o, dict):
            stack.extend(o.values())
    return best


def flat_scalars(obj, prefix="", out=None, depth=0):
    out = {} if out is None else out
    if isinstance(obj, dict) and depth < 3:
        for k, v in obj.items():
            flat_scalars(v, f"{prefix}{k}.", out, depth + 1)
    elif isinstance(obj, (int, float, str)) and not isinstance(obj, bool):
        out[prefix.rstrip(".")] = obj
    return out


def paraphase_rows(obj):
    rows = []
    for region, v in obj.items():
        if not isinstance(v, dict):
            continue
        row = {"region": region}
        for k in ("total_cn", "gene_cn", "genome_depth", "region_depth"):
            if k in v and isinstance(v[k], (int, float, str, type(None))):
                row[k] = v[k]
        for k, val in v.items():
            if k.endswith("_cn") and k not in row and isinstance(val, (int, float, str, type(None))):
                row[k] = val
        fh = v.get("final_haplotypes")
        if isinstance(fh, (list, dict)):
            row["haplotypes"] = len(fh)
        rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------------------
# Plot helpers
# --------------------------------------------------------------------------------------
PAGE = (11.69, 8.27)
thousands = FuncFormatter(lambda v, _: f"{int(v):,}")


def fmt_int(x):
    return f"{int(x):,}" if pd.notna(x) else "-"


def safe_page(label, fn, *a, **kw):
    try:
        fn(*a, **kw)
    except Exception as e:
        warn(f"could not draw {label} page ({type(e).__name__}: {e})")
        plt.close("all")


def lens_area(r1, r2, d):
    if d >= r1 + r2:
        return 0.0
    if d <= abs(r1 - r2):
        return math.pi * min(r1, r2) ** 2
    a = r1 * r1 * math.acos((d * d + r1 * r1 - r2 * r2) / (2 * d * r1))
    b = r2 * r2 * math.acos((d * d + r2 * r2 - r1 * r1) / (2 * d * r2))
    c = 0.5 * math.sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))
    return a + b - c


def venn2(ax, only_a, both, only_b, name_a, name_b, col_a, col_b, col_s):
    """Area-proportional two-set Venn diagram."""
    ax.axis("off")
    na, nb = only_a + both, only_b + both
    if na == 0 or nb == 0:
        ax.text(0.5, 0.5, "no data", ha="center", transform=ax.transAxes)
        return
    k = math.pi / max(na, nb)
    ra, rb = math.sqrt(na * k / math.pi), math.sqrt(nb * k / math.pi)
    lo, hi = abs(ra - rb), ra + rb
    for _ in range(60):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if lens_area(ra, rb, mid) > both * k else (lo, mid)
    d = (lo + hi) / 2
    ax.add_patch(Circle((0, 0), ra, fc=col_a, ec=col_a, alpha=0.55, lw=1.5))
    ax.add_patch(Circle((d, 0), rb, fc=col_b, ec=col_b, alpha=0.55, lw=1.5))
    fs = 9
    ax.text((-ra + min(d - rb, ra)) / 2, 0, f"{only_a:,}", ha="center", va="center", fontsize=fs)
    ax.text((max(ra, d - rb) + d + rb) / 2, 0, f"{only_b:,}", ha="center", va="center", fontsize=fs)
    ax.text(((d - rb) + ra) / 2 if d - rb < ra else d / 2, 0, f"{both:,}", ha="center", va="center",
            fontsize=fs, fontweight="bold", color="black",
            bbox=dict(boxstyle="round,pad=0.2", fc=col_s, ec="none", alpha=0.7))
    top = max(ra, rb) + 0.12
    ax.text(-ra, top, name_a, ha="left", va="bottom", fontsize=10, color=col_a, fontweight="bold")
    ax.text(d + rb, top, name_b, ha="right", va="bottom", fontsize=10, color=col_b, fontweight="bold")
    ax.set_xlim(-ra - 0.1, d + rb + 0.1)
    ax.set_ylim(-max(ra, rb) - 0.1, max(ra, rb) + 0.45)
    ax.set_aspect("equal")


def stacked_fraction_bars(ax, labels, only_a, both, only_b, col_a, col_s, col_b, title):
    y = np.arange(len(labels))[::-1]
    tot = np.maximum(np.array(only_a) + both + np.array(only_b), 1)
    fa, fs, fb = np.array(only_a) / tot, np.array(both) / tot, np.array(only_b) / tot
    ax.barh(y, fa, color=col_a)
    ax.barh(y, fs, left=fa, color=col_s)
    ax.barh(y, fb, left=fa + fs, color=col_b)
    for yi, t in zip(y, np.array(only_a) + both + np.array(only_b)):
        ax.text(1.02, yi, f"n={int(t):,}", va="center", fontsize=7)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction", fontsize=8)
    ax.set_title(title, fontsize=9, loc="left")
    ax.tick_params(axis="x", labelsize=7)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def genome_strip(ax, genome, mito_n):
    """Chromosome strip. chrM is not to scale and coloured in the mitorsaw colour by its call count."""
    for i, c in enumerate(genome.names):
        off, w = genome.off_arr[i], genome.disp[i]
        if c == "chrM":
            if mito_n:
                alpha = 0.2 + 0.8 * min(mito_n, 20) / 20
                ax.add_patch(Rectangle((off, 0.1), w, 0.8, fc=COLORS["mitorsaw"], ec="black", lw=0.5, alpha=alpha))
            else:
                ax.add_patch(Rectangle((off, 0.1), w, 0.8, fc="#e0e0e0", ec="black", lw=0.5))
            ax.text(off + w / 2, 0.5, "M", ha="center", va="center", fontsize=6.5)
            if mito_n is not None:
                ax.text(off + w / 2, -0.05, f"mitorsaw: {mito_n}", ha="center", va="top", fontsize=6,
                        color=COLORS["mitorsaw"])
        else:
            ax.add_patch(Rectangle((off, 0.1), w, 0.8, fc="#d5d5d5" if i % 2 == 0 else "#b3b3b3", ec="none"))
            narrow = w < 0.022 * genome.total and len(c) > 4
            ax.text(off + w / 2, 0.5, c[3:], ha="center", va="center", fontsize=6.5, rotation=90 if narrow else 0)
    ax.set_ylim(-0.4, 1)
    ax.set_xlim(0, genome.total)
    ax.axis("off")


def _track_fig(genome, mito_n, title):
    fig = plt.figure(figsize=PAGE)
    gs = fig.add_gridspec(3, 1, height_ratios=[5, 1.7, 0.6], hspace=0.07, left=0.07, right=0.985, top=0.9, bottom=0.07)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    ax3 = fig.add_subplot(gs[2], sharex=ax1)
    fig.suptitle(title, x=0.02, ha="left", fontsize=11.5, fontweight="bold")
    for ax in (ax1, ax2):
        for off in genome.off_arr[1:]:
            ax.axvline(off, color="#e6e6e6", lw=0.5, zorder=0)
        ax.tick_params(axis="x", labelbottom=False, length=0)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    genome_strip(ax3, genome, mito_n)
    ax1.set_xlim(0, genome.total)
    return fig, ax1, ax2


def _window_counts(ax, genome, x, shared, win, col_tool, col_shared, ylabel):
    edges = np.arange(0, genome.total + win, win)
    sh, _ = np.histogram(x[shared], bins=edges)
    pr, _ = np.histogram(x[~shared], bins=edges)
    ax.bar(edges[:-1], sh, width=win, align="edge", color=col_shared, linewidth=0)
    ax.bar(edges[:-1], pr, width=win, bottom=sh, align="edge", color=col_tool, linewidth=0)
    ax.set_ylim(0, max(int((sh + pr).max()), 1) * 1.5)
    done = []
    for i in np.argsort(pr)[::-1]:
        if len(done) == 5 or pr[i] <= 0:
            break
        if any(abs(i - j) <= 3 for j in done):
            continue  # one label per hotspot
        done.append(i)
        c, p = genome.locate(edges[i] + win / 2)
        ax.text(edges[i] + win / 2, sh[i] + pr[i], f" {c[3:]}:{p / 1e6:.0f}Mb ({pr[i]:,})",
                rotation=90, ha="center", va="bottom", fontsize=6)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.tick_params(axis="y", labelsize=7)


def track_page_snv(pdf, genome, sample, tool, other, T, shared, win, args, mito_n):
    col_t, col_s = tool_color(tool), COLORS["shared"]
    x = genome.x(T["ch"].to_numpy(), T["pos"].to_numpy())
    vaf = T["vaf"].to_numpy()
    use_vaf = np.isfinite(vaf).mean() > 0.5
    y = vaf if use_vaf else np.log10(np.nan_to_num(T["qual"].to_numpy(), nan=0.0) + 1)
    rng = np.random.default_rng(1)
    idx_sh, idx_pr = np.flatnonzero(shared), np.flatnonzero(~shared)
    down = len(idx_sh) > args.max_points or len(idx_pr) > args.max_points
    if len(idx_sh) > args.max_points:
        idx_sh = rng.choice(idx_sh, args.max_points, replace=False)
    if len(idx_pr) > args.max_points:
        idx_pr = rng.choice(idx_pr, args.max_points, replace=False)
    fig, ax1, ax2 = _track_fig(
        genome, mito_n,
        f"{sample}: {tool} small variants along the genome  "
        f"(n={len(T):,}; shared with {other}: {int(shared.sum()):,}; {tool}-only: {int((~shared).sum()):,})")
    ax1.scatter(x[idx_sh], y[idx_sh], s=0.3, c=col_s, linewidths=0, alpha=0.5, rasterized=True, zorder=2)
    ax1.scatter(x[idx_pr], y[idx_pr], s=0.9, c=col_t, linewidths=0, alpha=0.8, rasterized=True, zorder=3)
    ax1.set_ylabel("variant allele fraction" if use_vaf else "log10(QUAL+1)", fontsize=9)
    if use_vaf:
        ax1.set_ylim(0, 1.03)
    ax1.legend(handles=[Line2D([], [], marker="o", ls="", color=col_t, label=f"{tool} only"),
                        Line2D([], [], marker="o", ls="", color=col_s, label=f"shared with {other}")],
               loc="upper right", fontsize=8, ncol=2, frameon=False, bbox_to_anchor=(1, 1.08))
    _window_counts(ax2, genome, x, shared, win, col_t, col_s, f"calls / {win / 1e6:g} Mb")
    note = ("tool-only points are drawn on top of the shared ones; top-5 windows by tool-only calls are labelled"
            + (f"; points downsampled to <= {args.max_points:,} per class" if down else ""))
    fig.text(0.07, 0.015, note, fontsize=7, color="#555555")
    pdf.savefig(fig)
    plt.close(fig)


def track_page_sv(pdf, genome, sample, tool, other, R, win, args, mito_n):
    col_t, col_s = tool_color(tool), COLORS["shared"]
    x = genome.x(R["ch"].to_numpy(), R["start"].to_numpy())
    shared = (R["mask"].to_numpy() == 3)
    typ = R["type"].to_numpy()
    size = R["size"].to_numpy().astype(float)
    y_bnd = 7.45
    y = np.where(typ == "BND", y_bnd, np.log10(np.clip(size, 10, None)))
    fig, ax1, ax2 = _track_fig(
        genome, mito_n,
        f"{sample}: {tool} SVs along the genome  (n={len(R):,}; shared with {other}: {int(shared.sum()):,}; "
        f"{tool}-only: {int((~shared).sum()):,})")
    for flag, col, z in ((True, col_s, 2), (False, col_t, 3)):
        for t in SV_TYPES + ["other"]:
            m = (shared == flag) & ((typ == t) if t != "other" else ~np.isin(typ, SV_TYPES))
            if m.any():
                mk = SV_MARKERS.get(t, "o")
                ax1.scatter(x[m], y[m], s=9 if mk != "x" else 12, marker=mk, c=col, alpha=0.8,
                            linewidths=0.6 if mk == "x" else 0, rasterized=True, zorder=z)
    ticks = [1.7, 2, 3, 4, 5, 6, 7, y_bnd]
    ax1.set_yticks(ticks)
    ax1.set_yticklabels(["50", "100", "1k", "10k", "100k", "1M", "10M", "BND"], fontsize=7)
    ax1.set_ylim(1.5, 7.75)
    ax1.set_ylabel("SV size (bp)", fontsize=9)
    h = [Line2D([], [], marker="o", ls="", color=col_t, label=f"{tool} only"),
         Line2D([], [], marker="o", ls="", color=col_s, label=f"shared with {other}")]
    h += [Line2D([], [], marker=SV_MARKERS[t], ls="", color="#444444", label=t, markersize=5) for t in SV_TYPES]
    ax1.legend(handles=h, loc="upper right", fontsize=7.5, ncol=7, frameon=False, bbox_to_anchor=(1, 1.08))
    _window_counts(ax2, genome, x, shared, win, col_t, col_s, f"SVs / {win / 1e6:g} Mb")
    fig.text(0.07, 0.015, "tool-only SVs are drawn on top of the shared ones; top-5 windows by tool-only SVs are labelled",
             fontsize=7, color="#555555")
    pdf.savefig(fig)
    plt.close(fig)


def heatmap_page(pdf, genome, sample, entries, binsize, title):
    """entries: [(name, chrom codes, positions)] for exactly two tools -> A, B, log2(A/B)"""
    code_m = genome.code.get("chrM", -1)
    rows = [i for i, c in enumerate(genome.names) if c != "chrM"]
    ncol = int(genome.len[rows].max() // binsize) + 1
    mats = []
    for _, ch, pos in entries:
        ch, pos = np.asarray(ch, dtype=np.int64), np.asarray(pos, dtype=np.int64)
        keep = ch != code_m
        flat = ch[keep] * ncol + pos[keep] // binsize
        mats.append(np.bincount(flat, minlength=len(genome.names) * ncol).reshape(len(genome.names), ncol)[rows])
    valid = np.arange(ncol)[None, :] * binsize < genome.len[rows][:, None]
    fig, axes = plt.subplots(3, 1, figsize=PAGE, sharex=True)
    fig.suptitle(title, x=0.02, ha="left", fontsize=11.5, fontweight="bold")
    panels = [(entries[0][0], np.log10(mats[0] + 1), "viridis", None, "log10(count + 1)"),
              (entries[1][0], np.log10(mats[1] + 1), "viridis", None, "log10(count + 1)"),
              (f"log2({entries[0][0]} / {entries[1][0]})", np.clip(np.log2((mats[0] + 1) / (mats[1] + 1)), -3, 3),
               "coolwarm", (-3, 3), "log2 ratio")]
    vmax = max(panels[0][1][valid].max(), panels[1][1][valid].max(), 1e-9)
    for ax, (name, m, cmap, lim, lab) in zip(axes, panels):
        mm = np.ma.masked_where(~valid, m)
        kw = dict(vmin=lim[0], vmax=lim[1]) if lim else dict(vmin=0, vmax=vmax)
        im = ax.imshow(mm, aspect="auto", cmap=cmap, interpolation="nearest", **kw)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([genome.names[i][3:] for i in rows], fontsize=6)
        ax.set_title(name, fontsize=9, loc="left")
        fig.colorbar(im, ax=ax, pad=0.01, fraction=0.025).set_label(lab, fontsize=7)
    step = max(1, int(50_000_000 // binsize))
    axes[-1].set_xticks(np.arange(0, ncol, step))
    axes[-1].set_xticklabels([f"{int(i * binsize / 1e6)}" for i in np.arange(0, ncol, step)], fontsize=7)
    axes[-1].set_xlabel(f"position (Mb), {binsize / 1e6:g} Mb windows", fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    pdf.savefig(fig)
    plt.close(fig)


def table_page(pdf, title, tables):
    fig, axes = plt.subplots(len(tables), 1, figsize=PAGE,
                             gridspec_kw={"height_ratios": [len(df) + 3 for _, df in tables]})
    axes = np.atleast_1d(axes)
    fig.suptitle(title, x=0.02, ha="left", fontsize=14, fontweight="bold")
    for ax, (sub, df) in zip(axes, tables):
        ax.axis("off")
        ax.set_title(sub, loc="left", fontsize=10)
        tbl = ax.table(cellText=df.astype(str).values.tolist(), colLabels=list(df.columns),
                       loc="upper left", cellLoc="center")
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(9)
        tbl.scale(1, 1.4)
        tbl.auto_set_column_width(list(range(len(df.columns))))
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    pdf.savefig(fig)
    plt.close(fig)


def small_table(ax, df, fontsize=7.5, title=None):
    ax.axis("off")
    if title:
        ax.set_title(title, loc="left", fontsize=9)
    if df is None or df.empty:
        ax.text(0, 0.9, "nothing to show", transform=ax.transAxes, fontsize=8)
        return
    tbl = ax.table(cellText=df.astype(str).values.tolist(), colLabels=list(df.columns), loc="upper left", cellLoc="left")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(fontsize)
    tbl.scale(1, 1.25)
    tbl.auto_set_column_width(list(range(len(df.columns))))


# --------------------------------------------------------------------------------------
# Comparison dashboards
# --------------------------------------------------------------------------------------
def snv_dashboard(pdf, sample, nA, nB, A, B, c):
    cA, cB, cS = tool_color(nA), tool_color(nB), COLORS["shared"]
    sa, sb = c["sa"], c["sb"]
    n_sh = int(sa.sum())
    fig = plt.figure(figsize=PAGE)
    gs = fig.add_gridspec(2, 3, hspace=0.5, wspace=0.4, left=0.06, right=0.97, top=0.88, bottom=0.2)
    fig.suptitle(f"{sample}: small variants, {nA} vs {nB}", x=0.02, ha="left", fontsize=13, fontweight="bold")

    venn2(fig.add_subplot(gs[0, 0]), len(A) - n_sh, n_sh, len(B) - n_sh, nA, nB, cA, cB, cS)

    oa = [int(((A["vtype"] == i) & ~sa).sum()) for i in range(4)]
    both = np.array([int(((A["vtype"] == i) & sa).sum()) for i in range(4)])
    ob = [int(((B["vtype"] == i) & ~sb).sum()) for i in range(4)]
    stacked_fraction_bars(fig.add_subplot(gs[0, 1]), VTYPES, oa, both, ob, cA, cS, cB,
                          f"by type: {nA}-only | shared | {nB}-only")

    ax = fig.add_subplot(gs[0, 2])
    w = 0.27
    for k, (lab, sub, col) in enumerate(((f"{nA} only", A["sub"].to_numpy()[~sa & (A["vtype"].to_numpy() == 0)], cA),
                                         ("shared", A["sub"].to_numpy()[sa & (A["vtype"].to_numpy() == 0)], cS),
                                         (f"{nB} only", B["sub"].to_numpy()[~sb & (B["vtype"].to_numpy() == 0)], cB))):
        cnt = np.bincount(sub[sub >= 0], minlength=6)
        ax.bar(np.arange(6) + (k - 1) * w, cnt / max(cnt.sum(), 1), width=w, color=col, label=lab)
    ax.set_xticks(range(6))
    ax.set_xticklabels(SUBS, fontsize=8)
    ax.set_title("substitution spectrum (SNPs, fraction)", fontsize=9, loc="left")
    ax.legend(fontsize=7, frameon=False)
    ax.tick_params(axis="y", labelsize=7)

    ax = fig.add_subplot(gs[1, 0])
    gA, gB = A["gt"].to_numpy()[c["ia"]], B["gt"].to_numpy()[c["ib"]]
    M = np.zeros((3, 3), dtype=np.int64)
    np.add.at(M, (gA, gB), 1)
    ax.imshow(np.log10(M + 1), cmap="Purples", aspect="auto")
    for i in range(3):
        for j in range(3):
            ax.text(j, i, f"{M[i, j]:,}", ha="center", va="center", fontsize=7,
                    color="white" if np.log10(M[i, j] + 1) > 0.55 * np.log10(M.max() + 1) else "black")
    ax.set_xticks(range(3))
    ax.set_yticks(range(3))
    ax.set_xticklabels(GTCLS, fontsize=8)
    ax.set_yticklabels(GTCLS, fontsize=8)
    ax.set_xlabel(nB, fontsize=8)
    ax.set_ylabel(nA, fontsize=8)
    ax.set_title(f"genotype of shared sites ({100 * np.trace(M) / max(M.sum(), 1):.1f}% identical)", fontsize=9, loc="left")

    ax = fig.add_subplot(gs[1, 1])
    bins = np.arange(-20.5, 21.5)
    for lab, il, col in ((f"{nA} only", A["ilen"].to_numpy()[~sa & np.isin(A["vtype"], (1, 2))], cA),
                         ("shared", A["ilen"].to_numpy()[sa & np.isin(A["vtype"], (1, 2))], cS),
                         (f"{nB} only", B["ilen"].to_numpy()[~sb & np.isin(B["vtype"], (1, 2))], cB)):
        h, _ = np.histogram(il, bins=bins)
        ax.step(np.arange(-20, 21), h / max(h.sum(), 1), where="mid", color=col, label=lab)
    ax.set_title("indel length (bp, + insertion / - deletion), fraction", fontsize=9, loc="left")
    ax.legend(fontsize=7, frameon=False)
    ax.tick_params(labelsize=7)

    ax = fig.add_subplot(gs[1, 2])
    use_vaf = np.isfinite(A["vaf"].to_numpy()).mean() > 0.5 and np.isfinite(B["vaf"].to_numpy()).mean() > 0.5
    for T, s, nm, col in ((A, sa, nA, cA), (B, sb, nB, cB)):
        v = T["vaf"].to_numpy() if use_vaf else np.log10(np.nan_to_num(T["qual"].to_numpy(), nan=0) + 1)
        rng = (0, 1) if use_vaf else (0, max(np.nanmax(v), 1))
        for flag, ls, lab in ((False, "-", "only"), (True, "--", "shared")):
            h, e = np.histogram(v[s == flag], bins=50, range=rng)
            ax.plot((e[:-1] + e[1:]) / 2, h / max(h.sum(), 1), ls=ls, color=col, label=f"{nm} {lab}")
    ax.set_title(("variant allele fraction" if use_vaf else "log10(QUAL+1)") + " (fraction of calls)", fontsize=9, loc="left")
    ax.legend(fontsize=6.5, frameon=False)
    ax.tick_params(labelsize=7)

    sA, sB = A["sub"].to_numpy(), B["sub"].to_numpy()
    lines = [
        f"{nA}: {len(A):,} calls, {100 * n_sh / max(len(A), 1):.1f}% also in {nB}   |   "
        f"{nB}: {len(B):,} calls, {100 * n_sh / max(len(B), 1):.1f}% also in {nA}   |   shared: {n_sh:,}",
        f"ts/tv   {nA}-only {tstv(sA[~sa]):.2f}   shared {tstv(sA[sa]):.2f}   {nB}-only {tstv(sB[~sb]):.2f}   "
        "(low ts/tv in private calls often means false positives)"]
    fig.text(0.06, 0.03, "\n".join(lines), fontsize=8)
    pdf.savefig(fig)
    plt.close(fig)


def sv_dashboard(pdf, sample, nA, nB, allsv, cl, args):
    cA, cB, cS = tool_color(nA), tool_color(nB), COLORS["shared"]
    fig = plt.figure(figsize=PAGE)
    gs = fig.add_gridspec(2, 3, hspace=0.5, wspace=0.4, left=0.06, right=0.97, top=0.88, bottom=0.2)
    fig.suptitle(f"{sample}: structural variants, {nA} vs {nB} (clusters of matching calls)", x=0.02, ha="left",
                 fontsize=13, fontweight="bold")
    mk = cl["mask"].to_numpy()
    oa, both, ob = int((mk == 1).sum()), int((mk == 3).sum()), int((mk == 2).sum())
    venn2(fig.add_subplot(gs[0, 0]), oa, both, ob, nA, nB, cA, cB, cS)

    types = [t for t in SV_TYPES if (cl["type"] == t).any()] + (["other"] if (~cl["type"].isin(SV_TYPES)).any() else [])
    def cnt(m, t):
        sel = (cl["type"] == t) if t != "other" else ~cl["type"].isin(SV_TYPES)
        return int((sel.to_numpy() & (mk == m)).sum())
    stacked_fraction_bars(fig.add_subplot(gs[0, 1]), types, [cnt(1, t) for t in types],
                          np.array([cnt(3, t) for t in types]), [cnt(2, t) for t in types], cA, cS, cB,
                          f"by type: {nA}-only | shared | {nB}-only")

    ax = fig.add_subplot(gs[0, 2])
    edges = [50, 100, 500, 1000, 5000, 50000, np.inf]
    labels = ["50-100", "100-500", "0.5-1k", "1-5k", "5-50k", ">50k"]
    size_any = cl[["size_0", "size_1"]].mean(axis=1).to_numpy()
    nonbnd = (cl["type"] != "BND").to_numpy()
    w = 0.27
    for k, (m, col, lab) in enumerate(((1, cA, f"{nA} only"), (3, cS, "shared"), (2, cB, f"{nB} only"))):
        h, _ = np.histogram(size_any[(mk == m) & nonbnd], bins=edges)
        ax.bar(np.arange(6) + (k - 1) * w, h, width=w, color=col, label=lab)
    ax.set_xticks(range(6))
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_title("SV size (bp), clusters", fontsize=9, loc="left")
    ax.legend(fontsize=7, frameon=False)
    ax.tick_params(axis="y", labelsize=7)

    ax = fig.add_subplot(gs[1, 0])
    sh = cl[(cl["mask"] == 3) & (cl["type"] != "BND")]
    for t, mkr in SV_MARKERS.items():
        s = sh[sh["type"] == t]
        if len(s) and t != "BND":
            ax.scatter(np.clip(s["size_0"], 10, None), np.clip(s["size_1"], 10, None), s=6, marker=mkr, label=t, alpha=0.6)
    lim = [40, max(float(sh[["size_0", "size_1"]].max().max()) if len(sh) else 1e4, 1e3) * 1.5]
    ax.plot(lim, lim, color="#888888", lw=0.8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(f"{nA} size (bp)", fontsize=8)
    ax.set_ylabel(f"{nB} size (bp)", fontsize=8)
    ax.set_title("size of shared SVs", fontsize=9, loc="left")
    ax.legend(fontsize=6.5, frameon=False)
    ax.tick_params(labelsize=7)

    ax = fig.add_subplot(gs[1, 1])
    off = (sh["start_0"] - sh["start_1"]).dropna().to_numpy()
    ax.hist(np.clip(off, -args.max_dist, args.max_dist), bins=41, color=cS)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_title(f"start({nA}) - start({nB}) of shared SVs (bp)", fontsize=9, loc="left")
    ax.tick_params(labelsize=7)

    ax = fig.add_subplot(gs[1, 2])
    bins = np.linspace(1.5, 7, 45)
    for i, (nm, col) in enumerate(((nA, cA), (nB, cB))):
        s = allsv[(allsv["tool"] == i) & (allsv["type"] != "BND") & (allsv["size"] > 0)]["size"].to_numpy()
        h, _ = np.histogram(np.log10(s), bins=bins)
        ax.step((bins[:-1] + bins[1:]) / 2, h, where="mid", color=col, label=f"{nm} (n={len(s):,})")
    ax.set_xticks([2, 3, 4, 5, 6, 7])
    ax.set_xticklabels(["100", "1k", "10k", "100k", "1M", "10M"], fontsize=7)
    ax.set_yscale("log")
    ax.set_title("SV size distribution, all calls", fontsize=9, loc="left")
    ax.legend(fontsize=7, frameon=False)
    ax.tick_params(axis="y", labelsize=7)

    nrec = [int((allsv["tool"] == i).sum()) for i in (0, 1)]
    fig.text(0.06, 0.03,
             f"{nA}: {nrec[0]:,} calls    {nB}: {nrec[1]:,} calls    clusters: {nA}-only {oa:,}, shared {both:,}, {nB}-only {ob:,}\n"
             f"matching: same type and chromosome, breakpoints <= {args.max_dist} bp or reciprocal overlap >= "
             f"{args.reciprocal_overlap} (INS: position and size ratio); a BND pair counts once", fontsize=8)
    pdf.savefig(fig)
    plt.close(fig)


# --------------------------------------------------------------------------------------
# Special tool pages
# --------------------------------------------------------------------------------------
STATUS_NAMES = {-1: "no status", 0: "normal", 1: "pre_mutation", 2: "full_mutation"}
STATUS_COLORS = {-1: COLORS["trgt"], 0: "#A9B7D0", 1: "#F58518", 2: "#D62728"}


def trgt_page(pdf, genome, sample, T, mito_n, args, annotated):
    x = genome.x(T["ch"].to_numpy(), T["pos"].to_numpy())
    y = np.log10(np.clip(T["al_max"].to_numpy(), 1, None))
    st = T["status"].to_numpy()
    fig = plt.figure(figsize=PAGE)
    gs = fig.add_gridspec(3, 2, height_ratios=[3.4, 0.55, 3], hspace=0.35, wspace=0.12,
                          left=0.07, right=0.985, top=0.9, bottom=0.05)
    fig.suptitle(f"{sample}: TRGT repeat loci ({len(T):,}); " +
                 ("STR_STATUS from stranger" if annotated else "no stranger annotation found, loci not classified"),
                 x=0.02, ha="left", fontsize=11.5, fontweight="bold")
    ax1 = fig.add_subplot(gs[0, :])
    ax3 = fig.add_subplot(gs[1, :], sharex=ax1)
    rng = np.random.default_rng(1)
    for code in (-1, 0, 1, 2):
        idx = np.flatnonzero((st == code) & np.isfinite(y))
        if not len(idx):
            continue
        if len(idx) > args.max_points:
            idx = rng.choice(idx, args.max_points, replace=False)
        ax1.scatter(x[idx], y[idx], s=0.6 if code <= 0 else 14, c=STATUS_COLORS[code], linewidths=0,
                    alpha=0.5 if code <= 0 else 0.95, rasterized=True, zorder=2 + max(code, 0),
                    label=f"{STATUS_NAMES[code]} ({int((st == code).sum()):,})")
    for off in genome.off_arr[1:]:
        ax1.axvline(off, color="#e6e6e6", lw=0.5, zorder=0)
    ax1.set_xlim(0, genome.total)
    ax1.set_yticks([1, 2, 3, 4])
    ax1.set_yticklabels(["10", "100", "1k", "10k"], fontsize=7)
    ax1.set_ylabel("longest allele (bp)", fontsize=9)
    ax1.tick_params(axis="x", labelbottom=False, length=0)
    ax1.legend(fontsize=7.5, frameon=False, ncol=4, loc="upper right", markerscale=2, bbox_to_anchor=(1, 1.08))
    for sp in ("top", "right"):
        ax1.spines[sp].set_visible(False)
    genome_strip(ax3, genome, mito_n)

    ax = fig.add_subplot(gs[2, 0])
    called = int(T["called"].sum())
    het = int((T["called"] & T["het"]).sum())
    labels = ["called", "no call", "het", "hom"]
    vals = [called, len(T) - called, het, called - het]
    cols = [COLORS["trgt"], "#bbbbbb", "#7f7f7f", "#4d4d4d"]
    ax.barh(range(4)[::-1], vals, color=cols)
    ax.set_yticks(range(4)[::-1])
    ax.set_yticklabels(labels, fontsize=8)
    for i, v in zip(range(4)[::-1], vals):
        ax.text(v, i, f" {v:,} ({100 * v / max(len(T), 1):.1f}%)", va="center", fontsize=7)
    ax.set_xlim(0, max(vals) * 1.35)
    ax.set_title("genotyping success and zygosity", fontsize=9, loc="left")
    ax.tick_params(axis="x", labelsize=7)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    ax = fig.add_subplot(gs[2, 1])
    if annotated and (st >= 1).any():
        sel = T[st >= 1].sort_values(["status", "al_max"], ascending=False).head(14)
        title = "loci above the pre-mutation threshold (stranger)"
    else:
        sel = T[T["called"]].sort_values("al_max", ascending=False).head(14)
        title = "longest alleles" + ("" if annotated else " (no stranger annotation)")
    tab = pd.DataFrame({
        "locus": sel["trid"].str.slice(0, 22), "pos": [f"{genome.names[c][3:]}:{p:,}" for c, p in zip(sel["ch"], sel["pos"])],
        "longest (bp)": sel["al_max"].round(0).astype("Int64"), "units": sel["mc_max"].round(0).astype("Int64"),
        "status": [STATUS_NAMES[s] for s in sel["status"]], "disease": sel["disease"].str.slice(0, 22)})
    small_table(ax, tab, 6.5, title)
    pdf.savefig(fig)
    plt.close(fig)


def mito_gene(pos):
    for name, a, b in MT_GENES:
        if a <= pos <= b:
            return name
    return "tRNA/other"


def mito_page(pdf, sample, V, cov, hap_scalars, args):
    fig = plt.figure(figsize=PAGE)
    gs = fig.add_gridspec(4, 2, height_ratios=[1.6, 0.8, 2, 2.4], hspace=0.55, wspace=0.12,
                          left=0.07, right=0.985, top=0.9, bottom=0.05)
    fig.suptitle(f"{sample}: mitorsaw chrM variants ({len(V)}); VCF FORMAT tags: {V.attrs.get('format_tags', '?')}",
                 x=0.02, ha="left", fontsize=11, fontweight="bold")
    col = COLORS["mitorsaw"]
    axc = fig.add_subplot(gs[0, :])
    if cov is not None:
        axc.fill_between(np.arange(1, len(cov) + 1), cov, color="#9bd29b", lw=0)
        axc.set_ylabel("coverage", fontsize=8)
    else:
        axc.text(0.5, 0.5, "no per-position coverage found in coverage_stats.json", ha="center", transform=axc.transAxes, fontsize=8)
    axc.set_xlim(1, 16569)
    axc.tick_params(labelsize=7, labelbottom=False)
    axg = fig.add_subplot(gs[1, :], sharex=axc)
    for i, (name, a, b) in enumerate(MT_GENES):
        row = i % 2
        axg.add_patch(Rectangle((a, row * 0.5), b - a, 0.42, fc="#d9ead9" if name != "D-loop" else "#e0e0e0", ec="#555555", lw=0.4))
        axg.text((a + b) / 2, row * 0.5 + 0.21, name, ha="center", va="center", fontsize=6)
    axg.set_ylim(0, 1)
    axg.axis("off")
    axv = fig.add_subplot(gs[2, :], sharex=axc)
    if len(V):
        af = V["af"].to_numpy(dtype=float)
        has_af = np.isfinite(af).any()
        y = np.where(np.isfinite(af), af, 1.0)
        axv.vlines(V["pos"], 0, y, color=col, lw=1)
        axv.scatter(V["pos"], y, s=22, c=np.where(y < 0.95, "#ffffff", col), edgecolors=col, zorder=3)
        axv.set_ylim(0, 1.08)
        axv.set_ylabel("heteroplasmy fraction" if has_af else "variant (no AF field)", fontsize=8)
        axv.text(0.995, 0.95, "open = heteroplasmic (< 0.95)", transform=axv.transAxes, ha="right", fontsize=7)
    axv.set_xlabel("chrM position", fontsize=8)
    axv.tick_params(labelsize=7)
    for ax in (axc, axv):
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axt = fig.add_subplot(gs[3, 0])
    top = V.sort_values("pos").head(16)
    tab = pd.DataFrame({"pos": top["pos"].astype("Int64"), "change": top["ref"] + ">" + top["alt"].str.slice(0, 12),
                        "AF": top["af"].round(3), "DP": top["dp"].astype("Int64"), "gene": [mito_gene(p) for p in top["pos"]]})
    small_table(axt, tab, 7, "variants (first 16 by position)")
    axs = fig.add_subplot(gs[3, 1])
    items = list(hap_scalars.items())[:14]
    small_table(axs, pd.DataFrame({"hap_stats.json": [k[:34] for k, _ in items], "value": [str(v)[:24] for _, v in items]}) if items else None,
                7, "haplotype statistics (as found in the JSON)")
    pdf.savefig(fig)
    plt.close(fig)


def paraphase_page(pdf, sample, P):
    fig = plt.figure(figsize=PAGE)
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 2.2], wspace=0.15, left=0.06, right=0.985, top=0.86, bottom=0.08)
    cn = pd.to_numeric(P.get("total_cn"), errors="coerce") if "total_cn" in P else pd.Series(dtype=float)
    def nice(v):
        return int(v) if isinstance(v, float) and v.is_integer() else v
    smn = {c: nice(P[c].dropna().iloc[0]) for c in P.columns if c.startswith("smn") and c.endswith("_cn") and P[c].notna().any()}
    fig.suptitle(f"{sample}: paraphase, {len(P)} regions" + ("   |   " + ", ".join(f"{k}={v}" for k, v in smn.items()) if smn else ""),
                 x=0.02, ha="left", fontsize=11.5, fontweight="bold")
    ax = fig.add_subplot(gs[0])
    if cn.notna().any():
        counts = cn.dropna().astype(int).value_counts().sort_index()
        ax.bar(counts.index.astype(str), counts.values, color=[COLORS["paraphase"] if i != 2 else "#bbbbbb" for i in counts.index])
        for i, v in enumerate(counts.values):
            ax.text(i, v, str(v), ha="center", va="bottom", fontsize=7)
        ax.set_xlabel("total copy number of the region", fontsize=8)
        ax.set_ylabel("regions", fontsize=8)
        ax.set_title("copy number distribution (grey = 2 copies)", fontsize=9, loc="left")
    else:
        ax.axis("off")
        ax.text(0, 0.9, "no 'total_cn' field found in the JSON", transform=ax.transAxes, fontsize=8)
    ax.tick_params(labelsize=7)
    ax2 = fig.add_subplot(gs[1])
    sel = P[cn.reindex(P.index).ne(2)] if cn.notna().any() else P
    cols = [c for c in P.columns][:8]
    show = sel[cols].head(30).astype(object).map(lambda v: "-" if pd.isna(v) else nice(v))
    small_table(ax2, show, 7, "regions with total_cn != 2 (first 30)" if cn.notna().any() else "regions")
    pdf.savefig(fig)
    plt.close(fig)


def text_page(pdf, sample, found, genome_src, args):
    lines = [f"Variant comparison report: {sample}", f"generated {dt.datetime.now():%Y-%m-%d %H:%M}", "", "Input files:"]
    for name, (kind, path) in found.items():
        lines += textwrap.wrap(f"{name:13s} [{kind:7s}] {path}", 150, subsequent_indent=" " * 24)
    lines += ["", f"Genome: contig lengths from {genome_src}; chr1-22,X,Y,M only; chrM drawn not to scale",
              f"SNV pair: {args.snv_pair[0]} vs {args.snv_pair[1]}   SV pair: {args.sv_pair[0]} vs {args.sv_pair[1]}",
              f"PASS-only={not args.no_pass_filter}  SV min size={args.min_sv_size}  max breakpoint distance={args.max_dist}  "
              f"reciprocal overlap={args.reciprocal_overlap}", "", f"Warnings ({len(WARNINGS)}):"]
    for w in WARNINGS[:30]:
        lines += textwrap.wrap(w, 150, initial_indent="  - ", subsequent_indent="    ")
    if len(WARNINGS) > 30:
        lines.append(f"  ... {len(WARNINGS) - 30} more (see terminal)")
    fig = plt.figure(figsize=PAGE)
    fig.text(0.04, 0.96, "\n".join(lines), va="top", ha="left", family="monospace", fontsize=7.5)
    pdf.savefig(fig)
    plt.close(fig)


# --------------------------------------------------------------------------------------
# Per sample
# --------------------------------------------------------------------------------------
def process_sample(sample, args):
    registry = {n: (k, d, p) for n, k, d, p in REGISTRY}
    needed = set(args.snv_pair) | set(args.sv_pair) | {"trgt", "mitorsaw", "paraphase"}
    found = {}
    for name in needed:
        if name not in registry:
            continue
        kind, dirpat, pats = registry[name]
        path = find_file(sample, args.variants_dir, dirpat, pats)
        if path is None:
            warn(f"{sample}: {name}: not found ({dirpat.replace('{s}', sample)}/{'|'.join(pats).replace('{s}', sample)})")
        else:
            found[name] = (kind, path)
    for spec in args.add:
        m = re.match(r"^([^:=]+):(snv|sv)=(.+)$", spec)
        if not m:
            warn(f"ignoring --add '{spec}' (expected NAME:KIND=PATH with KIND = snv|sv)")
        elif not Path(m.group(3)).is_file():
            warn(f"{sample}: {m.group(1)}: file not found: {m.group(3)}")
        else:
            found[m.group(1)] = (m.group(2), Path(m.group(3)))
    annotated_trgt = None
    if "trgt" in found:
        annotated_trgt = find_file(sample, args.annotated_dir, STRANGER[0], STRANGER[1])
        if annotated_trgt is None:
            warn(f"{sample}: stranger-annotated TRGT VCF not found, TRGT loci will not be classified")
    if not found:
        warn(f"{sample}: nothing found, no report written")
        return

    vcfs = [p for n, (k, p) in found.items() if str(p).endswith((".vcf", ".vcf.gz"))]
    lengths, genome_src = load_lengths(args, vcfs)
    genome = Genome(lengths)
    info(f"[{sample}] genome: {len(genome.names)} contigs from {genome_src}")

    def worker(task):
        what, name, path = task
        info(f"[{sample}] {what}: {name}")
        if what == "snv":
            return snv_table(path, args, genome)
        if what == "sv":
            return read_svs(path, args, genome)
        if what == "trgt":
            return read_trgt(path, genome, args)
        if what == "mito":
            return read_mito(path, args)
        if what == "para":
            return paraphase_rows(load_json(path))
        return bcftools_stats(path, args)

    tasks = []
    for n in args.snv_pair:
        if n in found:
            tasks.append(("snv", n, found[n][1]))
    for n in args.sv_pair:
        if n in found:
            tasks.append(("sv", n, found[n][1]))
    if "trgt" in found:
        tasks.append(("trgt", "trgt", annotated_trgt or found["trgt"][1]))
    if "mitorsaw" in found:
        tasks.append(("mito", "mitorsaw", found["mitorsaw"][1]))
    if "paraphase" in found:
        tasks.append(("para", "paraphase", found["paraphase"][1]))
    if not args.no_stats:
        for n in list(args.snv_pair) + ["trgt", "mitorsaw"]:
            if n in found:
                tasks.append(("stats", n, found[n][1]))
    results = {}
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        futs = [(t, ex.submit(worker, t)) for t in tasks]
        for t, fut in futs:
            try:
                results[(t[0], t[1])] = fut.result()
            except Exception as e:
                warn(f"{sample}: {t[1]}: {t[0]} step failed, skipped ({e})")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    mito = results.get(("mito", "mitorsaw"))
    mito_n = len(mito) if mito is not None else None

    nA, nB = args.snv_pair
    A, B = results.get(("snv", nA)), results.get(("snv", nB))
    snv_cmp = compare_snv(A, B) if A is not None and B is not None else None
    if snv_cmp is None:
        warn(f"{sample}: SNV comparison {nA} vs {nB} skipped (both VCFs are needed)")
    sA, sB = args.sv_pair
    SA, SB = results.get(("sv", sA)), results.get(("sv", sB))
    allsv = cl = None
    if SA is not None and SB is not None:
        allsv = sv_cluster({sA: SA[0], sB: SB[0]}, args)
        cl = sv_clusters(allsv)
    else:
        warn(f"{sample}: SV comparison {sA} vs {sB} skipped (both VCFs are needed)")

    stat_rows = []
    for n in list(args.snv_pair) + ["trgt", "mitorsaw"]:
        s = results.get(("stats", n))
        if s is None:
            continue
        het, hom = s.get("het"), s.get("homalt")
        stat_rows.append({"tool": n, "records": s.get("number of records"), "SNPs": s.get("number of SNPs"),
                          "indels": s.get("number of indels"), "MNPs": s.get("number of MNPs"),
                          "other": s.get("number of others"), "multiallelic": s.get("number of multiallelic sites"),
                          "ts/tv": s.get("ts/tv", "-"), "het": het, "hom-alt": hom,
                          "het/hom": round(het / hom, 2) if het is not None and hom else "-"})
    stats_df = pd.DataFrame(stat_rows)
    sv_rows = []
    for n, r in ((sA, SA), (sB, SB)):
        if r is None:
            continue
        df, npass, nsmall = r
        row = {"tool": n, "SVs kept": len(df)}
        for t in SV_TYPES:
            row[t] = int((df["type"] == t).sum())
        row["other"] = int((~df["type"].isin(SV_TYPES)).sum())
        nb = df[(df["type"] != "BND") & (df["size"] > 0)]["size"]
        row["median size"] = int(nb.median()) if len(nb) else "-"
        row["dropped non-PASS"], row["dropped < min size"] = npass, nsmall
        sv_rows.append(row)
    sv_df = pd.DataFrame(sv_rows)

    pdf_path = outdir / f"{sample}_variant_overlap.pdf"
    with PdfPages(pdf_path) as pdf:
        text_page(pdf, sample, found, genome_src, args)
        tabs = []
        if not stats_df.empty:
            show = stats_df.copy()
            for c in ["records", "SNPs", "indels", "MNPs", "other", "multiallelic", "het", "hom-alt"]:
                show[c] = show[c].map(fmt_int)
            tabs.append(("bcftools stats (PASS only unless --no-pass-filter)", show))
            stats_df.to_csv(outdir / f"{sample}_small_variant_stats.tsv", sep="\t", index=False)
        if not sv_df.empty:
            show = sv_df.copy()
            show["SVs kept"] = show["SVs kept"].map(fmt_int)
            tabs.append((f"SV callers (PASS, size >= {args.min_sv_size} bp)", show))
            sv_df.to_csv(outdir / f"{sample}_sv_summary.tsv", sep="\t", index=False)
        if tabs:
            safe_page("overview tables", table_page, pdf, f"{sample}: overview", tabs)

        if snv_cmp is not None:
            safe_page("SNV dashboard", snv_dashboard, pdf, sample, nA, nB, A, B, snv_cmp)
            if not args.no_tracks:
                safe_page(f"{nA} track", track_page_snv, pdf, genome, sample, nA, nB, A, snv_cmp["sa"],
                          args.window_snv * 1e6, args, mito_n)
                safe_page(f"{nB} track", track_page_snv, pdf, genome, sample, nB, nA, B, snv_cmp["sb"],
                          args.window_snv * 1e6, args, mito_n)
            if not args.no_heatmaps:
                safe_page("SNV heatmaps", heatmap_page, pdf, genome, sample,
                          [(nA, A["ch"], A["pos"]), (nB, B["ch"], B["pos"])], HEAT_BIN_SNV,
                          f"{sample}: small variant density ({nA} vs {nB}), {HEAT_BIN_SNV / 1e6:g} Mb windows")
            pd.DataFrame({"combination": [f"{nA} only", "shared", f"{nB} only"],
                          "n_variants": [int((~snv_cmp['sa']).sum()), int(snv_cmp['sa'].sum()), int((~snv_cmp['sb']).sum())]}
                         ).to_csv(outdir / f"{sample}_snv_overlap.tsv", sep="\t", index=False)

        if cl is not None:
            safe_page("SV dashboard", sv_dashboard, pdf, sample, sA, sB, allsv, cl, args)
            if not args.no_tracks:
                for i, (nm, other) in enumerate(((sA, sB), (sB, sA))):
                    safe_page(f"{nm} SV track", track_page_sv, pdf, genome, sample, nm, other,
                              allsv[allsv["tool"] == i], args.window_sv * 1e6, args, mito_n)
            if not args.no_heatmaps:
                ents = [(nm, allsv[allsv["tool"] == i]["ch"], allsv[allsv["tool"] == i]["start"])
                        for i, nm in enumerate((sA, sB))]
                safe_page("SV heatmaps", heatmap_page, pdf, genome, sample, ents, HEAT_BIN_SV,
                          f"{sample}: SV density ({sA} vs {sB}), {HEAT_BIN_SV / 1e6:g} Mb windows")
            mk = cl["mask"].to_numpy()
            pd.DataFrame({"combination": [f"{sA} only", "shared", f"{sB} only"],
                          "n_clusters": [int((mk == 1).sum()), int((mk == 3).sum()), int((mk == 2).sum())]}
                         ).to_csv(outdir / f"{sample}_sv_overlap.tsv", sep="\t", index=False)

        T = results.get(("trgt", "trgt"))
        if T is not None and len(T):
            safe_page("TRGT", trgt_page, pdf, genome, sample, T, mito_n, args, annotated_trgt is not None)
        if mito is not None:
            cov = hap = None
            mpath = found["mitorsaw"][1]
            try:
                cj = sorted(mpath.parent.rglob("coverage_stats.json"))
                cov = longest_numeric_list(load_json(cj[0])) if cj else None
            except Exception as e:
                warn(f"{sample}: mitorsaw coverage_stats.json not readable ({e})")
            try:
                hj = sorted(mpath.parent.glob("*hap_stats.json"))
                hap = flat_scalars(load_json(hj[0])) if hj else {}
            except Exception as e:
                warn(f"{sample}: mitorsaw hap_stats.json not readable ({e})")
                hap = {}
            safe_page("mitorsaw", mito_page, pdf, sample, mito, cov, hap or {}, args)
        P = results.get(("para", "paraphase"))
        if P is not None and len(P):
            safe_page("paraphase", paraphase_page, pdf, sample, P)
            P.to_csv(outdir / f"{sample}_paraphase_regions.tsv", sep="\t", index=False)
        elif P is not None:
            warn(f"{sample}: paraphase JSON has no per-region dictionaries, nothing drawn")
    info(f"[{sample}] wrote {pdf_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-d", "--variants-dir", action="append", required=True, metavar="DIR",
                    help="directory with the <tool>_<sample> folders (repeatable)")
    ap.add_argument("-s", "--sample", nargs="+", required=True, help="sample name(s); one PDF per sample")
    ap.add_argument("-o", "--outdir", default=".", help="output directory (default .)")
    ap.add_argument("--annotated-dir", action="append", metavar="DIR",
                    help="directory with stranger_<sample>/ (default: ../annotated_variants next to --variants-dir)")
    ap.add_argument("--ref", help="reference FASTA (left-aligns indels; <ref>.fai gives the contig lengths)")
    ap.add_argument("--fai", help="contig lengths from this .fai (otherwise VCF headers or built-in GRCh38)")
    ap.add_argument("--snv-pair", nargs=2, default=["deepvariant", "nanocaller"], metavar=("A", "B"))
    ap.add_argument("--sv-pair", nargs=2, default=["sawfish", "sniffles"], metavar=("A", "B"))
    ap.add_argument("--add", action="append", default=[], metavar="NAME:KIND=PATH",
                    help="extra VCF usable in a pair, KIND = snv|sv (repeatable)")
    ap.add_argument("--no-pass-filter", action="store_true", help="also use records with FILTER != PASS/.")
    ap.add_argument("--min-sv-size", type=int, default=50)
    ap.add_argument("--reciprocal-overlap", type=float, default=0.5)
    ap.add_argument("--max-dist", type=int, default=500, help="SV breakpoint distance in bp (default 500)")
    ap.add_argument("--window-snv", type=float, default=2.0, help="window for the SNV count bars in Mb (default 2)")
    ap.add_argument("--window-sv", type=float, default=5.0, help="window for the SV count bars in Mb (default 5)")
    ap.add_argument("--max-points", type=int, default=1_500_000, help="max plotted points per class (default 1.5M)")
    ap.add_argument("--jobs", type=int, default=4, help="parallel bcftools jobs (default 4)")
    ap.add_argument("--no-stats", action="store_true", help="skip the bcftools stats table")
    ap.add_argument("--no-tracks", action="store_true", help="skip the genome track pages")
    ap.add_argument("--no-heatmaps", action="store_true", help="skip the density heatmaps")
    args = ap.parse_args()
    if not args.annotated_dir:
        args.annotated_dir = [str(Path(d).resolve().parent / "annotated_variants") for d in args.variants_dir]
    if shutil.which("bcftools") is None:
        sys.exit("bcftools not found in PATH (activate the conda env that has it)")
    for d in args.variants_dir:
        if not Path(d).is_dir():
            warn(f"variants dir does not exist: {d}")
    for sample in args.sample:
        process_sample(sample, args)


if __name__ == "__main__":
    main()