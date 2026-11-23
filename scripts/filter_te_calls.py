#!/usr/bin/env python3
"""Filter sTELLeR non-reference TE calls into a high-confidence set + an annotated table.

usage:
  filter_te_calls.py <steller_repeats.vcf[.gz]> <out_prefix> <bam|none> <sex: male|female|unknown|auto> [options]

Writes:
  <out_prefix>.high_confidence.vcf   calls that pass every hard filter
  <out_prefix>.annotated.tsv         ALL calls with evidence, pass/fail and the reasons
                                     (result is column 20, sv_match is column 17)

Evidence layers (each optional except the VCF):
- sTELLeR call itself (family, subfamily, length, genotype)
- BAM (haplotagged or not): reads with an insertion CIGAR or soft/hard clip near the site,
    local depth, allele fraction, median insertion length, supporting reads per HP tag,
    HP-tag balance for genotyping (het calls concentrate on one haplotype)
- --sv-vcf: Sniffles/Sawfish insertions near the site
- --exclude-bed: regions to distrust (e.g. blacklist, centromeres)
- sex (male/female/auto) for hemizygous-chromosome genotype checks (GRCh38 PAR coordinates)

Hard filters (-> not in the high-confidence set): oversize / short, old Alu subfamily,
no or low BAM support, low allele fraction, excluded region, het call on a hemizygous
chromosome, chrY call in a female, duplicate of a better call.
Soft flags (only reported): length discordant between sTELLeR and BAM, genotype
disagreement (either direction), low depth, no SV-caller match, mostly-untagged or
ambiguous HP tags, HP genotype conflicting with allele fraction.
"""
import argparse
import bisect
import gzip
import re
import statistics
import sys
from collections import Counter, defaultdict

PAR_GRCH38 = {"X": [(10001, 2781479), (155701383, 156030895)],
              "Y": [(10001, 2781479), (56887903, 57217415)]}
DEFAULT_MAX = {"Alu": 700, "L1": 7000, "SVA": 4500, "ERV": 10000}
DEFAULT_MIN = {"Alu": 200, "L1": 200, "SVA": 300, "ERV": 500}
DEFAULT_SUB = {"Alu": r"^AluY"}
PRIMARY = re.compile(r"^(chr)?([0-9]{1,2}|X|Y|M|MT)$")
FETCH_MARGIN = 1000   # extra slack when merging BAM fetch windows (covers long reads/clips)
HP_HOM_MIN_AF = 0.4   # min AF to believe an HP-based 1/1 when tags disagree with dosage

def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def norm_chr(c):
    return c[3:] if c.lower().startswith("chr") else c

def chr_key(c):
    n = norm_chr(c)
    if n.isdigit():
        return (0, int(n), "")
    return (1, {"X": 1, "Y": 2, "M": 3, "MT": 3}.get(n, 9), n)

def kv_list(items, default, cast):
    d = dict(default)
    for item in items or []:
        k, v = item.split("=", 1)
        d[k] = cast(v)
    return d

# ----------------------------------------------------------------- input parsing
def read_steller(path):
    calls, sample, skipped = [], "sample", 0
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                sample = f[9] if len(f) > 9 else sample
                continue
            try:
                info = dict(x.partition("=")[::2] for x in f[7].split(";"))
                parts = info["MEINFO"].split(",")
                fam, _, sub = parts[0].partition(":")
                start, end = int(parts[1]), int(parts[2])
                keys = f[8].split(":")
                vals = f[9].split(":")
                # sTELLeR writes TE as "Family:Subfamily", so the sample column has one more
                # ':'-separated value than FORMAT has keys; merge them (-> Family_Subfamily)
                if len(vals) > len(keys):
                    k = len(vals) - len(keys) + 1
                    vals = [":".join(vals[:k]).replace(":", "_")] + vals[k:]
                sd = dict(zip(keys, vals))
                gt = sd.get("GT", "./.")
                # the VCF spec wants GT first in FORMAT (sTELLeR writes it last)
                order = (["GT"] if "GT" in sd else []) + [k for k in keys if k != "GT"]
                fmt_out, sample_out = ":".join(order), ":".join(sd[k] for k in order)
                calls.append(dict(chrom=f[0], pos=int(f[1]), ref=f[3].upper(), alt=f[4],
                                  meinfo=info["MEINFO"], fam=fam, sub=sub or fam,
                                  tlen=end - start, gt=gt, fmt=fmt_out, sample_field=sample_out,
                                  sd=sd, order=order))
            except (IndexError, KeyError, ValueError):
                skipped += 1
                sys.stderr.write(f"skipping malformed line: {line[:60].rstrip()}\n")
    return calls, sample, skipped

def read_sv_ins(path, min_len=50):
    """Insertions from an SV VCF -> {chrom: sorted list of (pos, length, id)}."""
    out = defaultdict(list)
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            info = dict(x.partition("=")[::2] for x in f[7].split(";"))
            length = None
            if "SVLEN" in info:
                try:
                    length = abs(int(float(info["SVLEN"].split(",")[0])))
                except ValueError:
                    pass
            alt = f[4]
            if length is None and not alt.startswith("<") and len(alt) > len(f[3]):
                length = len(alt) - len(f[3])
            is_ins = info.get("SVTYPE") == "INS" or alt in ("<INS>",) or \
                (not alt.startswith("<") and len(alt) > len(f[3]))
            if is_ins and length and length >= min_len:
                out[norm_chr(f[0])].append((int(f[1]), length, f[2]))
    for c in out:
        out[c].sort()
    return out

def read_bed(path):
    """BED regions -> {chrom: sorted, *merged* list of (start, end)}."""
    raw = defaultdict(list)
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            f = line.split("\t")
            raw[norm_chr(f[0])].append((int(f[1]), int(f[2])))
    out = defaultdict(list)
    for c, iv in raw.items():
        iv.sort()
        merged = []
        for s, e in iv:
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        out[c] = [(s, e) for s, e in merged]
    return out

def in_bed(bed, chrom, pos):
    iv = bed.get(norm_chr(chrom), [])
    i = bisect.bisect_right(iv, (pos, float("inf"))) - 1
    return i >= 0 and iv[i][0] < pos <= iv[i][1]

# ----------------------------------------------------------------- BAM evidence
def plan_bam_fetches(calls, window):
    """Merge overlapping fetch windows so clustered calls share one BAM fetch.

    Sets c["region"] on every call; returns a list of [chrom, start, end, reads, failed].
    Extra reads pulled in by the margin are harmless: bam_evidence still applies
    position-specific spanning/insertion checks per call. Calls must be sorted;
    main() releases each region's reads once the loop moves past it.
    """
    regions = []
    for c in calls:
        s = max(0, c["pos"] - 1 - window)
        e = c["pos"] + window
        if regions and regions[-1][0] == c["chrom"] and s <= regions[-1][2] + FETCH_MARGIN:
            regions[-1][2] = max(regions[-1][2], e)
            c["region"] = len(regions) - 1
        else:
            regions.append([c["chrom"], s, e, None, False])
            c["region"] = len(regions) - 1
    return regions

def region_reads(bam, regions, idx):
    r = regions[idx]
    if r[3] is None and not r[4]:
        try:
            r[3] = list(bam.fetch(r[0], r[1], r[2]))
        except ValueError:
            r[3], r[4] = [], True
    return None if r[4] else r[3]

def bam_evidence(call, a, reads):
    pos = call["pos"]
    support, spanning, ins_lens, hp = set(), set(), [], Counter()
    for r in reads:
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < a.mapq:
            continue
        if r.reference_start <= pos - 1 < r.reference_end:
            spanning.add(r.query_name)
        cig = r.cigartuples or []
        refpos, hit, ilen = r.reference_start, False, None
        for op, ln in cig:
            # breakpoint of a clip is the current reference position: for a leading
            # clip that's reference_start, for an internal clip the anchor next to it
            if op == 1 and ln >= a.min_ins and abs(refpos - (pos - 1)) <= a.window:
                hit, ilen = True, ln if ilen is None else max(ilen, ln)
            elif op in (4, 5) and ln >= a.min_clip and abs(refpos - (pos - 1)) <= a.window:
                hit = True
            if op in (0, 2, 3, 7, 8):
                refpos += ln
        if hit:
            support.add(r.query_name)
            if ilen:
                ins_lens.append(ilen)
            hp[r.get_tag("HP") if r.has_tag("HP") else 0] += 1
    depth = len(spanning | support)  # reads clipped at the site don't span it but are real coverage
    return dict(depth=depth, support=len(support), ins_n=len(ins_lens),
                ins_len=int(statistics.median(ins_lens)) if ins_lens else None, hp=hp)

# ----------------------------------------------------------------- sex
def infer_sex(bam):
    """male/female/unknown from mapped-read density of chrX vs autosomes (from the BAM index)."""
    try:
        mapped = {s.contig: s.mapped for s in bam.get_index_statistics()}
    except (ValueError, AttributeError):
        return "unknown"
    lens = dict(zip(bam.references, bam.lengths))

    def density(names):
        length = sum(lens[c] for c in names)
        return sum(mapped.get(c, 0) for c in names) / length if length else None

    auto = [c for c in bam.references if norm_chr(c).isdigit() and PRIMARY.match(c)]
    xs = [c for c in bam.references if norm_chr(c) == "X"]
    ys = [c for c in bam.references if norm_chr(c) == "Y"]
    da, dx, dy = density(auto), density(xs), density(ys)
    if not da or dx is None:
        return "unknown"
    rx = dx / da
    sex = "male" if rx < 0.75 else ("female" if rx > 0.85 else "unknown")
    ry = f"{dy / da:.2f}" if dy is not None else "NA"
    sys.stderr.write(f"sex inference: chrX/autosome read density {rx:.2f} (male ~0.5, female ~1.0), "
                     f"chrY/autosome {ry} -> {sex}\n")
    return sex

# ----------------------------------------------------------------- rules
def in_par(chrom, pos):
    return any(s <= pos <= e for s, e in PAR_GRCH38.get(norm_chr(chrom), []))

def evaluate(call, ev, sv, a, bed, maxlen, minlen, subre):
    fails, flags = [], []
    fam = call["fam"]

    # effective length: BAM insertion length when enough reads carry an insertion CIGAR
    size, source = call["tlen"], "steller"
    if ev and ev["ins_len"] and ev["ins_n"] >= a.min_ins_reads and a.size_source != "steller":
        size, source = ev["ins_len"], "bam"
        ratio = call["tlen"] / ev["ins_len"]
        if ratio > 2 or ratio < 0.5:
            flags.append(f"len_discordant(steller={call['tlen']},bam={ev['ins_len']})")
    call["size"], call["size_source"] = size, source

    if fam in maxlen and size > maxlen[fam]:
        fails.append(f"oversize({size}>{maxlen[fam]},{source})")
    if fam in minlen and size < minlen[fam]:
        fails.append(f"short({size}<{minlen[fam]},{source})")
    if fam in subre and not subre[fam].search(call["sub"]):
        fails.append(f"subfamily({call['sub']})")
    if bed and in_bed(bed, call["chrom"], call["pos"]):
        fails.append("excluded_region")
    if not a.keep_alt and not PRIMARY.match(call["chrom"]):
        fails.append(f"non_primary_contig({call['chrom']})")

    af, bam_gt, gt_method = None, None, None
    if ev is not None:
        af = ev["support"] / ev["depth"] if ev["depth"] else 0.0
        if ev["support"] < a.min_support:
            fails.append(f"low_bam_support({ev['support']}<{a.min_support})")
        elif af < a.min_af:
            fails.append(f"low_af({af:.2f}<{a.min_af})")
        if ev["depth"] < a.min_depth:
            flags.append(f"low_depth({ev['depth']})")

        # HP-tag-balance genotyping: a heterozygous insertion concentrates its reads
        # on one haplotype, a homozygous one appears on both; falls back to allele
        # fraction when tags are missing, rare, or ambiguous. Het is tested first
        # (mis-tagged reads scale with depth), hom needs a substantial minor fraction.
        hp_gt = None
        if ev["support"] > 0 and a.use_hp_gt:
            h1, h2, h0 = ev["hp"].get(1, 0), ev["hp"].get(2, 0), ev["hp"].get(0, 0)
            tagged = h1 + h2
            if tagged == 0:
                pass  # untagged BAM -> AF-based genotyping only
            elif tagged / ev["support"] < a.hp_min_tagged:
                flags.append(f"mostly_untagged({h0}/{ev['support']})")
            else:
                hi, lo = max(h1, h2), min(h1, h2)
                minor = lo / (hi + lo)
                if minor <= a.hp_het_max_minor:
                    hp_gt = "0/1"
                elif lo >= a.hp_hom_min and minor >= a.hp_hom_min_minor:
                    hp_gt = "1/1"
                else:
                    flags.append(f"hp_ambiguous({h1}/{h2})")
        # HP tags can be partial or wrong; only trust them when they roughly
        # agree with the allele fraction
        if hp_gt and ((hp_gt == "0/1" and af >= a.hom_af) or
                      (hp_gt == "1/1" and af < HP_HOM_MIN_AF)):
            flags.append(f"hp_af_conflict(hp={hp_gt},af={af:.2f})")
            hp_gt = None
        if hp_gt:
            bam_gt, gt_method = hp_gt, "hp"
        else:
            bam_gt = "1/1" if af >= a.hom_af else ("0/1" if af >= a.min_af else "0/0")
            gt_method = "af"

        steller_gt = call["gt"].replace("|", "/")
        steller_hom = steller_gt in ("1/1",)
        steller_ref = steller_gt in ("0/0", ".", "./.")
        if bam_gt in ("1/1", "0/1") and (bam_gt == "1/1") != steller_hom:
            flags.append(f"gt_disagree(steller={call['gt']},bam={bam_gt})")
        elif bam_gt == "0/0" and not steller_ref:
            flags.append(f"gt_disagree(steller={call['gt']},bam=0/0)")
    call["af"], call["bam_gt"], call["gt_method"] = af, bam_gt, gt_method

    c = norm_chr(call["chrom"])
    if bam_gt is not None:
        het = bam_gt == "0/1"
    elif af is not None:
        het = af < a.hom_af
    else:
        het = len(set(call["gt"].replace("|", "/").split("/"))) > 1
    if a.sex == "female" and c == "Y":
        fails.append("chrY_in_female")
    elif a.sex == "male" and c in ("X", "Y") and not in_par(call["chrom"], call["pos"]) and het:
        fails.append("het_on_hemizygous")

    if sv is not None:
        if not call.get("sv_match"):
            flags.append("no_sv_caller_match")
            if a.require_sv_match or fam in a.require_fams:
                fails.append("no_sv_caller_match")
    return fails, flags

def find_sv(sv, call, window):
    lst = sv.get(norm_chr(call["chrom"]), [])
    i = bisect.bisect_left(lst, (call["pos"] - window, 0, ""))
    best = None
    while i < len(lst) and lst[i][0] <= call["pos"] + window:
        d = abs(lst[i][0] - call["pos"])
        if best is None or d < best[0]:
            best = (d, lst[i])
        i += 1
    return best[1] if best else None

# ----------------------------------------------------------------- output
VCF_HEADER = """##fileformat=VCFv4.2
##source=filter_te_calls.py
##ALT=<ID=INS:ME,Description="Mobile element insertion">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Insertion length used by the filter (BAM-based when available)">
##INFO=<ID=MEINFO,Number=.,Type=String,Description="sTELLeR mobile element info: family:subfamily,start,end">
##INFO=<ID=TEFAM,Number=1,Type=String,Description="TE family">
##INFO=<ID=TESUB,Number=1,Type=String,Description="TE subfamily">
##INFO=<ID=BAM_SUPPORT,Number=1,Type=Integer,Description="Reads with insertion CIGAR or soft clip at the site">
##INFO=<ID=BAM_AF,Number=1,Type=Float,Description="BAM support / local depth">
##INFO=<ID=BAM_HP,Number=3,Type=Integer,Description="Supporting reads on HP1,HP2,untagged">
##INFO=<ID=SV_MATCH,Number=0,Type=Flag,Description="Insertion in the SV-caller VCF within the match window">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype from sTELLeR">
##FORMAT=<ID=TE,Number=1,Type=String,Description="TE type as in the TE fasta">
##FORMAT=<ID=TS,Number=1,Type=Integer,Description="Transposon start position">
##FORMAT=<ID=HT,Number=1,Type=String,Description="Haplotag (0 if not phased)">
"""

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("vcf")
    p.add_argument("out_prefix")
    p.add_argument("bam", help="BAM (haplotagged or not) or 'none'")
    p.add_argument("sex", choices=["male", "female", "unknown", "auto"],
                   help="auto = infer from chrX read density in the BAM (needs a BAM index)")
    p.add_argument("--sv-vcf")
    p.add_argument("--sv-window", type=int, default=200)
    p.add_argument("--require-sv-match", action="store_true", help="fail every call without a matching SV-caller insertion")
    p.add_argument("--require-sv-match-fam", metavar="FAM,FAM", help="same, but only for these families, e.g. L1,SVA,ERV")
    p.add_argument("--exclude-bed")
    p.add_argument("--fai", help="reference .fai, adds ##contig lines to the VCF")
    p.add_argument("--max-len", action="append", metavar="FAM=N", help="default Alu=700 L1=7000 SVA=4500 ERV=10000")
    p.add_argument("--min-len", action="append", metavar="FAM=N", help="default Alu=200 L1=200 SVA=300 ERV=500")
    p.add_argument("--subfamily", action="append", metavar="FAM=REGEX", help="default Alu=^AluY")
    p.add_argument("--size-source", choices=["auto", "steller"], default="auto")
    p.add_argument("--window", type=int, default=150)
    p.add_argument("--min-ins", type=int, default=50, help="min CIGAR insertion length counted as support")
    p.add_argument("--min-clip", type=int, default=100, help="min soft/hard clip length counted as support")
    p.add_argument("--min-ins-reads", type=int, default=2, help="reads needed to trust the BAM insertion length")
    p.add_argument("--mapq", type=int, default=20)
    p.add_argument("--min-support", type=int, default=3)
    p.add_argument("--min-af", type=float, default=0.15)
    p.add_argument("--hom-af", type=float, default=0.75)
    p.add_argument("--min-depth", type=int, default=8)
    p.add_argument("--merge-dist", type=int, default=1000)
    p.add_argument("--keep-alt", action="store_true", help="do not fail calls on alt/decoy/unplaced contigs")
    p.add_argument("--gt-source", choices=["steller", "bam"], default="steller",
                   help="genotype written to the VCF; 'bam' = from HP balance or allele fraction")
    p.add_argument("--no-hp-gt", action="store_true",
                   help="disable HP-tag-balance genotyping; genotype from allele fraction only")
    p.add_argument("--hp-min-tagged", type=float, default=0.5,
                   help="min fraction of HP-tagged supporting reads to use HP genotyping")
    p.add_argument("--hp-hom-min", type=int, default=3,
                   help="min supporting reads on the minor haplotype to call 1/1 from HP tags")
    p.add_argument("--hp-het-max-minor", type=float, default=0.25,
                   help="max minor-haplotype fraction among tagged reads to call 0/1 from HP tags")
    p.add_argument("--hp-hom-min-minor", type=float, default=0.35,
                   help="min minor-haplotype fraction among tagged reads to call 1/1 from HP tags")
    a = p.parse_args()
    a.use_hp_gt = not a.no_hp_gt

    maxlen = kv_list(a.max_len, DEFAULT_MAX, int)
    minlen = kv_list(a.min_len, DEFAULT_MIN, int)
    subre = {k: re.compile(v) for k, v in kv_list(a.subfamily, DEFAULT_SUB, str).items()}
    a.require_fams = {x for x in (a.require_sv_match_fam or "").split(",") if x}

    calls, sample, n_skipped = read_steller(a.vcf)
    calls.sort(key=lambda c: (chr_key(c["chrom"]), c["pos"]))
    bed = read_bed(a.exclude_bed) if a.exclude_bed else None
    sv = read_sv_ins(a.sv_vcf) if a.sv_vcf else None
    bam = None
    if a.bam.lower() != "none":
        import pysam
        bam = pysam.AlignmentFile(a.bam, "rb")
        # an unindexed BAM raises the same ValueError as a missing contig, which
        # would silently disable every BAM filter -> fail loudly instead
        if not bam.has_index():
            sys.exit(f"error: no index for {a.bam} - run: samtools index {a.bam}")

    if a.sex == "auto":
        a.sex = infer_sex(bam) if bam else "unknown"
    if a.require_fams and sv is None:
        sys.stderr.write("warning: --require-sv-match-fam needs --sv-vcf; ignoring it\n")

    regions = plan_bam_fetches(calls, a.window) if bam else []
    last_region = None
    for c in calls:
        ev = None
        if bam:
            if last_region is not None and c["region"] != last_region:
                # done with that region: drop its reads ([] also prevents a refetch)
                regions[last_region][3] = []
            last_region = c["region"]
            reads = region_reads(bam, regions, c["region"])
            if reads is None:
                c["no_contig"] = True
            else:
                ev = bam_evidence(c, a, reads)
        c["ev"] = ev
        if sv is not None:
            hit = find_sv(sv, c, a.sv_window)
            c["sv_match"], c["sv"] = bool(hit), hit
        c["fails"], c["flags"] = evaluate(c, ev, sv, a, bed, maxlen, minlen, subre)
        if c.get("no_contig"):
            c["flags"].append("contig_not_in_bam")

    # duplicates: same family, calls within merge-dist of the cluster anchor -> keep
    # the best call of each cluster (anchored, so distant calls stay separate)
    groups = defaultdict(list)
    for c in calls:
        groups[(c["chrom"], c["fam"])].append(c)
    for lst in groups.values():
        lst.sort(key=lambda x: x["pos"])
        cluster = [lst[0]]
        for c in lst[1:] + [None]:
            if c is not None and c["pos"] - cluster[0]["pos"] <= a.merge_dist:
                cluster.append(c)
                continue
            if len(cluster) > 1:
                ok = [x for x in cluster if not x["fails"]] or cluster
                best = max(ok, key=lambda x: (x["ev"]["support"] if x["ev"] else 0, -x["pos"]))
                for x in cluster:
                    if x is not best:
                        x["fails"].append(f"duplicate_of({best['chrom']}:{best['pos']})")
            cluster = [c] if c is not None else []

    n_nocontig = sum(1 for c in calls if c.get("no_contig"))
    if n_nocontig:
        sys.stderr.write(f"warning: {n_nocontig}/{len(calls)} calls are on contigs missing from the BAM header; "
                         "BAM filters were skipped for them\n")

    # outputs
    kept = [c for c in calls if not c["fails"]]
    with open(f"{a.out_prefix}.annotated.tsv", "w") as o:
        o.write("\t".join(["chrom", "pos", "family", "subfamily", "steller_gt", "steller_len",
                           "bam_depth", "bam_support", "bam_af", "bam_ins_reads", "bam_ins_len",
                           "size_used", "size_source", "hp_support(HP1/HP2/none)", "bam_gt",
                           "bam_gt_method", "sv_match", "sv_pos", "sv_len",
                           "result", "fail_reasons", "flags"]) + "\n")
        for c in calls:
            ev = c["ev"]
            hp = ev["hp"] if ev else {}
            sv_hit = c.get("sv")
            row = [c["chrom"], c["pos"], c["fam"], c["sub"], c["gt"], c["tlen"],
                   ev["depth"] if ev else "NA", ev["support"] if ev else "NA",
                   f"{c['af']:.2f}" if c["af"] is not None else "NA",
                   ev["ins_n"] if ev else "NA", (ev["ins_len"] if ev and ev["ins_len"] else "NA"),
                   c["size"], c["size_source"],
                   f"{hp.get(1, 0)}/{hp.get(2, 0)}/{hp.get(0, 0)}" if ev else "NA",
                   c["bam_gt"] or "NA", c["gt_method"] or "NA",
                   ("yes" if sv_hit else "no") if sv is not None else "NA",
                   sv_hit[0] if sv_hit else "NA", sv_hit[1] if sv_hit else "NA",
                   "PASS" if not c["fails"] else "FAIL",
                   ";".join(c["fails"]) or "-", ";".join(c["flags"]) or "-"]
            o.write("\t".join(str(x) for x in row) + "\n")

    contigs = ""
    if a.fai:
        with open(a.fai) as fh:
            contigs = "".join(f"##contig=<ID={l.split()[0]},length={l.split()[1]}>\n" for l in fh)
    elif bam is not None:  # no .fai needed: take the contigs used by the output from the BAM header
        used = {c["chrom"] for c in kept}
        contigs = "".join(f"##contig=<ID={n},length={ln}>\n"
                          for n, ln in zip(bam.references, bam.lengths) if n in used)
    with open(f"{a.out_prefix}.high_confidence.vcf", "w") as o:
        o.write(VCF_HEADER + contigs)
        o.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample]) + "\n")
        for c in kept:
            info = [f"SVTYPE=INS", f"SVLEN={c['size']}", f"MEINFO={c['meinfo']}",
                    f"TEFAM={c['fam']}", f"TESUB={c['sub']}"]
            if c["ev"]:
                hp = c["ev"]["hp"]
                info += [f"BAM_SUPPORT={c['ev']['support']}", f"BAM_AF={c['af']:.2f}",
                         f"BAM_HP={hp.get(1, 0)},{hp.get(2, 0)},{hp.get(0, 0)}"]
            if c.get("sv_match"):
                info.append("SV_MATCH")
            sample_out = c["sample_field"]
            if a.gt_source == "bam" and c["bam_gt"] in ("1/1", "0/1"):
                sd = dict(c["sd"], GT=c["bam_gt"])
                sample_out = ":".join(sd[k] for k in c["order"])
            o.write("\t".join([c["chrom"], str(c["pos"]), ".", c["ref"], c["alt"], ".", "PASS",
                               ";".join(info), c["fmt"], sample_out]) + "\n")

    reasons = Counter(r.split("(")[0] for c in calls for r in c["fails"])
    summary = f"{len(calls)} calls, {len(kept)} high-confidence"
    if n_skipped:
        summary += f", {n_skipped} malformed VCF lines skipped"
    sys.stderr.write(summary + "\n")
    for k, v in reasons.most_common():
        sys.stderr.write(f"  fail {k}: {v}\n")

if __name__ == "__main__":
    main()