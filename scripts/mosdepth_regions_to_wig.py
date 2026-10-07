#!/usr/bin/env python3
"""Convert mosdepth <prefix>.regions.bed.gz (fixed windows) to a fixedStep wig.

usage: mosdepth_regions_to_wig.py <regions.bed.gz> <out.wig> <step> <strip_chr 0|1>

Keeps chr1-22, X, Y only. strip_chr=1 writes 1,2,..,X instead of chr1,..
(check which naming your chromograph version expects if plots come out empty).
"""
import gzip
import sys

KEEP = {str(i) for i in range(1, 23)} | {"X", "Y"}


def main(bed, out, step, strip):
    step = int(step)
    strip = strip == "1"
    cur = None
    with gzip.open(bed, "rt") as fh, open(out, "w") as o:
        for line in fh:
            c, s, e, v = line.split("\t")[:4]
            base = c[3:] if c.startswith("chr") else c
            if base not in KEEP:
                continue
            if c != cur:
                name = base if strip else c
                o.write(f"fixedStep chrom={name} start=1 step={step} span={step}\n")
                cur = c
            o.write(f"{float(v):.2f}\n")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])