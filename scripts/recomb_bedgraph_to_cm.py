#!/usr/bin/env python3
"""
Convert recombination rate bedGraph (cM/Mb) to cumulative cM "knots" TSV.

Default output columns: position <tab> cM
Optional: include a chromosome column and strip 'chr' prefix.

Assumptions:
- bedGraph columns: chrom  start  end  rate_cM_per_Mb
- Positions are 0-based half-open [start, end) like standard bedGraph.
- Cumulative cM is reset when chromosome changes (if multi-chr input).
"""
import sys
import argparse

def open_in(path):
    return sys.stdin if path == "-" else open(path, "rt")

def open_out(path):
    return sys.stdout if path == "-" else open(path, "wt")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", default="-", help="Input bedGraph file (default: stdin)")
    ap.add_argument("--out", dest="out", default="-", help="Output TSV file (default: stdout)")
    ap.add_argument("--include-chrom", action="store_true",
                    help="Include chromosome column as the first column")
    ap.add_argument("--strip-chr", action="store_true",
                    help="Strip leading 'chr' from chromosome names in output")
    args = ap.parse_args()

    cm = 0.0
    last_chrom = None

    with open_in(args.inp) as f, open_out(args.out) as g:
        # Header
        if args.include_chrom:
            g.write("chr\tposition\tcM\n")
        else:
            g.write("position\tcM\n")

        for line in f:
            if not line.strip() or line.startswith(("track", "browser", "#")):
                continue
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            chrom, s, e, rate = parts[0], int(parts[1]), int(parts[2]), float(parts[3])

            # Reset cumulative when chromosome changes (handles multi-chr inputs)
            if chrom != last_chrom:
                cm = 0.0
                last_chrom = chrom

            # Integrate rate over [s, e): delta_cM = rate(cM/Mb) * bp / 1e6
            cm += rate * (e - s) / 1e6

            pos_end_1based = e  # 1-based knot at interval end
            out_chrom = chrom[3:] if (args.strip_chr and chrom.startswith("chr")) else chrom

            if args.include_chrom:
                g.write(f"{out_chrom}\t{pos_end_1based}\t{cm:.6f}\n")
            else:
                g.write(f"{pos_end_1based}\t{cm:.6f}\n")

if __name__ == "__main__":
    main()
