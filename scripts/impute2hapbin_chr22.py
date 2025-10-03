#!/usr/bin/env python3
"""
IMPUTE hap/legend -> hapbin .hap/.map
- Polarizes alleles using ancestral FASTA (ancestral -> 0)
- Interpolates cM from a TSV (position, cM)
Legend header can be 'ID position a0 a1' or lowercase; matching is case-insensitive.
"""
import sys, gzip, argparse, bisect

def auto_open(path, mode='rt'):
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)

def load_recomb_map(tsv_path, pos_col, cm_col, sep='\t', comment='#'):
    bp, cm = [], []
    with auto_open(tsv_path, 'rt') as f:
        header = f.readline().rstrip('\n').split(sep)
        # case-insensitive columns
        hdr = [h.strip() for h in header]
        low = {h.lower(): i for i, h in enumerate(hdr)}
        if pos_col.lower() not in low or cm_col.lower() not in low:
            raise ValueError(f"Recomb map missing columns. Have: {header}. Need: {pos_col}, {cm_col}")
        ipos, icm = low[pos_col.lower()], low[cm_col.lower()]
        for line in f:
            if not line.strip() or line.startswith(comment): 
                continue
            parts = line.rstrip('\n').split(sep)
            try:
                bp_i = int(parts[ipos]); cm_i = float(parts[icm])
            except Exception:
                continue
            bp.append(bp_i); cm.append(cm_i)
    if not bp:
        raise ValueError("Recomb map loaded 0 rows")
    return bp, cm

def interp_cm(bp_knots, cm_knots, x):
    i = bisect.bisect_left(bp_knots, x)
    if i == 0: return cm_knots[0]
    if i >= len(bp_knots): return cm_knots[-1]
    x0, x1 = bp_knots[i-1], bp_knots[i]
    y0, y1 = cm_knots[i-1], cm_knots[i]
    if x1 == x0: return y0
    return y0 + (y1 - y0) * ((x - x0) / (x1 - x0))

def flip_01_tokens(tokens):
    return ['1' if t == '0' else ('0' if t == '1' else t) for t in tokens]

def main():
    ap = argparse.ArgumentParser(description="IMPUTE->hapbin with ancestral polarization + cM interpolation")
    ap.add_argument("--hap", required=True)
    ap.add_argument("--legend", required=True)
    ap.add_argument("--anc-fasta", required=True)
    ap.add_argument("--chr", required=True)
    ap.add_argument("--recomb-tsv", required=True)
    ap.add_argument("--pos-col", default="position", help="Recomb map bp column (default: position)")
    ap.add_argument("--cm-col",  default="cM",       help="Recomb map cM column (default: cM)")
    ap.add_argument("--sep",     default="\t",       help="Recomb map separator (default: tab)")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    # pysam import (fail nicely if missing)
    try:
        import pysam
    except ImportError:
        sys.stderr.write("ERROR: pysam not installed. Use 'pixi add pysam' or 'python -m pip install --user pysam'.\n")
        sys.exit(1)

    bp_knots, cm_knots = load_recomb_map(args.recomb_tsv, args.pos_col, args.cm_col, sep=args.sep)
    fa = pysam.FastaFile(args.anc_fasta)

    hap_in = auto_open(args.hap, 'rt')
    leg_in = auto_open(args.legend, 'rt')
    hap_out = open(args.out_prefix + ".hap", "wt")
    map_out = open(args.out_prefix + ".map", "wt")

    # Legend: accept case/alias variants, typically "id position a0 a1"
    header = leg_in.readline().rstrip('\n').split()
    low = {h.lower(): i for i, h in enumerate(header)}
    required = ("id","position","a0","a1")
    missing = [r for r in required if r not in low]
    if missing:
        sys.stderr.write(f"Legend header looks unexpected. Got: {header}\n"
                         f"Need fields (case-insensitive): {required}\n")
        sys.exit(1)
    i_id, i_pos, i_a0, i_a1 = low["id"], low["position"], low["a0"], low["a1"]

    kept = skipped = 0
    line_no = 0

    while True:
        leg_line = leg_in.readline()
        if not leg_line:
            break
        leg_line = leg_line.rstrip("\n")
        if not leg_line:
            continue
        parts = leg_line.split()
        pos = int(parts[i_pos])
        a0 = parts[i_a0].upper()
        a1 = parts[i_a1].upper()
        locus_id = parts[i_id]

        hap_line = hap_in.readline()
        if not hap_line:
            sys.stderr.write("ERROR: hap and legend out of sync (legend has more rows).\n")
            break
        tokens = hap_line.rstrip("\n").split()

        anc = fa.fetch(args.chr, pos-1, pos).upper()
        if anc not in ("A","C","G","T"):
            skipped += 1
            continue

        if anc == a0:
            polarized = tokens
        elif anc == a1:
            polarized = flip_01_tokens(tokens)
        else:
            skipped += 1
            continue

        cm_val = interp_cm(bp_knots, cm_knots, pos)
        map_out.write(f"{args.chr} {locus_id} {cm_val:.6f} {pos}\n")
        hap_out.write(" ".join(polarized) + "\n")
        kept += 1
        line_no += 1
        if args.verbose and line_no % 100000 == 0:
            sys.stderr.write(f"...processed {line_no} SNPs (kept {kept}, skipped {skipped})\n")

    hap_in.close(); leg_in.close(); hap_out.close(); map_out.close()
    sys.stderr.write(f"Done. Kept {kept} sites; skipped {skipped}.\n")

if __name__ == "__main__":
    main()
