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

def load_recomb_map(tsv_path, pos_col, cm_col, sep='\t', comment='#', no_header=False):
    bp, cm = [], []
    with auto_open(tsv_path, 'rt') as f:
        if not no_header:
            header = f.readline().rstrip('\n').split(sep)
            # case-insensitive columns
            hdr = [h.strip() for h in header]
            low = {h.lower(): i for i, h in enumerate(hdr)}
            if pos_col.lower() not in low or cm_col.lower() not in low:
                raise ValueError(f"Recomb map missing columns. Have: {header}. Need: {pos_col}, {cm_col}")
            ipos, icm = low[pos_col.lower()], low[cm_col.lower()]
        else:
            # Use column indices directly (0-indexed)
            ipos, icm = int(pos_col), int(cm_col)
            
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
    sys.stderr.write(f"Loaded recomb map: {len(bp)} positions, range {bp[0]}-{bp[-1]} bp\n")
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
    ap.add_argument("--chr", required=True, help="Chromosome name (must match FASTA, e.g. 'chr21' or '21')")
    ap.add_argument("--recomb-tsv", required=True)
    ap.add_argument("--pos-col", default="position", help="Recomb map bp column (default: position) or 0-indexed column number if --no-header")
    ap.add_argument("--cm-col",  default="cM",       help="Recomb map cM column (default: cM) or 0-indexed column number if --no-header")
    ap.add_argument("--sep",     default="\t",       help="Recomb map separator (default: tab)")
    ap.add_argument("--no-header", action="store_true", help="Recomb map has no header; use --pos-col and --cm-col as 0-indexed column numbers")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    # pysam import (fail nicely if missing)
    try:
        import pysam
    except ImportError:
        sys.stderr.write("ERROR: pysam not installed. Use 'pixi add pysam' or 'python -m pip install --user pysam'.\n")
        sys.exit(1)

    bp_knots, cm_knots = load_recomb_map(args.recomb_tsv, args.pos_col, args.cm_col, sep=args.sep, no_header=args.no_header)
    
    sys.stderr.write(f"Opening ancestral FASTA: {args.anc_fasta}\n")
    fa = pysam.FastaFile(args.anc_fasta)
    
    # Check what chromosomes are available
    sys.stderr.write(f"Available chromosomes in FASTA: {fa.references[:10]}...\n")
    if args.chr not in fa.references:
        sys.stderr.write(f"WARNING: Chromosome '{args.chr}' not found in FASTA!\n")
        sys.stderr.write(f"Available: {', '.join(fa.references)}\n")
        sys.stderr.write("Trying to continue anyway...\n")

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
    skip_reasons = {"no_anc": 0, "anc_mismatch": 0}
    line_no = 0
    first_skip_examples = []

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

        # Fetch ancestral allele
        try:
            anc = fa.fetch(args.chr, pos-1, pos).upper()
        except Exception as e:
            if len(first_skip_examples) < 5:
                first_skip_examples.append(f"pos {pos}: fetch error - {e}")
            skip_reasons["no_anc"] += 1
            skipped += 1
            continue

        if anc not in ("A","C","G","T"):
            if len(first_skip_examples) < 5:
                first_skip_examples.append(f"pos {pos}: anc='{anc}' (not ACGT)")
            skip_reasons["no_anc"] += 1
            skipped += 1
            continue

        if anc == a0:
            polarized = tokens
        elif anc == a1:
            polarized = flip_01_tokens(tokens)
        else:
            if len(first_skip_examples) < 5:
                first_skip_examples.append(f"pos {pos}: a0={a0}, a1={a1}, anc={anc} (mismatch)")
            skip_reasons["anc_mismatch"] += 1
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
    
    sys.stderr.write(f"\n{'='*60}\n")
    sys.stderr.write(f"SUMMARY:\n")
    sys.stderr.write(f"  Total processed: {line_no}\n")
    sys.stderr.write(f"  Kept: {kept}\n")
    sys.stderr.write(f"  Skipped: {skipped}\n")
    sys.stderr.write(f"    - No/invalid ancestral: {skip_reasons['no_anc']}\n")
    sys.stderr.write(f"    - Ancestral mismatch: {skip_reasons['anc_mismatch']}\n")
    
    if first_skip_examples:
        sys.stderr.write(f"\nFirst few skip examples:\n")
        for ex in first_skip_examples:
            sys.stderr.write(f"  {ex}\n")
    
    if kept == 0:
        sys.stderr.write(f"\n⚠️  WARNING: ZERO SNPs were kept!\n")
        sys.stderr.write(f"Common issues:\n")
        sys.stderr.write(f"  1. Wrong --chr name (check if FASTA uses 'chr21' vs '21')\n")
        sys.stderr.write(f"  2. Ancestral FASTA doesn't cover these positions\n")
        sys.stderr.write(f"  3. All ancestral alleles differ from a0/a1 in legend\n")
    sys.stderr.write(f"{'='*60}\n")

if __name__ == "__main__":
    main()