#!/usr/bin/env python3
"""
Create a .map file for hapbin from an IMPUTE legend file and recombination map.
Map format: chromosome locus_id genetic_position physical_position
"""

import gzip
import argparse
import bisect
import sys

def auto_open(path, mode='rt'):
    """Open gzipped or regular files automatically."""
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)

def load_recomb_map(tsv_path):
    """
    Load recombination map with format: chromosome position cM
    No header expected.
    Returns: (bp_positions, cM_values)
    """
    bp, cm = [], []
    with auto_open(tsv_path, 'rt') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                # Format: chromosome position cM
                bp_val = int(parts[1])
                cm_val = float(parts[2])
                bp.append(bp_val)
                cm.append(cm_val)
            except (ValueError, IndexError):
                continue

    if not bp:
        raise ValueError("Recombination map is empty or could not be parsed")

    print(f"Loaded {len(bp)} recombination map positions", file=sys.stderr)
    print(f"Range: {bp[0]:,} - {bp[-1]:,} bp", file=sys.stderr)
    print(f"cM range: {cm[0]:.6f} - {cm[-1]:.6f}", file=sys.stderr)
    return bp, cm

def interpolate_cm(bp_knots, cm_knots, position):
    """Linearly interpolate cM value for a given position."""
    i = bisect.bisect_left(bp_knots, position)

    # Before first knot
    if i == 0:
        return cm_knots[0]

    # After last knot
    if i >= len(bp_knots):
        return cm_knots[-1]

    # Linear interpolation between i-1 and i
    x0, x1 = bp_knots[i-1], bp_knots[i]
    y0, y1 = cm_knots[i-1], cm_knots[i]

    if x1 == x0:
        return y0

    return y0 + (y1 - y0) * ((position - x0) / (x1 - x0))

def main():
    parser = argparse.ArgumentParser(
        description="Create .map file for hapbin from IMPUTE legend and recombination map"
    )
    parser.add_argument("--legend", required=True, help="IMPUTE legend file (.legend.gz)")
    parser.add_argument("--recomb", required=True, help="Recombination map TSV (chr pos cM)")
    parser.add_argument("--chr", required=True, help="Chromosome name (e.g., 'X' or 'chrX')")
    parser.add_argument("--out", required=True, help="Output .map file")

    args = parser.parse_args()

    print(f"Creating map file for chromosome {args.chr}", file=sys.stderr)
    print(f"Legend file: {args.legend}", file=sys.stderr)
    print(f"Recomb map: {args.recomb}", file=sys.stderr)
    print(f"Output: {args.out}", file=sys.stderr)
    print("", file=sys.stderr)

    # Load recombination map
    bp_knots, cm_knots = load_recomb_map(args.recomb)

    # Process legend file
    variant_count = 0
    with auto_open(args.legend, 'rt') as legend_in, open(args.out, 'w') as map_out:
        # Read header
        header = legend_in.readline().strip().split()
        print(f"Legend header: {' '.join(header)}", file=sys.stderr)

        # Find column indices (case-insensitive)
        header_lower = {h.lower(): i for i, h in enumerate(header)}

        if 'id' not in header_lower or 'position' not in header_lower:
            print(f"ERROR: Legend file must have 'id' and 'position' columns", file=sys.stderr)
            print(f"Found columns: {header}", file=sys.stderr)
            sys.exit(1)

        id_col = header_lower['id']
        pos_col = header_lower['position']

        # Process each variant
        for line_num, line in enumerate(legend_in, start=2):
            line = line.strip()
            if not line:
                continue

            parts = line.split()
            locus_id = parts[id_col]
            position = int(parts[pos_col])

            # Interpolate genetic position
            cm_val = interpolate_cm(bp_knots, cm_knots, position)

            # Write map line: chromosome locus_id genetic_position physical_position
            map_out.write(f"{args.chr} {locus_id} {cm_val:.6f} {position}\n")

            variant_count += 1

            if variant_count % 500000 == 0:
                print(f"Processed {variant_count:,} variants...", file=sys.stderr)

    print("", file=sys.stderr)
    print(f"Complete! Wrote {variant_count:,} variants to {args.out}", file=sys.stderr)

if __name__ == "__main__":
    main()
