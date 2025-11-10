#!/usr/bin/env python3
"""
Interpolate genetic map to have an entry for every SNP position
"""

import sys

def linear_interp(x, xp, yp):
    """Fast linear interpolation using binary search"""
    import bisect
    result = []
    for xi in x:
        # Find the two points to interpolate between
        if xi <= xp[0]:
            result.append(yp[0])
        elif xi >= xp[-1]:
            result.append(yp[-1])
        else:
            # Use binary search to find insertion point
            i = bisect.bisect_left(xp, xi)
            if i > 0:
                i -= 1
            # Linear interpolation
            t = (xi - xp[i]) / (xp[i+1] - xp[i])
            yi = yp[i] + t * (yp[i+1] - yp[i])
            result.append(yi)
    return result

# Read SNP positions from HAPS file
print("Reading SNP positions from HAPS file...", file=sys.stderr)
snp_positions = []
with open("/home/vanbruggenmit/mit-ihh-pib/data/relate/chrX_nonPAR/chrX_nonPAR_biallelic.haps") as f:
    for line in f:
        fields = line.strip().split()
        snp_positions.append(int(fields[2]))  # Position is 3rd column

print(f"Read {len(snp_positions)} SNP positions", file=sys.stderr)

# Read genetic map
print("Reading genetic map...", file=sys.stderr)
map_positions = []
map_cm = []
with open("/home/vanbruggenmit/mit-ihh-pib/data/relate/chrX_nonPAR/chrX_nonPAR.relate.map") as f:
    for line in f:
        fields = line.strip().split()
        map_positions.append(int(fields[1]))  # Position is 2nd column
        map_cm.append(float(fields[2]))  # cM is 3rd column

print(f"Read {len(map_positions)} map positions", file=sys.stderr)

# Interpolate genetic map for all SNP positions
print("Interpolating genetic map...", file=sys.stderr)
interpolated_cm = linear_interp(snp_positions, map_positions, map_cm)

# Write interpolated map
print("Writing interpolated map...", file=sys.stderr)
with open("/home/vanbruggenmit/mit-ihh-pib/data/relate/chrX_nonPAR/chrX_nonPAR_biallelic_interpolated.relate.map", "w") as f:
    for pos, cm in zip(snp_positions, interpolated_cm):
        f.write(f"chrX:{pos} {pos} {cm:.10g}\n")

print(f"Done! Wrote {len(snp_positions)} positions to interpolated map", file=sys.stderr)
