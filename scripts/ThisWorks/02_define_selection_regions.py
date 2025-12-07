#!/usr/bin/env python3
"""
Script: Define selection regions from iHS candidates
Purpose: Cluster nearby significant SNPs into selection regions

Approach (based on literature):
- Villegas-Miron 2021: Multiple contiguous SNPs needed, not single outliers
- Approach: 20kb windows requiring ≥20 significant SNPs
- We'll use 20kb windows with ≥20 SNPs for initial scan, more stringent with ≥50 
- Then merge overlapping windows

Key considerations for X chromosome:
- Use empirical 99th percentile (as calculated in 01 script) instead of fixed 2.5 threshold
- Account for lower candidate density
- May require adjusted window parameters

Author: MIT IHH PIB Project
Date: 2025-10-23
Updated: 2025-12-07
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Setup paths
base_dir = Path("/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib")
candidates_dir = base_dir / "results/ihs_EAS/analysis/candidates/candidates"
regions_dir = base_dir / "results/ihs_EAS/analysis/regions"
regions_dir.mkdir(parents=True, exist_ok=True)

print("=" * 70)
print("iHS Selection Region Definition")
print("=" * 70)
print()

# Configuration
WINDOW_SIZE = 20000  # windows
MIN_SNPS = 20          # Minimum SNPs per region 
MIN_SNPS_STRINGENT = 50  # More stringent filter

# Threshold configuration
# Use X-specific empirical threshold 
AUTOSOME_THRESHOLD = 2.5
X_THRESHOLD = 3.131  # X chromosome 99th percentile for that population

print(f"Configuration:")
print(f"  Window size: {WINDOW_SIZE:,} bp ({WINDOW_SIZE/1000:.0f} kb)")
print(f"  Minimum SNPs per region: {MIN_SNPS}")
print(f"  Autosome threshold: |Std iHS| >= {AUTOSOME_THRESHOLD}")
print(f"  X chromosome threshold: |Std iHS| >= {X_THRESHOLD}")
print()

def load_candidates(threshold_level="moderate"):
    """Load candidate variants for specified threshold"""
    file = candidates_dir / f"ALL_candidates_{threshold_level}.tsv"
    print(f"Loading candidates from: {file.name}")
    df = pd.read_csv(file, sep='\t')
    print(f"  Total candidates: {len(df):,}")
    return df

def apply_chromosome_specific_threshold(df, autosome_thresh, x_thresh):
    """Apply different thresholds for autosomes vs X chromosome"""
    autosome_mask = (df['chr_type'] == 'Autosome') & (np.abs(df['Std iHS']) >= autosome_thresh)
    x_mask = (df['chr_type'] == 'X chromosome') & (np.abs(df['Std iHS']) >= x_thresh)

    filtered = df[autosome_mask | x_mask].copy()

    print(f"  After chromosome-specific thresholds:")
    print(f"    Autosomes (|Std iHS| >= {autosome_thresh}): {autosome_mask.sum():,}")
    print(f"    X chromosome (|Std iHS| >= {x_thresh}): {x_mask.sum():,}")
    print(f"    Total retained: {len(filtered):,}")

    return filtered

def cluster_candidates_to_regions(candidates, window_size, min_snps):
    """
    Cluster candidates into selection regions using sliding windows

    Algorithm:
    1. Sort candidates by chromosome and position
    2. For each chromosome, scan with sliding window
    3. Identify windows with >= min_snps candidates
    4. Merge overlapping windows
    5. Calculate region statistics
    """
    regions = []

    # Process each chromosome separately
    # Sort: autosomes numerically, then X
    def sort_chr(x):
        if x == 'X':
            return (1, 0)  # X comes last
        else:
            # Handle both "21" and "chr21" formats, and integer types
            chr_str = str(x).replace('chr', '')
            return (0, int(chr_str))  # Autosomes sorted numerically

    for chr_name in sorted(candidates['chr'].unique(), key=sort_chr):
        chr_data = candidates[candidates['chr'] == chr_name].sort_values('pos').copy()
        chr_type = chr_data['chr_type'].iloc[0]

        if len(chr_data) < min_snps:
            print(f"    chr{chr_name}: {len(chr_data)} candidates (< {min_snps}, skipping)")
            continue

        print(f"    chr{chr_name}: {len(chr_data)} candidates", end=" ")

        # Sort by position
        chr_data = chr_data.sort_values('pos')
        positions = chr_data['pos'].values

        # Identify clusters using window approach
        clusters = []
        i = 0
        while i < len(positions):
            window_start = positions[i]
            window_end = window_start + window_size

            # Find all SNPs in this window
            in_window = (positions >= window_start) & (positions < window_end)
            window_indices = np.where(in_window)[0]

            if len(window_indices) >= min_snps:
                # Create a region from this window
                region_snps = chr_data.iloc[window_indices]
                region_start = region_snps['pos'].min()
                region_end = region_snps['pos'].max()

                # Check if this overlaps with previous region
                if clusters and clusters[-1]['end'] >= region_start:
                    # Merge with previous region
                    clusters[-1]['indices'].extend(window_indices.tolist())
                    clusters[-1]['indices'] = sorted(list(set(clusters[-1]['indices'])))
                    clusters[-1]['end'] = max(clusters[-1]['end'], region_end)
                else:
                    # New region
                    clusters.append({
                        'start': region_start,
                        'end': region_end,
                        'indices': window_indices.tolist()
                    })

                # Move to next window starting from last SNP in current window
                i = window_indices[-1] + 1
            else:
                i += 1

        # Calculate statistics for each region
        for cluster in clusters:
            region_snps = chr_data.iloc[cluster['indices']]

            # Calculate region statistics
            n_snps = len(region_snps)
            region_length = cluster['end'] - cluster['start']
            mean_std_ihs = region_snps['Std iHS'].mean()
            max_std_ihs = region_snps['Std iHS'].abs().max()
            mean_pvalue = region_snps['p_value'].mean()
            min_pvalue = region_snps['p_value'].min()

            # Get top SNP
            top_snp_idx = region_snps['Std iHS'].abs().idxmax()
            top_snp = region_snps.loc[top_snp_idx]

            regions.append({
                'chr': chr_name,
                'start': cluster['start'],
                'end': cluster['end'],
                'length': region_length,
                'n_snps': n_snps,
                'mean_std_ihs': mean_std_ihs,
                'max_std_ihs': max_std_ihs,
                'mean_pvalue': mean_pvalue,
                'min_pvalue': min_pvalue,
                'top_snp_location': top_snp['Location'],
                'top_snp_std_ihs': top_snp['Std iHS'],
                'chr_type': chr_type
            })

        print(f"-> {len(clusters)} regions")

    return pd.DataFrame(regions)

# Step 1: Load moderate threshold candidates
print("Step 1: Loading candidates...")
candidates = load_candidates("moderate")
print()

# Step 2: Apply chromosome-specific thresholds
print("Step 2: Applying chromosome-specific thresholds...")
filtered_candidates = apply_chromosome_specific_threshold(
    candidates,
    AUTOSOME_THRESHOLD,
    X_THRESHOLD
)
print()

# Step 3: Cluster into regions
print(f"Step 3: Clustering candidates into regions (min {MIN_SNPS} SNPs)...")
regions = cluster_candidates_to_regions(filtered_candidates, WINDOW_SIZE, MIN_SNPS)
print()

# Step 4: Summary statistics
print("=" * 70)
print("Region Summary")
print("=" * 70)
print()

print(f"Total regions identified: {len(regions)}")
print()

# By chromosome type
for chr_type in ['Autosome', 'X chromosome']:
    type_regions = regions[regions['chr_type'] == chr_type]
    if len(type_regions) == 0:
        continue

    print(f"{chr_type}:")
    print(f"  Regions: {len(type_regions)}")
    print(f"  Mean SNPs per region: {type_regions['n_snps'].mean():.1f}")
    print(f"  Median SNPs per region: {type_regions['n_snps'].median():.0f}")
    print(f"  Mean region length: {type_regions['length'].mean()/1000:.1f} kb")
    print(f"  Mean max |Std iHS|: {type_regions['max_std_ihs'].mean():.2f}")
    print()

# Top regions
print("Top 10 regions by max |Std iHS|:")
top_regions = regions.nlargest(10, 'max_std_ihs')
for idx, row in top_regions.iterrows():
    print(f"  {row['chr']}:{row['start']:,}-{row['end']:,} "
          f"({row['length']/1000:.0f} kb, {row['n_snps']} SNPs, "
          f"max |Std iHS|={row['max_std_ihs']:.2f})")
print()

# X chromosome specific
x_regions = regions[regions['chr_type'] == 'X chromosome']
if len(x_regions) > 0:
    print(f"X chromosome regions: {len(x_regions)}")
    print("Top 5 X chromosome regions:")
    for idx, row in x_regions.nlargest(5, 'max_std_ihs').iterrows():
        print(f"  {row['chr']}:{row['start']:,}-{row['end']:,} "
              f"({row['length']/1000:.0f} kb, {row['n_snps']} SNPs, "
              f"max |Std iHS|={row['max_std_ihs']:.2f})")
    print()

# Step 5: Save results
print("Step 5: Saving results...")

# Save all regions
output_file = regions_dir / "selection_regions_250kb_min2snps.tsv"
regions.to_csv(output_file, sep='\t', index=False)
print(f"  Saved: {output_file.name}")

# Save regions by chromosome type
autosome_regions = regions[regions['chr_type'] == 'Autosome']
autosome_file = regions_dir / "selection_regions_autosomes.tsv"
autosome_regions.to_csv(autosome_file, sep='\t', index=False)
print(f"  Saved: {autosome_file.name} ({len(autosome_regions)} regions)")

x_regions = regions[regions['chr_type'] == 'X chromosome']
if len(x_regions) > 0:
    x_file = regions_dir / "selection_regions_X_chromosome.tsv"
    x_regions.to_csv(x_file, sep='\t', index=False)
    print(f"  Saved: {x_file.name} ({len(x_regions)} regions)")

# Save stringent regions (min 3 SNPs)
print()
print(f"Creating stringent regions (min {MIN_SNPS_STRINGENT} SNPs)...")
regions_stringent = cluster_candidates_to_regions(
    filtered_candidates, WINDOW_SIZE, MIN_SNPS_STRINGENT
)
stringent_file = regions_dir / "selection_regions_250kb_min3snps.tsv"
regions_stringent.to_csv(stringent_file, sep='\t', index=False)
print(f"  Saved: {stringent_file.name} ({len(regions_stringent)} regions)")



print("=" * 70)
print("Region definition complete!")
print("=" * 70)
