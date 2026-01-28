#!/usr/bin/env python3
"""
Script: Define selection regions from iHS candidates
Purpose: Cluster nearby significant SNPs into selection regions

Approach (based on literature):
- Villegas-Miron 2021: Multiple contiguous SNPs needed, not single outliers
- Approach: 20kb windows requiring ≥20 significant SNPs
- We use 20kb windows with ≥20 SNPs for initial scan, more stringent with ≥50 
- Then merge overlapping windows

Author: Mit Van Bruggen
Date: 2025-10-23
Updated: 2025-12-07
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Setup paths (inefficient because you need to do it for 26 populations in the end and manually put in their Q99, which is why I made script 11 in the end)
base_dir = Path("/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib")
candidates_dir = base_dir / "results/ihs_EUR/analysis/candidates/candidates"
regions_dir = base_dir / "results/ihs_EUR/analysis/regions"
regions_dir.mkdir(parents=True, exist_ok=True)

print("=" * 70)
print("iHS Selection Region Definition")
print("=" * 70)
print()

# Configuration
WINDOW_SIZE = 20000  # 20kb windows
STEP_SIZE = 4000      # 4kb steps (80% overlap - Villegas-Mirón approach)
MIN_SNPS = 20          # Minimum SNPs per region
MIN_SNPS_STRINGENT = 50  # More stringent filter

# Threshold configuration
# Use X-specific empirical threshold 
AUTOSOME_THRESHOLD = 2.5
X_THRESHOLD = 3.024  # X chromosome 99th percentile for that population


def load_candidates(threshold_level="moderate"):
    """Load candidate variants for specified threshold"""
    file = candidates_dir / f"ALL_candidates_{threshold_level}.tsv"
    print(f"Loading candidates from: {file.name}")
    df = pd.read_csv(file, sep='\t')
    print(f"  Total candidates: {len(df):,}")
    return df

def apply_chromosome_specific_threshold(df, autosome_thresh, x_thresh):
    """Apply different thresholds for autosomes vs X chromosome""" #in the end we didn't analyze autosomes, refer to script 11 
    autosome_mask = (df['chr_type'] == 'Autosome') & (np.abs(df['Std iHS']) >= autosome_thresh)
    x_mask = (df['chr_type'] == 'X chromosome') & (np.abs(df['Std iHS']) >= x_thresh)

    filtered = df[autosome_mask | x_mask].copy()

    return filtered

def cluster_candidates_to_regions(candidates, window_size, min_snps):
    """
    Cluster candidates into selection regions using sliding windows
    (Villegas-Mirón et al. 2021 approach)

    Algorithm:
    1. Sort candidates by chromosome and position
    2. For each chromosome, scan with overlapping sliding windows (80% overlap)
    3. Identify windows with >= min_snps candidates
    4. Merge consecutive/overlapping windows into regions
    5. Calculate region statistics
    """
    regions = []

    # Process each chromosome separately
    # Sort: autosomes numerically, then X
    def sort_chr(x):
        if x == 'X':
            return (1, 0)  # X comes last
        else:
            # Handle both "21" and "chr21" formats, and integer types, had some issues with this initially so instead of changing the data I added this
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
        chr_start = int(positions.min())
        chr_end = int(positions.max())

        # SLIDING WINDOW SCAN with fixed step size
        outlier_windows = []
        for win_start in range(chr_start, chr_end, STEP_SIZE):
            win_end = win_start + window_size

            # Find all SNPs in this window
            in_window = (positions >= win_start) & (positions < win_end)
            window_indices = np.where(in_window)[0]

            # Keep windows with sufficient SNPs
            if len(window_indices) >= min_snps:
                outlier_windows.append({
                    'start': win_start,
                    'end': win_end,
                    'indices': window_indices.tolist()
                })

        # MERGE overlapping/consecutive windows into regions
        if len(outlier_windows) == 0:
            print(f"-> 0 regions")
            continue

        clusters = []
        current_cluster = {
            'start': outlier_windows[0]['start'],
            'end': outlier_windows[0]['end'],
            'indices': set(outlier_windows[0]['indices'])
        }

        for window in outlier_windows[1:]:
            # Check if window overlaps with current cluster
            if window['start'] <= current_cluster['end']:
                # Extend and merge
                current_cluster['end'] = max(current_cluster['end'], window['end'])
                current_cluster['indices'].update(window['indices'])
            else:
                # Save current cluster and start new one
                clusters.append(current_cluster)
                current_cluster = {
                    'start': window['start'],
                    'end': window['end'],
                    'indices': set(window['indices'])
                }

        # Don't forget last cluster
        clusters.append(current_cluster)

        # Calculate statistics for each region
        for cluster in clusters:
            # Convert indices back to list and get SNPs
            indices_list = sorted(list(cluster['indices']))
            region_snps = chr_data.iloc[indices_list]

            # Use actual SNP positions for region boundaries (not window boundaries)
            region_start = region_snps['pos'].min()
            region_end = region_snps['pos'].max()

            # Calculate region statistics
            n_snps = len(region_snps)
            region_length = region_end - region_start
            mean_std_ihs = region_snps['Std iHS'].mean()
            max_std_ihs = region_snps['Std iHS'].abs().max()
            mean_pvalue = region_snps['p_value'].mean()
            min_pvalue = region_snps['p_value'].min()

            # Get top SNP
            top_snp_idx = region_snps['Std iHS'].abs().idxmax()
            top_snp = region_snps.loc[top_snp_idx]

            regions.append({
                'chr': chr_name,
                'start': region_start,
                'end': region_end,
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
output_file = regions_dir / "selection_regions_20kb_min20snps.tsv"
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

# Save stringent regions (min 50 SNPs)
print()
print(f"Creating stringent regions (min {MIN_SNPS_STRINGENT} SNPs)...")
regions_stringent = cluster_candidates_to_regions(
    filtered_candidates, WINDOW_SIZE, MIN_SNPS_STRINGENT
)
stringent_file = regions_dir / "selection_regions_20kb_min50snps.tsv"
regions_stringent.to_csv(stringent_file, sep='\t', index=False)
print(f"  Saved: {stringent_file.name} ({len(regions_stringent)} regions)")



print("=" * 70)
print("Region definition complete!")
print("=" * 70)
