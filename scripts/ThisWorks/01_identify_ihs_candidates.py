#!/usr/bin/env python3
"""
Script: Identify iHS selection candidates across all chromosomes
Purpose: Extract significant iHS signals and compare X vs autosomes

Based on literature recommendations:
- Voight et al. (2006): Standardized iHS ~ N(0,1) under neutrality
- Paul et al. (2024): |iHS| >= 2.5 threshold
- Salazar-Tortosa et al. (2023): Check empirical distributions

Author: MIT IHH PIB Project
Date: 2025-01-23
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Setup paths
base_dir = Path("/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib")
ihs_dir = base_dir / "results/ihs_fixed_par1"
output_dir = base_dir / "results/ihs_fixed_par1/analysis/candidates"

# Create output directories
(output_dir / "candidates/per_chromosome").mkdir(parents=True, exist_ok=True)
(output_dir / "plots").mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("iHS Selection Candidate Identification")
print("=" * 60)
print()

# Define thresholds based on literature recommendations
# Voight et al. 2006: standardized iHS ~ N(0,1) under neutrality
# Common thresholds:
#   |iHS| > 2.0  ≈ top ~2.3% (2 std dev from mean)
#   |iHS| > 2.5  ≈ top ~0.6% (used in Paul et al. 2024)
#   |iHS| > 3.0  ≈ top ~0.13% (very stringent)
thresholds = {
    "liberal": 2.0,     # Voight et al. 2006 suggestion
    "moderate": 2.5,    # Paul et al. 2024, more stringent
    "stringent": 3.0    # Very stringent, strongest signals
}

# Will also compute p-values from standard normal distribution
# and empirical percentiles for comparison

# Define chromosomes
chroms = [str(i) for i in range(1, 23)] + ['X']

print("Step 1: Processing each chromosome separately...")
print()

# Process each chromosome and extract candidates
chr_stats_list = []
all_candidates = {}

for chr_name in chroms:
    file = ihs_dir / f"ALL.chr{chr_name}.ihs.tsv"

    if not file.exists():
        print(f"  WARNING: Missing file for chr{chr_name}")
        continue

    print(f"  Processing chr{chr_name}...", end=" ")

    # Read chromosome data
    chr_data = pd.read_csv(file, sep='\t')

    # Parse Location column - handle two formats:
    # Format 1: chr:pos:ref:alt (most chromosomes)
    # Format 2: chr21:pos_ref_alt (chr21 only)
    if chr_data['Location'].iloc[0].count(':') == 3:
        # Standard format: chr:pos:ref:alt
        location_parts = chr_data['Location'].str.split(':', expand=True)
        chr_data['chr'] = location_parts[0]
        chr_data['pos'] = location_parts[1].astype(int)
        chr_data['ref'] = location_parts[2]
        chr_data['alt'] = location_parts[3]
    else:
        # Alternative format: chr21:pos_ref_alt (underscore separators)
        # First split on first colon to get chr and rest
        first_split = chr_data['Location'].str.split(':', n=1, expand=True)
        chr_data['chr'] = first_split[0]
        # Then split the rest on underscores
        rest_parts = first_split[1].str.split('_', expand=True)
        chr_data['pos'] = rest_parts[0].astype(int)
        chr_data['ref'] = rest_parts[1]
        chr_data['alt'] = rest_parts[2]

    print(f"{len(chr_data):,} variants", end="")

    # Add p-values based on standard normal assumption
    # Two-tailed p-value: P(|Z| > |iHS|) where Z ~ N(0,1)
    chr_data['p_value'] = 2 * stats.norm.sf(np.abs(chr_data['Std iHS']))
    chr_data['neg_log10_p'] = -np.log10(chr_data['p_value'])

    # Calculate statistics for this chromosome
    chr_stats = {
        'chr': chr_name,
        'n_variants': len(chr_data),
        'mean_iHS': chr_data['iHS'].mean(),
        'median_iHS': chr_data['iHS'].median(),
        'sd_iHS': chr_data['iHS'].std(),
        'mean_std_iHS': chr_data['Std iHS'].mean(),
        'median_std_iHS': chr_data['Std iHS'].median(),
        'sd_std_iHS': chr_data['Std iHS'].std(),
        'max_abs_std_iHS': np.abs(chr_data['Std iHS']).max(),
        'q99': np.abs(chr_data['Std iHS']).quantile(0.99),
        'q995': np.abs(chr_data['Std iHS']).quantile(0.995),
        'q999': np.abs(chr_data['Std iHS']).quantile(0.999)
    }
    chr_stats_list.append(chr_stats)

    # Extract candidates at each threshold for this chromosome
    for level, thresh in thresholds.items():
        candidates = chr_data[np.abs(chr_data['Std iHS']) >= thresh].copy()
        candidates = candidates.sort_values('Std iHS', key=abs, ascending=False)

        if len(candidates) > 0:
            # Save per-chromosome candidates
            outfile = output_dir / f"candidates/per_chromosome/chr{chr_name}_candidates_{level}.tsv"
            candidates.to_csv(outfile, sep='\t', index=False)

            # Store for combined file
            key = f"{level}_{chr_name}"
            all_candidates[key] = candidates

    print(" - Done")

print()

# Combine statistics across chromosomes
chr_stats = pd.DataFrame(chr_stats_list)
chr_stats['chr'] = pd.Categorical(chr_stats['chr'],
                                   categories=[str(i) for i in range(1, 23)] + ['X'],
                                   ordered=True)

print("Step 2: Calculating summary statistics...")
print()

# Add chromosome type to stats
chr_stats['chr_type'] = chr_stats['chr'].apply(
    lambda x: 'X chromosome' if x == 'X' else 'Autosome'
)

# Overall statistics by chromosome type
overall_stats = chr_stats.groupby('chr_type').agg({
    'chr': 'count',
    'n_variants': 'sum',
    'mean_std_iHS': 'mean',
    'sd_std_iHS': 'mean',
    'q99': 'mean',
    'q995': 'mean'
}).rename(columns={'chr': 'n_chromosomes', 'n_variants': 'total_variants'})
overall_stats['mean_variants_per_chr'] = overall_stats['total_variants'] / overall_stats['n_chromosomes']

print("=== Distribution Check (Voight et al. 2006) ===")
print("Standardized iHS should have mean ~0, variance ~1 under neutrality")
print()
print("Overall Statistics by Chromosome Type:")
print(overall_stats.to_string())
print()

print("=== Empirical Percentiles (|Std iHS|) ===")
print("These show the empirical thresholds for top 1%, 0.5%, 0.1% of SNPs")
chr_percentiles = chr_stats[['chr', 'chr_type', 'q99', 'q995', 'q999']].sort_values('chr')
print(chr_percentiles.to_string(index=False))
print()

print("=== Full Per-chromosome Statistics ===")
print(chr_stats.to_string(index=False))
print()

chr_stats.to_csv(output_dir / "candidates/chromosome_statistics.tsv",
                 sep='\t', index=False)
print("  Saved: candidates/chromosome_statistics.tsv")

# Save empirical percentiles summary
percentile_summary = chr_stats[['chr', 'chr_type', 'n_variants', 'q99', 'q995', 'q999']].sort_values('chr')
percentile_summary.to_csv(output_dir / "candidates/empirical_percentiles.tsv",
                          sep='\t', index=False)
print("  Saved: candidates/empirical_percentiles.tsv")
print()

# Check if distributions are approximately normal
print("=== Normality Check ===")
print("For standard normal, mean should be ~0, SD should be ~1")
normality_check = chr_stats[['chr', 'chr_type', 'mean_std_iHS', 'sd_std_iHS']].copy()
normality_check['mean_deviation'] = np.abs(normality_check['mean_std_iHS'])
normality_check['sd_deviation'] = np.abs(normality_check['sd_std_iHS'] - 1)
print(normality_check.to_string(index=False))
print()
print("Note: Large deviations suggest non-normal distribution (Salazar-Tortosa et al. 2023)")
print("      Consider using empirical percentiles rather than fixed thresholds")
print()

print("Step 3: Combining candidates across chromosomes...")
print()

# Combine candidates at each threshold level
candidates_list = {}

for level, thresh in thresholds.items():
    # Get all candidates for this level across all chromosomes
    level_candidates = []
    for chr_name in chroms:
        key = f"{level}_{chr_name}"
        if key in all_candidates:
            level_candidates.append(all_candidates[key])

    if level_candidates:
        combined = pd.concat(level_candidates, ignore_index=True)
        combined['chr_type'] = combined['chr'].apply(
            lambda x: 'X chromosome' if x == 'X' else 'Autosome'
        )

        print(f"  {level.capitalize()} threshold (|Std iHS| >= {thresh}): {len(combined):,} variants")

        # Count by chromosome type
        by_type = combined.groupby('chr_type').size()
        pct = (by_type / len(combined) * 100)

        if 'Autosome' in by_type.index:
            print(f"    - Autosomes: {by_type['Autosome']:,} ({pct['Autosome']:.1f}%)")
        if 'X chromosome' in by_type.index:
            print(f"    - X chromosome: {by_type['X chromosome']:,} ({pct['X chromosome']:.1f}%)")

        # Save combined candidates
        outfile = output_dir / f"candidates/ALL_candidates_{level}.tsv"
        combined.to_csv(outfile, sep='\t', index=False)
        print(f"    Saved: candidates/ALL_candidates_{level}.tsv")
        print()

        candidates_list[level] = combined

print("Step 4: Testing for X chromosome enrichment...")
print()

# Test if X chromosome is enriched for selection signals
# Using moderate threshold (|Std iHS| >= 2.5)
candidates_moderate = candidates_list['moderate']

# Contingency table
n_candidates_x = (candidates_moderate['chr'] == 'X').sum()
n_candidates_auto = (candidates_moderate['chr'] != 'X').sum()
n_total_x = chr_stats[chr_stats['chr'] == 'X']['n_variants'].sum()
n_total_auto = chr_stats[chr_stats['chr'] != 'X']['n_variants'].sum()

contingency = np.array([
    [n_candidates_x, n_total_x - n_candidates_x],
    [n_candidates_auto, n_total_auto - n_candidates_auto]
])

# Fisher's exact test
odds_ratio, p_value = stats.fisher_exact(contingency)

print(f"  X chromosome candidates: {n_candidates_x:,} / {n_total_x:,} ({n_candidates_x/n_total_x*100:.2f}%)")
print(f"  Autosome candidates: {n_candidates_auto:,} / {n_total_auto:,} ({n_candidates_auto/n_total_auto*100:.2f}%)")
print(f"  Odds ratio: {odds_ratio:.3f}")
print(f"  P-value: {p_value:.2e}")
print()

# Save test results
enrichment_results = pd.DataFrame([{
    'test': "Fisher's exact test",
    'threshold': 'moderate (|Std iHS| >= 2.5)',
    'x_candidates': n_candidates_x,
    'x_total': n_total_x,
    'x_pct': n_candidates_x / n_total_x * 100,
    'auto_candidates': n_candidates_auto,
    'auto_total': n_total_auto,
    'auto_pct': n_candidates_auto / n_total_auto * 100,
    'odds_ratio': odds_ratio,
    'p_value': p_value
}])

enrichment_results.to_csv(output_dir / "candidates/x_enrichment_test.tsv",
                          sep='\t', index=False)
print("  Saved: candidates/x_enrichment_test.tsv")
print()

print("Step 5: Creating visualizations...")
print()

# Set plot style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300

# Plot 1: Candidate counts by chromosome
print("  Plotting candidate counts by chromosome...")

candidates_by_chr = candidates_moderate.groupby(['chr', 'chr_type']).size().reset_index(name='n')
candidates_by_chr['chr'] = pd.Categorical(candidates_by_chr['chr'],
                                          categories=[str(i) for i in range(1, 23)] + ['X'],
                                          ordered=True)
candidates_by_chr = candidates_by_chr.sort_values('chr')

fig, ax = plt.subplots(figsize=(12, 6))
colors = {'Autosome': 'gray', 'X chromosome': '#E41A1C'}
for chr_type in ['Autosome', 'X chromosome']:
    data = candidates_by_chr[candidates_by_chr['chr_type'] == chr_type]
    ax.bar(data['chr'], data['n'], label=chr_type, color=colors[chr_type])

ax.set_xlabel('Chromosome')
ax.set_ylabel('Number of Candidates')
ax.set_title('Number of Selection Candidates by Chromosome\nModerate threshold (|Std iHS| >= 2.5)')
ax.legend(title='Chromosome Type')
plt.tight_layout()
plt.savefig(output_dir / "plots/candidates_by_chromosome.png")
plt.close()
print("    Saved: plots/candidates_by_chromosome.png")

# Plot 2: Candidate rate by chromosome
print("  Plotting candidate rate by chromosome...")

candidate_rates = candidates_by_chr.merge(chr_stats[['chr', 'n_variants']], on='chr')
candidate_rates['rate'] = candidate_rates['n'] / candidate_rates['n_variants'] * 100

fig, ax = plt.subplots(figsize=(12, 6))
for chr_type in ['Autosome', 'X chromosome']:
    data = candidate_rates[candidate_rates['chr_type'] == chr_type]
    ax.bar(data['chr'], data['rate'], label=chr_type, color=colors[chr_type])

ax.set_xlabel('Chromosome')
ax.set_ylabel('Candidate Rate (%)')
ax.set_title('Selection Candidate Rate by Chromosome\nModerate threshold (|Std iHS| >= 2.5)')
ax.legend(title='Chromosome Type')
plt.tight_layout()
plt.savefig(output_dir / "plots/candidate_rate_by_chromosome.png")
plt.close()
print("    Saved: plots/candidate_rate_by_chromosome.png")

# Plot 3: Distribution of Std iHS in candidates (X vs autosomes)
print("  Plotting Std iHS distribution in candidates...")

fig, ax = plt.subplots(figsize=(10, 6))
for chr_type in ['Autosome', 'X chromosome']:
    data = candidates_moderate[candidates_moderate['chr_type'] == chr_type]['Std iHS']
    ax.hist(data, bins=50, alpha=0.6, label=chr_type, color=colors[chr_type])

ax.set_xlabel('Standardized iHS')
ax.set_ylabel('Count')
ax.set_title('Distribution of Standardized iHS in Selection Candidates\nModerate threshold (|Std iHS| >= 2.5)')
ax.legend(title='Chromosome Type')
plt.tight_layout()
plt.savefig(output_dir / "plots/candidates_std_ihs_distribution.png")
plt.close()
print("    Saved: plots/candidates_std_ihs_distribution.png")
print()

print("=" * 60)
print("Analysis Complete!")
print("=" * 60)
print()

print("Output files:")
print("\nPer-chromosome candidates:")
print("  - candidates/per_chromosome/chr*_candidates_*.tsv (69 files)")
print("    Each file includes: Location, chr, pos, ref, alt, iHH_0, iHH_1,")
print("                        iHS, Std_iHS, p_value, neg_log10_p")
print("\nCombined files:")
print("  - candidates/chromosome_statistics.tsv (includes empirical percentiles)")
print("  - candidates/empirical_percentiles.tsv (99th, 99.5th, 99.9th percentiles)")
print("  - candidates/ALL_candidates_liberal.tsv (|Std iHS| >= 2.0)")
print("  - candidates/ALL_candidates_moderate.tsv (|Std iHS| >= 2.5)")
print("  - candidates/ALL_candidates_stringent.tsv (|Std iHS| >= 3.0)")
print("  - candidates/x_enrichment_test.tsv")
print("\nPlots:")
print("  - plots/candidates_by_chromosome.png")
print("  - plots/candidate_rate_by_chromosome.png")
print("  - plots/candidates_std_ihs_distribution.png")
print()

print("=== Threshold Interpretation (Literature-based) ===")
print("Based on Voight et al. 2006, Paul et al. 2024, and others:")
print()
print("Fixed thresholds (assuming Std iHS ~ N(0,1)):")
print("  |Std iHS| >= 2.0: ~2.3% of genome (liberal, Voight et al. 2006)")
print("  |Std iHS| >= 2.5: ~0.6% of genome (moderate, Paul et al. 2024)")
print("  |Std iHS| >= 3.0: ~0.13% of genome (stringent)")
print()
print("Equivalent p-value thresholds:")
print("  |Std iHS| >= 2.0: p < 0.046, -log10(p) > 1.34")
print("  |Std iHS| >= 2.5: p < 0.012, -log10(p) > 1.91")
print("  |Std iHS| >= 3.0: p < 0.0027, -log10(p) > 2.57")
print("  Paul et al. 2024 used: -log10(p) > 4 (p < 0.0001)")
print()
print("IMPORTANT: Check empirical percentiles in your data!")
print("  If distributions deviate from normality, use empirical percentiles")
print("  rather than fixed thresholds (Salazar-Tortosa et al. 2023)")
print()

print("Next steps:")
print("  1. Review the normality check and empirical percentiles")
print("  2. Check if X chromosome distribution differs from autosomes")
print("  3. Decide on threshold: fixed (2.0/2.5/3.0) or empirical (top 1%/0.5%)")
print("  4. Review plots to understand selection patterns")
print("  5. Run 02_define_selection_regions.py to cluster candidates into peaks")
print("      (Following literature: require multiple SNPs per region, not single outliers)")
print("  6. Annotate peaks with genes")
print("  7. Perform enrichment analysis for reproductive/immune genes")
print()
