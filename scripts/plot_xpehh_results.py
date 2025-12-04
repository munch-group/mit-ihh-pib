#!/usr/bin/env python3
"""
Plot XP-EHH results for chromosome X
Generates Manhattan plots, distribution plots, and identifies top candidate regions
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10

# Directories
RESULTS_DIR = Path("/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/xp_ehh")
PLOTS_DIR = Path("/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/plots/xp_ehh")
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

# Comparisons
comparisons = [
    "chrX_EUR_vs_EAS",
    "chrX_EUR_vs_AFR",
    "chrX_AFR_vs_EAS"
]

def parse_position(variant_id):
    """Extract position from variant ID (format: X:position:ref:alt)"""
    try:
        parts = variant_id.split(':')
        return int(parts[1])
    except:
        return None

def load_xpehh_data(filepath):
    """Load and process XP-EHH data"""
    print(f"Loading {filepath.name}...")
    df = pd.read_csv(filepath, sep='\t')

    # Parse position from ID
    df['Position'] = df['ID'].apply(parse_position)

    # Remove rows with invalid positions or nan standardized scores
    df = df.dropna(subset=['Position', 'std XPEHH'])
    df = df[np.isfinite(df['std XPEHH'])]

    # Convert position to Mb for plotting
    df['Position_Mb'] = df['Position'] / 1e6

    print(f"  Loaded {len(df):,} variants")
    print(f"  Position range: {df['Position'].min():,} - {df['Position'].max():,}")
    print(f"  std XPEHH range: {df['std XPEHH'].min():.2f} - {df['std XPEHH'].max():.2f}")

    return df

def plot_manhattan(df, comparison_name, output_dir):
    """Create Manhattan plot for XP-EHH scores"""
    fig, ax = plt.subplots(figsize=(14, 6))

    # Plot all points
    scatter = ax.scatter(df['Position_Mb'], df['std XPEHH'],
                        c=df['std XPEHH'], cmap='RdBu_r',
                        s=1, alpha=0.6, vmin=-5, vmax=5)

    # Add significance thresholds (commonly |XP-EHH| > 2 or 3)
    ax.axhline(y=2, color='red', linestyle='--', linewidth=0.8, alpha=0.5, label='|XP-EHH| = 2')
    ax.axhline(y=-2, color='red', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.axhline(y=3, color='darkred', linestyle='--', linewidth=0.8, alpha=0.5, label='|XP-EHH| = 3')
    ax.axhline(y=-3, color='darkred', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5, alpha=0.3)

    ax.set_xlabel('Position on Chromosome X (Mb)', fontsize=12)
    ax.set_ylabel('Standardized XP-EHH', fontsize=12)
    ax.set_title(f'XP-EHH Manhattan Plot: {comparison_name}', fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, label='std XP-EHH')

    plt.tight_layout()
    output_file = output_dir / f'{comparison_name}_manhattan.png'
    plt.savefig(output_file, bbox_inches='tight')
    print(f"  Saved: {output_file}")
    plt.close()

def plot_distribution(df, comparison_name, output_dir):
    """Plot distribution of XP-EHH scores"""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Histogram
    ax = axes[0]
    ax.hist(df['std XPEHH'], bins=100, color='steelblue', alpha=0.7, edgecolor='black')
    ax.axvline(x=2, color='red', linestyle='--', linewidth=1, label='|XP-EHH| = 2')
    ax.axvline(x=-2, color='red', linestyle='--', linewidth=1)
    ax.axvline(x=3, color='darkred', linestyle='--', linewidth=1, label='|XP-EHH| = 3')
    ax.axvline(x=-3, color='darkred', linestyle='--', linewidth=1)
    ax.set_xlabel('Standardized XP-EHH', fontsize=11)
    ax.set_ylabel('Frequency', fontsize=11)
    ax.set_title('Distribution of XP-EHH Scores', fontsize=12)
    ax.legend()

    # QQ plot
    ax = axes[1]
    from scipy import stats
    stats.probplot(df['std XPEHH'], dist="norm", plot=ax)
    ax.set_title('Q-Q Plot', fontsize=12)
    ax.grid(True, alpha=0.3)

    fig.suptitle(f'{comparison_name}', fontsize=14, fontweight='bold')
    plt.tight_layout()

    output_file = output_dir / f'{comparison_name}_distribution.png'
    plt.savefig(output_file, bbox_inches='tight')
    print(f"  Saved: {output_file}")
    plt.close()

def identify_candidates(df, comparison_name, threshold=2.0):
    """Identify candidate regions under selection"""
    # Positive XP-EHH: selection in pop A
    # Negative XP-EHH: selection in pop B

    pop_a = comparison_name.split('_vs_')[0].replace('chrX_', '')
    pop_b = comparison_name.split('_vs_')[1]

    candidates_a = df[df['std XPEHH'] > threshold].copy()
    candidates_b = df[df['std XPEHH'] < -threshold].copy()

    print(f"\n{comparison_name}:")
    print(f"  Variants with evidence of selection in {pop_a} (XP-EHH > {threshold}): {len(candidates_a):,}")
    print(f"  Variants with evidence of selection in {pop_b} (XP-EHH < -{threshold}): {len(candidates_b):,}")

    return candidates_a, candidates_b

def find_peak_regions(df, threshold=2.0, window_size=500000):
    """Find genomic regions with clusters of high XP-EHH scores"""
    significant = df[np.abs(df['std XPEHH']) > threshold].copy()

    if len(significant) == 0:
        return pd.DataFrame()

    significant = significant.sort_values('Position')

    # Group nearby variants into regions
    regions = []
    current_region = {
        'start': significant.iloc[0]['Position'],
        'end': significant.iloc[0]['Position'],
        'variants': [significant.iloc[0]],
        'max_score': significant.iloc[0]['std XPEHH']
    }

    for idx in range(1, len(significant)):
        pos = significant.iloc[idx]['Position']
        score = significant.iloc[idx]['std XPEHH']

        if pos - current_region['end'] <= window_size:
            # Extend current region
            current_region['end'] = pos
            current_region['variants'].append(significant.iloc[idx])
            if abs(score) > abs(current_region['max_score']):
                current_region['max_score'] = score
        else:
            # Save current region and start new one
            regions.append(current_region)
            current_region = {
                'start': pos,
                'end': pos,
                'variants': [significant.iloc[idx]],
                'max_score': score
            }

    # Add last region
    regions.append(current_region)

    # Convert to dataframe
    region_df = pd.DataFrame([
        {
            'chr': 'X',
            'start': r['start'],
            'end': r['end'],
            'n_variants': len(r['variants']),
            'max_xpehh': r['max_score'],
            'length_kb': (r['end'] - r['start']) / 1000
        }
        for r in regions
    ])

    region_df = region_df.sort_values('max_xpehh', key=abs, ascending=False)

    return region_df

def plot_comparison_overview(all_data, output_dir):
    """Create overview plot comparing all three population comparisons"""
    fig, axes = plt.subplots(3, 1, figsize=(14, 12), sharex=True)

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

    for idx, (comparison, df) in enumerate(all_data.items()):
        ax = axes[idx]

        # Scatter plot
        scatter = ax.scatter(df['Position_Mb'], df['std XPEHH'],
                           c=df['std XPEHH'], cmap='RdBu_r',
                           s=1, alpha=0.5, vmin=-5, vmax=5)

        # Thresholds
        ax.axhline(y=2, color='red', linestyle='--', linewidth=0.5, alpha=0.4)
        ax.axhline(y=-2, color='red', linestyle='--', linewidth=0.5, alpha=0.4)
        ax.axhline(y=0, color='black', linestyle='-', linewidth=0.3, alpha=0.3)

        ax.set_ylabel('std XP-EHH', fontsize=11)
        ax.set_title(comparison.replace('_', ' '), fontsize=12, fontweight='bold')
        ax.set_ylim(-6, 6)
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel('Position on Chromosome X (Mb)', fontsize=12)
    fig.suptitle('XP-EHH Analysis: All Population Comparisons',
                 fontsize=14, fontweight='bold', y=0.995)

    plt.tight_layout()
    output_file = output_dir / 'all_comparisons_overview.png'
    plt.savefig(output_file, bbox_inches='tight')
    print(f"  Saved: {output_file}")
    plt.close()

def main():
    print("=" * 60)
    print("XP-EHH Analysis Plotting")
    print("=" * 60)

    all_data = {}
    all_candidates = {}

    for comparison in comparisons:
        print(f"\nProcessing {comparison}...")
        filepath = RESULTS_DIR / comparison

        if not filepath.exists():
            print(f"  WARNING: {filepath} not found, skipping...")
            continue

        # Load data
        df = load_xpehh_data(filepath)
        all_data[comparison] = df

        # Manhattan plot
        print("  Creating Manhattan plot...")
        plot_manhattan(df, comparison, PLOTS_DIR)

        # Distribution plot
        print("  Creating distribution plot...")
        plot_distribution(df, comparison, PLOTS_DIR)

        # Identify candidates
        print("  Identifying candidate regions...")
        cand_a, cand_b = identify_candidates(df, comparison, threshold=2.0)
        all_candidates[comparison] = {'pop_a': cand_a, 'pop_b': cand_b}

        # Find peak regions
        print("  Finding peak regions...")
        regions = find_peak_regions(df, threshold=2.0, window_size=500000)

        if len(regions) > 0:
            print(f"\n  Top 10 candidate regions (threshold |XP-EHH| > 2):")
            print(regions.head(10).to_string(index=False))

            # Save to file
            regions_file = PLOTS_DIR / f'{comparison}_candidate_regions.tsv'
            regions.to_csv(regions_file, sep='\t', index=False)
            print(f"\n  Saved all candidate regions to: {regions_file}")

    # Overview comparison plot
    if len(all_data) > 0:
        print("\nCreating overview comparison plot...")
        plot_comparison_overview(all_data, PLOTS_DIR)

    print("\n" + "=" * 60)
    print("Analysis complete!")
    print(f"All plots saved to: {PLOTS_DIR}")
    print("=" * 60)

if __name__ == "__main__":
    main()
