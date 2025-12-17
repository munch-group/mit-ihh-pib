#!/usr/bin/env python3
"""
Pairwise comparison of RELATE and iHS selection signals by individual populations.

This script uses a percentile-based approach for iHS filtering and bidirectional matching:
1. Calculates population-specific 99.5th percentile threshold for each population
2. Ensures comparable stringency across populations with different distributions
3. Uses BIDIRECTIONAL matching: finds iHS signals with nearby RELATE AND RELATE signals with nearby iHS
4. Compares each individual population's iHS vs RELATE results pairwise
5. Creates visualizations showing separate signals and overlaps for each population
6. Aggregates the pairwise results into continental groups (AFR, AMR, EAS, EUR, SAS)

The bidirectional approach ensures that both iHS and RELATE peaks are properly highlighted,
addressing the issue where multiple signals in a peak might match the same counterpart signal.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Tuple, Optional

sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300

# Population mapping to continental groups (1000 Genomes populations)
POPULATION_MAPPING = {
    'AFR': ['ACB', 'ASW', 'ESN', 'GWD', 'LWK', 'MSL', 'YRI'],
    'AMR': ['CLM', 'MXL', 'PEL', 'PUR'],  # Admixed American populations
    'EAS': ['CDX', 'CHB', 'CHS', 'JPT', 'KHV'],
    'EUR': ['CEU', 'FIN', 'GBR', 'IBS', 'TSI'],
    'SAS': ['BEB', 'GIH', 'ITU', 'PJL', 'STU']  # South Asian populations
}

# Reverse mapping: population -> continental group
POP_TO_CONTINENT = {}
for continent, pops in POPULATION_MAPPING.items():
    for pop in pops:
        POP_TO_CONTINENT[pop] = continent


def calculate_ihs_threshold(pop_code: str, percentile: float = 99.5) -> Optional[float]:
    """
    Calculate the percentile-based threshold for iHS scores for a population.

    Args:
        pop_code: Population code (e.g., 'CEU', 'YRI')
        percentile: Percentile to use for threshold (default: 99.0)

    Returns:
        Threshold value or None if file not found (claude addition)
    """
    ihs_file = f'/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/ihs_{pop_code}/{pop_code}.chrX.ihs.tsv'

    if not Path(ihs_file).exists():
        return None

    # Read the iHS results
    df = pd.read_csv(ihs_file, sep='\t')

    # Calculate absolute standardized iHS scores, removing NaN values (needed because ASW has Nan Std.iHS scores)
    abs_std_ihs = df['Std iHS'].abs()
    abs_std_ihs = abs_std_ihs[~abs_std_ihs.isna()]  # Remove NaN values

    if len(abs_std_ihs) == 0:
        return None

    # Calculate the percentile threshold
    threshold = np.percentile(abs_std_ihs, percentile)

    return threshold


def read_ihs_results(pop_code: str,
                     percentile_threshold: float = 99.5) -> Optional[pd.DataFrame]:
    """
    Read iHS results for a single population and filter by percentile-based threshold.

    This approach uses the raw iHS data and calculates a population-specific threshold
    based on the empirical distribution, ensuring comparable stringency across populations.

    Args:
        pop_code: Population code (e.g., 'CEU', 'YRI')
        percentile_threshold: Percentile to use for threshold (default: 99.0, but trying 99.5 to make it more stringent)

    Returns:
        DataFrame with significant iHS signals or None if file not found
    """
    ihs_file = f'/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/ihs_{pop_code}/{pop_code}.chrX.ihs.tsv'

    if not Path(ihs_file).exists():
        print(f"  Warning: iHS file not found for {pop_code}")
        return None

    # Read the raw iHS results
    df = pd.read_csv(ihs_file, sep='\t')

    # Parse position from Location column
    df['pos'] = df['Location'].str.split(':').str[1].astype(int)

    # Calculate absolute standardized iHS scores, removing NaN values
    df['abs_std_ihs'] = df['Std iHS'].abs()

    # Remove rows with NaN Std iHS values (e.g., ASW has some)
    df_clean = df[~df['abs_std_ihs'].isna()].copy()

    if len(df_clean) == 0:
        print(f"  Warning: No valid iHS scores for {pop_code}")
        return None

    # Calculate population-specific threshold at the given percentile
    threshold = np.percentile(df_clean['abs_std_ihs'], percentile_threshold)

    # Filter to signals above threshold
    significant = df_clean[df_clean['abs_std_ihs'] >= threshold].copy()

    print(f"  {pop_code}: {len(df_clean):,} total variants, threshold={threshold:.3f} (p{percentile_threshold})")
    print(f"  {pop_code}: {len(significant):,} signals >= {threshold:.3f}")

    # Add population label
    significant['population'] = pop_code

    return significant


def read_relate_results(pop_code: str,
                       relate_threshold: float = 5.0) -> Optional[pd.DataFrame]:
    """
    Read RELATE results for a single population and filter by significance threshold.

    Args:
        pop_code: Population code (e.g., 'CEU', 'YRI')
        relate_threshold: -log10(p) threshold for selection

    Returns:
        DataFrame with significant RELATE signals or None if file not found
    """
    relate_file = f'/home/vanbruggenmit/mit-ihh-pib/data/relate_chrX/{pop_code}/haplotypes_demog_sele.sele'

    if not Path(relate_file).exists():
        print(f"  Warning: RELATE file not found for {pop_code}")
        return None

    # Read RELATE results (space-delimited)
    df = pd.read_csv(relate_file, sep=' ', low_memory=False)

    # First column is position, last column is the selection score
    positions = df.iloc[:, 0].values
    scores = df.iloc[:, -1].values

    # Create clean dataframe
    relate_df = pd.DataFrame({
        'pos': positions,
        'relate_score': np.abs(-1 * scores),  # Convert to positive
        'population': pop_code
    })

    # Filter by significance threshold
    significant = relate_df[relate_df['relate_score'] >= relate_threshold].copy()

    print(f"  {pop_code}: {len(significant):,} significant RELATE signals (score >= {relate_threshold})")

    return significant


def find_regional_overlaps(ihs_df: pd.DataFrame,
                          relate_df: pd.DataFrame,
                          window_size: int = 50000) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Find overlapping selection signals between iHS and RELATE within a genomic window.
    Uses bidirectional matching and identifies mutual overlaps.

    Args:
        ihs_df: DataFrame with iHS signals
        relate_df: DataFrame with RELATE signals
        window_size: Window size in base pairs for considering signals as overlapping

    Returns:
        Tuple of (overlapping_ihs_snps, overlapping_relate_snps, mutual_ihs, mutual_relate, overlap_details)
        - overlapping_ihs_snps: iHS signals with nearby RELATE (one-directional)
        - overlapping_relate_snps: RELATE signals with nearby iHS (one-directional)
        - mutual_ihs: iHS signals that are in mutual overlap
        - mutual_relate: RELATE signals that are in mutual overlap
    """
    overlapping_ihs = []
    overlapping_relate = []
    overlap_details = []

    # Direction 1: For each iHS signal, find nearby RELATE signals
    for _, ihs_snp in ihs_df.iterrows():
        ihs_pos = ihs_snp['pos']

        # Find RELATE signals within window
        nearby_relate = relate_df[
            (relate_df['pos'] >= ihs_pos - window_size) &
            (relate_df['pos'] <= ihs_pos + window_size)
        ]

        if len(nearby_relate) > 0:
            overlapping_ihs.append(ihs_snp.to_dict())

            # Get the closest/highest scoring RELATE match
            best_match = nearby_relate.loc[nearby_relate['relate_score'].idxmax()]

            overlap_details.append({
                'ihs_pos': ihs_pos,
                'ihs_std_score': ihs_snp['Std iHS'],
                'ihs_abs_score': ihs_snp['abs_std_ihs'],
                'relate_pos': best_match['pos'],
                'relate_score': best_match['relate_score'],
                'distance': abs(best_match['pos'] - ihs_pos),
                'n_relate_nearby': len(nearby_relate)
            })

    # Direction 2: For each RELATE signal, check if it has nearby iHS signals
    for _, relate_snp in relate_df.iterrows():
        relate_pos = relate_snp['pos']

        # Find iHS signals within window
        nearby_ihs = ihs_df[
            (ihs_df['pos'] >= relate_pos - window_size) &
            (ihs_df['pos'] <= relate_pos + window_size)
        ]

        if len(nearby_ihs) > 0:
            overlapping_relate.append(relate_snp.to_dict())

    overlap_ihs_df = pd.DataFrame(overlapping_ihs) if overlapping_ihs else pd.DataFrame()
    overlap_relate_df = pd.DataFrame(overlapping_relate) if overlapping_relate else pd.DataFrame()
    overlap_details_df = pd.DataFrame(overlap_details) if overlap_details else pd.DataFrame()

    # Identify mutual overlaps
    # Mutual iHS: iHS signals that have nearby RELATE, AND those RELATE signals also have nearby iHS
    mutual_ihs = []
    if len(overlap_ihs_df) > 0 and len(overlap_relate_df) > 0:
        overlap_relate_positions = set(overlap_relate_df['pos'].values)

        for _, ihs_snp in overlap_ihs_df.iterrows():
            ihs_pos = ihs_snp['pos']
            # Find RELATE signals within window
            nearby_relate = relate_df[
                (relate_df['pos'] >= ihs_pos - window_size) &
                (relate_df['pos'] <= ihs_pos + window_size)
            ]
            # Check if any of these nearby RELATE signals are also in overlap_relate_df
            if any(pos in overlap_relate_positions for pos in nearby_relate['pos'].values):
                mutual_ihs.append(ihs_snp.to_dict())

    # Mutual RELATE: RELATE signals that have nearby iHS, AND those iHS signals also have nearby RELATE
    mutual_relate = []
    if len(overlap_ihs_df) > 0 and len(overlap_relate_df) > 0:
        overlap_ihs_positions = set(overlap_ihs_df['pos'].values)

        for _, relate_snp in overlap_relate_df.iterrows():
            relate_pos = relate_snp['pos']
            # Find iHS signals within window
            nearby_ihs = ihs_df[
                (ihs_df['pos'] >= relate_pos - window_size) &
                (ihs_df['pos'] <= relate_pos + window_size)
            ]
            # Check if any of these nearby iHS signals are also in overlap_ihs_df
            if any(pos in overlap_ihs_positions for pos in nearby_ihs['pos'].values):
                mutual_relate.append(relate_snp.to_dict())

    mutual_ihs_df = pd.DataFrame(mutual_ihs) if mutual_ihs else pd.DataFrame()
    mutual_relate_df = pd.DataFrame(mutual_relate) if mutual_relate else pd.DataFrame()

    return overlap_ihs_df, overlap_relate_df, mutual_ihs_df, mutual_relate_df, overlap_details_df


def create_pairwise_comparison_plot(ihs_df: pd.DataFrame,
                                   relate_df: pd.DataFrame,
                                   overlap_ihs_df: pd.DataFrame,
                                   overlap_relate_df: pd.DataFrame,
                                   mutual_ihs_df: pd.DataFrame,
                                   mutual_relate_df: pd.DataFrame,
                                   overlap_details_df: pd.DataFrame,
                                   pop_code: str,
                                   output_dir: Path):
    """
    Create a 3-panel plot showing iHS signals, RELATE signals, and their overlaps.
    Uses three-color scheme:
    - Orange: One-directional overlap (iHS→RELATE or RELATE→iHS)
    - Red: Mutual/bidirectional overlap (strong reciprocal evidence)
    - Blue/Purple: No overlap
    """
    fig, axes = plt.subplots(3, 1, figsize=(18, 14), sharex=True)

    # Panel 1: iHS signals
    ax = axes[0]
    y_col = 'abs_std_ihs'
    y_label = '|Standardized iHS|'

    # Get positions for categorization
    mutual_ihs_positions = set(mutual_ihs_df['pos'].values) if len(mutual_ihs_df) > 0 else set()
    overlap_ihs_positions = set(overlap_ihs_df['pos'].values) if len(overlap_ihs_df) > 0 else set()

    # One-directional only: in overlap but not in mutual
    onedir_ihs_positions = overlap_ihs_positions - mutual_ihs_positions

    # Plot all signals first
    ax.scatter(ihs_df['pos']/1e6, ihs_df[y_col],
              c='dodgerblue', s=15, alpha=0.5, label='No overlap', rasterized=True)

    # Plot one-directional overlaps (orange)
    if len(onedir_ihs_positions) > 0:
        onedir_ihs = ihs_df[ihs_df['pos'].isin(onedir_ihs_positions)]
        ax.scatter(onedir_ihs['pos']/1e6, onedir_ihs[y_col],
                  c='orange', s=30, alpha=0.7, label='One-directional (iHS→RELATE)',
                  zorder=4, edgecolors='darkorange', linewidths=0.5)

    # Plot mutual overlaps (red) - plot last so they're on top
    if len(mutual_ihs_df) > 0:
        ax.scatter(mutual_ihs_df['pos']/1e6, mutual_ihs_df[y_col],
                  c='red', s=35, alpha=0.9, label='Mutual overlap',
                  zorder=5, edgecolors='darkred', linewidths=0.6)

    ax.set_ylabel(y_label, fontsize=12, fontweight='bold')
    ax.set_title(f'{pop_code}: iHS Selection Signals (n={len(ihs_df):,})',
                fontsize=13, fontweight='bold')
    ax.legend(loc='upper right', fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.3, linestyle=':')
    ax.axvspan(0, 2.7, alpha=0.08, color='green', zorder=0)
    ax.axvspan(58.1, 63.8, alpha=0.08, color='gray', zorder=0)

    # Panel 2: RELATE signals
    ax = axes[1]

    # Get positions for categorization
    mutual_relate_positions = set(mutual_relate_df['pos'].values) if len(mutual_relate_df) > 0 else set()
    overlap_relate_positions = set(overlap_relate_df['pos'].values) if len(overlap_relate_df) > 0 else set()

    # One-directional only: in overlap but not in mutual
    onedir_relate_positions = overlap_relate_positions - mutual_relate_positions

    # Plot all signals first
    ax.scatter(relate_df['pos']/1e6, relate_df['relate_score'],
              c='mediumpurple', s=15, alpha=0.5, label='No overlap', rasterized=True)

    # Plot one-directional overlaps (orange/gold)
    if len(onedir_relate_positions) > 0:
        onedir_relate = relate_df[relate_df['pos'].isin(onedir_relate_positions)]
        ax.scatter(onedir_relate['pos']/1e6, onedir_relate['relate_score'],
                  c='orange', s=35, alpha=0.7, label='One-directional (RELATE→iHS)',
                  zorder=4, edgecolors='darkorange', linewidths=0.6)

    # Plot mutual overlaps (red) - plot last so they're on top
    if len(mutual_relate_df) > 0:
        ax.scatter(mutual_relate_df['pos']/1e6, mutual_relate_df['relate_score'],
                  c='red', s=40, alpha=0.9, label='Mutual overlap',
                  zorder=5, edgecolors='darkred', linewidths=0.7)

    ax.set_ylabel('RELATE -log₁₀(P)', fontsize=12, fontweight='bold')
    ax.set_title(f'{pop_code}: RELATE Selection Signals (n={len(relate_df):,})',
                fontsize=13, fontweight='bold')
    ax.legend(loc='upper right', fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.3, linestyle=':')
    ax.axvspan(0, 2.7, alpha=0.08, color='green', zorder=0)
    ax.axvspan(58.1, 63.8, alpha=0.08, color='gray', zorder=0)

    # Panel 3: Regional density comparison
    ax = axes[2]
    bins = np.linspace(0, 156, 80)

    ax.hist(ihs_df['pos']/1e6, bins=bins, alpha=0.4, color='dodgerblue',
           label=f'iHS signals (n={len(ihs_df):,})',
           edgecolor='blue', linewidth=0.5)

    ax.hist(relate_df['pos']/1e6, bins=bins, alpha=0.4, color='mediumpurple',
           label=f'RELATE signals (n={len(relate_df):,})',
           edgecolor='purple', linewidth=0.5)

    if len(overlap_ihs_df) > 0:
        ax.hist(overlap_ihs_df['pos']/1e6, bins=bins, alpha=0.7, color='red',
               label=f'Overlapping (n={len(overlap_ihs_df):,})',
               edgecolor='darkred', linewidth=0.7)

    ax.set_xlabel('Position on X Chromosome (Mb)', fontsize=12, fontweight='bold')
    ax.set_ylabel('SNP Count per bin', fontsize=12, fontweight='bold')
    ax.set_title('Regional Distribution of Selection Signals', fontsize=13, fontweight='bold')
    ax.legend(loc='upper right', fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.3, linestyle=':', axis='y')
    ax.axvspan(0, 2.7, alpha=0.08, color='green', zorder=0)
    ax.axvspan(58.1, 63.8, alpha=0.08, color='gray', zorder=0)

    plt.tight_layout()
    output_file = output_dir / f'{pop_code}_pairwise_comparison.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"    Saved: {output_file.name}")


def create_venn_diagram(n_ihs: int, n_relate: int, n_overlap: int,
                       pop_code: str, output_dir: Path):
    """Create Venn diagram showing overlap between methods."""
    try:
        from matplotlib_venn import venn2

        fig, ax = plt.subplots(figsize=(10, 8))

        v = venn2(subsets=(n_ihs - n_overlap, n_relate - n_overlap, n_overlap),
                 set_labels=('iHS Signals', 'RELATE Signals'),
                 ax=ax)

        # Set colors for patches that exist (need ifs because for example for CHB there is no overlap found and then it crashes)
        if v.get_patch_by_id('10') is not None:
            v.get_patch_by_id('10').set_color('dodgerblue')
            v.get_patch_by_id('10').set_alpha(0.5)
        if v.get_patch_by_id('01') is not None:
            v.get_patch_by_id('01').set_color('mediumpurple')
            v.get_patch_by_id('01').set_alpha(0.5)
        if v.get_patch_by_id('11') is not None:
            v.get_patch_by_id('11').set_color('red')
            v.get_patch_by_id('11').set_alpha(0.7)

        overlap_pct_ihs = n_overlap / n_ihs * 100 if n_ihs > 0 else 0
        overlap_pct_relate = n_overlap / n_relate * 100 if n_relate > 0 else 0
        jaccard = n_overlap / (n_ihs + n_relate - n_overlap) if (n_ihs + n_relate - n_overlap) > 0 else 0

        ax.set_title(f'{pop_code}: Method Concordance\n' +
                    f'{overlap_pct_ihs:.1f}% of iHS confirmed | ' +
                    f'{overlap_pct_relate:.1f}% of RELATE confirmed',
                    fontsize=14, fontweight='bold', pad=20)

        textstr = f'iHS signals: {n_ihs:,}\n'
        textstr += f'RELATE signals: {n_relate:,}\n'
        textstr += f'Overlapping: {n_overlap:,}\n'
        textstr += f'Jaccard index: {jaccard:.3f}'

        ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=11,
               verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

        plt.tight_layout()
        output_file = output_dir / f'{pop_code}_venn_diagram.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"    Saved: {output_file.name}")

    except ImportError:
        print("    Note: matplotlib_venn not available, skipping Venn diagram")


def process_single_population(pop_code: str,
                              ihs_percentile: float,
                              relate_threshold: float,
                              window_size: int,
                              output_dir: Path) -> Optional[dict]:
    """
    Process a single population: read data, find overlaps, create plots.

    Args:
        pop_code: Population code
        ihs_percentile: Percentile threshold for iHS (e.g., 99.0 for top 1%)
        relate_threshold: Fixed RELATE threshold
        window_size: Window size for overlaps
        output_dir: Output directory

    Returns:
        Dictionary with summary statistics or None if data unavailable
    """
    print(f"\n  Processing {pop_code}:")

    # Read iHS and RELATE data (iHS filtered by percentile)
    ihs_df = read_ihs_results(pop_code, percentile_threshold=ihs_percentile)
    relate_df = read_relate_results(pop_code, relate_threshold)

    if ihs_df is None or relate_df is None:
        print(f"    Skipping {pop_code} due to missing data")
        return None

    if len(ihs_df) == 0 or len(relate_df) == 0:
        print(f"    Skipping {pop_code} due to no significant signals")
        return None

    # Find overlaps (bidirectional with mutual identification)
    overlap_ihs_df, overlap_relate_df, mutual_ihs_df, mutual_relate_df, overlap_details_df = find_regional_overlaps(
        ihs_df, relate_df, window_size
    )

    n_overlap_ihs = len(overlap_ihs_df)  # Number of iHS signals with nearby RELATE
    n_overlap_relate = len(overlap_relate_df)  # Number of RELATE signals with nearby iHS
    n_mutual_ihs = len(mutual_ihs_df)  # Number of iHS in mutual overlap
    n_mutual_relate = len(mutual_relate_df)  # Number of RELATE in mutual overlap

    overlap_pct_ihs = n_overlap_ihs / len(ihs_df) * 100 if len(ihs_df) > 0 else 0
    overlap_pct_relate = n_overlap_relate / len(relate_df) * 100 if len(relate_df) > 0 else 0
    mutual_pct_ihs = n_mutual_ihs / len(ihs_df) * 100 if len(ihs_df) > 0 else 0
    mutual_pct_relate = n_mutual_relate / len(relate_df) * 100 if len(relate_df) > 0 else 0

    print(f"    Overlaps (±{window_size/1000:.0f} kb):")
    print(f"      {n_overlap_ihs:,} iHS signals have nearby RELATE ({overlap_pct_ihs:.1f}%)")
    print(f"      {n_overlap_relate:,} RELATE signals have nearby iHS ({overlap_pct_relate:.1f}%)")
    print(f"      {n_mutual_ihs:,} iHS in mutual overlap ({mutual_pct_ihs:.1f}%)")
    print(f"      {n_mutual_relate:,} RELATE in mutual overlap ({mutual_pct_relate:.1f}%)")

    # Save overlap tables
    pop_output_dir = output_dir / 'pairwise_bidirectional' / pop_code
    pop_output_dir.mkdir(parents=True, exist_ok=True)

    if n_overlap_ihs > 0:
        overlap_ihs_df.to_csv(pop_output_dir / f'{pop_code}_ihs_overlapping.tsv',
                            sep='\t', index=False)
        overlap_details_df.to_csv(pop_output_dir / f'{pop_code}_overlap_details.tsv',
                                 sep='\t', index=False)

    if n_overlap_relate > 0:
        overlap_relate_df.to_csv(pop_output_dir / f'{pop_code}_relate_overlapping.tsv',
                                sep='\t', index=False)

    if n_mutual_ihs > 0:
        mutual_ihs_df.to_csv(pop_output_dir / f'{pop_code}_ihs_mutual.tsv',
                            sep='\t', index=False)

    if n_mutual_relate > 0:
        mutual_relate_df.to_csv(pop_output_dir / f'{pop_code}_relate_mutual.tsv',
                               sep='\t', index=False)

    # Create visualizations
    create_pairwise_comparison_plot(ihs_df, relate_df, overlap_ihs_df, overlap_relate_df,
                                   mutual_ihs_df, mutual_relate_df,
                                   overlap_details_df, pop_code, pop_output_dir)
    create_venn_diagram(len(ihs_df), len(relate_df), n_overlap_ihs,
                       pop_code, pop_output_dir)

    # Return summary statistics
    # For Jaccard index, use n_overlap_ihs as the primary metric (iHS signals confirmed by RELATE)
    jaccard = n_overlap_ihs / (len(ihs_df) + len(relate_df) - n_overlap_ihs) if (len(ihs_df) + len(relate_df) - n_overlap_ihs) > 0 else 0

    # Get the actual threshold used (percentile-based)
    ihs_threshold = ihs_df['abs_std_ihs'].min() if len(ihs_df) > 0 else None

    return {
        'population': pop_code,
        'continent': POP_TO_CONTINENT.get(pop_code, 'Unknown'),
        'ihs_threshold': ihs_threshold,  # Minimum |Std iHS| in filtered set (= percentile threshold)
        'n_ihs': len(ihs_df),
        'n_relate': len(relate_df),
        'n_overlap_ihs': n_overlap_ihs,  # iHS signals with nearby RELATE
        'n_overlap_relate': n_overlap_relate,  # RELATE signals with nearby iHS
        'n_mutual_ihs': n_mutual_ihs,  # iHS in mutual overlap
        'n_mutual_relate': n_mutual_relate,  # RELATE in mutual overlap
        'pct_ihs_confirmed': overlap_pct_ihs,
        'pct_relate_confirmed': overlap_pct_relate,
        'pct_mutual_ihs': mutual_pct_ihs,
        'pct_mutual_relate': mutual_pct_relate,
        'jaccard_index': jaccard
    }


def aggregate_to_continents(all_results: list, output_dir: Path):
    """
    Aggregate the pairwise results into continental groups for comparison.
    """
    print(f"\n{'='*70}")
    print("AGGREGATING TO CONTINENTAL GROUPS")
    print(f"{'='*70}")

    results_df = pd.DataFrame(all_results)

    continental_summary = []

    for continent in ['AFR', 'EAS', 'EUR']:
        continent_pops = results_df[results_df['continent'] == continent]

        if len(continent_pops) == 0:
            continue

        print(f"\n{continent}:")
        print(f"  Populations: {', '.join(continent_pops['population'].values)}")
        print(f"  Mean overlap: {continent_pops['pct_ihs_confirmed'].mean():.1f}% of iHS confirmed")
        print(f"  Range: {continent_pops['pct_ihs_confirmed'].min():.1f}% - {continent_pops['pct_ihs_confirmed'].max():.1f}%")

        continental_summary.append({
            'continent': continent,
            'n_populations': len(continent_pops),
            'populations': ','.join(continent_pops['population'].values),
            'mean_ihs_signals': continent_pops['n_ihs'].mean(),
            'mean_relate_signals': continent_pops['n_relate'].mean(),
            'mean_overlap_ihs': continent_pops['n_overlap_ihs'].mean(),
            'mean_overlap_relate': continent_pops['n_overlap_relate'].mean(),
            'mean_pct_ihs_confirmed': continent_pops['pct_ihs_confirmed'].mean(),
            'mean_pct_relate_confirmed': continent_pops['pct_relate_confirmed'].mean(),
            'min_pct_ihs_confirmed': continent_pops['pct_ihs_confirmed'].min(),
            'max_pct_ihs_confirmed': continent_pops['pct_ihs_confirmed'].max(),
            'min_pct_relate_confirmed': continent_pops['pct_relate_confirmed'].min(),
            'max_pct_relate_confirmed': continent_pops['pct_relate_confirmed'].max(),
            'mean_jaccard': continent_pops['jaccard_index'].mean()
        })

    # Save continental summary
    continental_df = pd.DataFrame(continental_summary)
    continental_df.to_csv(output_dir / 'continental_summary.tsv', sep='\t', index=False)

    print(f"\n{'='*70}")
    print("CONTINENTAL SUMMARY:")
    print(f"{'='*70}\n")
    print(continental_df.to_string(index=False))

    return continental_df


def create_comparison_figure(results_df: pd.DataFrame, output_dir: Path):
    """
    Create a figure comparing overlap percentages across all populations.
    """
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Sort by continent and overlap percentage
    results_df = results_df.sort_values(['continent', 'pct_ihs_confirmed'])

    # Color by continent
    continent_colors = {
        'AFR': 'orange',
        'AMR': 'red',
        'EAS': 'green',
        'EUR': 'blue',
        'SAS': 'purple'
    }
    colors = [continent_colors.get(c, 'gray') for c in results_df['continent']]

    # Panel 1: Percentage of iHS confirmed by RELATE
    ax = axes[0]
    y_pos = np.arange(len(results_df))
    ax.barh(y_pos, results_df['pct_ihs_confirmed'], color=colors, alpha=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(results_df['population'])
    ax.set_xlabel('% of iHS signals confirmed by RELATE', fontsize=11, fontweight='bold')
    ax.set_title('iHS Confirmation by RELATE (per population)', fontsize=12, fontweight='bold')
    ax.axvline(20, color='red', linestyle='--', alpha=0.5, label='Expected range (20-40%)')
    ax.axvline(40, color='red', linestyle='--', alpha=0.5)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='x')

    # Panel 2: Percentage of RELATE confirmed by iHS
    ax = axes[1]
    ax.barh(y_pos, results_df['pct_relate_confirmed'], color=colors, alpha=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(results_df['population'])
    ax.set_xlabel('% of RELATE signals confirmed by iHS', fontsize=11, fontweight='bold')
    ax.set_title('RELATE Confirmation by iHS (per population)', fontsize=12, fontweight='bold')
    ax.axvline(20, color='red', linestyle='--', alpha=0.5, label='Expected range (20-40%)')
    ax.axvline(40, color='red', linestyle='--', alpha=0.5)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='x')

    plt.tight_layout()
    output_file = output_dir / 'all_populations_overlap_comparison.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\nSaved comparison figure: {output_file.name}")


def main():
    """Main analysis pipeline."""
    print(f"{'='*70}")
    print("PAIRWISE RELATE vs iHS COMPARISON (BIDIRECTIONAL)")
    print(f"{'='*70}")

    # Setup
    output_dir = Path('/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/pairwise_relate_ihs_comparison_bidirectional')
    output_dir.mkdir(parents=True, exist_ok=True)

    # Parameters
    ihs_percentile = 99.5  # Use 99th percentile for each population / now 99.5
    relate_threshold = 5.0  # -log10(p) >= 5.0
    window_size = 50000  # 50kb window for overlaps

    print(f"\nParameters:")
    print(f"  iHS threshold: {ihs_percentile}th percentile (population-specific, ensures comparable stringency)")
    print(f"  RELATE threshold: -log10(p) >= {relate_threshold}")
    print(f"  Overlap window: ±{window_size/1000:.0f} kb")

    # Get all populations
    all_populations = []
    for pops in POPULATION_MAPPING.values():
        all_populations.extend(pops)

    print(f"\nProcessing {len(all_populations)} populations:")
    for continent, pops in POPULATION_MAPPING.items():
        print(f"  {continent}: {', '.join(pops)}")

    # Process each population individually
    print(f"\n{'='*70}")
    print("PAIRWISE COMPARISONS:")
    print(f"{'='*70}")

    all_results = []
    for pop_code in all_populations:
        result = process_single_population(
            pop_code, ihs_percentile, relate_threshold, window_size, output_dir
        )
        if result is not None:
            all_results.append(result)

    # Save pairwise results
    if all_results:
        results_df = pd.DataFrame(all_results)
        results_df.to_csv(output_dir / 'pairwise_population_results.tsv',
                         sep='\t', index=False)

        print(f"\n{'='*70}")
        print("PAIRWISE RESULTS SUMMARY:")
        print(f"{'='*70}\n")
        print(results_df.to_string(index=False))

        # Show threshold summary
        print(f"\n{'='*70}")
        print(f"iHS FILTERING: {ihs_percentile}th percentile (population-specific)")
        print(f"{'='*70}\n")
        if 'ihs_threshold' in results_df.columns:
            threshold_summary = results_df[['population', 'continent', 'ihs_threshold', 'n_ihs']].copy()
            threshold_summary = threshold_summary.sort_values(['continent', 'population'])
            print(threshold_summary.to_string(index=False))
        else:
            print(results_df[['population', 'continent', 'n_ihs']].to_string(index=False))

        # Create comparison figure
        create_comparison_figure(results_df, output_dir)

        # Aggregate to continental groups
        continental_df = aggregate_to_continents(all_results, output_dir)

        # Final interpretation
        print(f"\n{'='*70}")
        print("INTERPRETATION:")
        print(f"{'='*70}")

        for _, row in continental_df.iterrows():
            continent = row['continent']
            mean_pct = row['mean_pct_ihs_confirmed']

            print(f"\n{continent}:")
            print(f"  Populations: {row['n_populations']}")
            print(f"  Mean overlap: {mean_pct:.1f}% (range: {row['min_pct_ihs_confirmed']:.1f}% - {row['max_pct_ihs_confirmed']:.1f}%)")

            if mean_pct >= 20 and mean_pct <= 40:
                print(f"  EXCELLENT: Within expected range (20-40%) from Relate paper")
            elif mean_pct > 40:
                print(f"  STRONG: Above expected range - high concordance")
            elif mean_pct >= 10:
                print(f"  MODERATE: Below expected but reasonable")
            else:
                print(f"  WEAK: Below expected range")

        print(f"\nAll results saved to: {output_dir}")

    else:
        print("\nNo results generated - check data availability")


if __name__ == '__main__':
    main()
