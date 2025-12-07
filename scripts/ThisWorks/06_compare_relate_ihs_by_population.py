#!/usr/bin/env python3
"""
Compare RELATE and iHS selection signals by population
Print statements and fancy Warnings by Claude Code 
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300

# Which populations belong to which continental group
POPULATION_MAPPING = {
    'AFR': ['ACB', 'ASW', 'ESN', 'GWD', 'LWK', 'MSL', 'YRI'],
    'EAS': ['CDX', 'CHB', 'CHS', 'JPT', 'KHV'],
    'EUR': ['CEU', 'FIN', 'GBR', 'IBS', 'TSI']
}

def read_ihs_candidates(pop_group):
    # Read the iHS candidates for this population
    ihs_file = f'/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/ihs_{pop_group}/analysis/candidates/candidates/per_chromosome/chrX_candidates_stringent.tsv'

    if not Path(ihs_file).exists():
        print(f"Warning: iHS file not found: {ihs_file}")
        return None

    print(f"\n  Reading iHS candidates from {pop_group}")
    df = pd.read_csv(ihs_file, sep='\t')
    print(f"  Found {len(df)} SNPs")

    return df

def read_relate_for_population(pop_code, relate_dir, threshold=5.0):
    # Read RELATE file for one population
    relate_file = Path(relate_dir) / pop_code / 'haplotypes_demog_sele.sele'

    if not relate_file.exists():
        return None

    # Read the file - it's space delimited
    df = pd.read_csv(relate_file, sep=' ', low_memory=False)

    # Get positions (first column) and scores (last column)
    positions = df.iloc[:, 0].values
    scores = df.iloc[:, -1].values

    # Make dataframe with position and score
    relate_df = pd.DataFrame({
        'pos': positions,
        'relate_score': np.abs(-1 * scores),  # convert to positive
        'population': pop_code
    })

    # Only keep SNPs above threshold
    significant = relate_df[relate_df['relate_score'] >= threshold].copy()

    return significant

def pool_relate_signals(pop_group, relate_dir='/home/vanbruggenmit/mit-ihh-pib/data/relate_chrX',
                        threshold=5.0):
    # Combine RELATE signals from all populations in this group
    print(f"\n  Getting RELATE signals for {pop_group}:")

    pop_codes = POPULATION_MAPPING[pop_group]
    all_signals = []

    # Loop through each population and read their RELATE data
    for pop_code in pop_codes:
        signals = read_relate_for_population(pop_code, relate_dir, threshold)
        if signals is not None and len(signals) > 0:
            all_signals.append(signals)
            print(f"    {pop_code}: {len(signals)} SNPs")

    if len(all_signals) == 0:
        print(f"    No RELATE signals found")
        return None

    # Combine all the dataframes
    pooled = pd.concat(all_signals, ignore_index=True)

    # Remove duplicates - keep highest scoring one if SNP found in multiple pops
    pooled = pooled.sort_values('relate_score', ascending=False)
    pooled = pooled.drop_duplicates('pos', keep='first')

    print(f"  Total: {len(pooled)} unique SNPs")

    return pooled

def find_regional_overlaps(ihs_candidates, relate_signals, window_size=50000):
    # Find which iHS SNPs have a RELATE signal nearby
    overlapping_ihs = []
    overlap_details = []

    # Check each iHS SNP
    for _, ihs_snp in ihs_candidates.iterrows():
        ihs_pos = ihs_snp['pos']

        # Look for RELATE signals within +/- window_size bp
        nearby_relate = relate_signals[
            (relate_signals['pos'] >= ihs_pos - window_size) &
            (relate_signals['pos'] <= ihs_pos + window_size)
        ]

        # If we found any nearby RELATE signals, record it
        if len(nearby_relate) > 0:
            overlapping_ihs.append(ihs_snp.to_dict())

            # Get the best matching RELATE SNP (highest score)
            best_match = nearby_relate.loc[nearby_relate['relate_score'].idxmax()]

            overlap_details.append({
                'ihs_pos': ihs_pos,
                'ihs_neg_log10_p': ihs_snp['neg_log10_p'],
                'relate_pos': best_match['pos'],
                'relate_score': best_match['relate_score'],
                'relate_population': best_match['population'],
                'distance': abs(best_match['pos'] - ihs_pos),
                'n_relate_nearby': len(nearby_relate)
            })

    # Convert to dataframes
    if len(overlapping_ihs) > 0:
        overlap_df = pd.DataFrame(overlapping_ihs)
        details_df = pd.DataFrame(overlap_details)
    else:
        overlap_df = pd.DataFrame()
        details_df = pd.DataFrame()

    return overlap_df, details_df

def create_comparison_plot(ihs_candidates, relate_signals, overlap_df, details_df, pop_group, output_dir):
    # Make a 3-panel plot showing iHS, RELATE, and overlaps

    fig, axes = plt.subplots(3, 1, figsize=(18, 14), sharex=True)

    # Panel 1: iHS signals
    ax = axes[0]
    # Plot all iHS candidates
    ax.scatter(ihs_candidates['pos']/1e6, ihs_candidates['neg_log10_p'],
              c='dodgerblue', s=15, alpha=0.5, label='All iHS candidates', rasterized=True)

    # Highlight overlapping ones in red
    if len(overlap_df) > 0:
        ax.scatter(overlap_df['pos']/1e6, overlap_df['neg_log10_p'],
                  c='red', s=30, alpha=0.8, label='Also in RELATE',
                  zorder=5, edgecolors='darkred', linewidths=0.5)

    ax.set_ylabel('-log₁₀(p) iHS', fontsize=12, fontweight='bold')
    ax.set_title(f'{pop_group}: iHS Selection Candidates (n={len(ihs_candidates):,})',
                fontsize=13, fontweight='bold')
    ax.legend(loc='upper right', fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.3, linestyle=':')

    # Add genomic features
    ax.axvspan(0, 2.7, alpha=0.08, color='green', label='PAR1', zorder=0)
    ax.axvspan(58.1, 63.8, alpha=0.08, color='gray', label='Centromere', zorder=0)

    # Panel 2: RELATE signals
    ax = axes[1]
    if relate_signals is not None and len(relate_signals) > 0:
        # Color code by which population the signal came from
        if 'population' in relate_signals.columns:
            populations = relate_signals['population'].unique()
            colors = plt.cm.tab10(np.linspace(0, 1, len(populations)))
            color_map = dict(zip(populations, colors))

            for pop in populations:
                pop_signals = relate_signals[relate_signals['population'] == pop]
                ax.scatter(pop_signals['pos']/1e6, pop_signals['relate_score'],
                          c=[color_map[pop]], s=15, alpha=0.5,
                          label=pop, rasterized=True)
        else:
            # If no population info, just plot in purple
            ax.scatter(relate_signals['pos']/1e6, relate_signals['relate_score'],
                      c='mediumpurple', s=15, alpha=0.5, rasterized=True)

        # Mark the RELATE positions that overlap with iHS in red
        if len(details_df) > 0:
            ax.scatter(details_df['relate_pos']/1e6, details_df['relate_score'],
                      c='red', s=40, alpha=0.9, label='Also in iHS',
                      zorder=5, edgecolors='darkred', linewidths=0.7)

    ax.set_ylabel('RELATE -log₁₀(P)', fontsize=12, fontweight='bold')
    ax.set_title(f'{pop_group}: RELATE Selection Signals (n={len(relate_signals):,})',
                fontsize=13, fontweight='bold')
    ax.legend(loc='upper right', fontsize=9, framealpha=0.9, ncol=2)
    ax.grid(True, alpha=0.3, linestyle=':')
    ax.axvspan(0, 2.7, alpha=0.08, color='green', zorder=0)
    ax.axvspan(58.1, 63.8, alpha=0.08, color='gray', zorder=0)

    # Plot 3: Density comparison
    ax = axes[2]
    bins = np.linspace(0, 156, 80)

    # iHS density
    ax.hist(ihs_candidates['pos']/1e6, bins=bins, alpha=0.4, color='dodgerblue',
           label=f'iHS candidates (n={len(ihs_candidates):,})',
           edgecolor='blue', linewidth=0.5)

    # RELATE density
    if relate_signals is not None and len(relate_signals) > 0:
        ax.hist(relate_signals['pos']/1e6, bins=bins, alpha=0.4, color='mediumpurple',
               label=f'RELATE signals (n={len(relate_signals):,})',
               edgecolor='purple', linewidth=0.5)

    # Overlap density
    if len(overlap_df) > 0:
        ax.hist(overlap_df['pos']/1e6, bins=bins, alpha=0.7, color='red',
               label=f'Overlapping (n={len(overlap_df):,})',
               edgecolor='darkred', linewidth=0.7)

    ax.set_xlabel('Position on X Chromosome (Mb)', fontsize=12, fontweight='bold')
    ax.set_ylabel('SNP Count per bin', fontsize=12, fontweight='bold')
    ax.set_title('Regional Distribution of Selection Signals', fontsize=13, fontweight='bold')
    ax.legend(loc='upper right', fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.3, linestyle=':', axis='y')
    ax.axvspan(0, 2.7, alpha=0.08, color='green', zorder=0)
    ax.axvspan(58.1, 63.8, alpha=0.08, color='gray', zorder=0)

    plt.tight_layout()
    output_file = Path(output_dir) / f'{pop_group}_ihs_relate_comparison.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  ✓ Saved: {output_file}")

def create_venn_diagram(n_ihs, n_relate, n_overlap, pop_group, output_dir):
    """Create Venn diagram"""
    try:
        from matplotlib_venn import venn2

        fig, ax = plt.subplots(figsize=(10, 8))

        v = venn2(subsets=(n_ihs - n_overlap, n_relate - n_overlap, n_overlap),
                 set_labels=('iHS Candidates', 'RELATE Signals'),
                 ax=ax)

        # Colors
        v.get_patch_by_id('10').set_color('dodgerblue')
        v.get_patch_by_id('10').set_alpha(0.5)
        v.get_patch_by_id('01').set_color('mediumpurple')
        v.get_patch_by_id('01').set_alpha(0.5)
        v.get_patch_by_id('11').set_color('red')
        v.get_patch_by_id('11').set_alpha(0.7)

        # Statistics
        overlap_pct_ihs = n_overlap / n_ihs * 100 if n_ihs > 0 else 0
        overlap_pct_relate = n_overlap / n_relate * 100 if n_relate > 0 else 0
        jaccard = n_overlap / (n_ihs + n_relate - n_overlap) if (n_ihs + n_relate - n_overlap) > 0 else 0

        ax.set_title(f'{pop_group}: Method Concordance\n' +
                    f'{overlap_pct_ihs:.1f}% of iHS confirmed | ' +
                    f'{overlap_pct_relate:.1f}% of RELATE confirmed',
                    fontsize=14, fontweight='bold', pad=20)

        textstr = f'iHS candidates: {n_ihs:,}\n'
        textstr += f'RELATE signals: {n_relate:,}\n'
        textstr += f'Overlapping: {n_overlap:,}\n'
        textstr += f'Jaccard index: {jaccard:.3f}'

        ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=11,
               verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

        plt.tight_layout()
        output_file = Path(output_dir) / f'{pop_group}_venn_diagram.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"  ✓ Saved: {output_file}")

    except ImportError:
        print("  Note: matplotlib_venn not available, skipping Venn diagram")

def main():
    # Setup
    output_dir = Path('/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/relate_ihs_comparison')
    output_dir.mkdir(parents=True, exist_ok=True)

    relate_threshold = 5.0  # RELATE score threshold
    window_size = 50000  # Look for overlaps within 50kb

    results_summary = []

    # Process each continental population
    for pop_group in ['AFR', 'EAS', 'EUR']:
        print(f"\n{'='*70}")
        print(f"Analyzing {pop_group}")
        print(f"{'='*70}")

        # Read iHS candidates
        ihs_candidates = read_ihs_candidates(pop_group)
        if ihs_candidates is None:
            continue

        # Pool RELATE signals from all populations in this group
        relate_signals = pool_relate_signals(pop_group, threshold=relate_threshold)
        if relate_signals is None:
            print(f"  Warning: No RELATE signals for {pop_group}")
            continue

        # Find overlaps
        overlap_df, details_df = find_regional_overlaps(
            ihs_candidates, relate_signals, window_size
        )

        n_overlap = len(overlap_df)
        print(f"\n  Overlaps (±{window_size/1000:.0f} kb window): {n_overlap:,}")

        if n_overlap > 0:
            overlap_pct_ihs = n_overlap / len(ihs_candidates) * 100
            overlap_pct_relate = n_overlap / len(relate_signals) * 100

            print(f"  {overlap_pct_ihs:.1f}% of iHS candidates confirmed by RELATE")
            print(f"  {overlap_pct_relate:.1f}% of RELATE signals confirmed by iHS")

            # Save overlap details
            overlap_df.to_csv(output_dir / f'{pop_group}_ihs_overlapping.tsv',
                            sep='\t', index=False)
            details_df.to_csv(output_dir / f'{pop_group}_overlap_details.tsv',
                            sep='\t', index=False)

            print(f"  ✓ Saved overlap tables")

        # Create visualizations
        create_comparison_plot(ihs_candidates, relate_signals, overlap_df, details_df,
                             pop_group, output_dir)

        create_venn_diagram(len(ihs_candidates), len(relate_signals),
                          n_overlap, pop_group, output_dir)

        # Summary statistics
        jaccard = n_overlap / (len(ihs_candidates) + len(relate_signals) - n_overlap) if (len(ihs_candidates) + len(relate_signals) - n_overlap) > 0 else 0

        results_summary.append({
            'Population': pop_group,
            'RELATE_populations': ','.join(POPULATION_MAPPING[pop_group]),
            'iHS_candidates': len(ihs_candidates),
            'RELATE_signals': len(relate_signals),
            'Overlapping_regions': n_overlap,
            'iHS_overlap_pct': overlap_pct_ihs if n_overlap > 0 else 0,
            'RELATE_overlap_pct': overlap_pct_relate if n_overlap > 0 else 0,
            'Jaccard_index': jaccard
        })

    # Save summary
    if results_summary:
        summary_df = pd.DataFrame(results_summary)
        summary_df.to_csv(output_dir / 'population_comparison_summary.tsv',
                         sep='\t', index=False)

        print(f"\n{'='*70}")
        print("SUMMARY OF RESULTS:")
        print(f"{'='*70}\n")
        print(summary_df.to_string(index=False))

        # Interpretation
        print(f"\n{'='*70}")
        print("INTERPRETATION:")
        print(f"{'='*70}")

        for _, row in summary_df.iterrows():
            pop = row['Population']
            ihs_pct = row['iHS_overlap_pct']
            relate_pct = row['RELATE_overlap_pct']

            print(f"\n{pop}:")
            print(f"  iHS candidates: {row['iHS_candidates']:,}")
            print(f"  RELATE signals (pooled): {row['RELATE_signals']:,}")
            print(f"  Overlapping: {row['Overlapping_regions']:,}")

            if ihs_pct > 40:
                print(f"  ✓ STRONG: {ihs_pct:.1f}% of iHS confirmed by RELATE")
            elif ihs_pct > 20:
                print(f"  ≈ MODERATE: {ihs_pct:.1f}% of iHS confirmed by RELATE")
            else:
                print(f"  ✗ WEAK: {ihs_pct:.1f}% of iHS confirmed by RELATE")

            if relate_pct > 40:
                print(f"  ✓ STRONG: {relate_pct:.1f}% of RELATE confirmed by iHS")
            elif relate_pct > 20:
                print(f"  ≈ MODERATE: {relate_pct:.1f}% of RELATE confirmed by iHS")
            else:
                print(f"  → RELATE detecting additional signals (older/polygenic?)")

        print(f"\n✓ All results saved to: {output_dir}")

if __name__ == '__main__':
    main()
