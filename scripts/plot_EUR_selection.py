#!/usr/bin/env python3
"""
Create comprehensive plots for EUR-specific selection analysis
Visualizes the integration of EUR iHS + XP-EHH results
"""

import sys
import math
from pathlib import Path

# Directories
RESULTS_DIR = Path("/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/plots/EUR_selection_analysis")
OUTPUT_DIR = RESULTS_DIR

def load_variants(filepath):
    """Load variant-level data"""
    data = {'positions': [], 'eur_ihs': [], 'eur_vs_eas': [], 'eur_vs_afr': []}

    with open(filepath, 'r') as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                try:
                    data['positions'].append(int(parts[1]) / 1e6)  # Convert to Mb
                    data['eur_ihs'].append(float(parts[2]))
                    data['eur_vs_eas'].append(float(parts[3]))
                    data['eur_vs_afr'].append(float(parts[4]))
                except:
                    continue

    return data

def load_regions(filepath):
    """Load region-level data"""
    regions = []

    with open(filepath, 'r') as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 10:
                try:
                    regions.append({
                        'start': int(parts[1]) / 1e6,
                        'end': int(parts[2]) / 1e6,
                        'n_variants': int(parts[3]),
                        'combined_score': float(parts[4]),
                        'length_kb': float(parts[5]),
                        'peak_pos': int(parts[6]) / 1e6,
                        'peak_ihs': float(parts[7]),
                        'peak_eur_vs_eas': float(parts[8]),
                        'peak_eur_vs_afr': float(parts[9])
                    })
                except:
                    continue

    return regions

def create_matplotlib_plots():
    """Create plots using matplotlib"""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
        from matplotlib.gridspec import GridSpec

        print("Creating plots with matplotlib...")

        # Set style
        plt.rcParams['figure.dpi'] = 300
        plt.rcParams['savefig.dpi'] = 300
        plt.rcParams['font.size'] = 10

        # Load data
        print("Loading data...")
        high_conf_vars = load_variants(RESULTS_DIR / "HIGH_CONFIDENCE_EUR_selection_variants.tsv")
        high_conf_regions = load_regions(RESULTS_DIR / "HIGH_CONFIDENCE_EUR_selection_regions.tsv")
        xpehh_only_regions = load_regions(RESULTS_DIR / "XPEHH_ONLY_EUR_selection_regions.tsv")
        ihs_only_regions = load_regions(RESULTS_DIR / "iHS_ONLY_EUR_selection_regions.tsv")

        # Plot 1: Genome-wide overview with all signal categories
        print("Creating genome-wide overview plot...")
        fig = plt.figure(figsize=(16, 10))
        gs = GridSpec(4, 1, height_ratios=[1, 1, 1, 0.5], hspace=0.3)

        # Panel A: EUR iHS (use absolute value for visualization since we filter on |iHS|)
        ax1 = fig.add_subplot(gs[0])
        if len(high_conf_vars['positions']) > 0:
            ihs_abs = [abs(x) for x in high_conf_vars['eur_ihs']]
            ax1.scatter(high_conf_vars['positions'], ihs_abs,
                       c='red', s=10, alpha=0.6, label='|iHS|>2 + both XP-EHH>3')
        ax1.axhline(y=2, color='gray', linestyle='--', linewidth=0.8, alpha=0.5, label='|iHS|=2 (suggestive)')
        ax1.axhline(y=2.9, color='darkgray', linestyle='--', linewidth=0.8, alpha=0.7, label='|iHS|=2.9 (strong)')
        ax1.set_ylabel('|EUR iHS|', fontsize=11, fontweight='bold')
        ax1.set_title('EUR-Specific Selection Signals on Chromosome X\n(Note: All concordant variants show selection on ancestral allele)',
                     fontsize=13, fontweight='bold')
        ax1.legend(loc='upper right', fontsize=8)
        ax1.grid(True, alpha=0.3)

        # Panel B: EUR vs EAS XP-EHH
        ax2 = fig.add_subplot(gs[1], sharex=ax1)
        if len(high_conf_vars['positions']) > 0:
            ax2.scatter(high_conf_vars['positions'], high_conf_vars['eur_vs_eas'],
                       c='blue', s=10, alpha=0.6, label='Multi-method concordant')
        ax2.axhline(y=3, color='gray', linestyle='--', linewidth=0.8, alpha=0.5, label='XP-EHH=3 (strong)')
        ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.5, alpha=0.3)
        ax2.set_ylabel('EUR vs EAS\nXP-EHH', fontsize=11, fontweight='bold')
        ax2.legend(loc='upper right', fontsize=8)
        ax2.grid(True, alpha=0.3)

        # Panel C: EUR vs AFR XP-EHH
        ax3 = fig.add_subplot(gs[2], sharex=ax1)
        if len(high_conf_vars['positions']) > 0:
            ax3.scatter(high_conf_vars['positions'], high_conf_vars['eur_vs_afr'],
                       c='green', s=10, alpha=0.6, label='Multi-method concordant')
        ax3.axhline(y=3, color='gray', linestyle='--', linewidth=0.8, alpha=0.5, label='XP-EHH=3 (strong)')
        ax3.axhline(y=0, color='black', linestyle='-', linewidth=0.5, alpha=0.3)
        ax3.set_ylabel('EUR vs AFR\nXP-EHH', fontsize=11, fontweight='bold')
        ax3.legend(loc='upper right', fontsize=8)
        ax3.grid(True, alpha=0.3)

        # Panel D: Signal category regions
        ax4 = fig.add_subplot(gs[3], sharex=ax1)
        y_pos = 0

        # Plot high confidence regions
        for region in high_conf_regions:
            width = region['end'] - region['start']
            rect = patches.Rectangle((region['start'], y_pos - 0.3), width, 0.6,
                                    linewidth=0, facecolor='red', alpha=0.7)
            ax4.add_patch(rect)

        # Plot XP-EHH only regions
        y_pos = 1
        for region in xpehh_only_regions[:50]:  # Limit to top 50
            width = region['end'] - region['start']
            rect = patches.Rectangle((region['start'], y_pos - 0.3), width, 0.6,
                                    linewidth=0, facecolor='purple', alpha=0.5)
            ax4.add_patch(rect)

        # Plot iHS only regions
        y_pos = 2
        for region in ihs_only_regions[:50]:  # Limit to top 50
            width = region['end'] - region['start']
            rect = patches.Rectangle((region['start'], y_pos - 0.3), width, 0.6,
                                    linewidth=0, facecolor='orange', alpha=0.5)
            ax4.add_patch(rect)

        ax4.set_ylim(-0.5, 2.5)
        ax4.set_yticks([0, 1, 2])
        ax4.set_yticklabels(['All 3 Methods\nConcordant', 'XP-EHH Only\n(Near Fix.)', 'iHS Only\n(Incomplete)'],
                           fontsize=9)
        ax4.set_xlabel('Position on Chromosome X (Mb)', fontsize=12, fontweight='bold')
        ax4.set_xlim(0, 160)
        ax4.set_xticks(range(0, 161, 20))  # Add x-axis ticks every 20 Mb
        ax4.grid(True, alpha=0.3, axis='x')

        plt.tight_layout()
        outfile = OUTPUT_DIR / 'EUR_selection_genome_wide_overview.png'
        plt.savefig(outfile, bbox_inches='tight')
        print(f"  Saved: {outfile}")
        plt.close()

        # Plot 2: Method agreement scatter plots
        print("Creating method agreement plot...")
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))

        # EUR iHS vs EUR vs EAS XP-EHH
        ax = axes[0, 0]
        if len(high_conf_vars['positions']) > 0:
            scatter = ax.scatter(high_conf_vars['eur_ihs'], high_conf_vars['eur_vs_eas'],
                               c=high_conf_vars['eur_vs_afr'], cmap='RdYlGn',
                               s=20, alpha=0.6, vmin=0, vmax=10)
            plt.colorbar(scatter, ax=ax, label='EUR vs AFR XP-EHH')
        ax.axhline(y=3, color='gray', linestyle='--', linewidth=0.8, alpha=0.5, label='XP-EHH=3')
        ax.axvline(x=2, color='gray', linestyle='--', linewidth=0.8, alpha=0.5, label='|iHS|=2')
        ax.axvline(x=-2, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
        ax.axvline(x=2.9, color='darkgray', linestyle='--', linewidth=0.8, alpha=0.7, label='|iHS|=2.9')
        ax.axvline(x=-2.9, color='darkgray', linestyle='--', linewidth=0.8, alpha=0.7)
        ax.set_xlabel('EUR iHS (neg=ancestral, pos=derived)', fontsize=11)
        ax.set_ylabel('EUR vs EAS XP-EHH', fontsize=11)
        ax.set_title('Multi-method Concordant:\niHS vs EUR vs EAS XP-EHH', fontsize=12, fontweight='bold')
        ax.legend(loc='upper left', fontsize=7)
        ax.grid(True, alpha=0.3)

        # EUR iHS vs EUR vs AFR XP-EHH
        ax = axes[0, 1]
        if len(high_conf_vars['positions']) > 0:
            scatter = ax.scatter(high_conf_vars['eur_ihs'], high_conf_vars['eur_vs_afr'],
                               c=high_conf_vars['eur_vs_eas'], cmap='RdYlBu',
                               s=20, alpha=0.6, vmin=0, vmax=10)
            plt.colorbar(scatter, ax=ax, label='EUR vs EAS XP-EHH')
        ax.axhline(y=3, color='gray', linestyle='--', linewidth=0.8, alpha=0.5, label='XP-EHH=3')
        ax.axvline(x=2, color='gray', linestyle='--', linewidth=0.8, alpha=0.5, label='|iHS|=2')
        ax.axvline(x=-2, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
        ax.axvline(x=2.9, color='darkgray', linestyle='--', linewidth=0.8, alpha=0.7, label='|iHS|=2.9')
        ax.axvline(x=-2.9, color='darkgray', linestyle='--', linewidth=0.8, alpha=0.7)
        ax.set_xlabel('EUR iHS (neg=ancestral, pos=derived)', fontsize=11)
        ax.set_ylabel('EUR vs AFR XP-EHH', fontsize=11)
        ax.set_title('Multi-method Concordant:\niHS vs EUR vs AFR XP-EHH', fontsize=12, fontweight='bold')
        ax.legend(loc='upper left', fontsize=7)
        ax.grid(True, alpha=0.3)

        # EUR vs EAS vs EUR vs AFR XP-EHH
        ax = axes[1, 0]
        if len(high_conf_vars['positions']) > 0:
            scatter = ax.scatter(high_conf_vars['eur_vs_eas'], high_conf_vars['eur_vs_afr'],
                               c=high_conf_vars['eur_ihs'], cmap='plasma',
                               s=20, alpha=0.6, vmin=2, vmax=6)
            plt.colorbar(scatter, ax=ax, label='EUR iHS')
        ax.axhline(y=3, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
        ax.axvline(x=3, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
        ax.plot([0, 15], [0, 15], 'k--', linewidth=0.8, alpha=0.3, label='y=x')
        ax.set_xlabel('EUR vs EAS XP-EHH', fontsize=11)
        ax.set_ylabel('EUR vs AFR XP-EHH', fontsize=11)
        ax.set_title('High Confidence Variants:\nXP-EHH Comparisons', fontsize=12, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Summary statistics panel
        ax = axes[1, 1]
        ax.axis('off')

        # Count statistics
        n_high = len(high_conf_vars['positions'])
        n_high_regions = len(high_conf_regions)
        n_xpehh_regions = len(xpehh_only_regions)
        n_ihs_regions = len(ihs_only_regions)

        summary_text = f"""
EUR-SPECIFIC SELECTION SUMMARY

Thresholds Used:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  |iHS| > 2.0: Suggestive signal
  |iHS| ≥ 2.9: Strong signal (chrX)
  XP-EHH > 3.0: Strong signal

Signal Categories:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
MULTI-METHOD CONCORDANT
  • {n_high:,} variants
  • {n_high_regions} genomic regions
  • |iHS| > 2.0 AND both XP-EHH > 3.0
  • Strong EUR-specific selection

XP-EHH ONLY (near fixation)
  • {n_xpehh_regions} genomic regions
  • Selection completed
  • Allele at high frequency
  • iHS loses power

iHS ONLY (incomplete sweep)
  • {n_ihs_regions} genomic regions
  • Ongoing selection
  • Allele not yet fixed

Interpretation:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Points along y=x diagonal:
  → Equal signal vs both populations
  → True EUR-specific (not drift)
        """

        ax.text(0.1, 0.95, summary_text, transform=ax.transAxes,
               fontsize=10, verticalalignment='top', family='monospace')

        fig.suptitle('EUR-Specific Selection: Method Agreement Analysis',
                    fontsize=14, fontweight='bold', y=0.995)
        plt.tight_layout()

        outfile = OUTPUT_DIR / 'EUR_selection_method_agreement.png'
        plt.savefig(outfile, bbox_inches='tight')
        print(f"  Saved: {outfile}")
        plt.close()

        # Plot 3: Top regions detail
        print("Creating top regions detail plot...")
        fig, axes = plt.subplots(5, 4, figsize=(16, 14))
        axes = axes.flatten()

        for i, region in enumerate(high_conf_regions[:20]):
            ax = axes[i]

            # Create a simple representation
            start = region['start']
            end = region['end']
            peak = region['peak_pos']

            # Bar showing the region
            ax.barh(1, end - start, left=start, height=0.5, color='red', alpha=0.6)
            ax.plot([peak, peak], [0.5, 1.5], 'k-', linewidth=2)

            # Add scores as text
            ax.text(start + (end - start) / 2, 2.2,
                   f"iHS: {region['peak_ihs']:.2f}\n"
                   f"EAS: {region['peak_eur_vs_eas']:.2f}\n"
                   f"AFR: {region['peak_eur_vs_afr']:.2f}",
                   ha='center', va='top', fontsize=8)

            ax.set_ylim(0, 3)
            ax.set_xlim(start - (end - start) * 0.1, end + (end - start) * 0.1)
            ax.set_title(f"#{i+1}: {start:.2f}-{end:.2f} Mb\n"
                        f"{region['n_variants']} vars, {region['length_kb']:.1f} kb",
                        fontsize=9)
            ax.set_yticks([])
            ax.set_xlabel('Position (Mb)', fontsize=8)
            ax.grid(True, alpha=0.3, axis='x')

        # Hide unused subplots
        for i in range(20, 20):
            axes[i].axis('off')

        fig.suptitle('Top 20 High Confidence EUR-Specific Selection Regions',
                    fontsize=14, fontweight='bold')
        plt.tight_layout()

        outfile = OUTPUT_DIR / 'EUR_selection_top20_regions.png'
        plt.savefig(outfile, bbox_inches='tight')
        print(f"  Saved: {outfile}")
        plt.close()

        # Plot 4: Signal strength distributions
        print("Creating signal strength distribution plot...")
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        # EUR iHS distribution (showing sign to indicate ancestral vs derived)
        ax = axes[0, 0]
        if len(high_conf_vars['eur_ihs']) > 0:
            ax.hist(high_conf_vars['eur_ihs'], bins=50, color='red', alpha=0.7, edgecolor='black')
            ax.axvline(x=2, color='gray', linestyle='--', linewidth=1.5, label='|iHS|=2')
            ax.axvline(x=-2, color='gray', linestyle='--', linewidth=1.5)
            ax.axvline(x=2.9, color='darkgray', linestyle='--', linewidth=1.5, label='|iHS|=2.9')
            ax.axvline(x=-2.9, color='darkgray', linestyle='--', linewidth=1.5)
            mean_ihs = sum(high_conf_vars['eur_ihs']) / len(high_conf_vars['eur_ihs'])
            ax.axvline(x=mean_ihs, color='darkred', linestyle='-', linewidth=1.5, label=f'Mean: {mean_ihs:.2f}')
        ax.set_xlabel('EUR iHS (negative=ancestral allele selection)', fontsize=11)
        ax.set_ylabel('Frequency', fontsize=11)
        ax.set_title('Multi-method Concordant Variants:\niHS Distribution', fontsize=12, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # EUR vs EAS XP-EHH distribution
        ax = axes[0, 1]
        if len(high_conf_vars['eur_vs_eas']) > 0:
            ax.hist(high_conf_vars['eur_vs_eas'], bins=50, color='blue', alpha=0.7, edgecolor='black')
            ax.axvline(x=3, color='gray', linestyle='--', linewidth=1.5, label='XP-EHH=3 (strong)')
            mean_eas = sum(high_conf_vars['eur_vs_eas']) / len(high_conf_vars['eur_vs_eas'])
            ax.axvline(x=mean_eas, color='darkblue', linestyle='-', linewidth=1.5, label=f'Mean: {mean_eas:.2f}')
        ax.set_xlabel('EUR vs EAS XP-EHH', fontsize=11)
        ax.set_ylabel('Frequency', fontsize=11)
        ax.set_title('Multi-method Concordant Variants:\nEUR vs EAS Distribution', fontsize=12, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # EUR vs AFR XP-EHH distribution
        ax = axes[1, 0]
        if len(high_conf_vars['eur_vs_afr']) > 0:
            ax.hist(high_conf_vars['eur_vs_afr'], bins=50, color='green', alpha=0.7, edgecolor='black')
            ax.axvline(x=3, color='gray', linestyle='--', linewidth=1.5, label='XP-EHH=3 (strong)')
            mean_afr = sum(high_conf_vars['eur_vs_afr']) / len(high_conf_vars['eur_vs_afr'])
            ax.axvline(x=mean_afr, color='darkgreen', linestyle='-', linewidth=1.5, label=f'Mean: {mean_afr:.2f}')
        ax.set_xlabel('EUR vs AFR XP-EHH', fontsize=11)
        ax.set_ylabel('Frequency', fontsize=11)
        ax.set_title('Multi-method Concordant Variants:\nEUR vs AFR Distribution', fontsize=12, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Combined score distribution
        ax = axes[1, 1]
        if len(high_conf_vars['eur_ihs']) > 0:
            combined = [abs(high_conf_vars['eur_ihs'][i]) +
                       high_conf_vars['eur_vs_eas'][i] +
                       high_conf_vars['eur_vs_afr'][i]
                       for i in range(len(high_conf_vars['eur_ihs']))]
            ax.hist(combined, bins=50, color='purple', alpha=0.7, edgecolor='black')
            mean_combined = sum(combined) / len(combined)
            ax.axvline(x=mean_combined, color='darkviolet', linestyle='-', linewidth=1.5,
                      label=f'Mean: {mean_combined:.2f}')
        ax.set_xlabel('Combined Score (|iHS| + XP-EHH_EAS + XP-EHH_AFR)', fontsize=11)
        ax.set_ylabel('Frequency', fontsize=11)
        ax.set_title('High Confidence Variants:\nCombined Score Distribution', fontsize=12, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)

        fig.suptitle('EUR-Specific Selection: Signal Strength Distributions',
                    fontsize=14, fontweight='bold', y=0.995)
        plt.tight_layout()

        outfile = OUTPUT_DIR / 'EUR_selection_signal_distributions.png'
        plt.savefig(outfile, bbox_inches='tight')
        print(f"  Saved: {outfile}")
        plt.close()

        return True

    except ImportError:
        print("matplotlib not available")
        return False

def create_summary_statistics():
    """Create a summary statistics file"""
    print("\nCreating summary statistics...")

    # Load all data
    high_conf_regions = load_regions(RESULTS_DIR / "HIGH_CONFIDENCE_EUR_selection_regions.tsv")
    xpehh_only_regions = load_regions(RESULTS_DIR / "XPEHH_ONLY_EUR_selection_regions.tsv")
    ihs_only_regions = load_regions(RESULTS_DIR / "iHS_ONLY_EUR_selection_regions.tsv")

    outfile = OUTPUT_DIR / "EUR_selection_summary_statistics.txt"
    with open(outfile, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("EUR-SPECIFIC SELECTION ANALYSIS - SUMMARY STATISTICS\n")
        f.write("=" * 80 + "\n\n")

        f.write("GENOMIC REGIONS BY CATEGORY:\n")
        f.write("-" * 80 + "\n")
        f.write(f"HIGH CONFIDENCE (all 3 methods): {len(high_conf_regions)} regions\n")
        f.write(f"XP-EHH ONLY (near fixation): {len(xpehh_only_regions)} regions\n")
        f.write(f"iHS ONLY (incomplete sweep): {len(ihs_only_regions)} regions\n\n")

        if len(high_conf_regions) > 0:
            f.write("HIGH CONFIDENCE REGIONS - DETAILED LIST:\n")
            f.write("-" * 80 + "\n")
            f.write(f"{'Rank':<6} {'Position (Mb)':<20} {'Variants':<10} {'Combined':<12} "
                   f"{'iHS':<8} {'EAS':<8} {'AFR':<8}\n")
            f.write("-" * 80 + "\n")

            for i, region in enumerate(high_conf_regions, 1):
                f.write(f"{i:<6} {region['start']:.2f}-{region['end']:.2f} ({region['length_kb']:.1f}kb)"
                       f"{'':<5} {region['n_variants']:<10} {region['combined_score']:<12.2f} "
                       f"{region['peak_ihs']:<8.2f} {region['peak_eur_vs_eas']:<8.2f} "
                       f"{region['peak_eur_vs_afr']:<8.2f}\n")

        f.write("\n" + "=" * 80 + "\n")
        f.write("INTERPRETATION:\n")
        f.write("=" * 80 + "\n\n")
        f.write("HIGH CONFIDENCE regions show EUR-specific selection supported by:\n")
        f.write("  1. EUR iHS > 2 (selection within EUR population)\n")
        f.write("  2. EUR vs EAS XP-EHH > 3 (EUR different from East Asian)\n")
        f.write("  3. EUR vs AFR XP-EHH > 3 (EUR different from African)\n")
        f.write("  → These are NOT due to drift, but true EUR-specific selection\n\n")

        f.write("XP-EHH ONLY regions (near fixation):\n")
        f.write("  - Strong cross-population signal but weak iHS\n")
        f.write("  - Likely completed sweeps with allele at high frequency\n")
        f.write("  - iHS loses power when allele is near fixation\n\n")

        f.write("iHS ONLY regions (incomplete sweep):\n")
        f.write("  - Strong within-population signal but weak XP-EHH\n")
        f.write("  - Likely ongoing or incomplete selective sweeps\n")
        f.write("  - Allele frequency may still be rising\n\n")

    print(f"  Saved: {outfile}")

def main():
    print("=" * 80)
    print("EUR-SPECIFIC SELECTION - PLOTTING")
    print("=" * 80)

    # Create summary statistics
    create_summary_statistics()

    # Try to create matplotlib plots
    print("\n" + "=" * 80)
    print("Creating Visualizations")
    print("=" * 80)

    success = create_matplotlib_plots()

    if success:
        print("\n" + "=" * 80)
        print("All plots created successfully!")
        print(f"Output directory: {OUTPUT_DIR}")
        print("=" * 80)
    else:
        print("\nNote: matplotlib not available for plotting")
        print("Summary statistics have been created.")

if __name__ == "__main__":
    main()
