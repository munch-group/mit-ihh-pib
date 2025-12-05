#!/usr/bin/env python3
"""
Plot population-specific iHS results for chromosome X
Creates 4-panel plot:
- Panel A: EUR iHS
- Panel B: AFR iHS
- Panel C: EAS iHS
- Panel D: All three populations overlaid
"""

import sys
from pathlib import Path

# Directories
RESULTS_DIR = Path("/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results")
PLOTS_DIR = Path("/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/plots")
OUTPUT_DIR = PLOTS_DIR / "population_ihs_comparison"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def parse_position(variant_id):
    """Extract position from variant ID (format: X:position:ref:alt)"""
    try:
        parts = variant_id.split(':')
        return int(parts[1])
    except:
        return None

def load_ihs_data(filepath, population):
    """Load iHS results for a population"""
    print(f"Loading {population} iHS data from {filepath.name}...")

    positions = []
    ihs_scores = []

    with open(filepath, 'r') as f:
        header = f.readline()  # Skip header

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                variant_id = parts[0]
                pos = parse_position(variant_id)

                if pos is None:
                    continue

                try:
                    std_ihs = float(parts[4])
                    # Check for valid values
                    if abs(std_ihs) < 100:  # Filter extreme outliers
                        positions.append(pos / 1e6)  # Convert to Mb
                        ihs_scores.append(std_ihs)
                except:
                    continue

    print(f"  Loaded {len(positions):,} variants")
    if len(ihs_scores) > 0:
        print(f"  iHS range: {min(ihs_scores):.2f} - {max(ihs_scores):.2f}")

    return positions, ihs_scores

def create_matplotlib_plots():
    """Create plots using matplotlib"""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec

        print("\nCreating 4-panel population comparison plot...")

        # Load all three populations
        eur_pos, eur_ihs = load_ihs_data(RESULTS_DIR / "ihs_EUR" / "EUR.chrX.ihs.tsv", "EUR")
        afr_pos, afr_ihs = load_ihs_data(RESULTS_DIR / "ihs_AFR" / "AFR.chrX.ihs.tsv", "AFR")
        eas_pos, eas_ihs = load_ihs_data(RESULTS_DIR / "ihs_EAS" / "EAS.chrX.ihs.tsv", "EAS")

        # Create figure with 4 panels
        fig = plt.figure(figsize=(16, 12))
        gs = GridSpec(4, 1, height_ratios=[1, 1, 1, 1], hspace=0.35)

        # Set style
        plt.rcParams['figure.dpi'] = 300
        plt.rcParams['savefig.dpi'] = 300
        plt.rcParams['font.size'] = 10

        # Panel A: EUR
        ax1 = fig.add_subplot(gs[0])
        ax1.scatter(eur_pos, eur_ihs, c='#1f77b4', s=0.5, alpha=0.6, label='EUR')
        ax1.axhline(y=2, color='gray', linestyle='--', linewidth=0.8, alpha=0.5, label='|iHS|=2')
        ax1.axhline(y=-2, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
        ax1.axhline(y=2.9, color='darkgray', linestyle='--', linewidth=0.8, alpha=0.7, label='|iHS|=2.9')
        ax1.axhline(y=-2.9, color='darkgray', linestyle='--', linewidth=0.8, alpha=0.7)
        ax1.axhline(y=0, color='black', linestyle='-', linewidth=0.5, alpha=0.3)
        ax1.set_ylabel('EUR iHS', fontsize=11, fontweight='bold')
        ax1.set_ylim(-7, 7)
        ax1.set_xlim(0, 160)
        ax1.set_xticks(range(0, 161, 20))
        ax1.legend(loc='upper right', fontsize=8)
        ax1.grid(True, alpha=0.3)

        # Panel B: AFR
        ax2 = fig.add_subplot(gs[1], sharex=ax1)
        ax2.scatter(afr_pos, afr_ihs, c='#ff7f0e', s=0.5, alpha=0.6, label='AFR')
        ax2.axhline(y=2, color='gray', linestyle='--', linewidth=0.8, alpha=0.5, label='|iHS|=2')
        ax2.axhline(y=-2, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
        ax2.axhline(y=2.9, color='darkgray', linestyle='--', linewidth=0.8, alpha=0.7, label='|iHS|=2.9')
        ax2.axhline(y=-2.9, color='darkgray', linestyle='--', linewidth=0.8, alpha=0.7)
        ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.5, alpha=0.3)
        ax2.set_ylabel('AFR iHS', fontsize=11, fontweight='bold')
        ax2.set_ylim(-7, 7)
        ax2.set_xticks(range(0, 161, 20))
        ax2.legend(loc='upper right', fontsize=8)
        ax2.grid(True, alpha=0.3)

        # Panel C: EAS
        ax3 = fig.add_subplot(gs[2], sharex=ax1)
        ax3.scatter(eas_pos, eas_ihs, c='#2ca02c', s=0.5, alpha=0.6, label='EAS')
        ax3.axhline(y=2, color='gray', linestyle='--', linewidth=0.8, alpha=0.5, label='|iHS|=2')
        ax3.axhline(y=-2, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
        ax3.axhline(y=2.9, color='darkgray', linestyle='--', linewidth=0.8, alpha=0.7, label='|iHS|=2.9')
        ax3.axhline(y=-2.9, color='darkgray', linestyle='--', linewidth=0.8, alpha=0.7)
        ax3.axhline(y=0, color='black', linestyle='-', linewidth=0.5, alpha=0.3)
        ax3.set_ylabel('EAS iHS', fontsize=11, fontweight='bold')
        ax3.set_ylim(-7, 7)
        ax3.set_xticks(range(0, 161, 20))
        ax3.legend(loc='upper right', fontsize=8)
        ax3.grid(True, alpha=0.3)

        # Panel D: All populations overlaid
        ax4 = fig.add_subplot(gs[3], sharex=ax1)
        ax4.scatter(eur_pos, eur_ihs, c='#1f77b4', s=0.3, alpha=0.4, label='EUR', rasterized=True)
        ax4.scatter(afr_pos, afr_ihs, c='#ff7f0e', s=0.3, alpha=0.4, label='AFR', rasterized=True)
        ax4.scatter(eas_pos, eas_ihs, c='#2ca02c', s=0.3, alpha=0.4, label='EAS', rasterized=True)
        ax4.axhline(y=2, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
        ax4.axhline(y=-2, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
        ax4.axhline(y=2.9, color='darkgray', linestyle='--', linewidth=0.8, alpha=0.7)
        ax4.axhline(y=-2.9, color='darkgray', linestyle='--', linewidth=0.8, alpha=0.7)
        ax4.axhline(y=0, color='black', linestyle='-', linewidth=0.5, alpha=0.3)
        ax4.set_ylabel('iHS (All Populations)', fontsize=11, fontweight='bold')
        ax4.set_xlabel('Position on Chromosome X (Mb)', fontsize=12, fontweight='bold')
        ax4.set_ylim(-7, 7)
        ax4.set_xlim(0, 160)
        ax4.set_xticks(range(0, 161, 20))
        ax4.legend(loc='upper right', fontsize=9)
        ax4.grid(True, alpha=0.3)

        # Overall title
        fig.suptitle('Population-Specific iHS Signals on Chromosome X',
                    fontsize=14, fontweight='bold', y=0.995)

        plt.tight_layout()

        # Save plot
        outfile = OUTPUT_DIR / 'population_ihs_chrX_comparison.png'
        plt.savefig(outfile, bbox_inches='tight', dpi=300)
        print(f"\n✓ Saved: {outfile}")
        plt.close()

        # Create a zoomed version focusing on strong signals
        print("\nCreating zoomed plot (strong signals only)...")
        fig = plt.figure(figsize=(16, 12))
        gs = GridSpec(4, 1, height_ratios=[1, 1, 1, 1], hspace=0.35)

        # Filter for strong signals (|iHS| > 2)
        eur_strong_pos = [p for p, s in zip(eur_pos, eur_ihs) if abs(s) > 2]
        eur_strong_ihs = [s for s in eur_ihs if abs(s) > 2]

        afr_strong_pos = [p for p, s in zip(afr_pos, afr_ihs) if abs(s) > 2]
        afr_strong_ihs = [s for s in afr_ihs if abs(s) > 2]

        eas_strong_pos = [p for p, s in zip(eas_pos, eas_ihs) if abs(s) > 2]
        eas_strong_ihs = [s for s in eas_ihs if abs(s) > 2]

        # Panel A: EUR strong signals
        ax1 = fig.add_subplot(gs[0])
        if len(eur_strong_pos) > 0:
            ax1.scatter(eur_strong_pos, eur_strong_ihs, c='#1f77b4', s=2, alpha=0.7, label=f'EUR (n={len(eur_strong_pos):,})')
        ax1.axhline(y=2.9, color='darkgray', linestyle='--', linewidth=1, alpha=0.7, label='|iHS|=2.9 (strong)')
        ax1.axhline(y=-2.9, color='darkgray', linestyle='--', linewidth=1, alpha=0.7)
        ax1.axhline(y=0, color='black', linestyle='-', linewidth=0.5, alpha=0.3)
        ax1.set_ylabel('EUR iHS', fontsize=11, fontweight='bold')
        ax1.set_ylim(-7, 7)
        ax1.set_xlim(0, 160)
        ax1.set_xticks(range(0, 161, 20))
        ax1.legend(loc='upper right', fontsize=8)
        ax1.grid(True, alpha=0.3)

        # Panel B: AFR strong signals
        ax2 = fig.add_subplot(gs[1], sharex=ax1)
        if len(afr_strong_pos) > 0:
            ax2.scatter(afr_strong_pos, afr_strong_ihs, c='#ff7f0e', s=2, alpha=0.7, label=f'AFR (n={len(afr_strong_pos):,})')
        ax2.axhline(y=2.9, color='darkgray', linestyle='--', linewidth=1, alpha=0.7, label='|iHS|=2.9 (strong)')
        ax2.axhline(y=-2.9, color='darkgray', linestyle='--', linewidth=1, alpha=0.7)
        ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.5, alpha=0.3)
        ax2.set_ylabel('AFR iHS', fontsize=11, fontweight='bold')
        ax2.set_ylim(-7, 7)
        ax2.set_xticks(range(0, 161, 20))
        ax2.legend(loc='upper right', fontsize=8)
        ax2.grid(True, alpha=0.3)

        # Panel C: EAS strong signals
        ax3 = fig.add_subplot(gs[2], sharex=ax1)
        if len(eas_strong_pos) > 0:
            ax3.scatter(eas_strong_pos, eas_strong_ihs, c='#2ca02c', s=2, alpha=0.7, label=f'EAS (n={len(eas_strong_pos):,})')
        ax3.axhline(y=2.9, color='darkgray', linestyle='--', linewidth=1, alpha=0.7, label='|iHS|=2.9 (strong)')
        ax3.axhline(y=-2.9, color='darkgray', linestyle='--', linewidth=1, alpha=0.7)
        ax3.axhline(y=0, color='black', linestyle='-', linewidth=0.5, alpha=0.3)
        ax3.set_ylabel('EAS iHS', fontsize=11, fontweight='bold')
        ax3.set_ylim(-7, 7)
        ax3.set_xticks(range(0, 161, 20))
        ax3.legend(loc='upper right', fontsize=8)
        ax3.grid(True, alpha=0.3)

        # Panel D: All strong signals overlaid
        ax4 = fig.add_subplot(gs[3], sharex=ax1)
        if len(eur_strong_pos) > 0:
            ax4.scatter(eur_strong_pos, eur_strong_ihs, c='#1f77b4', s=2, alpha=0.5, label=f'EUR (n={len(eur_strong_pos):,})', rasterized=True)
        if len(afr_strong_pos) > 0:
            ax4.scatter(afr_strong_pos, afr_strong_ihs, c='#ff7f0e', s=2, alpha=0.5, label=f'AFR (n={len(afr_strong_pos):,})', rasterized=True)
        if len(eas_strong_pos) > 0:
            ax4.scatter(eas_strong_pos, eas_strong_ihs, c='#2ca02c', s=2, alpha=0.5, label=f'EAS (n={len(eas_strong_pos):,})', rasterized=True)
        ax4.axhline(y=2.9, color='darkgray', linestyle='--', linewidth=1, alpha=0.7)
        ax4.axhline(y=-2.9, color='darkgray', linestyle='--', linewidth=1, alpha=0.7)
        ax4.axhline(y=0, color='black', linestyle='-', linewidth=0.5, alpha=0.3)
        ax4.set_ylabel('iHS (All Populations)', fontsize=11, fontweight='bold')
        ax4.set_xlabel('Position on Chromosome X (Mb)', fontsize=12, fontweight='bold')
        ax4.set_ylim(-7, 7)
        ax4.set_xlim(0, 160)
        ax4.set_xticks(range(0, 161, 20))
        ax4.legend(loc='upper right', fontsize=9)
        ax4.grid(True, alpha=0.3)

        fig.suptitle('Population-Specific iHS Signals on Chromosome X (|iHS| > 2 only)',
                    fontsize=14, fontweight='bold', y=0.995)

        plt.tight_layout()

        outfile = OUTPUT_DIR / 'population_ihs_chrX_comparison_strong_signals.png'
        plt.savefig(outfile, bbox_inches='tight', dpi=300)
        print(f"✓ Saved: {outfile}")
        plt.close()

        return True

    except ImportError as e:
        print(f"matplotlib not available: {e}")
        return False

def print_summary_statistics():
    """Print summary statistics for all populations"""
    print("\n" + "=" * 80)
    print("SUMMARY STATISTICS")
    print("=" * 80)

    for pop in ["EUR", "AFR", "EAS"]:
        filepath = RESULTS_DIR / f"ihs_{pop}" / f"{pop}.chrX.ihs.tsv"

        if not filepath.exists():
            print(f"\n{pop}: File not found")
            continue

        positions, ihs_scores = load_ihs_data(filepath, pop)

        if len(ihs_scores) == 0:
            continue

        # Count by threshold
        n_total = len(ihs_scores)
        n_above_2 = sum(1 for s in ihs_scores if abs(s) > 2.0)
        n_above_2p9 = sum(1 for s in ihs_scores if abs(s) >= 2.9)
        n_above_3 = sum(1 for s in ihs_scores if abs(s) > 3.0)

        print(f"\n{pop}:")
        print(f"  Total variants: {n_total:,}")
        print(f"  |iHS| > 2.0:    {n_above_2:,} ({n_above_2/n_total*100:.1f}%)")
        print(f"  |iHS| ≥ 2.9:    {n_above_2p9:,} ({n_above_2p9/n_total*100:.1f}%)")
        print(f"  |iHS| > 3.0:    {n_above_3:,} ({n_above_3/n_total*100:.1f}%)")
        print(f"  Range:          {min(ihs_scores):.2f} to {max(ihs_scores):.2f}")

def main():
    print("=" * 80)
    print("Population-Specific iHS Comparison Plot")
    print("=" * 80)

    # Print summary statistics
    print_summary_statistics()

    # Create plots
    print("\n" + "=" * 80)
    print("Creating Plots")
    print("=" * 80)

    success = create_matplotlib_plots()

    if success:
        print("\n" + "=" * 80)
        print("Plots created successfully!")
        print(f"Output directory: {OUTPUT_DIR}")
        print("=" * 80)
        print("\nFiles created:")
        print("  1. population_ihs_chrX_comparison.png - All variants")
        print("  2. population_ihs_chrX_comparison_strong_signals.png - Strong signals only (|iHS| > 2)")
    else:
        print("\nNote: matplotlib not available for plotting")

if __name__ == "__main__":
    main()
