#!/usr/bin/env python3
"""
Plot only EUR iHS results for chromosome X
Simple single-panel plot for comparison with EUR selection analysis
"""

from pathlib import Path

# Directories
RESULTS_DIR = Path("/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results")
PLOTS_DIR = Path("/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/plots")
OUTPUT_DIR = PLOTS_DIR / "EUR_ihs_only"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def parse_position(variant_id):
    """Extract position from variant ID (format: X:position:ref:alt)"""
    try:
        parts = variant_id.split(':')
        return int(parts[1])
    except:
        return None

def load_eur_ihs_data():
    """Load EUR iHS results"""
    filepath = RESULTS_DIR / "ihs_EUR" / "EUR.chrX.ihs.tsv"
    print(f"Loading EUR iHS data from {filepath.name}...")

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

def create_matplotlib_plot():
    """Create simple EUR iHS plot"""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        print("\nCreating EUR iHS plot...")

        # Load EUR data
        eur_pos, eur_ihs = load_eur_ihs_data()

        # Create figure
        fig, ax = plt.subplots(figsize=(16, 6))

        # Set style
        plt.rcParams['figure.dpi'] = 300
        plt.rcParams['savefig.dpi'] = 300
        plt.rcParams['font.size'] = 11

        # Plot EUR iHS
        ax.scatter(eur_pos, eur_ihs, c='#1f77b4', s=0.5, alpha=0.6, rasterized=True)

        # Add threshold lines
        ax.axhline(y=2, color='gray', linestyle='--', linewidth=1, alpha=0.6, label='|iHS|=2 (suggestive)')
        ax.axhline(y=-2, color='gray', linestyle='--', linewidth=1, alpha=0.6)
        ax.axhline(y=2.9, color='darkgray', linestyle='--', linewidth=1, alpha=0.8, label='|iHS|=2.9 (strong)')
        ax.axhline(y=-2.9, color='darkgray', linestyle='--', linewidth=1, alpha=0.8)
        ax.axhline(y=0, color='black', linestyle='-', linewidth=0.6, alpha=0.4)

        # Labels and formatting
        ax.set_xlabel('Position on Chromosome X (Mb)', fontsize=13, fontweight='bold')
        ax.set_ylabel('EUR iHS', fontsize=13, fontweight='bold')
        ax.set_title('EUR Population-Specific iHS on Chromosome X', fontsize=15, fontweight='bold')
        ax.set_xlim(0, 160)
        ax.set_ylim(-7, 7)
        ax.set_xticks(range(0, 161, 20))
        ax.legend(loc='upper right', fontsize=10)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        # Save plot
        outfile = OUTPUT_DIR / 'EUR_ihs_chrX.png'
        plt.savefig(outfile, bbox_inches='tight', dpi=300)
        print(f"\n✓ Saved: {outfile}")
        plt.close()

        # Also create a version showing only strong signals (|iHS| > 2)
        print("\nCreating EUR iHS plot (strong signals only)...")

        eur_strong_pos = [p for p, s in zip(eur_pos, eur_ihs) if abs(s) > 2]
        eur_strong_ihs = [s for s in eur_ihs if abs(s) > 2]

        fig, ax = plt.subplots(figsize=(16, 6))

        # Plot strong signals only
        ax.scatter(eur_strong_pos, eur_strong_ihs, c='#1f77b4', s=2, alpha=0.7, rasterized=True,
                  label=f'EUR (n={len(eur_strong_pos):,} variants with |iHS|>2)')

        # Add threshold lines
        ax.axhline(y=2.9, color='darkgray', linestyle='--', linewidth=1.2, alpha=0.8, label='|iHS|=2.9 (strong)')
        ax.axhline(y=-2.9, color='darkgray', linestyle='--', linewidth=1.2, alpha=0.8)
        ax.axhline(y=0, color='black', linestyle='-', linewidth=0.6, alpha=0.4)

        # Labels and formatting
        ax.set_xlabel('Position on Chromosome X (Mb)', fontsize=13, fontweight='bold')
        ax.set_ylabel('EUR iHS', fontsize=13, fontweight='bold')
        ax.set_title('EUR Population-Specific iHS on Chromosome X (|iHS| > 2 only)',
                    fontsize=15, fontweight='bold')
        ax.set_xlim(0, 160)
        ax.set_ylim(-7, 7)
        ax.set_xticks(range(0, 161, 20))
        ax.legend(loc='upper right', fontsize=10)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        outfile = OUTPUT_DIR / 'EUR_ihs_chrX_strong_signals.png'
        plt.savefig(outfile, bbox_inches='tight', dpi=300)
        print(f"✓ Saved: {outfile}")
        plt.close()

        return True

    except ImportError as e:
        print(f"matplotlib not available: {e}")
        return False

def print_summary_statistics():
    """Print summary statistics for EUR"""
    print("\n" + "=" * 80)
    print("EUR iHS SUMMARY STATISTICS")
    print("=" * 80)

    positions, ihs_scores = load_eur_ihs_data()

    if len(ihs_scores) == 0:
        print("No data loaded")
        return

    # Count by threshold
    n_total = len(ihs_scores)
    n_above_2 = sum(1 for s in ihs_scores if abs(s) > 2.0)
    n_above_2p9 = sum(1 for s in ihs_scores if abs(s) >= 2.9)
    n_above_3 = sum(1 for s in ihs_scores if abs(s) > 3.0)

    # Count positive vs negative
    n_positive = sum(1 for s in ihs_scores if s > 2.0)
    n_negative = sum(1 for s in ihs_scores if s < -2.0)

    print(f"\nTotal variants: {n_total:,}")
    print(f"\nVariants by threshold:")
    print(f"  |iHS| > 2.0:    {n_above_2:,} ({n_above_2/n_total*100:.2f}%)")
    print(f"  |iHS| ≥ 2.9:    {n_above_2p9:,} ({n_above_2p9/n_total*100:.2f}%)")
    print(f"  |iHS| > 3.0:    {n_above_3:,} ({n_above_3/n_total*100:.2f}%)")

    print(f"\nDirection of selection (|iHS| > 2):")
    print(f"  Positive (derived allele):   {n_positive:,} ({n_positive/n_above_2*100:.1f}% of significant)")
    print(f"  Negative (ancestral allele): {n_negative:,} ({n_negative/n_above_2*100:.1f}% of significant)")

    print(f"\niHS score range: {min(ihs_scores):.2f} to {max(ihs_scores):.2f}")

def main():
    print("=" * 80)
    print("EUR-Only iHS Plot")
    print("=" * 80)

    # Print summary statistics
    print_summary_statistics()

    # Create plots
    print("\n" + "=" * 80)
    print("Creating Plots")
    print("=" * 80)

    success = create_matplotlib_plot()

    if success:
        print("\n" + "=" * 80)
        print("Plots created successfully!")
        print(f"Output directory: {OUTPUT_DIR}")
        print("=" * 80)
        print("\nFiles created:")
        print("  1. EUR_ihs_chrX.png - All variants")
        print("  2. EUR_ihs_chrX_strong_signals.png - Strong signals only (|iHS| > 2)")
        print("\nCompare with:")
        print("  plots/EUR_selection_analysis/EUR_selection_genome_wide_overview.png")
    else:
        print("\nNote: matplotlib not available for plotting")

if __name__ == "__main__":
    main()
