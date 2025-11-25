#!/usr/bin/env python3
"""
Visualization of iHS selection signals along the X chromosome
Shows distribution of significant SNPs from iHS analysis
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('white')

# High-quality figures
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

def load_ihs_data():
    """Load iHS selection candidates"""
    data_file = '/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/analysis/candidates/per_chromosome/chrX_candidates_liberal.tsv'
    df = pd.read_csv(data_file, sep='\t')

    # Ensure position is numeric
    df['pos'] = pd.to_numeric(df['pos'], errors='coerce')

    return df

def create_simple_manhattan():
    """
    Simple Manhattan plot - overview of iHS signals
    """
    from geneinfo.plot import ChromIdeogram

    df = load_ihs_data()

    print(f"Total iHS candidates: {len(df)}")
    print(f"Position range: {df['pos'].min()/1e6:.2f} - {df['pos'].max()/1e6:.2f} Mb")
    print(f"-log10(p) range: {df['neg_log10_p'].min():.2f} - {df['neg_log10_p'].max():.2f}")

    # Temporarily reduce DPI to avoid huge figure with ChromIdeogram
    original_dpi = plt.rcParams['figure.dpi']
    plt.rcParams['figure.dpi'] = 100

    # Create ChromIdeogram and draw chromosome
    g = ChromIdeogram('chrX', assembly='hg38')
    g.draw_chromosomes(base=9, height=0.5, facecolor='lightgray')

    # Add axis for the Manhattan plot
    ax = g.add_axes(1, height_ratio=0.85)

    # Restore original DPI
    plt.rcParams['figure.dpi'] = original_dpi

    # Color by significance
    colors = plt.cm.Blues(plt.Normalize(
        vmin=df['neg_log10_p'].min(),
        vmax=df['neg_log10_p'].max()
    )(df['neg_log10_p']))

    # Main scatter
    scatter = ax.scatter(
        df['pos'] / 1e6,  # Convert to Mb
        df['neg_log10_p'],
        c=colors,
        s=df['neg_log10_p'] * 12,
        alpha=0.7,
        edgecolors='black',
        linewidths=0.5
    )

    # Threshold line (if there's a clear threshold)
    threshold = df['neg_log10_p'].min()
    ax.axhline(y=threshold, color='blue', linestyle='--', linewidth=1.5, alpha=0.5,
               label=f'Threshold: -log₁₀(p) = {threshold:.2f}')

    # Centromere region (approximate for hg38)
    ax.axvspan(58.1, 63.8, alpha=0.1, color='gray', label='Centromere')

    # PAR1 region
    ax.axvspan(0, 2.7, alpha=0.15, color='green', label='PAR1')

    ax.set_xlabel('Position on X Chromosome (Mb)', fontsize=12, weight='bold')
    ax.set_ylabel('-log₁₀(p-value)', fontsize=12, weight='bold')

    title = f'Distribution of {len(df)} iHS Selection Candidates on X Chromosome'
    ax.set_title(title, fontsize=13, weight='bold', pad=15)

    ax.set_xlim(0, 156)
    ax.set_ylim(df['neg_log10_p'].min() - 0.2, df['neg_log10_p'].max() + 0.5)

    ax.legend(loc='upper right', frameon=True, fontsize=10)
    ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)

    # Add colorbar for significance
    sm = plt.cm.ScalarMappable(
        cmap=plt.cm.Blues,
        norm=plt.Normalize(vmin=df['neg_log10_p'].min(), vmax=df['neg_log10_p'].max())
    )
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, pad=0.02, aspect=30, shrink=0.8)
    cbar.set_label('-log₁₀(p-value)', rotation=270, labelpad=20,
                  weight='bold', fontsize=10)

    sns.despine(ax=ax)
    plt.tight_layout()

    output_path = '/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/analysis/ihs_relate_comparison/relate_filtering/ihs_manhattan_chrX.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"✓ Saved Manhattan plot: {output_path}")
    plt.close()

def create_full_chromosome_view():
    """
    Create full X chromosome view with all SNPs and density plot
    """
    from geneinfo.plot import ChromIdeogram

    df = load_ihs_data()

    # Create ideogram
    g = ChromIdeogram('chrX', assembly='hg38', font_size=8)

    # Draw chromosome at specific position
    g.draw_chromosomes(base=3, height=1, facecolor='lightgray')

    # Create scatter plot axes above chromosome
    ax1, ax2 = g.add_axes(2, height_ratio=0.8, hspace=0.3)

    # Plot 1: All SNPs colored by significance
    colors = plt.cm.Blues(plt.Normalize(
        vmin=df['neg_log10_p'].min(),
        vmax=df['neg_log10_p'].max()
    )(df['neg_log10_p']))

    ax1.scatter(
        df['pos'],
        df['neg_log10_p'],
        c=colors,
        s=df['neg_log10_p'] * 8,
        alpha=0.7,
        edgecolors='black',
        linewidths=0.3
    )
    threshold = df['neg_log10_p'].min()
    ax1.axhline(y=threshold, color='blue', linestyle='--', linewidth=1, alpha=0.5,
               label=f'Threshold: {threshold:.2f}')
    ax1.set_ylabel('-log₁₀(p-value)', fontsize=10, weight='bold')
    ax1.set_title('iHS Selection Signals on X Chromosome', fontsize=12, weight='bold')
    ax1.legend(loc='upper right', frameon=True, fontsize=9)
    ax1.grid(True, alpha=0.3, linestyle=':')

    # Plot 2: Density/coverage
    bins = 100
    ax2.hist(df['pos'], bins=bins, color='darkblue', alpha=0.6, edgecolor='black')
    ax2.set_ylabel('SNP Count', fontsize=10, weight='bold')
    ax2.set_title('Selection Signal Distribution Density', fontsize=11, weight='bold')
    ax2.grid(True, alpha=0.3, linestyle=':', axis='y')

    # Clean up
    sns.despine(ax=ax1)
    sns.despine(ax=ax2)

    plt.tight_layout()

    output_path = '/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/analysis/ihs_relate_comparison/relate_filtering/ihs_full_chromosome.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"✓ Saved full chromosome view: {output_path}")
    plt.close()

def create_comparison_binned():
    """
    Create binned comparison showing signal density across chromosome
    """
    df = load_ihs_data()

    # Create bins along the chromosome (10 Mb bins)
    bin_size = 10e6  # 10 Mb
    max_pos = 156e6
    bins = np.arange(0, max_pos + bin_size, bin_size)
    bin_centers = (bins[:-1] + bins[1:]) / 2

    # Count signals per bin
    counts, _ = np.histogram(df['pos'], bins=bins)

    # Average significance per bin
    bin_labels = pd.cut(df['pos'], bins=bins, labels=bin_centers/1e6)
    avg_sig = df.groupby(bin_labels)['neg_log10_p'].mean()

    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

    # Plot 1: Signal count per bin
    ax = axes[0]
    ax.bar(bin_centers/1e6, counts, width=8, color='steelblue', alpha=0.7, edgecolor='black')
    ax.axvspan(0, 2.7, alpha=0.15, color='green', label='PAR1')
    ax.axvspan(58.1, 63.8, alpha=0.1, color='gray', label='Centromere')
    ax.set_ylabel('Number of Candidates', fontsize=11, weight='bold')
    ax.set_title('iHS Selection Candidates Distribution (10 Mb bins)', fontsize=12, weight='bold')
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3, linestyle=':', axis='y')
    sns.despine(ax=ax)

    # Plot 2: Average significance per bin
    ax = axes[1]
    ax.bar(bin_centers/1e6, avg_sig.values, width=8, color='darkblue', alpha=0.7, edgecolor='black')
    ax.axvspan(0, 2.7, alpha=0.15, color='green')
    ax.axvspan(58.1, 63.8, alpha=0.1, color='gray')
    ax.set_xlabel('Position on X Chromosome (Mb)', fontsize=11, weight='bold')
    ax.set_ylabel('Mean -log₁₀(p-value)', fontsize=11, weight='bold')
    ax.set_title('Average Selection Signal Strength per Region', fontsize=12, weight='bold')
    ax.grid(True, alpha=0.3, linestyle=':', axis='y')
    ax.set_xlim(0, 156)
    sns.despine(ax=ax)

    plt.tight_layout()

    output_path = '/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/analysis/ihs_relate_comparison/relate_filtering/ihs_binned_distribution.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"✓ Saved binned distribution: {output_path}")
    plt.close()

def create_summary_stats():
    """
    Create summary statistics figure
    """
    df = load_ihs_data()

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # 1. Distribution of -log10(p-values)
    ax = axes[0, 0]
    ax.hist(df['neg_log10_p'], bins=50, color='darkblue', alpha=0.7, edgecolor='black')
    ax.axvline(x=df['neg_log10_p'].mean(), color='red', linestyle='--',
               linewidth=2, label=f'Mean: {df["neg_log10_p"].mean():.2f}')
    ax.axvline(x=df['neg_log10_p'].median(), color='green', linestyle='--',
               linewidth=2, label=f'Median: {df["neg_log10_p"].median():.2f}')
    ax.set_xlabel('-log₁₀(p-value)', fontsize=10, weight='bold')
    ax.set_ylabel('Count', fontsize=10, weight='bold')
    ax.set_title('Distribution of Signal Strengths', fontsize=11, weight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3, linestyle=':')
    sns.despine(ax=ax)

    # 2. Distribution of |iHS| scores
    ax = axes[0, 1]
    ax.hist(df['Std iHS'].abs(), bins=50, color='purple', alpha=0.7, edgecolor='black')
    ax.axvline(x=df['Std iHS'].abs().mean(), color='red', linestyle='--',
               linewidth=2, label=f'Mean: {df["Std iHS"].abs().mean():.2f}')
    ax.set_xlabel('|Standardized iHS|', fontsize=10, weight='bold')
    ax.set_ylabel('Count', fontsize=10, weight='bold')
    ax.set_title('Distribution of |iHS| Scores', fontsize=11, weight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3, linestyle=':')
    sns.despine(ax=ax)

    # 3. Position distribution along chromosome
    ax = axes[1, 0]
    ax.hist(df['pos'] / 1e6, bins=50, color='steelblue', alpha=0.7, edgecolor='black')
    ax.axvspan(0, 2.7, alpha=0.2, color='green', label='PAR1')
    ax.axvspan(58.1, 63.8, alpha=0.2, color='gray', label='Centromere')
    ax.set_xlabel('Position (Mb)', fontsize=10, weight='bold')
    ax.set_ylabel('Count', fontsize=10, weight='bold')
    ax.set_title('Positional Distribution', fontsize=11, weight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3, linestyle=':')
    sns.despine(ax=ax)

    # 4. Top signals table
    ax = axes[1, 1]
    ax.axis('off')
    top10 = df.nlargest(10, 'neg_log10_p')[['chr', 'pos', 'Std iHS', 'neg_log10_p']]
    top10['pos'] = (top10['pos'] / 1e6).round(2)
    top10['Std iHS'] = top10['Std iHS'].round(2)
    top10['neg_log10_p'] = top10['neg_log10_p'].round(2)
    top10.columns = ['Chr', 'Pos (Mb)', 'Std iHS', '-log₁₀(p)']

    table = ax.table(cellText=top10.values, colLabels=top10.columns,
                    cellLoc='center', loc='center',
                    colWidths=[0.2, 0.3, 0.25, 0.25])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2)

    # Style header
    for i in range(len(top10.columns)):
        cell = table[(0, i)]
        cell.set_facecolor('darkblue')
        cell.set_text_props(weight='bold', color='white')

    ax.set_title('Top 10 Strongest Signals', fontsize=11, weight='bold', pad=20)

    plt.suptitle(f'iHS Selection Candidates Summary (n={len(df)})',
                fontsize=14, weight='bold')
    plt.tight_layout()

    output_path = '/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/analysis/ihs_relate_comparison/relate_filtering/ihs_summary_stats.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"✓ Saved summary statistics: {output_path}")
    plt.close()

def main():
    """Main execution"""
    print("\n" + "="*60)
    print("Creating iHS Signal Visualizations")
    print("="*60 + "\n")

    # Create output directory
    import os
    output_dir = '/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/analysis/ihs_relate_comparison/relate_filtering'
    os.makedirs(output_dir, exist_ok=True)

    print("1. Creating simple Manhattan plot...")
    try:
        create_simple_manhattan()
    except Exception as e:
        print(f"  Error: {e}")

    print("\n2. Creating full chromosome view...")
    try:
        create_full_chromosome_view()
    except Exception as e:
        print(f"  Error: {e}")

    print("\n3. Creating binned distribution...")
    try:
        create_comparison_binned()
    except Exception as e:
        print(f"  Error: {e}")

    print("\n4. Creating summary statistics...")
    try:
        create_summary_stats()
    except Exception as e:
        print(f"  Error: {e}")

    print("\n" + "="*60)
    print("✓ Visualization complete!")
    print("="*60 + "\n")

if __name__ == '__main__':
    main()
