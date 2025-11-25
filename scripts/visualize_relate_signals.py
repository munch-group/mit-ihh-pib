#!/usr/bin/env python3
"""
Visualization of RELATE selection signals along the X chromosome
Shows distribution of significant SNPs similar to visualize_highsig_geneinfo.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

sns.set_style('white')

# High-quality figures
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

def load_relate_data():
    """Load filtered RELATE data - filtered for negLog10p-values ≥ 6"""
    data_file = '/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/analysis/ihs_relate_comparison/relate_filtering/chrX_signals_log10p_ge_6.csv'
    df = pd.read_csv(data_file)

    # Ensure position is numeric
    df['Position'] = pd.to_numeric(df['Position'], errors='coerce')

    return df

def create_simple_manhattan():
    """
    Simple Manhattan plot without gene context - overview of RELATE signals
    """
    from geneinfo.plot import ChromIdeogram

    df = load_relate_data()

    # Get unique populations
    populations = df['Population'].unique()
    print(f"Populations in dataset: {', '.join(populations)}")
    print(f"Total SNPs: {len(df)}")

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
    colors = plt.cm.Reds(plt.Normalize(
        vmin=6,  # Start from threshold
        vmax=df['minus_log10_p'].max()
    )(df['minus_log10_p']))

    # Main scatter - color by population if multiple, otherwise by significance
    if len(populations) > 1:
        # Create color map for populations
        pop_colors = plt.cm.tab10(np.linspace(0, 1, len(populations)))
        pop_color_map = {pop: pop_colors[i] for i, pop in enumerate(populations)}
        point_colors = [pop_color_map[pop] for pop in df['Population']]
    else:
        point_colors = colors

    scatter = ax.scatter(
        df['Position'] / 1e6,  # Convert to Mb
        df['minus_log10_p'],
        c=point_colors,
        s=df['minus_log10_p'] * 12,
        alpha=0.7,
        edgecolors='black',
        linewidths=0.5
    )

    # Threshold line
    ax.axhline(y=6, color='red', linestyle='--', linewidth=1.5, alpha=0.5,
               label='Threshold: p = 10⁻⁶')

    # Centromere region (approximate for hg38)
    ax.axvspan(58.1, 63.8, alpha=0.1, color='gray', label='Centromere')

    # Add population legend if multiple populations
    if len(populations) > 1:
        # Create custom legend for populations
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor=pop_color_map[pop],
                                edgecolor='black',
                                label=pop)
                          for pop in populations]
        legend_elements.append(plt.Line2D([0], [0], color='red', linestyle='--',
                                         linewidth=1.5, label='Threshold: p = 10⁻⁶'))
        legend_elements.append(Patch(facecolor='gray', alpha=0.1,
                                    edgecolor='none', label='Centromere'))
        ax.legend(handles=legend_elements, loc='upper right', frameon=True,
                 fontsize=9, ncol=2 if len(populations) > 5 else 1)

    ax.set_xlabel('Position on X Chromosome (Mb)', fontsize=12, weight='bold')
    ax.set_ylabel('-log₁₀(p-value)', fontsize=12, weight='bold')

    title = f'Distribution of {len(df)} RELATE Selection Signals on X Chromosome'
    if len(populations) > 1:
        title += f'\n({len(populations)} populations)'
    ax.set_title(title, fontsize=13, weight='bold', pad=15)

    ax.set_xlim(0, 156)
    ax.set_ylim(5.8, df['minus_log10_p'].max() + 0.5)

    if len(populations) == 1:
        ax.legend(loc='upper right', frameon=True, fontsize=10)

    ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)

    # Add colorbar for significance (if single population)
    if len(populations) == 1:
        sm = plt.cm.ScalarMappable(
            cmap=plt.cm.Reds,
            norm=plt.Normalize(vmin=6, vmax=df['minus_log10_p'].max())
        )
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, pad=0.02, aspect=30, shrink=0.8)
        cbar.set_label('-log₁₀(p-value)', rotation=270, labelpad=20,
                      weight='bold', fontsize=10)

    sns.despine(ax=ax)
    plt.tight_layout()

    output_path = '/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/analysis/ihs_relate_comparison/relate_filtering/relate_manhattan_chrX.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"✓ Saved Manhattan plot: {output_path}")
    plt.close()

def create_full_chromosome_view():
    """
    Create full X chromosome view with all SNPs and density plot
    """
    from geneinfo.plot import ChromIdeogram

    df = load_relate_data()

    # Create ideogram
    g = ChromIdeogram('chrX', assembly='hg38', font_size=8)

    # Draw chromosome at specific position
    g.draw_chromosomes(base=3, height=1, facecolor='lightgray')

    # Create scatter plot axes above chromosome
    ax1, ax2 = g.add_axes(2, height_ratio=0.8, hspace=0.3)

    # Plot 1: All SNPs colored by significance
    colors = plt.cm.Reds(plt.Normalize(
        vmin=df['minus_log10_p'].min(),
        vmax=df['minus_log10_p'].max()
    )(df['minus_log10_p']))

    ax1.scatter(
        df['Position'],
        df['minus_log10_p'],
        c=colors,
        s=df['minus_log10_p'] * 8,
        alpha=0.7,
        edgecolors='black',
        linewidths=0.3
    )
    ax1.axhline(y=6, color='red', linestyle='--', linewidth=1, alpha=0.5, label='p = 10⁻⁶')
    ax1.set_ylabel('-log₁₀(p-value)', fontsize=10, weight='bold')
    ax1.set_title('RELATE Selection Signals on X Chromosome', fontsize=12, weight='bold')
    ax1.legend(loc='upper right', frameon=True, fontsize=9)
    ax1.grid(True, alpha=0.3, linestyle=':')

    # Plot 2: Density/coverage
    bins = 100
    ax2.hist(df['Position'], bins=bins, color='darkred', alpha=0.6, edgecolor='black')
    ax2.set_ylabel('SNP Count', fontsize=10, weight='bold')
    ax2.set_title('Selection Signal Distribution Density', fontsize=11, weight='bold')
    ax2.grid(True, alpha=0.3, linestyle=':', axis='y')

    # Clean up
    sns.despine(ax=ax1)
    sns.despine(ax=ax2)

    plt.tight_layout()

    output_path = '/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/analysis/ihs_relate_comparison/relate_filtering/relate_full_chromosome.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"✓ Saved full chromosome view: {output_path}")
    plt.close()

def create_population_comparison():
    """
    Create plots comparing signals across populations
    """
    df = load_relate_data()
    populations = df['Population'].unique()

    if len(populations) == 1:
        print("Only one population - skipping population comparison plot")
        return

    # Create figure with subplots for each population
    n_pops = len(populations)
    fig, axes = plt.subplots(n_pops, 1, figsize=(14, 4*n_pops), sharex=True)

    if n_pops == 1:
        axes = [axes]

    for idx, pop in enumerate(sorted(populations)):
        ax = axes[idx]
        pop_df = df[df['Population'] == pop]

        # Color by significance
        colors = plt.cm.Reds(plt.Normalize(
            vmin=6,
            vmax=pop_df['minus_log10_p'].max()
        )(pop_df['minus_log10_p']))

        ax.scatter(
            pop_df['Position'] / 1e6,
            pop_df['minus_log10_p'],
            c=colors,
            s=pop_df['minus_log10_p'] * 10,
            alpha=0.7,
            edgecolors='black',
            linewidths=0.4
        )

        ax.axhline(y=6, color='red', linestyle='--', linewidth=1, alpha=0.5)
        ax.axvspan(58.1, 63.8, alpha=0.1, color='gray')  # Centromere

        ax.set_ylabel('-log₁₀(p-value)', fontsize=10, weight='bold')
        ax.set_title(f'{pop} ({len(pop_df)} signals)', fontsize=11, weight='bold')
        ax.grid(True, alpha=0.3, linestyle=':')
        ax.set_xlim(0, 156)

        sns.despine(ax=ax)

    axes[-1].set_xlabel('Position on X Chromosome (Mb)', fontsize=11, weight='bold')

    plt.suptitle('RELATE Selection Signals by Population',
                fontsize=14, weight='bold', y=0.995)
    plt.tight_layout()

    output_path = '/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/analysis/ihs_relate_comparison/relate_filtering/relate_by_population.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"✓ Saved population comparison: {output_path}")
    plt.close()

def create_summary_stats():
    """
    Create summary statistics figure
    """
    df = load_relate_data()

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # 1. Distribution of -log10(p-values)
    ax = axes[0, 0]
    ax.hist(df['minus_log10_p'], bins=50, color='darkred', alpha=0.7, edgecolor='black')
    ax.axvline(x=df['minus_log10_p'].mean(), color='blue', linestyle='--',
               linewidth=2, label=f'Mean: {df["minus_log10_p"].mean():.2f}')
    ax.axvline(x=df['minus_log10_p'].median(), color='green', linestyle='--',
               linewidth=2, label=f'Median: {df["minus_log10_p"].median():.2f}')
    ax.set_xlabel('-log₁₀(p-value)', fontsize=10, weight='bold')
    ax.set_ylabel('Count', fontsize=10, weight='bold')
    ax.set_title('Distribution of Signal Strengths', fontsize=11, weight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3, linestyle=':')
    sns.despine(ax=ax)

    # 2. Signals per population
    ax = axes[0, 1]
    pop_counts = df['Population'].value_counts().sort_index()
    ax.bar(range(len(pop_counts)), pop_counts.values, color='steelblue',
           alpha=0.7, edgecolor='black')
    ax.set_xticks(range(len(pop_counts)))
    ax.set_xticklabels(pop_counts.index, rotation=45, ha='right')
    ax.set_ylabel('Number of Signals', fontsize=10, weight='bold')
    ax.set_title('Signals per Population', fontsize=11, weight='bold')
    ax.grid(True, alpha=0.3, linestyle=':', axis='y')
    sns.despine(ax=ax)

    # 3. Position distribution along chromosome
    ax = axes[1, 0]
    ax.hist(df['Position'] / 1e6, bins=50, color='purple', alpha=0.7, edgecolor='black')
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
    top10 = df.nlargest(10, 'minus_log10_p')[['Population', 'Position', 'minus_log10_p']]
    top10['Position'] = (top10['Position'] / 1e6).round(2)
    top10['minus_log10_p'] = top10['minus_log10_p'].round(2)
    top10.columns = ['Pop', 'Pos (Mb)', '-log₁₀(p)']

    table = ax.table(cellText=top10.values, colLabels=top10.columns,
                    cellLoc='center', loc='center',
                    colWidths=[0.3, 0.35, 0.35])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2)

    # Style header
    for i in range(len(top10.columns)):
        cell = table[(0, i)]
        cell.set_facecolor('darkred')
        cell.set_text_props(weight='bold', color='white')

    ax.set_title('Top 10 Strongest Signals', fontsize=11, weight='bold', pad=20)

    plt.suptitle(f'RELATE Selection Signals Summary (n={len(df)})',
                fontsize=14, weight='bold')
    plt.tight_layout()

    output_path = '/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/analysis/ihs_relate_comparison/relate_filtering/relate_summary_stats.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"✓ Saved summary statistics: {output_path}")
    plt.close()

def main():
    """Main execution"""
    print("\n" + "="*60)
    print("Creating RELATE Signal Visualizations")
    print("="*60 + "\n")

    # Create output directory (already exists, but just in case)
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

    print("\n3. Creating population comparison...")
    try:
        create_population_comparison()
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
