#!/usr/bin/env python3
"""
Create X chromosome ideogram visualizations for selection analysis

Creates three visualizations:
1. X chromosome with 43 high-significance SNPs
2. X chromosome with Phase 3 selection regions
3. Combined view with SNPs, genes, and functional categories

Author: MIT van Bruggen
Date: 2025-10-30
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from pathlib import Path

# Import geneinfo
import geneinfo.plot as gplt
import geneinfo.information as gi
import geneinfo.utils as utils

# Configure plotting
sns.set_style('white')
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

# Setup paths
PROJECT_DIR = Path('/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib')
RESULTS_DIR = PROJECT_DIR / 'results/analysis'
FIGURES_DIR = PROJECT_DIR / 'results/figures'
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Input files
HIGH_SIG_SNPS = RESULTS_DIR / 'candidates/chrX_high_significance.tsv'
SELECTION_REGIONS = RESULTS_DIR / 'regions/selection_regions_500kb_min2snps.tsv'
GENE_DETAILS = RESULTS_DIR / 'high_significance_genes/gene_details_complete.tsv'
SNP_TO_GENE = RESULTS_DIR / 'high_significance_genes/snp_to_nearest_gene_only.tsv'

print("="*70)
print("X Chromosome Ideogram Visualizations")
print("="*70)
print()

#==============================================================================
# 1. X Chromosome Ideogram with 43 High-Significance SNPs
#==============================================================================

def create_snp_ideogram():
    """Create ideogram showing high-significance SNPs on X chromosome"""
    print("Creating visualization 1: High-Significance SNPs")

    # Load SNP data
    snps = pd.read_csv(HIGH_SIG_SNPS, sep='\t')
    snp_mapping = pd.read_csv(SNP_TO_GENE, sep='\t')

    print(f"  Loaded {len(snps)} high-significance SNPs")

    # Label the 11 overlap genes (highest confidence)
    overlap_genes = ['DACH2', 'DIAPH2', 'FAAH2', 'GPC3', 'HTR2C',
                     'IL1RAPL2', 'KDM6A', 'NAP1L2', 'NCBP2L', 'TRPC5', 'ZMAT1']

    # Get gene positions for labeling on chromosome
    gene_labels_data = []
    snps_sorted = snps.sort_values('neg_log10_p', ascending=False)
    labeled_genes = set()

    for _, snp in snps_sorted.iterrows():
        gene = snp_mapping[snp_mapping['snp_location'] == snp['Location']]['nearest_gene'].values
        gene_name = gene[0] if len(gene) > 0 else ''

        if gene_name in overlap_genes and gene_name not in labeled_genes:
            gene_labels_data.append(('chrX', snp['pos'], gene_name, 'darkred'))
            labeled_genes.add(gene_name)

    # Create ideogram with labels on chromosome
    g = gplt.ChromIdeogram('chrX', assembly='hg38', font_size=9)
    g.draw_chromosomes(base=2, facecolor='lightgray', height=0.8)

    # Add gene labels on the chromosome
    if gene_labels_data:
        g.add_labels(gene_labels_data, base=g.ideogram_base,
                    min_height=g.ideogram_height*2.5, bold=overlap_genes)

    # Add three axes below the chromosome
    new_ax, new_ax2, new_ax3 = g.add_axes(3, height_ratio=0.5, hspace=0.5)

    import matplotlib.ticker as ticker

    #--------------------------------------------------------------------------
    # First axis: Manhattan-style scatter plot
    #--------------------------------------------------------------------------
    # Create color scale
    cmap = plt.get_cmap('Reds')
    norm = mcolors.Normalize(vmin=6, vmax=snps['neg_log10_p'].max())
    sizes = ((snps['neg_log10_p'] - 6) / (snps['neg_log10_p'].max() - 6)) * 100 + 20
    colors = [cmap(norm(p)) for p in snps['neg_log10_p']]

    new_ax.scatter(snps['pos'], snps['neg_log10_p'],
                  s=sizes, c=colors, alpha=0.8, edgecolors='black', linewidth=0.5,
                  label='High-significance SNPs', zorder=3)
    new_ax.axhline(y=6, color='red', linestyle='--', linewidth=1.5, alpha=0.7,
                   label='Threshold: p = 10⁻⁶', zorder=1)

    new_ax.set_ylabel('-log₁₀(p-value)', fontsize=11, fontweight='bold')
    new_ax.set_xlabel('')
    new_ax.set_ylim(5.5, snps['neg_log10_p'].max() + 0.5)
    new_ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5),
                 fontsize=9, frameon=False)
    new_ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f'{int(x/1e6)}'))
    new_ax.tick_params(axis='both', which='major', labelsize=10)
    new_ax.grid(True, alpha=0.2)

    #--------------------------------------------------------------------------
    # Second axis: Smoothed trend line
    #--------------------------------------------------------------------------
    new_ax2.scatter(snps['pos'], snps['neg_log10_p'],
                   s=15, c='black', alpha=0.5, label='Individual SNPs')

    # Use lowess smoothing from utils if available, otherwise use simple smoothing
    try:
        fit = utils.fit_lowess(snps.assign(x=snps['pos'], y=snps['neg_log10_p']),
                              'x', 'y', frac=0.2)
        new_ax2.plot(fit.x, fit.y, c='red', linewidth=2.5,
                    label='LOWESS trend', zorder=4)
    except:
        # Fallback to simple smoothing
        from scipy.signal import savgol_filter
        snps_sorted_pos = snps.sort_values('pos')
        if len(snps_sorted_pos) > 5:
            window = min(11, len(snps_sorted_pos) if len(snps_sorted_pos) % 2 == 1
                        else len(snps_sorted_pos) - 1)
            smoothed = savgol_filter(snps_sorted_pos['neg_log10_p'], window, 3)
            new_ax2.plot(snps_sorted_pos['pos'], smoothed,
                        c='red', linewidth=2.5, label='Smoothed trend', zorder=4)

    new_ax2.axhline(y=6, color='red', linestyle='--', linewidth=1.5, alpha=0.5)
    new_ax2.set_ylabel('-log₁₀(p-value)', fontsize=11, fontweight='bold')
    new_ax2.set_xlabel('')
    new_ax2.set_ylim(5.5, snps['neg_log10_p'].max() + 0.5)
    new_ax2.legend(loc='center left', bbox_to_anchor=(1.02, 0.5),
                  fontsize=9, frameon=False)
    new_ax2.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f'{int(x/1e6)}'))
    new_ax2.tick_params(axis='both', which='major', labelsize=10)
    new_ax2.grid(True, alpha=0.2)

    #--------------------------------------------------------------------------
    # Third axis: Histogram of SNP density
    #--------------------------------------------------------------------------
    x_length = utils.chrom_lengths['hg38']['chrX']
    bins = np.arange(0, x_length, 5_000_000)  # 5 Mb bins for cleaner visualization
    counts, bin_edges = np.histogram(snps['pos'], bins=bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    new_ax3.bar(bin_centers, counts, width=5_000_000*0.9,
               color='steelblue', alpha=0.7, edgecolor='black', linewidth=0.8,
               label='SNP density')

    new_ax3.set_ylabel('Number of SNPs\n(5 Mb bins)', fontsize=11, fontweight='bold')
    new_ax3.set_xlabel('Position on X chromosome (Mb)', fontsize=11, fontweight='bold')
    new_ax3.set_ylim(0, max(counts) + 1)
    new_ax3.legend(loc='center left', bbox_to_anchor=(1.02, 0.5),
                  fontsize=9, frameon=False)
    new_ax3.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f'{int(x/1e6)}'))
    new_ax3.tick_params(axis='both', which='major', labelsize=10)
    new_ax3.grid(True, alpha=0.2, axis='y')

    # Add title
    plt.suptitle('Distribution of 43 High-Significance SNPs on X Chromosome',
                fontsize=14, fontweight='bold', y=0.97)

    # Add subtitle
    plt.gcf().text(0.5, 0.945,
                  f'{len(overlap_genes)} high-confidence genes labeled on chromosome (overlap between high-sig SNPs and Phase 3 regions)',
                  ha='center', fontsize=9, style='italic', color='#555555')

    # Save
    output_file = FIGURES_DIR / 'X_chromosome_ideogram_high_significance_snps.png'
    plt.tight_layout(rect=[0, 0, 0.98, 0.94])
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {output_file.name}")
    plt.close()

#==============================================================================
# 2. X Chromosome Ideogram with Selection Regions
#==============================================================================

def create_regions_ideogram():
    """Create ideogram showing Phase 3 selection regions"""
    print("\nCreating visualization 2: Selection Regions")

    # Load selection regions
    regions = pd.read_csv(SELECTION_REGIONS, sep='\t')
    x_regions = regions[regions['chr'] == 'X'].copy()

    print(f"  Loaded {len(x_regions)} X chromosome selection regions")

    # Load high-sig SNPs to highlight regions containing them
    snps = pd.read_csv(HIGH_SIG_SNPS, sep='\t')

    # Create figure
    fig, ax = plt.subplots(figsize=(16, 8))

    # Create ideogram
    g = gplt.ChromIdeogram('chrX', assembly='hg38', font_size=8)
    g.draw_chromosomes(base=3, height=1)

    # Prepare region segments
    # Classify regions: those containing high-sig SNPs vs others
    regions_with_high_sig = []
    regions_without_high_sig = []

    for idx, region in x_regions.iterrows():
        # Check if any high-sig SNP falls in this region
        has_high_sig = any((snps['pos'] >= region['start']) &
                          (snps['pos'] <= region['end']))

        segment = ('chrX', region['start'], region['end'], region['max_std_ihs'])

        if has_high_sig:
            regions_with_high_sig.append(segment)
        else:
            regions_without_high_sig.append(segment)

    print(f"  Regions with high-sig SNPs: {len(regions_with_high_sig)}")
    print(f"  Regions without high-sig SNPs: {len(regions_without_high_sig)}")

    # Add regions without high-sig SNPs (blue, lower opacity)
    if regions_without_high_sig:
        g.add_segments(regions_without_high_sig, facecolor='tab:blue',
                      alpha=0.3, base=4, height=1.5, label='Selection regions')

    # Add regions with high-sig SNPs (red, higher opacity)
    if regions_with_high_sig:
        g.add_segments(regions_with_high_sig, facecolor='tab:red',
                      alpha=0.6, base=4, height=1.5,
                      label='Regions with high-sig SNPs')

    # Add histogram of region coverage
    new_ax = g.add_axes(1, height_ratio=0.5, hspace=0.3)

    # Create coverage track
    x_length = utils.chrom_lengths['hg38']['chrX']
    bins = np.arange(0, x_length, 500000)  # 500kb bins
    coverage = np.zeros(len(bins))

    for idx, region in x_regions.iterrows():
        start_bin = int(region['start'] / 500000)
        end_bin = int(region['end'] / 500000) + 1
        coverage[start_bin:end_bin] = np.maximum(
            coverage[start_bin:end_bin],
            region['max_std_ihs']
        )

    # Plot coverage
    new_ax.fill_between(bins, 0, coverage, alpha=0.5, color='tab:blue')
    new_ax.plot(bins, coverage, color='tab:blue', linewidth=1)
    new_ax.set_ylabel('Max |Std iHS| in 500kb window', fontsize=11, fontweight='bold')
    new_ax.set_xlabel('Position on X chromosome (Mb)', fontsize=11, fontweight='bold')
    new_ax.axhline(y=2.9, color='gray', linestyle='--', linewidth=1.5, alpha=0.7,
                  label='Phase 3 threshold (|Std iHS|≥2.9)')
    new_ax.grid(True, alpha=0.3)
    new_ax.legend(loc='upper right', fontsize=9)

    # Format x-axis to show Mb
    import matplotlib.ticker as ticker
    new_ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{int(x/1e6)}'))

    # Add legend
    g.add_legend(loc='upper left', fontsize=10)

    # Title
    fig.suptitle(f'X Chromosome: Selection Regions (n={len(x_regions)}, |Std iHS| ≥ 2.9)',
                fontsize=14, fontweight='bold', y=0.98)

    # Save
    output_file = FIGURES_DIR / 'X_chromosome_ideogram_selection_regions.png'
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_file.name}")
    plt.close()

#==============================================================================
# 3. Combined View: SNPs + Genes + Functional Categories
#==============================================================================

def create_combined_ideogram():
    """Create combined ideogram with SNPs, genes, and functional categories"""
    print("\nCreating visualization 3: Combined SNPs + Genes + Functions")

    # Load data
    snps = pd.read_csv(HIGH_SIG_SNPS, sep='\t')
    genes = pd.read_csv(GENE_DETAILS, sep='\t')
    snp_mapping = pd.read_csv(SNP_TO_GENE, sep='\t')

    # Filter for protein-coding genes
    protein_coding = genes[genes['gene_type'] == 'protein-coding'].copy()

    print(f"  {len(snps)} SNPs, {len(protein_coding)} protein-coding genes")

    # Create figure
    fig, ax = plt.subplots(figsize=(16, 10))

    # Create ideogram
    g = gplt.ChromIdeogram('chrX', assembly='hg38', font_size=7)
    g.draw_chromosomes(base=5, height=1)

    # Define functional category colors
    func_colors = {
        'immune': 'red',
        'neural': 'blue',
        'immune; neural': 'purple',
        'developmental': 'green',
        'other': 'gray'
    }

    # Prepare gene segments by functional category
    gene_segments_by_func = {func: [] for func in func_colors.keys()}

    for idx, gene in protein_coding.iterrows():
        func = gene['functional_category']
        if pd.notna(func):
            color = func_colors.get(func, 'gray')
            segment = ('chrX', gene['start'], gene['end'], color)
            gene_segments_by_func[func].append(segment)

    # Add gene segments in layers by function
    base_pos = 6
    for func, color in func_colors.items():
        if gene_segments_by_func[func]:
            g.add_segments(gene_segments_by_func[func], facecolor=color,
                          alpha=0.7, base=base_pos, height=0.5,
                          label=f'{func} genes')
            base_pos += 0.6

    # Add gene labels for overlap genes (those in both high-sig and Phase 3)
    overlap_genes = ['DACH2', 'DIAPH2', 'FAAH2', 'GPC3', 'HTR2C',
                     'IL1RAPL2', 'KDM6A', 'NAP1L2', 'NCBP2L', 'TRPC5', 'ZMAT1']

    gene_labels = []
    for gene_name in overlap_genes:
        gene_info = protein_coding[protein_coding['gene_symbol'] == gene_name]
        if not gene_info.empty:
            gene_row = gene_info.iloc[0]
            pos = (gene_row['start'] + gene_row['end']) / 2
            func = gene_row['functional_category']
            color = func_colors.get(func, 'gray')
            gene_labels.append(('chrX', pos, gene_name, color))

    if gene_labels:
        g.add_labels(gene_labels, base=base_pos, min_height=base_pos+1,
                    bold=overlap_genes)

    # Add SNPs as scatter in axis below
    new_ax = g.add_axes(1, height_ratio=0.4, hspace=0.3)

    # Color SNPs by their nearest gene's function
    snp_colors = []
    for idx, snp in snps.iterrows():
        mapping = snp_mapping[snp_mapping['snp_location'] == snp['Location']]
        if not mapping.empty:
            gene_name = mapping.iloc[0]['nearest_gene']
            gene_info = protein_coding[protein_coding['gene_symbol'] == gene_name]
            if not gene_info.empty:
                func = gene_info.iloc[0]['functional_category']
                snp_colors.append(func_colors.get(func, 'gray'))
            else:
                snp_colors.append('gray')
        else:
            snp_colors.append('gray')

    # Plot SNPs
    sizes = ((snps['neg_log10_p'] - 6) / (snps['neg_log10_p'].max() - 6)) * 100 + 30
    new_ax.scatter(snps['pos'], snps['neg_log10_p'],
                  s=sizes, c=snp_colors, alpha=0.8, edgecolors='black', linewidth=0.5)
    new_ax.set_ylabel('-log10(p-value) of SNP', fontsize=11, fontweight='bold')
    new_ax.set_xlabel('Position on X chromosome (Mb)', fontsize=11, fontweight='bold')
    new_ax.axhline(y=6, color='gray', linestyle='--', linewidth=1.5, alpha=0.7,
                   label='Threshold (p=10⁻⁶)')
    new_ax.grid(True, alpha=0.3)
    new_ax.legend(loc='upper right', fontsize=9)

    # Format x-axis to show Mb
    import matplotlib.ticker as ticker
    new_ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{int(x/1e6)}'))

    # Add legend
    g.add_legend(loc='upper left', fontsize=8, ncol=2)

    # Title
    fig.suptitle('X Chromosome: SNPs, Genes, and Functional Categories',
                fontsize=14, fontweight='bold', y=0.98)

    # Save
    output_file = FIGURES_DIR / 'X_chromosome_ideogram_functional_categories.png'
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_file.name}")
    plt.close()

#==============================================================================
# Main execution
#==============================================================================

def create_comprehensive_snp_table():
    """Create comprehensive table of SNPs with all relevant information"""
    print("\nCreating comprehensive SNP table...")

    # Load data
    snps = pd.read_csv(HIGH_SIG_SNPS, sep='\t')
    snp_mapping = pd.read_csv(SNP_TO_GENE, sep='\t')

    # Merge SNP data with gene mapping
    comprehensive = snps.merge(
        snp_mapping[['snp_location', 'nearest_gene', 'distance', 'overlap_type', 'position_relative']],
        left_on='Location',
        right_on='snp_location',
        how='left'
    )

    # Select and rename columns
    output_table = comprehensive[[
        'chr', 'pos', 'Location', 'neg_log10_p', 'Std iHS',
        'nearest_gene', 'distance', 'overlap_type', 'position_relative'
    ]].copy()

    # Rename columns for clarity
    output_table.columns = [
        'Chromosome', 'Position', 'SNP_ID', 'neg_log10_pvalue',
        'Std_iHS', 'Nearest_Gene', 'Distance_to_Gene_bp',
        'Location_Type', 'Position_Relative_to_Gene'
    ]

    # Add indicator for overlap genes
    overlap_genes = ['DACH2', 'DIAPH2', 'FAAH2', 'GPC3', 'HTR2C',
                     'IL1RAPL2', 'KDM6A', 'NAP1L2', 'NCBP2L', 'TRPC5', 'ZMAT1']
    output_table['Is_Overlap_Gene'] = output_table['Nearest_Gene'].isin(overlap_genes)

    # Sort by position
    output_table = output_table.sort_values('Position')

    # Save to file
    output_file = RESULTS_DIR / 'high_significance_genes/comprehensive_snp_table.tsv'
    output_table.to_csv(output_file, sep='\t', index=False)

    print(f"  Saved: {output_file.name}")
    print(f"  Total SNPs: {len(output_table)}")
    print(f"  SNPs near overlap genes: {output_table['Is_Overlap_Gene'].sum()}")

    # Create summary by gene
    gene_summary = output_table.groupby('Nearest_Gene').agg({
        'SNP_ID': 'count',
        'neg_log10_pvalue': ['min', 'max', 'mean'],
        'Std_iHS': ['min', 'max', 'mean'],
        'Is_Overlap_Gene': 'first'
    }).round(2)

    gene_summary.columns = ['N_SNPs', 'Min_log10p', 'Max_log10p', 'Mean_log10p',
                           'Min_StdIHS', 'Max_StdIHS', 'Mean_StdIHS', 'Is_Overlap_Gene']
    gene_summary = gene_summary.sort_values('Max_log10p', ascending=False)

    summary_file = RESULTS_DIR / 'high_significance_genes/gene_summary_table.tsv'
    gene_summary.to_csv(summary_file, sep='\t')

    print(f"  Saved: {summary_file.name}")
    print(f"  Unique genes: {len(gene_summary)}")

    return output_table, gene_summary


def main():
    """Run all visualizations"""
    try:
        # Create visualization 1
        create_snp_ideogram()

        # Create visualization 2
        create_regions_ideogram()

        # Create visualization 3
        create_combined_ideogram()

        # Create comprehensive tables
        snp_table, gene_summary = create_comprehensive_snp_table()

        print("\n" + "="*70)
        print("All visualizations and tables created successfully!")
        print("="*70)
        print(f"\nOutput directory: {FIGURES_DIR}")
        print("\nVisualization files created:")
        print("  1. X_chromosome_ideogram_high_significance_snps.png")
        print("  2. X_chromosome_ideogram_selection_regions.png")
        print("  3. X_chromosome_ideogram_functional_categories.png")
        print(f"\nData table files created:")
        print("  4. comprehensive_snp_table.tsv - All SNPs with details")
        print("  5. gene_summary_table.tsv - Summary statistics by gene")

    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1

    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main())
