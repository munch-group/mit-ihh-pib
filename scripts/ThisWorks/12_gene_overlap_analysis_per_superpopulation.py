#!/usr/bin/env python3
"""
Gene Overlap Analysis Per Superpopulation

Purpose:
    Analyze protein-coding gene overlaps within each superpopulation (AFR, EUR, EAS, AMR, SAS)
    across their constituent subpopulations.

Outputs per superpopulation:
    1. Gene presence/absence table (rows=genes, columns=subpopulations)
    2. Venn diagrams showing gene overlap patterns
    3. Gene lists by overlap category (shared vs unique)
    4. Summary statistics

Based on: Population_Gene_Overlap_Analysis.py (but in this script I was still pooling populations which doesn't work with iHS)

Author: Mit Van Bruggen 
Date: 2025-12-18
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib_venn import venn2, venn3, venn2_circles, venn3_circles
import seaborn as sns
from pathlib import Path
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

# Superpopulation structure
SUPERPOPULATIONS = {
    'AFR': ['ACB', 'ASW', 'ESN', 'GWD', 'LWK', 'MSL', 'YRI'],
    'AMR': ['CLM', 'MXL', 'PEL', 'PUR'],
    'EAS': ['CDX', 'CHB', 'CHS', 'JPT', 'KHV'],
    'EUR': ['CEU', 'FIN', 'GBR', 'IBS', 'TSI'],
    'SAS': ['BEB', 'GIH', 'ITU', 'PJL', 'STU']
}

# Paths
base_dir = Path("/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib")
output_base = base_dir / "results/gene_overlap_per_superpopulation"
output_base.mkdir(parents=True, exist_ok=True)

# Set style
sns.set_style("whitegrid")

print("=" * 80)
print("Gene Overlap Analysis Per Superpopulation")
print("=" * 80)
print()
print("Analyzing protein-coding genes under selection within each superpopulation")
print("Comparing subpopulation overlap patterns")
print()


# ============================================================================
# Data Loading Functions
# ============================================================================

def load_gene_set(pop_code: str):
    """
    Load protein-coding genes for a population
    Returns: set of gene names
    """
    genes_file = base_dir / f"results/ihs_{pop_code}/analysis/gene_annotation/genes_in_regions_protein_coding.tsv"

    if not genes_file.exists():
        return None

    df = pd.read_csv(genes_file, sep='\t')

    # Get unique gene names
    genes = set(df['gene_name'].unique())

    return genes


# ============================================================================
# Analysis Functions
# ============================================================================

def create_gene_presence_table(gene_sets: dict, superpop: str):
    """
    Create a binary table showing gene presence/absence across populations
    Rows = genes, Columns = populations
    """
    # Get all unique genes
    all_genes = set()
    for genes in gene_sets.values():
        all_genes.update(genes)

    # Build presence/absence matrix
    presence_data = []
    for gene in sorted(all_genes):
        gene_row = {'gene_name': gene}
        for pop, pop_genes in gene_sets.items():
            gene_row[pop] = 1 if gene in pop_genes else 0

        # Add summary statistics
        n_pops = sum([gene in pop_genes for pop_genes in gene_sets.values()])
        gene_row['n_populations'] = n_pops
        gene_row['populations'] = ','.join([pop for pop, pop_genes in gene_sets.items() if gene in pop_genes])

        presence_data.append(gene_row)

    presence_df = pd.DataFrame(presence_data)

    # Sort by number of populations (descending), then by gene name
    presence_df = presence_df.sort_values(['n_populations', 'gene_name'], ascending=[False, True])

    return presence_df


def analyze_overlap_patterns(gene_sets: dict, superpop: str):
    """
    Analyze gene overlap patterns within a superpopulation
    Returns: dictionary with overlap categories
    """
    n_pops = len(gene_sets)
    all_genes = set()
    for genes in gene_sets.values():
        all_genes.update(genes)

    # Categorize genes by number of populations
    overlap_categories = {}

    # Genes in ALL populations
    if n_pops > 0:
        all_pops_genes = set.intersection(*gene_sets.values())
        overlap_categories['all_populations'] = all_pops_genes
    else:
        overlap_categories['all_populations'] = set()

    # Genes in exactly N populations (for N from n_pops-1 down to 1)
    for n in range(n_pops - 1, 0, -1):
        # Find genes in exactly n populations
        genes_in_n_pops = set()
        for gene in all_genes:
            n_occurrences = sum([gene in pop_genes for pop_genes in gene_sets.values()])
            if n_occurrences == n:
                genes_in_n_pops.add(gene)

        overlap_categories[f'exactly_{n}_populations'] = genes_in_n_pops

    # Population-specific genes (for detailed breakdown)
    unique_genes = {}
    for pop, pop_genes in gene_sets.items():
        # Genes only in this population
        other_genes = set()
        for other_pop, other_pop_genes in gene_sets.items():
            if other_pop != pop:
                other_genes.update(other_pop_genes)
        unique = pop_genes - other_genes
        unique_genes[pop] = unique

    overlap_categories['unique_per_population'] = unique_genes

    return overlap_categories

def create_upset_style_plot(gene_sets: dict, superpop: str, output_dir: Path):
    """
    Create UpSet-style plot for 4+ population overlaps (which is the case for all of mine)
    Shows intersection sizes in a bar chart with gene names listed
    """
    from itertools import combinations

    pop_names = list(gene_sets.keys())
    n_pops = len(pop_names)

    # Calculate all possible intersections
    intersection_data = []

    # For each possible combination size (from all populations down to 1)
    for size in range(n_pops, 0, -1):
        for combo in combinations(pop_names, size):
            # Calculate intersection: genes in ALL of these populations but NONE of the others
            combo_set = set(combo)
            other_pops = set(pop_names) - combo_set

            # Start with genes in first population of combo
            intersection = gene_sets[combo[0]].copy()

            # Intersect with all other populations in combo
            for pop in combo[1:]:
                intersection &= gene_sets[pop]

            # Remove genes that are in any other population
            for other_pop in other_pops:
                intersection -= gene_sets[other_pop]

            if len(intersection) > 0:
                intersection_data.append({
                    'populations': combo,
                    'n_populations': len(combo),
                    'size': len(intersection),
                    'genes': sorted(intersection)
                })

    # Sort by size
    intersection_data.sort(key=lambda x: x['size'], reverse=True)

    # Decide layout based on number of bars (if too many, use legend instead of text boxes)
    use_legend = len(intersection_data) > 15

    if use_legend:
        # For many bars: use legend on the side
        fig = plt.figure(figsize=(28, 16))
        gs = fig.add_gridspec(2, 2, height_ratios=[4, 1], width_ratios=[3, 1], hspace=0.3, wspace=0.1)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[1, 0])
        ax_legend = fig.add_subplot(gs[:, 1])
    else:
        # For few bars: use text boxes under bars
        fig = plt.figure(figsize=(24, 16))
        gs = fig.add_gridspec(2, 1, height_ratios=[3.5, 1], hspace=0.5)
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1])

    # Top panel: Bar chart
    x_pos = np.arange(len(intersection_data))
    sizes = [d['size'] for d in intersection_data]
    colors = plt.cm.tab20(np.linspace(0, 1, len(intersection_data)))  # More colors for many bars

    bars = ax1.bar(x_pos, sizes, color=colors, edgecolor='black', linewidth=1.5, width=0.8)
    ax1.set_ylabel('Number of Genes', fontsize=14, fontweight='bold')
    ax1.set_title(f'{superpop}: Protein-Coding Gene Overlaps Between Subpopulations',
                 fontsize=16, fontweight='bold', pad=20)
    ax1.grid(True, alpha=0.3, axis='y')

    # Add value labels on top of bars
    for i, size in enumerate(sizes):
        ax1.text(i, size + max(sizes)*0.02, str(size),
                ha='center', va='bottom', fontsize=11, fontweight='bold')

    # Set y-axis to start at 0 and add space below for text
    ax1.set_ylim(0, max(sizes) * 1.15)

    if use_legend:
        # Add Set labels on x-axis
        ax1.set_xticks(x_pos)
        ax1.set_xticklabels([f'Set {i+1}' for i in range(len(intersection_data))],
                           fontsize=9, rotation=45, ha='right')
        ax1.set_xlabel('Intersection Sets', fontsize=12, fontweight='bold')

        # Create legend panel with gene lists
        ax_legend.axis('off')
        legend_text = []
        for i, data in enumerate(intersection_data):
            pops_str = '+'.join(data['populations'])
            genes_display = data['genes'][:15]  # Show first 15 genes
            genes_str = ', '.join(genes_display)
            if len(data['genes']) > 15:
                genes_str += f' (+{len(data['genes']) - 15} more)'

            legend_text.append(f"Set {i+1} [{pops_str}]:\n{genes_str}\n")

        # Display legend text
        legend_str = '\n'.join(legend_text)
        ax_legend.text(0.05, 0.98, legend_str, transform=ax_legend.transAxes,
                      fontsize=8, verticalalignment='top', fontfamily='monospace',
                      bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.2))
        ax_legend.set_title('Gene Lists', fontsize=12, fontweight='bold', pad=10)

    else:
        # Add gene names under each bar as text boxes
        text_fontsize = max(9, min(12, 120 / len(intersection_data)))

        for i, data in enumerate(intersection_data):
            # Format gene list with max 3 genes per line (show first 12 genes)
            # BUT: if a gene name starts with ENSG (long identifier), put it on its own line
            genes_display = data['genes'][:15]  # Show up to 15 genes
            genes_lines = []
            current_line = []

            for gene in genes_display:
                # Check if gene is a long ENSG identifier
                if gene.startswith('ENSG'):
                    # Flush current line if it has genes
                    if current_line:
                        genes_lines.append(', '.join(current_line))
                        current_line = []
                    # Add ENSG gene on its own line
                    genes_lines.append(gene)
                else:
                    # Regular gene name - add to current line (max 2 per line)
                    current_line.append(gene)
                    if len(current_line) >= 2:
                        genes_lines.append(', '.join(current_line))
                        current_line = []

            # Flush any remaining genes in current line
            if current_line:
                genes_lines.append(', '.join(current_line))

            genes_str = '\n'.join(genes_lines)
            if len(data['genes']) > 15:
                genes_str += f'\n+{len(data['genes']) - 15} more'

            # Add text annotation below the axis with more vertical offset to avoid overlap with grid
            ax1.annotate(genes_str,
                        xy=(i, 0), xycoords='data',
                        xytext=(0, -30), textcoords='offset points',
                        ha='center', va='top',
                        fontsize=text_fontsize,
                        bbox=dict(boxstyle='round,pad=0.5',
                                 facecolor=colors[i],
                                 edgecolor='black',
                                 alpha=0.3, linewidth=1.5),
                        wrap=True)

        # Remove x-axis labels (replaced by text boxes)
        ax1.set_xticks([])
        ax1.set_xlabel('')
        ax1.margins(y=0.05)

    # Bottom panel: Population membership grid
    # Create binary matrix showing which populations are in each intersection
    matrix = np.zeros((len(pop_names), len(intersection_data)))
    for j, data in enumerate(intersection_data):
        for i, pop in enumerate(pop_names):
            if pop in data['populations']:
                matrix[i, j] = 1

    ax2.imshow(matrix, cmap='RdYlGn', aspect='auto', vmin=0, vmax=1)

    # Set labels
    ax2.set_yticks(range(len(pop_names)))
    ax2.set_yticklabels(pop_names, fontsize=11, fontweight='bold')
    ax2.set_xticks(range(len(intersection_data)))
    ax2.set_xticklabels([f'Set {i+1}' for i in range(len(intersection_data))],
                        fontsize=9, rotation=0)
    ax2.set_xlabel('Population Membership', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Populations', fontsize=12, fontweight='bold')

    # Add grid
    ax2.set_xticks(np.arange(len(intersection_data)) - 0.5, minor=True)
    ax2.set_yticks(np.arange(len(pop_names)) - 0.5, minor=True)
    ax2.grid(which='minor', color='black', linewidth=1.5)

    plt.tight_layout()
    output_file = output_dir / f"{superpop}_overlap_upset.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.savefig(output_file.with_suffix('.pdf'), bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file.name}")

    # Also save a detailed text file with ALL genes for each intersection
    text_file = output_dir / f"{superpop}_overlap_details.txt"
    with open(text_file, 'w') as f:
        f.write(f"{superpop}: Gene Overlap Details\n")
        f.write("=" * 80 + "\n\n")

        for i, data in enumerate(intersection_data):
            pops_str = ' + '.join(data['populations'])
            f.write(f"Set {i+1}: {pops_str}\n")
            f.write(f"  Populations: {len(data['populations'])}\n")
            f.write(f"  Number of genes: {len(data['genes'])}\n")
            f.write(f"  Genes: {', '.join(data['genes'])}\n")
            f.write("\n" + "-" * 80 + "\n\n")

    print(f"  Saved: {text_file.name} (full gene lists)")


# ============================================================================
# Main Processing Function
# ============================================================================

def process_superpopulation(superpop: str, populations: list):
    """
    Complete analysis for one superpopulation
    """
    print(f"\n{superpop}:")
    print("-" * 40)

    # Create output directory
    output_dir = output_base / superpop
    output_dir.mkdir(exist_ok=True)

    # Load gene sets
    print(f"  Loading gene sets for {len(populations)} populations...")
    gene_sets = {}
    for pop in populations:
        genes = load_gene_set(pop)
        if genes is not None and len(genes) > 0:
            gene_sets[pop] = genes
            print(f"    {pop}: {len(genes)} genes")
        else:
            print(f"    {pop}: No data")

    if len(gene_sets) == 0:
        print(f"  No data available for {superpop}")
        return None

    # 1. Create gene presence/absence table
    print(f"\n  Creating gene presence/absence table...")
    presence_table = create_gene_presence_table(gene_sets, superpop)
    presence_file = output_dir / f"{superpop}_gene_presence_table.tsv"
    presence_table.to_csv(presence_file, sep='\t', index=False)
    print(f"    Saved: {presence_file.name}")
    print(f"    Total unique genes: {len(presence_table)}")

    # 2. Analyze overlap patterns
    print(f"\n  Analyzing overlap patterns...")
    overlap_categories = analyze_overlap_patterns(gene_sets, superpop)

    # Save overlap categories
    n_all = len(overlap_categories['all_populations'])
    print(f"    Genes in ALL {len(gene_sets)} populations: {n_all}")

    if n_all > 0:
        all_file = output_dir / f"{superpop}_genes_all_populations.txt"
        with open(all_file, 'w') as f:
            f.write('\n'.join(sorted(overlap_categories['all_populations'])))
        print(f"      Saved: {all_file.name}")

    # Save population-specific genes
    unique_genes = overlap_categories['unique_per_population']
    for pop, genes in unique_genes.items():
        if len(genes) > 0:
            unique_file = output_dir / f"{superpop}_{pop}_unique_genes.txt"
            with open(unique_file, 'w') as f:
                f.write('\n'.join(sorted(genes)))
            print(f"    {pop} unique genes: {len(genes)} -> {unique_file.name}")

    # 3. Create Venn diagram or UpSet plot
    print(f"\n  Creating visualization...")
    create_venn_diagram_n_way(gene_sets, superpop, output_dir)

    # 4. Summary statistics
    all_genes = set()
    for genes in gene_sets.values():
        all_genes.update(genes)

    summary = {
        'superpopulation': superpop,
        'n_populations_with_data': len(gene_sets),
        'total_unique_genes': len(all_genes),
        'genes_in_all_populations': len(overlap_categories['all_populations']),
        'min_genes_per_population': min([len(genes) for genes in gene_sets.values()]),
        'max_genes_per_population': max([len(genes) for genes in gene_sets.values()]),
        'mean_genes_per_population': np.mean([len(genes) for genes in gene_sets.values()])
    }

    return summary


# ============================================================================
# Main Execution
# ============================================================================

def main():
    summary_stats = []

    for superpop, populations in SUPERPOPULATIONS.items():
        result = process_superpopulation(superpop, populations)
        if result:
            summary_stats.append(result)

    # Save overall summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80 + "\n")

    summary_df = pd.DataFrame(summary_stats)
    summary_file = output_base / "superpopulation_summary.tsv"
    summary_df.to_csv(summary_file, sep='\t', index=False)

    print(summary_df.to_string(index=False))
    print()
    print(f"Summary saved: {summary_file}")
    print(f"All results saved to: {output_base}")
    print()


if __name__ == "__main__":
    main()
