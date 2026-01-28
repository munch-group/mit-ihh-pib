#!/usr/bin/env python3
"""
Cross-Continental Gene Comparison with Provenance Tracking

Purpose:
    Compare genes under selection across continents (superpopulations) by analyzing
    the union of genes from all subpopulations within each continent, while tracking
    which specific subpopulations contribute each gene.

Analysis:
    1. Create gene unions for each superpopulation with provenance tracking
    2. Compare gene sets between continents (pairwise and overall)
    3. Classify overlap patterns: broad-broad, broad-narrow, narrow-narrow
    4. Identify continent-specific vs shared selection signals

Outputs:
    1. Superpopulation gene union tables with subpopulation provenance
    2. Pairwise continental comparison tables
    3. Provenance analysis for overlapping genes
    4. Venn diagrams for pairwise comparisons
    5. Heatmaps showing cross-continental overlap patterns
    6. Summary statistics and classifications

Author: Mit Van Bruggen (got help from Claude to automate comparisons across all)
Date: 2025-12-18
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2, venn2_circles, venn3, venn3_circles
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
output_base = base_dir / "results/cross_continental_gene_comparison"
output_base.mkdir(parents=True, exist_ok=True)

# Set style
sns.set_style("whitegrid")
plt.rcParams['font.size'] = 11

print("=" * 80)
print("Cross-Continental Gene Comparison with Provenance Tracking")
print("=" * 80)
print()
print("Comparing genes under selection across continental populations")
print("Tracking which subpopulations contribute each gene")
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


def create_superpopulation_unions():
    """
    Create gene unions for each superpopulation with provenance tracking

    Returns:
        dict: {superpop: {gene: [list of subpopulations]}}
    """
    print("Creating superpopulation gene unions with provenance...")
    print()

    superpop_unions = {}

    for superpop, subpopulations in SUPERPOPULATIONS.items():
        print(f"{superpop}:")
        gene_provenance = {}

        # Load genes for each subpopulation
        for subpop in subpopulations:
            genes = load_gene_set(subpop)
            if genes is not None and len(genes) > 0:
                print(f"  {subpop}: {len(genes)} genes")

                # Add to provenance
                for gene in genes:
                    if gene not in gene_provenance:
                        gene_provenance[gene] = []
                    gene_provenance[gene].append(subpop)
            else:
                print(f"  {subpop}: No data")

        if gene_provenance:
            superpop_unions[superpop] = gene_provenance
            total_genes = len(gene_provenance)
            print(f"  → Total unique genes in {superpop}: {total_genes}")
        else:
            print(f"  → No data available for {superpop}")

        print()

    return superpop_unions


# ============================================================================
# Analysis Functions
# ============================================================================

def classify_signal_breadth(n_subpops, total_subpops, threshold=0.5):
    """
    Classify selection signal as 'broad' or 'narrow'

    Args:
        n_subpops: number of subpopulations with the gene
        total_subpops: total number of subpopulations in superpopulation
        threshold: proportion threshold for 'broad' classification

    Returns:
        str: 'broad' or 'narrow'
    """
    proportion = n_subpops / total_subpops
    return 'broad' if proportion >= threshold else 'narrow'


def analyze_pairwise_overlap(superpop_unions, superpop1, superpop2):
    """
    Analyze gene overlap between two superpopulations

    Returns:
        DataFrame with provenance information for all genes
    """
    genes1 = superpop_unions[superpop1]
    genes2 = superpop_unions[superpop2]

    all_genes = set(genes1.keys()) | set(genes2.keys())

    overlap_data = []

    for gene in sorted(all_genes):
        in_sp1 = gene in genes1
        in_sp2 = gene in genes2

        # Get subpopulation lists
        subpops1 = genes1.get(gene, [])
        subpops2 = genes2.get(gene, [])

        n_subpops1 = len(subpops1)
        n_subpops2 = len(subpops2)

        # Total subpopulations in each superpopulation
        total_subpops1 = len(SUPERPOPULATIONS[superpop1])
        total_subpops2 = len(SUPERPOPULATIONS[superpop2])

        # Classify breadth
        breadth1 = classify_signal_breadth(n_subpops1, total_subpops1) if in_sp1 else 'absent'
        breadth2 = classify_signal_breadth(n_subpops2, total_subpops2) if in_sp2 else 'absent'

        # Overall category
        if in_sp1 and in_sp2:
            category = 'shared'
            pattern = f"{breadth1}-{breadth2}"
        elif in_sp1:
            category = f'{superpop1}_specific'
            pattern = breadth1
        else:
            category = f'{superpop2}_specific'
            pattern = breadth2

        overlap_data.append({
            'gene_name': gene,
            'category': category,
            'pattern': pattern,
            f'in_{superpop1}': int(in_sp1),
            f'in_{superpop2}': int(in_sp2),
            f'{superpop1}_n_subpops': n_subpops1,
            f'{superpop1}_total_subpops': total_subpops1,
            f'{superpop1}_proportion': n_subpops1 / total_subpops1 if in_sp1 else 0,
            f'{superpop1}_subpops': ','.join(sorted(subpops1)) if subpops1 else '',
            f'{superpop2}_n_subpops': n_subpops2,
            f'{superpop2}_total_subpops': total_subpops2,
            f'{superpop2}_proportion': n_subpops2 / total_subpops2 if in_sp2 else 0,
            f'{superpop2}_subpops': ','.join(sorted(subpops2)) if subpops2 else '',
        })

    df = pd.DataFrame(overlap_data)

    # Sort by category, then by pattern
    df = df.sort_values(['category', 'pattern', 'gene_name'])

    return df


def create_pairwise_venn(superpop_unions, superpop1, superpop2, output_dir):
    """
    Create Venn diagram comparing two superpopulations with gene annotations
    """
    genes1 = set(superpop_unions[superpop1].keys())
    genes2 = set(superpop_unions[superpop2].keys())

    # Calculate overlaps
    shared_genes = genes1 & genes2
    unique1 = genes1 - genes2
    unique2 = genes2 - genes1

    # Create larger figure to accommodate annotations
    fig, ax = plt.subplots(figsize=(18, 14))

    venn = venn2([genes1, genes2], set_labels=[superpop1, superpop2], ax=ax,
                 set_colors=('#3498db', '#e74c3c'), alpha=0.6)
    venn_circles = venn2_circles([genes1, genes2], ax=ax)

    # Style circles
    for circle in venn_circles:
        circle.set_linewidth(2)
        circle.set_edgecolor('black')

    # Add title with counts
    n_shared = len(shared_genes)
    n_unique1 = len(unique1)
    n_unique2 = len(unique2)

    title = f'Gene Overlap: {superpop1} vs {superpop2}\n'
    title += f'Shared: {n_shared} | {superpop1} only: {n_unique1} | {superpop2} only: {n_unique2}'
    ax.set_title(title, fontsize=16, fontweight='bold', pad=30)

    # Add gene names with text boxes and arrows
    # Annotations: (gene_set, label, text_position, arrow_target, horizontal_align, vertical_align, arrow_curve)
    annotations = [
        # Unique to superpop1 (left circle only)
        (unique1, f'{superpop1} only', (-0.75, 0.40), (-0.35, 0.0), 'right', 'center', 'arc3,rad=0.2'),
        # Unique to superpop2 (right circle only)
        (unique2, f'{superpop2} only', (0.75, 0.40), (0.35, 0.0), 'left', 'center', 'arc3,rad=-0.2'),
        # Shared (overlap region)
        (shared_genes, 'Shared', (0.0, -0.70), (0.0, 0.0), 'center', 'top', 'arc3,rad=0'),
    ]

    for gene_set, label, text_pos, arrow_pos, ha, va, connection_style in annotations:
        if gene_set:  # Only annotate if there are genes in this region
            # Limit number of genes shown (show first 30, indicate if more)
            genes_to_show = sorted(gene_set)[:30]
            gene_text = '\n'.join(genes_to_show)

            if len(gene_set) > 30:
                gene_text += f'\n... +{len(gene_set) - 30} more'

            full_text = f"{label} (n={len(gene_set)}):\n{gene_text}"

            # Add text box with arrow
            ax.annotate(
                full_text,
                xy=arrow_pos,           # Arrow points to this location (inside Venn region)
                xytext=text_pos,        # Text box location
                fontsize=10,
                family='monospace',
                ha=ha,
                va=va,
                bbox=dict(
                    boxstyle='round,pad=0.5',
                    facecolor='wheat',
                    alpha=0.9,
                    edgecolor='#8B4513',
                    linewidth=1.5
                ),
                arrowprops=dict(
                    arrowstyle='-|>',
                    connectionstyle=connection_style,
                    color='#555555',
                    lw=1.5,
                    mutation_scale=12
                )
            )

    # Adjust axis limits to accommodate all annotations
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1.0, 0.8)
    ax.axis('off')

    output_file = output_dir / f"venn_{superpop1}_vs_{superpop2}.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.savefig(output_file.with_suffix('.pdf'), bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()

    return output_file.name


def create_overlap_heatmap(superpop_unions, output_dir):  
    """
    Create heatmap showing gene overlap between all superpopulations
    """
    superpops = list(superpop_unions.keys())
    n = len(superpops)

    # Create overlap matrix
    overlap_matrix = np.zeros((n, n))
    jaccard_matrix = np.zeros((n, n))

    for i, sp1 in enumerate(superpops):
        for j, sp2 in enumerate(superpops):
            genes1 = set(superpop_unions[sp1].keys())
            genes2 = set(superpop_unions[sp2].keys())

            overlap = len(genes1 & genes2)
            union = len(genes1 | genes2)

            overlap_matrix[i, j] = overlap
            jaccard_matrix[i, j] = overlap / union if union > 0 else 0

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

    # Heatmap 1: Absolute overlap counts
    sns.heatmap(overlap_matrix, annot=True, fmt='.0f', cmap='YlOrRd',
                xticklabels=superpops, yticklabels=superpops,
                ax=ax1, cbar_kws={'label': 'Number of Shared Genes'},
                linewidths=1, linecolor='black')
    ax1.set_title('Gene Overlap Between Continents\n(Absolute Counts)',
                  fontsize=14, fontweight='bold', pad=15)
    ax1.set_xlabel('Superpopulation', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Superpopulation', fontsize=12, fontweight='bold')

    # Heatmap 2: Jaccard similarity
    sns.heatmap(jaccard_matrix, annot=True, fmt='.3f', cmap='viridis',
                xticklabels=superpops, yticklabels=superpops,
                ax=ax2, cbar_kws={'label': 'Jaccard Similarity'},
                linewidths=1, linecolor='black', vmin=0, vmax=1)
    ax2.set_title('Gene Overlap Between Continents\n(Jaccard Similarity)',
                  fontsize=14, fontweight='bold', pad=15)
    ax2.set_xlabel('Superpopulation', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Superpopulation', fontsize=12, fontweight='bold')

    plt.tight_layout()
    output_file = output_dir / "overlap_heatmap_all_continents.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.savefig(output_file.with_suffix('.pdf'), bbox_inches='tight')
    plt.close()

    return overlap_matrix, jaccard_matrix, output_file.name


def create_pattern_breakdown_plot(pairwise_results, superpop1, superpop2, output_dir):
    """
    Create bar plot showing breakdown of overlap patterns
    """
    # Count patterns
    shared_genes = pairwise_results[pairwise_results['category'] == 'shared']

    if len(shared_genes) == 0:
        print(f"    No shared genes between {superpop1} and {superpop2}")
        return None

    pattern_counts = shared_genes['pattern'].value_counts()

    # Create plot
    fig, ax = plt.subplots(figsize=(10, 6))

    colors = {'broad-broad': '#2ecc71', 'broad-narrow': '#f39c12',
              'narrow-broad': '#f39c12', 'narrow-narrow': '#e74c3c'}

    x_pos = range(len(pattern_counts))
    bars = ax.bar(x_pos, pattern_counts.values,
                  color=[colors.get(p, '#95a5a6') for p in pattern_counts.index],
                  edgecolor='black', linewidth=1.5)

    # Add value labels
    for i, v in enumerate(pattern_counts.values):
        ax.text(i, v + max(pattern_counts.values)*0.02, str(v),
                ha='center', va='bottom', fontsize=12, fontweight='bold')

    ax.set_xticks(x_pos)
    ax.set_xticklabels(pattern_counts.index, fontsize=11, fontweight='bold')
    ax.set_ylabel('Number of Genes', fontsize=12, fontweight='bold')
    ax.set_title(f'Shared Gene Patterns: {superpop1} vs {superpop2}\n' +
                 f'Total Shared: {len(shared_genes)} genes',
                 fontsize=14, fontweight='bold', pad=15)
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    output_file = output_dir / f"pattern_breakdown_{superpop1}_vs_{superpop2}.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.savefig(output_file.with_suffix('.pdf'), bbox_inches='tight')
    plt.close()

    return output_file.name


def create_provenance_summary_table(superpop_unions):
    """
    Create summary table showing gene distribution across subpopulations
    """
    summary_data = []

    for superpop, gene_provenance in superpop_unions.items():
        total_subpops = len(SUPERPOPULATIONS[superpop])

        # Count genes by number of subpopulations
        n_subpop_counts = {}
        for gene, subpops in gene_provenance.items():
            n = len(subpops)
            if n not in n_subpop_counts:
                n_subpop_counts[n] = 0
            n_subpop_counts[n] += 1

        # Calculate broad vs narrow
        n_broad = sum(count for n, count in n_subpop_counts.items()
                     if n / total_subpops >= 0.5)
        n_narrow = sum(count for n, count in n_subpop_counts.items()
                      if n / total_subpops < 0.5)

        summary_data.append({
            'superpopulation': superpop,
            'total_unique_genes': len(gene_provenance),
            'n_subpopulations': total_subpops,
            'genes_broad_signal': n_broad,
            'genes_narrow_signal': n_narrow,
            'proportion_broad': n_broad / len(gene_provenance) if len(gene_provenance) > 0 else 0,
            'mean_subpops_per_gene': np.mean([len(subpops) for subpops in gene_provenance.values()]),
            'median_subpops_per_gene': np.median([len(subpops) for subpops in gene_provenance.values()])
        })

    return pd.DataFrame(summary_data)


def analyze_multiway_overlap(superpop_unions):
    """
    Analyze gene overlap across ALL superpopulations

    Returns:
        DataFrame with genes categorized by number of continents
    """
    # Get all unique genes
    all_genes = set()
    for gene_provenance in superpop_unions.values():
        all_genes.update(gene_provenance.keys())

    # For each gene, determine which superpopulations contain it
    gene_data = []

    for gene in sorted(all_genes):
        superpops_with_gene = []
        superpop_details = {}

        for superpop, gene_provenance in superpop_unions.items():
            if gene in gene_provenance:
                superpops_with_gene.append(superpop)
                subpops = gene_provenance[gene]
                total_subpops = len(SUPERPOPULATIONS[superpop])
                breadth = classify_signal_breadth(len(subpops), total_subpops)

                superpop_details[superpop] = {
                    'n_subpops': len(subpops),
                    'total_subpops': total_subpops,
                    'proportion': len(subpops) / total_subpops,
                    'breadth': breadth,
                    'subpops': ','.join(sorted(subpops))
                }

        n_continents = len(superpops_with_gene)

        # Create row
        row = {
            'gene_name': gene,
            'n_continents': n_continents,
            'continents': ','.join(sorted(superpops_with_gene))
        }

        # Add per-continent details
        for superpop in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']:
            if superpop in superpop_details:
                details = superpop_details[superpop]
                row[f'{superpop}_present'] = 1
                row[f'{superpop}_n_subpops'] = details['n_subpops']
                row[f'{superpop}_proportion'] = details['proportion']
                row[f'{superpop}_breadth'] = details['breadth']
                row[f'{superpop}_subpops'] = details['subpops']
            else:
                row[f'{superpop}_present'] = 0
                row[f'{superpop}_n_subpops'] = 0
                row[f'{superpop}_proportion'] = 0
                row[f'{superpop}_breadth'] = 'absent'
                row[f'{superpop}_subpops'] = ''

        gene_data.append(row)

    df = pd.DataFrame(gene_data)
    df = df.sort_values(['n_continents', 'gene_name'], ascending=[False, True])

    return df


def create_multiway_upset_plot(multiway_df, output_dir):
    """
    Create UpSet-style plot showing gene overlap across all continents
    """
    # Count genes by number of continents
    continent_counts = multiway_df['n_continents'].value_counts().sort_index(ascending=False)

    # Create figure with two panels
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.3)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    # Top panel: Bar chart of gene counts by number of continents
    x_pos = range(len(continent_counts))
    colors = plt.cm.RdYlGn(np.linspace(0.3, 0.9, len(continent_counts)))

    bars = ax1.bar(x_pos, continent_counts.values, color=colors,
                   edgecolor='black', linewidth=1.5, width=0.6)

    ax1.set_xticks(x_pos)
    ax1.set_xticklabels([f'{n} continents' for n in continent_counts.index],
                        fontsize=12, fontweight='bold')
    ax1.set_ylabel('Number of Genes', fontsize=14, fontweight='bold')
    ax1.set_title('Gene Distribution Across Continents\n(How many continents share each gene?)',
                 fontsize=16, fontweight='bold', pad=20)
    ax1.grid(True, alpha=0.3, axis='y')

    # Add value labels on bars
    for i, (n_cont, count) in enumerate(continent_counts.items()):
        ax1.text(i, count + max(continent_counts.values)*0.02, str(count),
                ha='center', va='bottom', fontsize=12, fontweight='bold')

        # Add percentage
        total_genes = len(multiway_df)
        pct = count / total_genes * 100
        ax1.text(i, count/2, f'{pct:.1f}%',
                ha='center', va='center', fontsize=10, fontweight='bold',
                color='white' if pct > 20 else 'black')

    ax1.set_ylim(0, max(continent_counts.values) * 1.15)

    # Bottom panel: Summary text with example genes
    ax2.axis('off')

    summary_lines = []
    summary_lines.append("GENE EXAMPLES BY CONTINENTAL BREADTH:\n")

    for n_cont in sorted(continent_counts.index, reverse=True):
        genes_at_n = multiway_df[multiway_df['n_continents'] == n_cont]
        n_genes = len(genes_at_n)

        if n_cont == 5:
            label = f"Universal (All 5 continents): {n_genes} genes"
        elif n_cont == 1:
            label = f"Continent-specific (1 continent only): {n_genes} genes"
        else:
            label = f"In exactly {n_cont} continents: {n_genes} genes"

        # Show first 10 gene examples
        example_genes = sorted(genes_at_n['gene_name'].tolist())[:10]
        genes_str = ', '.join(example_genes)
        if len(genes_at_n) > 10:
            genes_str += f' (+{len(genes_at_n)-10} more)'

        summary_lines.append(f"{label}")
        summary_lines.append(f"  Examples: {genes_str}\n")

    summary_text = '\n'.join(summary_lines)
    ax2.text(0.05, 0.95, summary_text, transform=ax2.transAxes,
            fontsize=10, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

    plt.tight_layout()
    output_file = output_dir / "multiway_continent_overlap.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(output_file.with_suffix('.pdf'), bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"  Saved: {output_file.name}")

    return output_file.name


def create_continent_combination_analysis(multiway_df, output_dir):
    """
    Analyze which specific combinations of continents share genes
    Focus on genes in 2-4 continents to see common patterns
    """
    # Filter to genes in 2-4 continents (most interesting range)
    subset = multiway_df[(multiway_df['n_continents'] >= 2) &
                        (multiway_df['n_continents'] <= 4)].copy()

    if len(subset) == 0:
        print("  No genes found in 2-4 continents")
        return None

    # Count combinations
    combo_counts = subset['continents'].value_counts().head(20)

    # Create bar plot
    fig, ax = plt.subplots(figsize=(14, 8))

    y_pos = range(len(combo_counts))
    colors = plt.cm.viridis(np.linspace(0.3, 0.9, len(combo_counts)))

    bars = ax.barh(y_pos, combo_counts.values, color=colors,
                   edgecolor='black', linewidth=1.5)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(combo_counts.index, fontsize=10, fontweight='bold')
    ax.set_xlabel('Number of Genes', fontsize=12, fontweight='bold')
    ax.set_title('Top 20 Continent Combinations Sharing Genes\n(Genes in 2-4 continents)',
                fontsize=14, fontweight='bold', pad=15)
    ax.grid(True, alpha=0.3, axis='x')

    # Add value labels
    for i, v in enumerate(combo_counts.values):
        ax.text(v + max(combo_counts.values)*0.01, i, str(v),
               va='center', fontsize=10, fontweight='bold')

    ax.set_xlim(0, max(combo_counts.values) * 1.15)

    plt.tight_layout()
    output_file = output_dir / "continent_combinations_top20.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(output_file.with_suffix('.pdf'), bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"  Saved: {output_file.name}")

    # Also save detailed table of all combinations
    combo_details = []
    for combo_str, count in subset['continents'].value_counts().items():
        genes = subset[subset['continents'] == combo_str]['gene_name'].tolist()
        combo_details.append({
            'continent_combination': combo_str,
            'n_continents': len(combo_str.split(',')),
            'n_genes': count,
            'genes': ','.join(sorted(genes))
        })

    combo_df = pd.DataFrame(combo_details)
    combo_df = combo_df.sort_values(['n_continents', 'n_genes'], ascending=[False, False])

    combo_file = output_dir / "continent_combinations_detailed.tsv"
    combo_df.to_csv(combo_file, sep='\t', index=False)
    print(f"  Saved: {combo_file.name}")

    return output_file.name


def create_threeway_venn_afr_eas_eur(superpop_unions, output_dir):
    """
    Create 3-way Venn diagram for AFR, EAS, and EUR with gene annotations
    Similar to the original Population_Gene_Overlap_Analysis.py but using
    the union of all subpopulations within each continent
    """
    # Check if all three superpopulations are available
    if not all(sp in superpop_unions for sp in ['AFR', 'EAS', 'EUR']):
        print("  WARNING: Not all required superpopulations (AFR, EAS, EUR) available")
        return None

    # Get gene sets (union of all subpopulations)
    afr_genes = set(superpop_unions['AFR'].keys())
    eas_genes = set(superpop_unions['EAS'].keys())
    eur_genes = set(superpop_unions['EUR'].keys())

    # Calculate overlaps
    all_three = afr_genes & eas_genes & eur_genes
    eur_eas_only = (eur_genes & eas_genes) - afr_genes
    eur_afr_only = (eur_genes & afr_genes) - eas_genes
    eas_afr_only = (eas_genes & afr_genes) - eur_genes
    eur_unique = eur_genes - eas_genes - afr_genes
    eas_unique = eas_genes - eur_genes - afr_genes
    afr_unique = afr_genes - eur_genes - eas_genes

    # Create larger figure to accommodate annotations
    fig, ax = plt.subplots(figsize=(20, 16))

    venn = venn3(
        [eur_genes, eas_genes, afr_genes],
        set_labels=('EUR', 'EAS', 'AFR'),
        set_colors=('#3498db', '#e74c3c', '#2ecc71'),
        alpha=0.6,
        ax=ax
    )
    venn3_circles(
        [eur_genes, eas_genes, afr_genes],
        linewidth=2,
        ax=ax
    )

    plt.title('Overlap of Genes Under Selection\nAcross EUR, EAS, and AFR (Union of All Subpopulations)',
              fontsize=18, fontweight='bold', pad=30)

    # Add gene names with text boxes and arrows
    # Each tuple: (gene_set, label, text_position, arrow_target, horizontal_align, vertical_align, arrow_curve)
    annotations = [
        # EUR only - top left
        (eur_unique, 'EUR only', (-0.72, 0.30), (-0.28, 0.32), 'right', 'bottom', 'arc3,rad=0.2'),
        # EAS only - top right
        (eas_unique, 'EAS only', (0.72, 0.60), (0.28, 0.32), 'left', 'bottom', 'arc3,rad=-0.2'),
        # AFR only - bottom center (moved lower to avoid overlap)
        (afr_unique, 'AFR only', (-0.20, -0.80), (0.0, -0.35), 'center', 'top', 'arc3,rad=0'),
        # EUR & EAS only - top center
        (eur_eas_only, 'EUR & EAS', (0.0, 0.80), (-0.10, 0.38), 'center', 'bottom', 'arc3,rad=0'),
        # EUR & AFR only - far left
        (eur_afr_only, 'EUR & AFR', (-0.72, 0.10), (-0.30, 0.20), 'right', 'center', 'arc3,rad=-0.2'),
        # EAS & AFR only - far right
        (eas_afr_only, 'EAS & AFR', (0.82, -0.15), (0.15, -0.05), 'left', 'center', 'arc3,rad=0.2'),
        # All 3 populations - TOP CENTER, arrow pointing straight down to center
        (all_three, 'All 3 continents', (-0.50, 0.75), (-0.15, 0.24), 'center', 'bottom', 'arc3,rad=0'),
    ]

    for gene_set, label, text_pos, arrow_pos, ha, va, connection_style in annotations:
        if gene_set:  # Only annotate if there are genes in this region
            sorted_genes = sorted(gene_set)

            # For AFR only category, show all genes in multiple columns if needed
            if label == 'AFR only':
                # Display all genes in columns
                n_genes = len(sorted_genes)

                if n_genes <= 25:
                    # Single column if 25 or fewer genes
                    gene_text = '\n'.join(sorted_genes)
                else:
                    # Multiple columns for more than 25 genes
                    # Determine number of columns needed (aim for ~20-25 genes per column)
                    genes_per_col = 20
                    n_cols = (n_genes + genes_per_col - 1) // genes_per_col

                    # Create columns
                    columns = []
                    for col_idx in range(n_cols):
                        start_idx = col_idx * genes_per_col
                        end_idx = min(start_idx + genes_per_col, n_genes)
                        columns.append(sorted_genes[start_idx:end_idx])

                    # Format as side-by-side columns with padding
                    max_rows = max(len(col) for col in columns)
                    lines = []
                    for row_idx in range(max_rows):
                        row_parts = []
                        for col in columns:
                            if row_idx < len(col):
                                row_parts.append(col[row_idx].ljust(20))
                            else:
                                row_parts.append(' ' * 20)
                        lines.append('  '.join(row_parts))

                    gene_text = '\n'.join(lines)
            else:
                # For other categories, limit to first 25 genes
                genes_to_show = sorted_genes[:25]
                gene_text = '\n'.join(genes_to_show)

                if len(gene_set) > 25:
                    gene_text += f'\n... +{len(gene_set) - 25} more'

            full_text = f"{label} (n={len(gene_set)}):\n{gene_text}"

            # Add text box with arrow
            ax.annotate(
                full_text,
                xy=arrow_pos,           # Arrow points to this location (inside Venn region)
                xytext=text_pos,        # Text box location
                fontsize=11,
                family='monospace',
                ha=ha,
                va=va,
                bbox=dict(
                    boxstyle='round,pad=0.5',
                    facecolor='wheat',
                    alpha=0.9,
                    edgecolor='#8B4513',
                    linewidth=1.5
                ),
                arrowprops=dict(
                    arrowstyle='-|>',
                    connectionstyle=connection_style,
                    color='#555555',
                    lw=1.5,
                    mutation_scale=12
                )
            )

    # Adjust axis limits to accommodate all annotations
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1.15, 0.95)
    ax.axis('off')

    output_file = output_dir / "venn_3way_AFR_EAS_EUR.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.savefig(output_file.with_suffix('.pdf'), bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()

    print(f"  Saved: {output_file.name}")

    # Also save gene lists for each category
    categories = {
        'all_three_continents': sorted(all_three),
        'EUR_EAS_only': sorted(eur_eas_only),
        'EUR_AFR_only': sorted(eur_afr_only),
        'EAS_AFR_only': sorted(eas_afr_only),
        'EUR_unique': sorted(eur_unique),
        'EAS_unique': sorted(eas_unique),
        'AFR_unique': sorted(afr_unique)
    }

    for category_name, genes in categories.items():
        if genes:
            gene_list_file = output_dir / f"AFR_EAS_EUR_{category_name}.txt"
            with open(gene_list_file, 'w') as f:
                f.write('\n'.join(genes))
            print(f"    {category_name}: {len(genes)} genes -> {gene_list_file.name}")

    # Print summary
    print()
    print(f"  3-way Venn summary (AFR, EAS, EUR):")
    print(f"    All 3 continents: {len(all_three)}")
    print(f"    EUR & EAS only: {len(eur_eas_only)}")
    print(f"    EUR & AFR only: {len(eur_afr_only)}")
    print(f"    EAS & AFR only: {len(eas_afr_only)}")
    print(f"    EUR unique: {len(eur_unique)}")
    print(f"    EAS unique: {len(eas_unique)}")
    print(f"    AFR unique: {len(afr_unique)}")

    return output_file.name


# ============================================================================
# Main Processing
# ============================================================================

def main():
    # Step 1: Create superpopulation unions with provenance
    superpop_unions = create_superpopulation_unions()

    if len(superpop_unions) == 0:
        print("ERROR: No data available for any superpopulation")
        return

    print("=" * 80)
    print()

    # Step 2: Save superpopulation union tables
    print("Saving superpopulation gene union tables...")
    print()

    for superpop, gene_provenance in superpop_unions.items():
        union_data = []
        total_subpops = len(SUPERPOPULATIONS[superpop])

        for gene, subpops in sorted(gene_provenance.items()):
            n_subpops = len(subpops)
            breadth = classify_signal_breadth(n_subpops, total_subpops)

            union_data.append({
                'gene_name': gene,
                'n_subpopulations': n_subpops,
                'total_subpopulations': total_subpops,
                'proportion': n_subpops / total_subpops,
                'signal_breadth': breadth,
                'subpopulations': ','.join(sorted(subpops))
            })

        union_df = pd.DataFrame(union_data)
        union_df = union_df.sort_values(['n_subpopulations', 'gene_name'],
                                       ascending=[False, True])

        union_file = output_base / f"{superpop}_gene_union_with_provenance.tsv"
        union_df.to_csv(union_file, sep='\t', index=False)
        print(f"  {superpop}: {len(union_df)} genes → {union_file.name}")

    print()

    # Step 3: Create provenance summary table
    print("Creating provenance summary table...")
    summary_df = create_provenance_summary_table(superpop_unions)
    summary_file = output_base / "superpopulation_summary.tsv"
    summary_df.to_csv(summary_file, sep='\t', index=False)
    print(f"  Saved: {summary_file.name}")
    print()
    print(summary_df.to_string(index=False))
    print()
    print("=" * 80)
    print()

    # Step 4: Pairwise comparisons
    print("Performing pairwise continental comparisons...")
    print()

    superpops = list(superpop_unions.keys())
    pairwise_summary = []

    for superpop1, superpop2 in combinations(superpops, 2):
        print(f"{superpop1} vs {superpop2}:")

        # Analyze overlap
        overlap_df = analyze_pairwise_overlap(superpop_unions, superpop1, superpop2)

        # Save detailed results
        output_file = output_base / f"comparison_{superpop1}_vs_{superpop2}.tsv"
        overlap_df.to_csv(output_file, sep='\t', index=False)
        print(f"  Saved: {output_file.name}")

        # Count categories
        category_counts = overlap_df['category'].value_counts()
        n_shared = category_counts.get('shared', 0)

        print(f"  Shared genes: {n_shared}")

        if n_shared > 0:
            # Pattern breakdown
            shared_df = overlap_df[overlap_df['category'] == 'shared']
            pattern_counts = shared_df['pattern'].value_counts()
            for pattern, count in pattern_counts.items():
                print(f"    {pattern}: {count}")

            # Create pattern breakdown plot
            pattern_plot = create_pattern_breakdown_plot(overlap_df, superpop1,
                                                        superpop2, output_base)
            if pattern_plot:
                print(f"  Saved: {pattern_plot}")

        # Create Venn diagram
        venn_file = create_pairwise_venn(superpop_unions, superpop1, superpop2,
                                        output_base)
        print(f"  Saved: {venn_file}")

        # Summary statistics
        pairwise_summary.append({
            'comparison': f'{superpop1}_vs_{superpop2}',
            'superpop1': superpop1,
            'superpop2': superpop2,
            'genes_superpop1': len(superpop_unions[superpop1]),
            'genes_superpop2': len(superpop_unions[superpop2]),
            'shared_genes': n_shared,
            'superpop1_specific': category_counts.get(f'{superpop1}_specific', 0),
            'superpop2_specific': category_counts.get(f'{superpop2}_specific', 0),
            'jaccard_similarity': n_shared / len(overlap_df) if len(overlap_df) > 0 else 0
        })

        print()

    # Save pairwise summary
    pairwise_df = pd.DataFrame(pairwise_summary)
    pairwise_file = output_base / "pairwise_comparison_summary.tsv"
    pairwise_df.to_csv(pairwise_file, sep='\t', index=False)
    print(f"Pairwise summary saved: {pairwise_file.name}")
    print()

    # Step 5: Create overall heatmap
    print("Creating overall overlap heatmap...")
    overlap_matrix, jaccard_matrix, heatmap_file = create_overlap_heatmap(
        superpop_unions, output_base
    )
    print(f"  Saved: {heatmap_file}")
    print()

    # Save matrices as tables
    superpops = list(superpop_unions.keys())
    overlap_df = pd.DataFrame(overlap_matrix, index=superpops, columns=superpops)
    jaccard_df = pd.DataFrame(jaccard_matrix, index=superpops, columns=superpops)

    overlap_table_file = output_base / "overlap_matrix_counts.tsv"
    jaccard_table_file = output_base / "overlap_matrix_jaccard.tsv"

    overlap_df.to_csv(overlap_table_file, sep='\t')
    jaccard_df.to_csv(jaccard_table_file, sep='\t')

    print(f"  Saved: {overlap_table_file.name}")
    print(f"  Saved: {jaccard_table_file.name}")
    print()

    # Step 6: Multi-way overlap analysis (all continents)
    print("=" * 80)
    print()
    print("Analyzing multi-way overlap across ALL continents...")
    print()

    multiway_df = analyze_multiway_overlap(superpop_unions)

    # Save detailed table
    multiway_file = output_base / "multiway_continent_overlap_detailed.tsv"
    multiway_df.to_csv(multiway_file, sep='\t', index=False)
    print(f"  Saved: {multiway_file.name}")
    print()

    # Print summary statistics
    print("Multi-way overlap summary:")
    for n_cont in sorted(multiway_df['n_continents'].unique(), reverse=True):
        n_genes = len(multiway_df[multiway_df['n_continents'] == n_cont])
        pct = n_genes / len(multiway_df) * 100
        if n_cont == 5:
            print(f"  Genes in ALL 5 continents (universal): {n_genes} ({pct:.1f}%)")
        elif n_cont == 1:
            print(f"  Genes in exactly 1 continent (continent-specific): {n_genes} ({pct:.1f}%)")
        else:
            print(f"  Genes in exactly {n_cont} continents: {n_genes} ({pct:.1f}%)")
    print()

    # Save gene lists by number of continents
    for n_cont in sorted(multiway_df['n_continents'].unique(), reverse=True):
        genes_at_n = multiway_df[multiway_df['n_continents'] == n_cont]
        if len(genes_at_n) > 0:
            gene_list_file = output_base / f"genes_in_{n_cont}_continents.txt"
            with open(gene_list_file, 'w') as f:
                f.write('\n'.join(sorted(genes_at_n['gene_name'].tolist())))
            print(f"  Saved: {gene_list_file.name} ({len(genes_at_n)} genes)")
    print()

    # Create visualizations
    print("Creating multi-way overlap visualizations...")
    create_multiway_upset_plot(multiway_df, output_base)
    create_continent_combination_analysis(multiway_df, output_base)
    print()

    # Step 7: Create 3-way Venn diagram for AFR, EAS, EUR
    print("Creating 3-way Venn diagram for AFR, EAS, and EUR...")
    print()
    create_threeway_venn_afr_eas_eur(superpop_unions, output_base)
    print()

    # Final summary
    print("=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print()
    print(f"All results saved to: {output_base}")
    print()
    print("Key outputs:")
    print("  - Superpopulation gene unions with provenance (per continent)")
    print("  - Pairwise comparison tables (detailed gene-level data)")
    print("  - Pattern breakdown plots (broad-broad, broad-narrow, etc.)")
    print("  - Venn diagrams (pairwise comparisons with gene annotations)")
    print("  - 3-way Venn diagram (AFR, EAS, EUR with gene annotations)")
    print("  - Overlap heatmaps (all continents)")
    print("  - Multi-way overlap analysis (genes in 1-5 continents)")
    print("  - Continent combination analysis (which combinations are common)")
    print("  - Summary statistics tables")
    print()


if __name__ == "__main__":
    main()
