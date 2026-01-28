#!/usr/bin/env python3
"""
Comprehensive analysis of genes identified in selection regions across subpopulations.
This script was made to get comprehensive visualizations and statistics to include in my project report.
Got some help from Claude to be able to combine data from all subpopulations and do the various analyses.
Author: Mit Van Bruggen
This script:
1. Counts total unique genes across all subpopulations
2. Provides breakdown by gene biotype
3. Calculates proportion of X chromosome genes
4. Analyzes X-inactivation status (escapees vs subjects)
5. Identifies which genes are escapees
6. Tests for enrichment of escapee genes in selection regions
7. Identifies genes with strongest selection evidence based on |Std.iHS| scores
8. Visualizes top protein-coding genes along the X chromosome
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import fisher_exact
import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def find_genes_files(base_dir):
    """Find all genes_in_regions.tsv files."""
    base_path = Path(base_dir)
    files = list(base_path.glob("**/genes_in_regions.tsv"))
    return sorted(files)


def load_all_genes(files):
    """Load and combine all gene files."""
    dfs = []
    for file in files:
        # Extract population name from path
        pop = file.parent.parent.parent.name.replace("ihs_", "")
        df = pd.read_csv(file, sep="\t")
        df['population'] = pop
        dfs.append(df)

    combined = pd.concat(dfs, ignore_index=True)
    return combined


def calculate_basic_statistics(df):
    """Calculate basic gene statistics."""
    print("=" * 80)
    print("BASIC GENE STATISTICS")
    print("=" * 80)

    # Total unique genes
    total_genes = df['gene_id'].nunique()
    print(f"\nTotal unique genes identified: {total_genes:,}")

    # By gene name
    total_gene_names = df['gene_name'].nunique()
    print(f"Total unique gene names: {total_gene_names:,}")

    # Genes per population
    print(f"\nGenes per population:")
    genes_per_pop = df.groupby('population')['gene_id'].nunique().sort_values(ascending=False)
    for pop, count in genes_per_pop.items():
        print(f"  {pop}: {count:,}")

    return total_genes


def analyze_gene_biotypes(df):
    """Analyze gene biotype distribution."""
    print("\n" + "=" * 80)
    print("GENE BIOTYPE BREAKDOWN")
    print("=" * 80)

    # Get unique genes (some may appear multiple times across populations/regions)
    unique_genes = df.drop_duplicates(subset=['gene_id'])

    biotype_counts = unique_genes['gene_type'].value_counts()
    total = len(unique_genes)

    print(f"\nTotal unique genes: {total:,}\n")
    print(f"{'Biotype':<30} {'Count':>10} {'Percentage':>12}")
    print("-" * 54)

    for biotype, count in biotype_counts.items():
        pct = (count / total) * 100
        print(f"{biotype:<30} {count:>10,} {pct:>11.2f}%")

    return biotype_counts


def analyze_x_chromosome(df):
    """Analyze X chromosome gene representation."""
    print("\n" + "=" * 80)
    print("X CHROMOSOME ANALYSIS")
    print("=" * 80)

    # Get unique genes
    unique_genes = df.drop_duplicates(subset=['gene_id'])

    # Count X chromosome genes
    x_genes = unique_genes[unique_genes['chr'] == 'X']
    total_genes = len(unique_genes)
    x_count = len(x_genes)

    print(f"\nTotal unique genes: {total_genes:,}")
    print(f"X chromosome genes: {x_count:,}")
    print(f"Proportion on X: {(x_count/total_genes)*100:.2f}%")

    # Autosomal genes
    autosomal = unique_genes[unique_genes['chr'] != 'X']
    print(f"Autosomal genes: {len(autosomal):,}")

    # Breakdown by chromosome
    print("\nTop 10 chromosomes by gene count:")
    chr_counts = unique_genes['chr'].value_counts().head(10)
    for chrom, count in chr_counts.items():
        print(f"  Chr {chrom}: {count:,}")

    return x_genes, x_count, total_genes


def analyze_x_inactivation(df):
    """Analyze X-inactivation status of X chromosome genes."""
    print("\n" + "=" * 80)
    print("X-INACTIVATION STATUS ANALYSIS")
    print("=" * 80)

    # Filter for X chromosome genes only
    x_genes = df[df['chr'] == 'X'].copy()

    if len(x_genes) == 0:
        print("\nNo X chromosome genes found!")
        return None, None

    # Get unique X genes
    unique_x = x_genes.drop_duplicates(subset=['gene_id'])

    # Count by X-inactivation status
    xi_status = unique_x['x_inactivation_status'].value_counts()
    total_x = len(unique_x)

    print(f"\nTotal unique X chromosome genes: {total_x:,}\n")
    print(f"{'Status':<20} {'Count':>10} {'Percentage':>12}")
    print("-" * 44)

    for status, count in xi_status.items():
        pct = (count / total_x) * 100
        print(f"{status:<20} {count:>10,} {pct:>11.2f}%")

    # Get escapee genes
    escapees = unique_x[unique_x['x_inactivation_status'] == 'escapee']
    subjects = unique_x[unique_x['x_inactivation_status'] == 'subject']

    n_escapees = len(escapees)
    n_subjects = len(subjects)

    print(f"\nSummary:")
    print(f"  Escapee genes: {n_escapees:,}")
    print(f"  Subject genes: {n_subjects:,}")

    return escapees, subjects


def list_escapee_genes(escapees):
    """List all escapee genes identified."""
    print("\n" + "=" * 80)
    print("ESCAPEE GENES IDENTIFIED")
    print("=" * 80)

    if escapees is None or len(escapees) == 0:
        print("\nNo escapee genes found!")
        return

    # Sort by gene name
    escapee_list = escapees[['gene_name', 'gene_id', 'gene_type']].sort_values('gene_name')

    print(f"\nTotal escapee genes: {len(escapee_list):,}\n")
    print(f"{'Gene Name':<20} {'Gene ID':<20} {'Biotype':<30}")
    print("-" * 72)

    for _, row in escapee_list.iterrows():
        print(f"{row['gene_name']:<20} {row['gene_id']:<20} {row['gene_type']:<30}")

    return escapee_list


def test_escapee_enrichment(df):
    """Test for enrichment of escapee genes in selection regions using Fisher's exact test."""
    print("\n" + "=" * 80)
    print("ESCAPEE GENE ENRICHMENT ANALYSIS")
    print("=" * 80)

    # For this analysis, we need a background set of all X chromosome genes
    # We'll compare the proportion of escapees in selection regions vs. expected

    x_genes = df[df['chr'] == 'X'].drop_duplicates(subset=['gene_id'])

    if len(x_genes) == 0:
        print("\nNo X chromosome genes to analyze!")
        return

    # Count escapees and subjects in selection regions
    in_selection_escapees = len(x_genes[x_genes['x_inactivation_status'] == 'escapee'])
    in_selection_subjects = len(x_genes[x_genes['x_inactivation_status'] == 'subject'])
    in_selection_total = in_selection_escapees + in_selection_subjects

    print("\nGenes in selection regions (X chromosome only):")
    print(f"  Escapees: {in_selection_escapees:,}")
    print(f"  Subjects: {in_selection_subjects:,}")
    print(f"  Total: {in_selection_total:,}")
    print(f"  Escapee proportion: {(in_selection_escapees/in_selection_total)*100:.2f}%")

    # Note: For a proper enrichment test, we would need the background distribution
    # of all X chromosome genes (not just those in selection regions)
    # This would typically come from a genome annotation file

    print("\n" + "-" * 80)
    print("NOTE: For formal enrichment testing, we would need:")
    print("  1. Total number of escapee genes on X chromosome (genome-wide)")
    print("  2. Total number of subject genes on X chromosome (genome-wide)")
    print("\nWithout this background, we can only report the observed proportions")
    print("in selection regions, which is shown above.")
    print("-" * 80)

    return in_selection_escapees, in_selection_subjects


def identify_strongest_selection_genes(df, top_n=50):
    """Identify genes with strongest selection evidence based on |Std.iHS| scores."""
    print("\n" + "=" * 80)
    print(f"TOP {top_n} GENES WITH STRONGEST SELECTION EVIDENCE")
    print("=" * 80)

    # For each gene, get the maximum |Std.iHS| score from overlapping SNPs
    # This is stored in region_max_std_ihs column

    # Get unique gene-region combinations with their max iHS scores
    gene_scores = df[['gene_name', 'gene_id', 'gene_type', 'chr',
                      'x_inactivation_status', 'population',
                      'region_max_std_ihs']].copy()

    # Get the maximum score for each gene across all populations/regions
    max_scores = gene_scores.groupby('gene_id').agg({
        'gene_name': 'first',
        'gene_type': 'first',
        'chr': 'first',
        'x_inactivation_status': 'first',
        'region_max_std_ihs': 'max',
        'population': lambda x: ','.join(sorted(set(x)))
    }).reset_index()

    # Sort by absolute value of iHS score
    max_scores['abs_ihs'] = max_scores['region_max_std_ihs'].abs()
    top_genes = max_scores.nlargest(top_n, 'abs_ihs')

    print(f"\n{'Rank':<6} {'Gene':<15} {'Chr':<5} {'Biotype':<25} {'Max|iHS|':<10} {'XI Status':<10} {'Populations'}")
    print("-" * 130)

    for i, (_, row) in enumerate(top_genes.iterrows(), 1):
        xi_status = row['x_inactivation_status'] if pd.notna(row['x_inactivation_status']) else 'N/A'
        # Truncate population list if too long
        pops = row['population']
        if len(pops) > 40:
            pops = pops[:37] + "..."

        print(f"{i:<6} {row['gene_name']:<15} {row['chr']:<5} {row['gene_type']:<25} "
              f"{row['abs_ihs']:<10.3f} {xi_status:<10} {pops}")

    return top_genes


def plot_top_x_chromosome_genes(df, top_n=15, output_file=None):
    """
    Plot the top protein-coding genes along the X chromosome based on selection evidence.
    """
    print("\n" + "=" * 80)
    print(f"VISUALIZING TOP {top_n} PROTEIN-CODING GENES ON X CHROMOSOME")
    print("=" * 80)

    # Filter for X chromosome protein-coding genes
    x_protein = df[(df['chr'] == 'X') & (df['gene_type'] == 'protein_coding')].copy()

    if len(x_protein) == 0:
        print("\nNo X chromosome protein-coding genes found!")
        return None, None

    # Get the gene with max position info and max iHS score
    gene_info = x_protein.groupby('gene_id').agg({
        'gene_name': 'first',
        'start': 'first',
        'end': 'first',
        'x_inactivation_status': 'first',
        'region_max_std_ihs': 'max'
    }).reset_index()

    gene_info['abs_ihs'] = gene_info['region_max_std_ihs'].abs()
    gene_info['midpoint'] = (gene_info['start'] + gene_info['end']) / 2

    # Get top N genes
    top_genes = gene_info.nlargest(top_n, 'abs_ihs')

    print(f"\nTop {top_n} protein-coding genes on X chromosome:\n")
    print(f"{'Rank':<6} {'Gene':<15} {'Position (Mb)':<20} {'Max|iHS|':<12} {'XI Status':<12}")
    print("-" * 70)

    for i, (_, row) in enumerate(top_genes.iterrows(), 1):
        pos_mb = row['midpoint'] / 1e6
        xi_status = row['x_inactivation_status'] if pd.notna(row['x_inactivation_status']) else 'N/A'
        print(f"{i:<6} {row['gene_name']:<15} {pos_mb:<20.2f} {row['abs_ihs']:<12.3f} {xi_status:<12}")

    # Create visualization
    fig, ax = plt.subplots(figsize=(14, 8))

    # X chromosome length (approximate)
    x_chr_length = 155270560  # bp

    # Draw chromosome ideogram
    ax.add_patch(plt.Rectangle((0, 0), x_chr_length/1e6, 1,
                               facecolor='lightgray', edgecolor='black', linewidth=2))

    # Plot genes
    colors = {'escapee': '#E74C3C', 'subject': '#3498DB'}

    for i, (_, row) in enumerate(top_genes.iterrows()):
        xi_status = row['x_inactivation_status'] if pd.notna(row['x_inactivation_status']) else 'unknown'
        color = colors.get(xi_status, '#95A5A6')

        # Gene position in Mb
        start_mb = row['start'] / 1e6
        end_mb = row['end'] / 1e6
        midpoint_mb = row['midpoint'] / 1e6
        width_mb = end_mb - start_mb

        # Draw gene as a rectangle
        y_pos = 2 + (i * 0.8)
        rect_height = 0.4
        ax.add_patch(plt.Rectangle((start_mb, y_pos - rect_height/2), width_mb, rect_height,
                                    facecolor=color, edgecolor='black', linewidth=1.5, alpha=0.8))

        # Add gene label
        ax.text(midpoint_mb, y_pos + 0.3, row['gene_name'],
               fontsize=9, ha='center', fontweight='bold')

        # Add iHS score
        ax.text(x_chr_length/1e6 + 5, y_pos, f"|iHS|={row['abs_ihs']:.2f}",
               fontsize=8, va='center')

    # Formatting
    ax.set_xlim(-5, x_chr_length/1e6 + 20)
    ax.set_ylim(-0.5, 2 + (top_n * 0.8) + 1)
    ax.set_xlabel('Position on X Chromosome (Mb)', fontsize=12, fontweight='bold')
    ax.set_title(f'Top {top_n} Protein-Coding Genes Under Selection on X Chromosome',
                fontsize=14, fontweight='bold', pad=20)

    # Remove y-axis
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Add legend
    legend_elements = [
        mpatches.Patch(facecolor=colors['escapee'], edgecolor='black', label='X-inactivation Escapee'),
        mpatches.Patch(facecolor=colors['subject'], edgecolor='black', label='X-inactivation Subject')
    ]
    ax.legend(handles=legend_elements, loc='upper left', frameon=True, fontsize=10)

    # Add grid for x-axis
    ax.grid(axis='x', alpha=0.3, linestyle='--')

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"\nPlot saved to: {output_file}")

    plt.close()

    return fig, top_genes


def generate_summary_stats(df, output_file=None):
    """Generate summary statistics and optionally save to file."""
    stats = {
        'total_unique_genes': df['gene_id'].nunique(),
        'total_unique_gene_names': df['gene_name'].nunique(),
        'total_populations': df['population'].nunique(),
        'total_observations': len(df),
    }

    # X chromosome stats
    unique_genes = df.drop_duplicates(subset=['gene_id'])
    x_genes = unique_genes[unique_genes['chr'] == 'X']
    stats['x_chromosome_genes'] = len(x_genes)
    stats['x_chromosome_proportion'] = len(x_genes) / len(unique_genes)

    # X-inactivation stats
    if len(x_genes) > 0:
        stats['escapee_genes'] = len(x_genes[x_genes['x_inactivation_status'] == 'escapee'])
        stats['subject_genes'] = len(x_genes[x_genes['x_inactivation_status'] == 'subject'])
        stats['escapee_proportion'] = stats['escapee_genes'] / len(x_genes)

    if output_file:
        with open(output_file, 'w') as f:
            f.write("Summary Statistics\n")
            f.write("=" * 50 + "\n")
            for key, value in stats.items():
                if isinstance(value, float):
                    f.write(f"{key}: {value:.4f}\n")
                else:
                    f.write(f"{key}: {value}\n")

    return stats


def main():
    parser = argparse.ArgumentParser(
        description="Analyze genes identified in selection regions across subpopulations"
    )
    parser.add_argument(
        "--base-dir",
        default="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results",
        help="Base directory containing ihs_* subdirectories"
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=50,
        help="Number of top genes to report (default: 50)"
    )
    parser.add_argument(
        "--top-plot",
        type=int,
        default=15,
        help="Number of top X chromosome genes to plot (default: 15)"
    )
    parser.add_argument(
        "--output-dir",
        help="Directory to save output files (optional)"
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Skip generating plots"
    )

    args = parser.parse_args()

    # Find all gene files
    print("Searching for genes_in_regions.tsv files...")
    files = find_genes_files(args.base_dir)
    print(f"Found {len(files)} files\n")

    # Load all data
    print("Loading data...")
    df = load_all_genes(files)
    print(f"Loaded {len(df):,} total gene-region observations\n")

    # Run analyses
    calculate_basic_statistics(df)
    analyze_gene_biotypes(df)
    x_genes, x_count, total_genes = analyze_x_chromosome(df)
    escapees, subjects = analyze_x_inactivation(df)
    list_escapee_genes(escapees)
    test_escapee_enrichment(df)
    top_genes = identify_strongest_selection_genes(df, top_n=args.top_n)

    # Plot top X chromosome genes
    if not args.no_plot:
        plot_file = None
        if args.output_dir:
            output_path = Path(args.output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
            plot_file = output_path / f"top_{args.top_plot}_x_genes_plot.png"

        fig, top_x_genes = plot_top_x_chromosome_genes(df, top_n=args.top_plot, output_file=plot_file)

    # Generate summary stats
    if args.output_dir:
        output_path = Path(args.output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Save summary stats
        stats_file = output_path / "summary_statistics.txt"
        generate_summary_stats(df, stats_file)
        print(f"\nSummary statistics saved to: {stats_file}")

        # Save escapee genes list
        if escapees is not None and len(escapees) > 0:
            escapee_file = output_path / "escapee_genes.tsv"
            escapees[['gene_name', 'gene_id', 'gene_type', 'chr', 'start', 'end']].drop_duplicates().to_csv(
                escapee_file, sep='\t', index=False
            )
            print(f"Escapee genes saved to: {escapee_file}")

        # Save top genes
        top_genes_file = output_path / f"top_{args.top_n}_selection_genes.tsv"
        top_genes.to_csv(top_genes_file, sep='\t', index=False)
        print(f"Top genes saved to: {top_genes_file}")

    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
