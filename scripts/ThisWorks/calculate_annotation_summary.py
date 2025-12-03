#!/usr/bin/env python3
"""
Calculate and save summary statistics for gene annotation results
"""

import pandas as pd
from pathlib import Path

# Setup paths
base_dir = Path("/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib")
annotation_dir = base_dir / "results/analysis/gene_annotation"

# Load data
print("Loading annotation data...")
annotated_regions = pd.read_csv(annotation_dir / "selection_regions_annotated.tsv", sep='\t')
genes_in_regions = pd.read_csv(annotation_dir / "genes_in_selection_regions_detailed.tsv", sep='\t')

# Get unique genes
all_genes = genes_in_regions.drop_duplicates(subset=['gene_id'])

# Calculate statistics
print("Calculating statistics...")

summary_stats = {
    # Region statistics
    'total_regions': len(annotated_regions),
    'regions_with_genes': len(annotated_regions[annotated_regions['n_genes_total'] > 0]),
    'regions_with_genes_percent': len(annotated_regions[annotated_regions['n_genes_total'] > 0]) / len(annotated_regions) * 100,
    'mean_genes_per_region': annotated_regions['n_genes_total'].mean(),
    'median_genes_per_region': annotated_regions['n_genes_total'].median(),
    'mean_contained_genes': annotated_regions['n_genes_contained'].mean(),
    'mean_overlapping_genes': annotated_regions['n_genes_overlapping'].mean(),
    'mean_flanking_genes': annotated_regions['n_genes_flanking'].mean(),

    # Gene statistics
    'total_unique_genes': len(all_genes),
    'x_chromosome_genes': len(all_genes[all_genes['chr'] == 'X']),
    'autosomal_genes': len(all_genes[all_genes['chr'] != 'X']),
    'protein_coding_genes': len(all_genes[all_genes['gene_type'] == 'protein_coding']),
    'protein_coding_percent': len(all_genes[all_genes['gene_type'] == 'protein_coding']) / len(all_genes) * 100,
    'x_protein_coding': len(all_genes[(all_genes['chr'] == 'X') & (all_genes['gene_type'] == 'protein_coding')]),
}

# X-inactivation statistics
x_genes = all_genes[all_genes['chr'] == 'X']
if len(x_genes) > 0:
    escapees = x_genes[x_genes['x_inactivation_status'] == 'escape']
    summary_stats['x_inactivation_escapees'] = len(escapees)
    summary_stats['x_inactivation_escapees_percent'] = len(escapees) / len(x_genes) * 100
    summary_stats['x_inactivation_subject'] = len(x_genes) - len(escapees)

# Gene type breakdown
gene_type_counts = all_genes['gene_type'].value_counts()

# Chromosome breakdown
autosome_regions = annotated_regions[annotated_regions['chr_type'] == 'Autosome']
x_regions = annotated_regions[annotated_regions['chr_type'] == 'X chromosome']

summary_stats['autosome_regions'] = len(autosome_regions)
summary_stats['x_chromosome_regions'] = len(x_regions)
summary_stats['autosome_regions_with_genes'] = len(autosome_regions[autosome_regions['n_genes_total'] > 0])
summary_stats['x_regions_with_genes'] = len(x_regions[x_regions['n_genes_total'] > 0])

# Save summary statistics
print("Saving summary statistics...")

# Main summary file
summary_df = pd.DataFrame([summary_stats]).T
summary_df.columns = ['value']
summary_df.to_csv(annotation_dir / "annotation_summary_statistics.tsv", sep='\t', header=True)
print(f"  Saved: annotation_summary_statistics.tsv")

# Gene type breakdown
gene_type_df = pd.DataFrame({
    'gene_type': gene_type_counts.index,
    'count': gene_type_counts.values,
    'percent': (gene_type_counts.values / len(all_genes) * 100)
})
gene_type_df.to_csv(annotation_dir / "gene_type_breakdown.tsv", sep='\t', index=False)
print(f"  Saved: gene_type_breakdown.tsv")

# Top regions by gene count
top_regions = annotated_regions.nlargest(20, 'n_genes_total')[
    ['chr', 'start', 'end', 'length', 'n_snps', 'n_genes_total',
     'n_genes_contained', 'n_genes_overlapping', 'n_genes_flanking',
     'max_std_ihs', 'protein_coding_genes', 'chr_type']
]
top_regions.to_csv(annotation_dir / "top_regions_by_gene_count.tsv", sep='\t', index=False)
print(f"  Saved: top_regions_by_gene_count.tsv")

# Top regions by iHS score
top_regions_ihs = annotated_regions.nlargest(20, 'max_std_ihs')[
    ['chr', 'start', 'end', 'length', 'n_snps', 'n_genes_total',
     'n_genes_contained', 'max_std_ihs', 'protein_coding_genes', 'chr_type']
]
top_regions_ihs.to_csv(annotation_dir / "top_regions_by_ihs_score.tsv", sep='\t', index=False)
print(f"  Saved: top_regions_by_ihs_score.tsv")

# X chromosome specific summary
if len(x_regions) > 0:
    x_summary = {
        'x_regions_total': len(x_regions),
        'x_regions_with_genes': len(x_regions[x_regions['n_genes_total'] > 0]),
        'x_mean_genes_per_region': x_regions['n_genes_total'].mean(),
        'x_median_genes_per_region': x_regions['n_genes_total'].median(),
        'x_total_unique_genes': len(x_genes),
        'x_protein_coding_genes': len(x_genes[x_genes['gene_type'] == 'protein_coding']),
        'x_escapee_genes': len(escapees) if len(x_genes) > 0 else 0,
        'x_subject_genes': len(x_genes) - len(escapees) if len(x_genes) > 0 else 0,
    }

    x_summary_df = pd.DataFrame([x_summary]).T
    x_summary_df.columns = ['value']
    x_summary_df.to_csv(annotation_dir / "X_chromosome_summary.tsv", sep='\t', header=True)
    print(f"  Saved: X_chromosome_summary.tsv")

    # Top X regions
    top_x_regions = x_regions.nlargest(10, 'max_std_ihs')[
        ['chr', 'start', 'end', 'length', 'n_snps', 'n_genes_total',
         'max_std_ihs', 'protein_coding_genes']
    ]
    top_x_regions.to_csv(annotation_dir / "top_X_regions.tsv", sep='\t', index=False)
    print(f"  Saved: top_X_regions.tsv")

# Overlap type statistics
overlap_counts = genes_in_regions['overlap_type'].value_counts()
overlap_df = pd.DataFrame({
    'overlap_type': overlap_counts.index,
    'count': overlap_counts.values,
    'percent': (overlap_counts.values / len(genes_in_regions) * 100)
})
overlap_df.to_csv(annotation_dir / "overlap_type_distribution.tsv", sep='\t', index=False)
print(f"  Saved: overlap_type_distribution.tsv")

print("\nSummary statistics saved successfully!")
print(f"\nKey Results:")
print(f"  Total regions: {summary_stats['total_regions']:,}")
print(f"  Regions with genes: {summary_stats['regions_with_genes']:,} ({summary_stats['regions_with_genes_percent']:.1f}%)")
print(f"  Total unique genes: {summary_stats['total_unique_genes']:,}")
print(f"  Protein-coding genes: {summary_stats['protein_coding_genes']:,} ({summary_stats['protein_coding_percent']:.1f}%)")
print(f"  X chromosome genes: {summary_stats['x_chromosome_genes']:,}")
print(f"  X protein-coding: {summary_stats['x_protein_coding']:,}")
if 'x_inactivation_escapees' in summary_stats:
    print(f"  X-inactivation escapees: {summary_stats['x_inactivation_escapees']:,} ({summary_stats['x_inactivation_escapees_percent']:.1f}%)")
