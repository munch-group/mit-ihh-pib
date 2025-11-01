#!/usr/bin/env python3
"""
Retrieve GO terms for all X chromosome genes for enrichment analysis.

This script retrieves GO annotations for:
1. All X chromosome protein-coding genes (664 genes)
2. All Phase 3 X selection genes (241 genes)

Output is saved as a comprehensive GO annotation table that can be used
for fast enrichment analysis.

Features:
- Efficient bulk retrieval using gene_annotation_table()
- No need for batching or checkpointing (completes in seconds)
- Progress tracking and detailed logging
"""

import geneinfo.ontology as go
from pathlib import Path
import pandas as pd
from goatools.obo_parser import GODag
import geneinfo
import os
import sys
from datetime import datetime

# Set email for NCBI queries
go.email('au799024@uni.au.dk')

# Define paths
BASE_DIR = Path('/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib')
GENE_LISTS_DIR = BASE_DIR / 'results/analysis/gene_annotation/gene_name_lists'
OUTPUT_DIR = BASE_DIR / 'results/analysis/functional_enrichment'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 70)
print("GO Term Retrieval for X Chromosome Genes")
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 70)
print()

# ============================================================================
# Load gene lists
# ============================================================================

print("Step 1: Loading gene lists...")

def load_gene_list(filepath):
    """Load gene list from file, return as list."""
    with open(filepath) as f:
        genes = [line.strip() for line in f if line.strip()]
    return genes

# Load all X chromosome genes
all_x_genes = load_gene_list(GENE_LISTS_DIR / 'X_chromosome_genes.txt')
x_selection_genes = load_gene_list(GENE_LISTS_DIR / 'X_protein_coding_genes.txt')

# Combine into one unique list
all_genes = sorted(set(all_x_genes + x_selection_genes))

print(f"  Total X genes to annotate: {len(all_genes)}")
print(f"    All X chromosome genes: {len(all_x_genes)}")
print(f"    X genes under selection: {len(x_selection_genes)}")
print()

# ============================================================================
# Load GO database for term name lookup
# ============================================================================

print("Step 2: Loading GO database...")
geneinfo_dir = os.path.dirname(geneinfo.__file__)
obo_file = Path(geneinfo_dir) / 'cache' / 'go-basic.obo'

if not obo_file.exists():
    print(f"  ERROR: GO database not found at {obo_file}")
    sys.exit(1)

print(f"  Loading GO DAG from {obo_file}")
godag = GODag(str(obo_file), optional_attrs='relationship')
print(f"  Loaded {len(godag)} GO terms")
print()

# ============================================================================
# Retrieve GO annotations for all human genes at once
# ============================================================================

print("Step 3: Retrieving all human gene-GO annotations...")
print("  This will download the complete human GO annotation table")
print("  (this may take 1-2 minutes for the first download)")
print()

try:
    # Get all human gene-GO annotations in one call
    # This returns: taxid, GeneID, GO_ID, Evidence, Qualifier, GO_term, PubMed, Category
    all_annotations = go.go_annotation_table(taxid=9606)
    print(f"  Successfully retrieved {len(all_annotations)} gene-GO annotations")
    print(f"  Covering {all_annotations['GeneID'].nunique()} unique gene IDs")
    print()
except Exception as e:
    print(f"  ERROR: Failed to retrieve gene annotations: {e}")
    sys.exit(1)

# Get gene ID to symbol mapping
print("Step 4: Getting gene ID to symbol mapping...")
gene_info = go.gene_annotation_table(taxid=9606)
gene_id_to_symbol = dict(zip(gene_info['GeneID'], gene_info['Symbol']))
print(f"  Loaded mapping for {len(gene_id_to_symbol)} genes")
print()

# ============================================================================
# Filter annotations for X chromosome genes
# ============================================================================

print("Step 5: Filtering annotations for X chromosome genes...")

# Add gene symbols to annotations
all_annotations['Symbol'] = all_annotations['GeneID'].map(gene_id_to_symbol)

# Filter for our genes of interest
x_gene_annotations = all_annotations[all_annotations['Symbol'].isin(all_genes)].copy()

print(f"  Found annotations for {x_gene_annotations['Symbol'].nunique()} / {len(all_genes)} genes")
print(f"  Total GO annotations: {len(x_gene_annotations)}")
print()

# Check which genes have no annotations
genes_with_annotations = set(x_gene_annotations['Symbol'].unique())
genes_without_annotations = set(all_genes) - genes_with_annotations

if genes_without_annotations:
    print(f"  Genes without GO annotations: {len(genes_without_annotations)}")
    print(f"    Examples: {', '.join(list(genes_without_annotations)[:10])}")
    if len(genes_without_annotations) > 10:
        print(f"    ... and {len(genes_without_annotations) - 10} more")
    print()

# ============================================================================
# Process annotations and add metadata
# ============================================================================

print("Step 6: Processing annotations and adding metadata...")

# Add flags for gene lists
x_gene_annotations['In_X_Selection'] = x_gene_annotations['Symbol'].isin(x_selection_genes)
x_gene_annotations['In_All_X'] = x_gene_annotations['Symbol'].isin(all_x_genes)

# Rename and select columns
# Available columns: taxid, GeneID, GO_ID, Evidence, Qualifier, GO_term, PubMed, Category, Symbol
x_gene_annotations = x_gene_annotations.rename(columns={
    'Symbol': 'Gene',
    'GO_term': 'GO_Name',
    'Category': 'GO_Category'
})

# Select output columns
output_columns = ['Gene', 'GO_ID', 'GO_Name', 'GO_Category', 'Evidence', 'In_X_Selection', 'In_All_X']
gene_go_data = x_gene_annotations[output_columns].copy()

genes_with_terms = len(genes_with_annotations)
genes_without_terms = len(genes_without_annotations)

print(f"  Processing complete!")
print(f"    Genes with GO terms: {genes_with_terms}/{len(all_genes)} ({genes_with_terms/len(all_genes)*100:.1f}%)")
print(f"    Genes without GO terms: {genes_without_terms}")
print(f"    Total GO annotations: {len(gene_go_data)}")
print()

# ============================================================================
# Save results
# ============================================================================

print("Step 7: Saving results...")

# Save main GO annotations
if len(gene_go_data) > 0:
    output_file = OUTPUT_DIR / 'X_chromosome_GO_annotations_complete.tsv'
    gene_go_data.to_csv(output_file, sep='\t', index=False)
    print(f"  Saved GO annotations: {output_file}")

    # Summary statistics
    print()
    print("  Summary by GO category:")
    if 'GO_Category' in gene_go_data.columns:
        category_counts = gene_go_data.groupby('GO_Category').size().sort_values(ascending=False)
        for cat, count in category_counts.items():
            print(f"    {cat:25s}: {count:5d} annotations")
    else:
        print("    Warning: GO_Category column not found")

    print()
    print("  Genes with most GO terms:")
    gene_term_counts = gene_go_data.groupby('Gene').size().sort_values(ascending=False).head(10)
    for gene, count in gene_term_counts.items():
        print(f"    {gene:15s}: {count:3d} GO terms")

# Save genes without annotations
if genes_without_annotations:
    missing_file = OUTPUT_DIR / 'X_chromosome_genes_without_GO.txt'
    with open(missing_file, 'w') as f:
        for gene in sorted(genes_without_annotations):
            f.write(f"{gene}\n")
    print()
    print(f"  Saved genes without GO annotations: {missing_file}")

# Save summary statistics
unique_go_terms = len(gene_go_data['GO_ID'].unique()) if len(gene_go_data) > 0 else 0

summary = {
    'Total_Genes': len(all_genes),
    'Genes_With_GO_Terms': genes_with_terms,
    'Genes_Without_GO_Terms': genes_without_terms,
    'Total_GO_Annotations': len(gene_go_data),
    'Unique_GO_Terms': unique_go_terms,
    'Completion_Time': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
}

summary_file = OUTPUT_DIR / 'X_chromosome_GO_retrieval_summary.txt'
with open(summary_file, 'w') as f:
    f.write("GO Term Retrieval Summary\n")
    f.write("=" * 50 + "\n\n")
    for key, value in summary.items():
        f.write(f"{key:30s}: {value}\n")

print()
print(f"  Saved summary: {summary_file}")

print()
print("=" * 70)
print("GO term retrieval complete!")
print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 70)
