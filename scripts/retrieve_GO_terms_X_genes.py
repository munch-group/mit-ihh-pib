#!/usr/bin/env python3
"""
Retrieve GO terms for all X chromosome genes for enrichment analysis.

This script retrieves GO annotations for:
1. All X chromosome protein-coding genes (664 genes)
2. All Phase 3 X selection genes (241 genes)

Output is saved as a comprehensive GO annotation table that can be used
for fast enrichment analysis.
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
go.email('vanbruggenmit@birc.au.dk')

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
# Retrieve GO terms for all genes
# ============================================================================

print("Step 3: Retrieving GO terms for all genes...")
print(f"  This will take approximately {len(all_genes) * 3 / 3600:.1f} hours")
print(f"  (estimating ~3 seconds per gene)")
print()

gene_go_data = []
genes_with_terms = 0
genes_without_terms = 0
errors = []

for i, gene in enumerate(all_genes):
    if (i + 1) % 10 == 0:
        elapsed_pct = ((i + 1) / len(all_genes)) * 100
        print(f"  Progress: {i+1}/{len(all_genes)} genes ({elapsed_pct:.1f}%) - "
              f"Success: {genes_with_terms}, No terms: {genes_without_terms}, "
              f"Errors: {len(errors)}")

    try:
        # Get GO term IDs for this gene
        go_term_ids = go.get_go_terms_for_genes([gene])

        if go_term_ids:
            genes_with_terms += 1

            # Get details for each GO term
            for go_id in go_term_ids:
                if go_id in godag:
                    go_term = godag[go_id]
                    gene_go_data.append({
                        'Gene': gene,
                        'GO_ID': go_id,
                        'GO_Name': go_term.name,
                        'GO_Category': go_term.namespace,
                        'In_X_Selection': gene in x_selection_genes,
                        'In_All_X': gene in all_x_genes
                    })
                else:
                    # GO term not in database
                    gene_go_data.append({
                        'Gene': gene,
                        'GO_ID': go_id,
                        'GO_Name': 'Unknown',
                        'GO_Category': 'Unknown',
                        'In_X_Selection': gene in x_selection_genes,
                        'In_All_X': gene in all_x_genes
                    })
        else:
            genes_without_terms += 1

    except Exception as e:
        errors.append({'Gene': gene, 'Error': str(e)})
        print(f"  ERROR retrieving GO terms for {gene}: {e}")

print()
print(f"  Completed GO term retrieval!")
print(f"    Genes with GO terms: {genes_with_terms}/{len(all_genes)} ({genes_with_terms/len(all_genes)*100:.1f}%)")
print(f"    Genes without GO terms: {genes_without_terms}")
print(f"    Errors: {len(errors)}")
print(f"    Total GO annotations: {len(gene_go_data)}")
print()

# ============================================================================
# Save results
# ============================================================================

print("Step 4: Saving results...")

# Save main GO annotations
if gene_go_data:
    df = pd.DataFrame(gene_go_data)
    output_file = OUTPUT_DIR / 'X_chromosome_GO_annotations_complete.tsv'
    df.to_csv(output_file, sep='\t', index=False)
    print(f"  Saved GO annotations: {output_file}")

    # Summary statistics
    print()
    print("  Summary by GO category:")
    category_counts = df.groupby('GO_Category').size().sort_values(ascending=False)
    for cat, count in category_counts.items():
        print(f"    {cat:25s}: {count:5d} annotations")

    print()
    print("  Genes with most GO terms:")
    gene_term_counts = df.groupby('Gene').size().sort_values(ascending=False).head(10)
    for gene, count in gene_term_counts.items():
        print(f"    {gene:15s}: {count:3d} GO terms")

# Save error log
if errors:
    error_df = pd.DataFrame(errors)
    error_file = OUTPUT_DIR / 'X_chromosome_GO_retrieval_errors.tsv'
    error_df.to_csv(error_file, sep='\t', index=False)
    print()
    print(f"  Saved error log: {error_file}")

# Save summary statistics
summary = {
    'Total_Genes': len(all_genes),
    'Genes_With_GO_Terms': genes_with_terms,
    'Genes_Without_GO_Terms': genes_without_terms,
    'Genes_With_Errors': len(errors),
    'Total_GO_Annotations': len(gene_go_data),
    'Unique_GO_Terms': len(df['GO_ID'].unique()) if gene_go_data else 0,
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
