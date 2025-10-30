#!/usr/bin/env python3
"""
Fast GO Enrichment Analysis for X Chromosome Selection Study

Uses pre-computed GO annotations to perform efficient enrichment analysis.
Focuses on the 11 overlap genes with known GO terms.

Performs three enrichment analyses:
A. High-confidence genes (11 overlap genes) vs all X genes
B. High-confidence genes (11 overlap genes) vs X genes under selection
C. X genes under selection vs autosomal genes under selection
"""

import pandas as pd
from pathlib import Path
from scipy import stats
import numpy as np
from collections import defaultdict

# Define paths
BASE_DIR = Path('/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib')
GENE_LISTS_DIR = BASE_DIR / 'results/analysis/gene_annotation/gene_name_lists'
HIGH_SIG_DIR = BASE_DIR / 'results/analysis/high_significance_genes'
OUTPUT_DIR = BASE_DIR / 'results/analysis/functional_enrichment'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Load pre-computed GO annotations
GO_ANNOTATIONS_FILE = OUTPUT_DIR / 'gene_GO_annotations.tsv'

print("=" * 70)
print("Fast GO Enrichment Analysis for X Chromosome Selection")
print("=" * 70)
print()

# ============================================================================
# Load gene lists
# ============================================================================

print("Step 1: Loading gene lists...")
print()

def load_gene_list(filepath):
    """Load gene list from file, return as list."""
    with open(filepath) as f:
        genes = [line.strip() for line in f if line.strip()]
    return genes

# Test sets - use the 11 overlap genes we have GO terms for
overlap_genes = load_gene_list(HIGH_SIG_DIR / 'overlap_genes.txt')
x_selection_genes = load_gene_list(GENE_LISTS_DIR / 'X_protein_coding_genes.txt')

# Background sets
all_x_genes = load_gene_list(GENE_LISTS_DIR / 'X_chromosome_genes.txt')
autosome_selection_genes = load_gene_list(GENE_LISTS_DIR / 'autosome_protein_coding_in_selection.txt')

print(f"  Overlap genes (test set):              {len(overlap_genes):5d}")
print(f"  X genes under selection (test set C):  {len(x_selection_genes):5d}")
print(f"  All X genes (background A):            {len(all_x_genes):5d}")
print(f"  Autosomal genes (background C):        {len(autosome_selection_genes):5d}")
print()

# ============================================================================
# Load pre-computed GO annotations
# ============================================================================

print("Step 2: Loading pre-computed GO annotations...")
print()

go_annotations = pd.read_csv(GO_ANNOTATIONS_FILE, sep='\t')
print(f"  Loaded {len(go_annotations)} GO term annotations for {go_annotations['Gene'].nunique()} genes")
print()

# Create gene -> GO terms mapping
gene_to_go = defaultdict(set)
for _, row in go_annotations.iterrows():
    gene_to_go[row['Gene']].add(row['GO_ID'])

# Create GO term -> genes mapping
go_to_genes = defaultdict(set)
for _, row in go_annotations.iterrows():
    go_to_genes[row['GO_ID']].add(row['Gene'])

# GO term metadata
go_metadata = go_annotations[['GO_ID', 'GO_Name', 'GO_Category']].drop_duplicates().set_index('GO_ID')

print(f"  Total unique GO terms: {len(go_to_genes)}")
print(f"  Genes with GO annotations: {len(gene_to_go)}")
print()

# ============================================================================
# Helper functions for enrichment analysis
# ============================================================================

def perform_go_enrichment_fast(test_genes, background_genes, analysis_name):
    """
    Perform GO enrichment analysis using pre-computed annotations.

    Args:
        test_genes: List of genes in test set
        background_genes: List of genes in background set
        analysis_name: Name for this analysis

    Returns:
        DataFrame with enrichment results
    """
    print(f"  Running {analysis_name}...")
    print(f"    Test set: {len(test_genes)} genes")
    print(f"    Background: {len(background_genes)} genes")

    # Filter to genes we have GO annotations for
    test_genes_annotated = [g for g in test_genes if g in gene_to_go]
    background_genes_annotated = [g for g in background_genes if g in gene_to_go]

    print(f"    Test genes with GO annotations: {len(test_genes_annotated)}")
    print(f"    Background genes with GO annotations: {len(background_genes_annotated)}")

    if len(test_genes_annotated) == 0:
        print("    Warning: No test genes have GO annotations!")
        return pd.DataFrame()

    if len(background_genes_annotated) == 0:
        print("    Warning: No background genes have GO annotations!")
        return pd.DataFrame()

    # Get all GO terms present in either test or background
    all_go_terms = set()
    for gene in test_genes_annotated:
        all_go_terms.update(gene_to_go[gene])
    for gene in background_genes_annotated:
        all_go_terms.update(gene_to_go[gene])

    print(f"    Testing {len(all_go_terms)} GO terms")

    # Perform Fisher's exact test for each GO term
    results = []
    n_test = len(test_genes_annotated)
    n_background = len(background_genes_annotated)

    for go_id in all_go_terms:
        # Count genes with this GO term
        test_with_term = sum(1 for g in test_genes_annotated if go_id in gene_to_go[g])
        test_without_term = n_test - test_with_term
        bg_with_term = sum(1 for g in background_genes_annotated if go_id in gene_to_go[g])
        bg_without_term = n_background - bg_with_term

        # Fisher's exact test (one-sided, testing for enrichment)
        oddsratio, pvalue = stats.fisher_exact(
            [[test_with_term, test_without_term],
             [bg_with_term, bg_without_term]],
            alternative='greater'
        )

        # Calculate fold enrichment
        test_freq = test_with_term / n_test if n_test > 0 else 0
        bg_freq = bg_with_term / n_background if n_background > 0 else 0
        fold_enrichment = test_freq / bg_freq if bg_freq > 0 else float('inf')

        # Get gene names
        genes_with_term = [g for g in test_genes_annotated if go_id in gene_to_go[g]]

        results.append({
            'GO_ID': go_id,
            'Test_Count': test_with_term,
            'Test_Total': n_test,
            'Test_Frequency': test_freq,
            'Background_Count': bg_with_term,
            'Background_Total': n_background,
            'Background_Frequency': bg_freq,
            'Fold_Enrichment': fold_enrichment,
            'P_value': pvalue,
            'Odds_Ratio': oddsratio,
            'Genes_With_Term': ', '.join(sorted(genes_with_term))
        })

    # Create DataFrame and sort by p-value
    df = pd.DataFrame(results)
    df = df.sort_values('P_value')

    # Apply FDR correction (Benjamini-Hochberg)
    from scipy.stats import false_discovery_control
    if len(df) > 0:
        df['FDR'] = false_discovery_control(df['P_value'].values)

    # Add GO term metadata
    df['GO_Name'] = df['GO_ID'].map(lambda x: go_metadata.loc[x, 'GO_Name'] if x in go_metadata.index else 'Unknown')
    df['GO_Category'] = df['GO_ID'].map(lambda x: go_metadata.loc[x, 'GO_Category'] if x in go_metadata.index else 'Unknown')

    print(f"    Completed enrichment testing for {len(df)} GO terms")
    print()

    return df


# ============================================================================
# Analysis A: High-confidence genes (11 overlap) vs all X genes
# ============================================================================

print("=" * 70)
print("ANALYSIS A: High-confidence genes (11 overlap) vs all X genes")
print("=" * 70)
print()

results_a = perform_go_enrichment_fast(
    overlap_genes,
    all_x_genes,
    "Analysis A"
)

if len(results_a) > 0:
    # Save results
    output_file_a = OUTPUT_DIR / 'GO_enrichment_high_confidence_vs_all_X.tsv'
    results_a.to_csv(output_file_a, sep='\t', index=False)
    print(f"  Saved: {output_file_a}")
    print()

# ============================================================================
# Analysis B: High-confidence genes (11 overlap) vs X genes under selection
# ============================================================================

print("=" * 70)
print("ANALYSIS B: High-confidence genes (11 overlap) vs X genes under selection")
print("=" * 70)
print()

results_b = perform_go_enrichment_fast(
    overlap_genes,
    x_selection_genes,
    "Analysis B"
)

if len(results_b) > 0:
    # Save results
    output_file_b = OUTPUT_DIR / 'GO_enrichment_high_confidence_vs_X_selection.tsv'
    results_b.to_csv(output_file_b, sep='\t', index=False)
    print(f"  Saved: {output_file_b}")
    print()

# ============================================================================
# Analysis C: X genes vs autosomal genes (both under selection)
# ============================================================================

print("=" * 70)
print("ANALYSIS C: X genes vs autosomal genes (both under selection)")
print("=" * 70)
print()
print("  Note: This analysis requires GO annotations for autosomal genes,")
print("  which we don't have pre-computed. Skipping for now.")
print("  To run this analysis, you would need to:")
print("  1. Retrieve GO terms for 241 X genes under selection")
print("  2. Retrieve GO terms for 16,323 autosomal genes")
print("  This would take several hours with individual NCBI queries.")
print()

results_c = pd.DataFrame()  # Empty for now

# ============================================================================
# Generate summary tables
# ============================================================================

print("=" * 70)
print("SUMMARY: Significantly Enriched GO Terms (FDR < 0.05)")
print("=" * 70)
print()

# GO terms of interest for X chromosome evolution
FOCUS_GO_TERMS = {
    'GO:0000003': 'reproduction',
    'GO:0019953': 'sexual reproduction',
    'GO:0007276': 'gamete generation',
    'GO:0007283': 'spermatogenesis',
    'GO:0048477': 'oogenesis',
    'GO:0007292': 'female gamete generation',
    'GO:0051321': 'meiotic cell cycle',
    'GO:0007126': 'meiotic cell cycle process',
    'GO:0002376': 'immune system process',
    'GO:0006952': 'defense response',
    'GO:0006954': 'inflammatory response',
    'GO:0002250': 'adaptive immune response',
    'GO:0045087': 'innate immune response',
    'GO:0006955': 'immune response',
    'GO:0002682': 'regulation of immune system process',
}

def summarize_enrichment(df, analysis_name):
    """Print summary of significant enrichment results."""
    if len(df) == 0:
        print(f"{analysis_name}: No results (skipped or no data)")
        print()
        return pd.DataFrame()

    sig = df[df['FDR'] < 0.05]

    print(f"{analysis_name}:")
    print(f"  Total GO terms tested: {len(df)}")
    print(f"  Significant terms (FDR < 0.05): {len(sig)}")

    if len(sig) > 0:
        print()
        print("  Top 10 enriched terms:")
        print("  " + "-" * 66)
        for i, row in sig.head(10).iterrows():
            go_name_short = row['GO_Name'][:35] if len(row['GO_Name']) <= 35 else row['GO_Name'][:32] + '...'
            print(f"  {row['GO_ID']} | {go_name_short:35s} | FDR={row['FDR']:.2e} | Fold={row['Fold_Enrichment']:.2f}")
        print()

        # Check for focus GO terms
        focus_terms = sig[sig['GO_ID'].isin(FOCUS_GO_TERMS.keys())]
        if len(focus_terms) > 0:
            print("  Focus terms (reproduction/immunity):")
            print("  " + "-" * 66)
            for i, row in focus_terms.iterrows():
                go_name_short = row['GO_Name'][:35] if len(row['GO_Name']) <= 35 else row['GO_Name'][:32] + '...'
                print(f"  {row['GO_ID']} | {go_name_short:35s} | FDR={row['FDR']:.2e} | Fold={row['Fold_Enrichment']:.2f}")
                print(f"    Genes: {row['Genes_With_Term']}")
        else:
            print("  No focus terms (reproduction/immunity) significantly enriched.")
    else:
        print("  No significantly enriched terms found.")

    print()
    return sig

sig_a = summarize_enrichment(results_a, "Analysis A: High-confidence (11 overlap) vs all X genes")
sig_b = summarize_enrichment(results_b, "Analysis B: High-confidence (11 overlap) vs X selection genes")
sig_c = summarize_enrichment(results_c, "Analysis C: X genes vs autosomal genes")

# Create combined summary table
summary_data = []

for analysis_name, sig_df in [('A_high_conf_vs_all_X', sig_a),
                               ('B_high_conf_vs_X_selection', sig_b)]:
    for _, row in sig_df.iterrows():
        summary_data.append({
            'Analysis': analysis_name,
            'GO_ID': row['GO_ID'],
            'GO_Name': row['GO_Name'],
            'GO_Category': row['GO_Category'],
            'Fold_Enrichment': row['Fold_Enrichment'],
            'P_value': row['P_value'],
            'FDR': row['FDR'],
            'Test_Genes': f"{row['Test_Count']}/{row['Test_Total']}",
            'Background_Genes': f"{row['Background_Count']}/{row['Background_Total']}",
            'Genes_With_Term': row['Genes_With_Term']
        })

if summary_data:
    summary_df = pd.DataFrame(summary_data)
    summary_file = OUTPUT_DIR / 'GO_enrichment_summary_significant_terms.tsv'
    summary_df.to_csv(summary_file, sep='\t', index=False)

    print("=" * 70)
    print("Analysis complete!")
    print("=" * 70)
    print()
    print("Output files:")
    if len(results_a) > 0:
        print(f"  {output_file_a}")
    if len(results_b) > 0:
        print(f"  {output_file_b}")
    print(f"  {summary_file}")
    print()
else:
    print("=" * 70)
    print("Analysis complete! No significant enrichments found.")
    print("=" * 70)
    print()
