#!/usr/bin/env python3
"""
GO Enrichment Analysis for X Chromosome Selection Study

Performs three enrichment analyses:
A. High-confidence genes vs all X genes
B. High-confidence genes vs X genes under selection
C. X genes under selection vs autosomal genes under selection

Uses geneinfo package for GO enrichment testing with FDR correction.
"""

import geneinfo.ontology as go
from pathlib import Path
import pandas as pd
from scipy import stats
import numpy as np

# Set email for NCBI queries
go.email('vanbruggenmit@birc.au.dk')

# Define paths
BASE_DIR = Path('/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib')
GENE_LISTS_DIR = BASE_DIR / 'results/analysis/gene_annotation/gene_name_lists'
HIGH_SIG_DIR = BASE_DIR / 'results/analysis/high_significance_genes'
OUTPUT_DIR = BASE_DIR / 'results/analysis/functional_enrichment'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

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

print("=" * 70)
print("GO Enrichment Analysis for X Chromosome Selection")
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

# Test sets
high_confidence_genes = load_gene_list(HIGH_SIG_DIR / 'high_confidence_X_genes.txt')
x_selection_genes = load_gene_list(GENE_LISTS_DIR / 'X_protein_coding_genes.txt')

# Background sets
all_x_genes = load_gene_list(GENE_LISTS_DIR / 'X_chromosome_genes.txt')
autosome_selection_genes = load_gene_list(GENE_LISTS_DIR / 'autosome_protein_coding_in_selection.txt')

print(f"  High-confidence genes (test set):     {len(high_confidence_genes):5d}")
print(f"  X genes under selection (test set C): {len(x_selection_genes):5d}")
print(f"  All X genes (background A):           {len(all_x_genes):5d}")
print(f"  Autosomal genes (background C):       {len(autosome_selection_genes):5d}")
print()

# ============================================================================
# Helper functions for enrichment analysis
# ============================================================================

def perform_go_enrichment(test_genes, background_genes, analysis_name):
    """
    Perform GO enrichment analysis using Fisher's exact test.

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

    # Get GO terms for test genes
    test_go_terms = {}
    for gene in test_genes:
        try:
            terms = go.get_go_terms_for_genes([gene])
            if terms:
                test_go_terms[gene] = terms
        except Exception as e:
            print(f"    Warning: Could not retrieve GO terms for {gene}: {e}")

    # Get GO terms for background genes
    print(f"    Retrieving GO terms for background genes (this may take a while)...")
    background_go_terms = {}
    for i, gene in enumerate(background_genes):
        if (i + 1) % 100 == 0:
            print(f"      Progress: {i+1}/{len(background_genes)} genes...")
        try:
            terms = go.get_go_terms_for_genes([gene])
            if terms:
                background_go_terms[gene] = terms
        except Exception as e:
            pass  # Silent fail for background genes

    print(f"    Retrieved GO terms for {len(test_go_terms)}/{len(test_genes)} test genes")
    print(f"    Retrieved GO terms for {len(background_go_terms)}/{len(background_genes)} background genes")

    # Count GO term occurrences
    go_term_counts = {}

    # Count in test set
    for gene, terms in test_go_terms.items():
        for term_id in terms:
            if term_id not in go_term_counts:
                go_term_counts[term_id] = {'test': 0, 'background': 0}
            go_term_counts[term_id]['test'] += 1

    # Count in background set
    for gene, terms in background_go_terms.items():
        for term_id in terms:
            if term_id not in go_term_counts:
                go_term_counts[term_id] = {'test': 0, 'background': 0}
            go_term_counts[term_id]['background'] += 1

    print(f"    Found {len(go_term_counts)} unique GO terms")

    # Perform Fisher's exact test for each GO term
    results = []
    n_test = len(test_go_terms)
    n_background = len(background_go_terms)

    for go_id, counts in go_term_counts.items():
        # 2x2 contingency table:
        # [[genes_with_term_in_test, genes_without_term_in_test],
        #  [genes_with_term_in_background, genes_without_term_in_background]]

        a = counts['test']  # Test genes with this GO term
        b = n_test - a      # Test genes without this GO term
        c = counts['background']  # Background genes with this GO term
        d = n_background - c      # Background genes without this GO term

        # Fisher's exact test (one-sided, testing for enrichment)
        oddsratio, pvalue = stats.fisher_exact([[a, b], [c, d]], alternative='greater')

        # Calculate fold enrichment
        test_freq = a / n_test if n_test > 0 else 0
        bg_freq = c / n_background if n_background > 0 else 0
        fold_enrichment = test_freq / bg_freq if bg_freq > 0 else float('inf')

        results.append({
            'GO_ID': go_id,
            'Test_Count': a,
            'Test_Total': n_test,
            'Test_Frequency': test_freq,
            'Background_Count': c,
            'Background_Total': n_background,
            'Background_Frequency': bg_freq,
            'Fold_Enrichment': fold_enrichment,
            'P_value': pvalue,
            'Odds_Ratio': oddsratio
        })

    # Create DataFrame and sort by p-value
    df = pd.DataFrame(results)
    df = df.sort_values('P_value')

    # Apply FDR correction (Benjamini-Hochberg)
    from scipy.stats import false_discovery_control
    df['FDR'] = false_discovery_control(df['P_value'].values)

    print(f"    Completed enrichment testing for {len(df)} GO terms")
    print()

    return df, test_go_terms, background_go_terms


def annotate_go_terms(df):
    """Add GO term names and categories to results DataFrame."""
    from goatools.obo_parser import GODag
    import geneinfo
    import os

    # Load GO DAG
    geneinfo_dir = os.path.dirname(geneinfo.__file__)
    obo_file = Path(geneinfo_dir) / 'cache' / 'go-basic.obo'
    godag = GODag(str(obo_file), optional_attrs='relationship')

    # Add GO term information
    go_names = []
    go_categories = []

    for go_id in df['GO_ID']:
        if go_id in godag:
            go_term = godag[go_id]
            go_names.append(go_term.name)
            go_categories.append(go_term.namespace)
        else:
            go_names.append('Unknown')
            go_categories.append('Unknown')

    df['GO_Name'] = go_names
    df['GO_Category'] = go_categories

    return df


def identify_contributing_genes(go_id, test_go_terms):
    """Identify which genes have a specific GO term."""
    genes_with_term = [gene for gene, terms in test_go_terms.items() if go_id in terms]
    return ', '.join(sorted(genes_with_term))


# ============================================================================
# Analysis A: High-confidence genes vs all X genes
# ============================================================================

print("=" * 70)
print("ANALYSIS A: High-confidence genes vs all X genes")
print("=" * 70)
print()

results_a, test_go_a, bg_go_a = perform_go_enrichment(
    high_confidence_genes,
    all_x_genes,
    "Analysis A"
)

# Annotate with GO term names
print("  Annotating GO terms with names and categories...")
results_a = annotate_go_terms(results_a)

# Add contributing genes for significant terms
results_a['Genes_With_Term'] = results_a['GO_ID'].apply(
    lambda x: identify_contributing_genes(x, test_go_a)
)

# Save results
output_file_a = OUTPUT_DIR / 'GO_enrichment_high_confidence_vs_all_X.tsv'
results_a.to_csv(output_file_a, sep='\t', index=False)
print(f"  Saved: {output_file_a}")
print()

# ============================================================================
# Analysis B: High-confidence genes vs X genes under selection
# ============================================================================

print("=" * 70)
print("ANALYSIS B: High-confidence genes vs X genes under selection")
print("=" * 70)
print()

results_b, test_go_b, bg_go_b = perform_go_enrichment(
    high_confidence_genes,
    x_selection_genes,
    "Analysis B"
)

# Annotate with GO term names
print("  Annotating GO terms with names and categories...")
results_b = annotate_go_terms(results_b)

# Add contributing genes
results_b['Genes_With_Term'] = results_b['GO_ID'].apply(
    lambda x: identify_contributing_genes(x, test_go_b)
)

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

results_c, test_go_c, bg_go_c = perform_go_enrichment(
    x_selection_genes,
    autosome_selection_genes,
    "Analysis C"
)

# Annotate with GO term names
print("  Annotating GO terms with names and categories...")
results_c = annotate_go_terms(results_c)

# Add contributing genes
results_c['Genes_With_Term'] = results_c['GO_ID'].apply(
    lambda x: identify_contributing_genes(x, test_go_c)
)

# Save results
output_file_c = OUTPUT_DIR / 'GO_enrichment_X_vs_autosomes.tsv'
results_c.to_csv(output_file_c, sep='\t', index=False)
print(f"  Saved: {output_file_c}")
print()

# ============================================================================
# Generate summary tables
# ============================================================================

print("=" * 70)
print("SUMMARY: Significantly Enriched GO Terms (FDR < 0.05)")
print("=" * 70)
print()

def summarize_enrichment(df, analysis_name):
    """Print summary of significant enrichment results."""
    sig = df[df['FDR'] < 0.05]

    print(f"{analysis_name}:")
    print(f"  Total GO terms tested: {len(df)}")
    print(f"  Significant terms (FDR < 0.05): {len(sig)}")

    if len(sig) > 0:
        print()
        print("  Top 10 enriched terms:")
        print("  " + "-" * 66)
        for i, row in sig.head(10).iterrows():
            print(f"  {row['GO_ID']} | {row['GO_Name'][:40]:40s} | FDR={row['FDR']:.2e} | Fold={row['Fold_Enrichment']:.2f}")
        print()

        # Check for focus GO terms
        focus_terms = sig[sig['GO_ID'].isin(FOCUS_GO_TERMS.keys())]
        if len(focus_terms) > 0:
            print("  Focus terms (reproduction/immunity):")
            print("  " + "-" * 66)
            for i, row in focus_terms.iterrows():
                print(f"  {row['GO_ID']} | {row['GO_Name'][:40]:40s} | FDR={row['FDR']:.2e} | Fold={row['Fold_Enrichment']:.2f}")
                print(f"    Genes: {row['Genes_With_Term']}")
        else:
            print("  No focus terms (reproduction/immunity) significantly enriched.")
    else:
        print("  No significantly enriched terms found.")

    print()
    return sig

sig_a = summarize_enrichment(results_a, "Analysis A: High-confidence vs all X genes")
sig_b = summarize_enrichment(results_b, "Analysis B: High-confidence vs X selection genes")
sig_c = summarize_enrichment(results_c, "Analysis C: X genes vs autosomal genes")

# Create combined summary table
summary_data = []

for analysis_name, sig_df in [('A_high_conf_vs_all_X', sig_a),
                               ('B_high_conf_vs_X_selection', sig_b),
                               ('C_X_vs_autosomes', sig_c)]:
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

summary_df = pd.DataFrame(summary_data)
summary_file = OUTPUT_DIR / 'GO_enrichment_summary_significant_terms.tsv'
summary_df.to_csv(summary_file, sep='\t', index=False)

print("=" * 70)
print("Analysis complete!")
print("=" * 70)
print()
print("Output files:")
print(f"  {output_file_a}")
print(f"  {output_file_b}")
print(f"  {output_file_c}")
print(f"  {summary_file}")
print()
