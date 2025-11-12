#!/usr/bin/env python3
"""
GO enrichment analysis with functional categories and DAG visualization

This script:
1. Searches for GO terms matching functional categories
2. Performs enrichment analysis for high-confidence genes
3. Visualizes GO term relationships using DAGs

Author: MIT van Bruggen
Date: 2025-11-11
"""

import pandas as pd
from pathlib import Path
from geneinfo.genelist import GeneList as glist
import geneinfo.ontology as go
import shutil
from scipy.stats import fisher_exact
import numpy as np

# Set your email for GO database access
go.email('au799024@uni.au.dk')

# Setup paths
PROJECT_DIR = Path('/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib')
RESULTS_DIR = PROJECT_DIR / 'results/analysis'
OUTPUT_DIR = RESULTS_DIR / 'functional_enrichment/GO_with_terms'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# My study genes (high-confidence X chromosome genes)
study_genes = ['DACH2', 'DIAPH2', 'FAAH2', 'GPC3', 'HTR2C',
               'IL1RAPL2', 'KDM6A', 'NAP1L2', 'NCBP2L',
               'TRPC5', 'ZMAT1']

# Create a GeneList object
study_gl = glist(study_genes)

print("="*70) #nice output separator
print("GO Enrichment Analysis with Functional Categories")
print("="*70)
print(f"\nStudy genes: {', '.join(study_genes)}")
print(f"Total: {len(study_genes)} genes")
print()

#==============================================================================
# Step 1: Define functional categories (same as get_GO_annotations.py)
#==============================================================================

print("Step 1: Defining functional categories...") #had some help defining the keywords 
print()

KEYWORD_CATEGORIES = {
    'reproduction': [
        'meiosis', 'gamete', 'fertility', 'gonad', 'sperm', 'oocyte',
        'spermatogenesis', 'oogenesis', 'germ cell', 'reproductive',
        'fertilization', 'ovary', 'testis'
    ],
    'immunity': [
        'immune', 'inflammation', 'defense', 'cytokine', 'antibody',
        'inflammatory', 'innate immune', 'adaptive immune', 'lymphocyte',
        'T cell', 'B cell', 'interferon', 'interleukin', 'antigen'
    ],
    'neurodevelopment': [
        'neuron', 'synapse', 'cognition', 'brain', 'neural', 'axon',
        'dendrite', 'synaptic', 'neuronal', 'nervous system', 'learning',
        'memory', 'behavior', 'neurogenesis', 'neurotransmitter'
    ],
    'development': [
        'development', 'differentiation', 'morphogenesis', 'organogenesis',
        'embryo', 'pattern', 'specification', 'determination'
    ]
}

#==============================================================================
# Step 2: Search for GO terms matching each category
#==============================================================================

print("Step 2: Searching for GO terms matching each functional category...")
print()

category_terms = {}

for category, keywords in KEYWORD_CATEGORIES.items(): #so look at each category and its keywords
    print(f"  Searching for {category} terms...")

    # Create regex pattern from keywords (combine with OR)
    pattern = '|'.join(keywords)

    # Search for matching GO terms
    terms = go.get_terms_for_go_regex(pattern)

    category_terms[category] = terms
    print(f"    Found {len(terms)} GO terms matching {category}")

print()
print("Summary of GO terms found:")
for category, terms in category_terms.items():
    print(f"  {category:20s}: {len(terms):4d} terms")

#==============================================================================
# Step 3: Save GO term lists for each category
#==============================================================================

print("\nStep 3: Saving GO term lists...")
print()

for category, terms in category_terms.items():
    output_file = OUTPUT_DIR / f'GO_terms_{category}.tsv'

    # Convert to DataFrame
    terms_df = pd.DataFrame([
        {'GO_ID': term_id, 'Category': category}
        for term_id in terms
    ])

    terms_df.to_csv(output_file, sep='\t', index=False)
    print(f"  Saved: {output_file}")

print()

#==============================================================================
# Step 4: Get genes associated with each GO term
#==============================================================================

print("\nStep 4: Getting genes associated with each GO term...")
print()

# Get genes for specific GO terms (returns DataFrames)
neurodevelopment_df = go.get_genes_for_go_terms(category_terms['neurodevelopment'])
reproduction_df = go.get_genes_for_go_terms(category_terms['reproduction'])
immunity_df = go.get_genes_for_go_terms(category_terms['immunity'])
development_df = go.get_genes_for_go_terms(category_terms['development'])

# Extract gene symbols as lists
neurodevelopment_genes = neurodevelopment_df['symbol'].unique().tolist()
reproduction_genes = reproduction_df['symbol'].unique().tolist()
immunity_genes = immunity_df['symbol'].unique().tolist()
development_genes = development_df['symbol'].unique().tolist()

# Convert to gene lists
neurodevelopment_gl = glist(neurodevelopment_genes)
reproduction_gl = glist(reproduction_genes)
immunity_gl = glist(immunity_genes)
development_gl = glist(development_genes)

print(f"  Neurodevelopment: {len(neurodevelopment_genes)} genes")
print(f"  Reproduction:     {len(reproduction_genes)} genes")
print(f"  Immunity:         {len(immunity_genes)} genes")
print(f"  Development:      {len(development_genes)} genes")

#==============================================================================
# Step 5: Define background gene sets
#==============================================================================

print("\n\nStep 5: Defining background gene sets...")
print()

# Background 1: All X chromosome genes
gene_table = go.gene_annotation_table()
all_x_genes = gene_table[gene_table['chromosome'] == 'X']['Symbol'].tolist()
x_background = glist(all_x_genes)
print(f"  Background 1: All X chromosome genes - {len(x_background)} genes")

# Background 2: All protein-coding genes
all_protein_coding = go.all_protein_coding()
pc_background = glist(all_protein_coding)
print(f"  Background 2: All protein-coding genes - {len(pc_background)} genes")

# Background 3: X chromosome genes under selection (from Phase 3)
x_selection_file = PROJECT_DIR / 'results/analysis/gene_annotation/X_protein_coding_genes.tsv'
if x_selection_file.exists():
    x_selection_df = pd.read_csv(x_selection_file, sep='\t')
    x_selection_genes = x_selection_df['gene_name'].unique().tolist()
    x_sel_background = glist(x_selection_genes)
    print(f"  Background 3: X chromosome genes under selection - {len(x_sel_background)} genes")
else:
    print(f"  Background 3: Warning - File not found: {x_selection_file}")

# Background 4: Only genes in specific GO categories
# (e.g., neuron-related genes as background)
neurodevelopment_background = glist(neurodevelopment_genes)
reproduction_background = glist(reproduction_genes)
immunity_background = glist(immunity_genes)
development_background = glist(development_genes)

print(f"  Background 4: Neurodevelopment genes - {len(neurodevelopment_background)} genes")
print(f"  Background 4: Reproduction genes - {len(reproduction_background)} genes")
print(f"  Background 4: Immunity genes - {len(immunity_background)} genes")
print(f"  Background 4: Development genes - {len(development_background)} genes")

print()

#==============================================================================
# Step 6: Visualize GO term hierarchies using DAGs
#==============================================================================

print("\nStep 6: Visualizing GO term hierarchies as DAGs...")
print()

# Visualize the hierarchy for each functional category
for category, terms in category_terms.items():
    print(f"  Creating DAG for {category} ({len(terms)} terms)...")

    # Create the DAG visualization
    # Note: show_go_dag_for_terms() saves to cache_dir/plot.png
    go.show_go_dag_for_terms(terms)

    # Copy the cached plot to our output directory
    cache_plot = Path(go.cache_dir) / 'plot.png'
    output_fig = OUTPUT_DIR / f'GO_DAG_{category}.png'

    if cache_plot.exists():
        shutil.copy2(cache_plot, output_fig)
        print(f"    ✓ DAG saved to: {output_fig}")
    else:
        print(f"    ✗ Warning: Cache plot not found at {cache_plot}")
    print()

#==============================================================================
# Step 7: Perform Fisher's Exact Tests for enrichment
#==============================================================================

print("\nStep 7: Performing Fisher's exact tests for enrichment...")
print()

def custom_enrichment_test(study_genes, category_genes, background_genes):
    """
    Perform Fisher's exact test for enrichment

    Parameters:
    - study_genes: Your genes of interest (e.g., high-confidence X genes)
    - category_genes: Genes in a specific GO category (e.g., neuron genes)
    - background_genes: All possible genes for comparison

    Returns: odds ratio, p-value, contingency table
    """
    # Convert to sets for faster operations
    study_set = set(study_genes)
    category_set = set(category_genes)
    background_set = set(background_genes)

    # Build 2x2 contingency table
    a = len(study_set & category_set)  # In study AND in category
    b = len(study_set - category_set)  # In study, NOT in category
    c = len(category_set - study_set)  # NOT in study, IN category
    d = len(background_set - study_set - category_set)  # Neither

    contingency_table = [[a, b], [c, d]]

    # Fisher's exact test (one-sided: testing for enrichment)
    odds_ratio, p_value = fisher_exact(contingency_table, alternative='greater')

    return odds_ratio, p_value, contingency_table, a

# Store results for each test
enrichment_results = []

# Define the backgrounds to test against
backgrounds = {
    'X chromosome genes': (x_background, all_x_genes),
    'All protein-coding genes': (pc_background, all_protein_coding),
    'X genes under selection': (x_sel_background, x_selection_genes) if x_selection_file.exists() else None
}

# Define the functional categories to test
category_gene_sets = {
    'neurodevelopment': neurodevelopment_genes,
    'reproduction': reproduction_genes,
    'immunity': immunity_genes,
    'development': development_genes
}

# Perform enrichment tests for each category against each background
for bg_name, bg_data in backgrounds.items():
    if bg_data is None:
        continue

    bg_genelist, bg_genes = bg_data

    print(f"\nTesting against background: {bg_name}")
    print(f"  Background size: {len(bg_genes)} genes")
    print()

    for category, category_genes in category_gene_sets.items():
        or_val, p_val, table, overlap = custom_enrichment_test(
            study_genes=study_genes,
            category_genes=category_genes,
            background_genes=bg_genes
        )

        # Store results
        enrichment_results.append({
            'Background': bg_name,
            'Category': category,
            'Study_genes': len(study_genes),
            'Category_genes': len(category_genes),
            'Overlap': overlap,
            'Odds_Ratio': or_val,
            'P_value': p_val,
            'In_study_in_category': table[0][0],
            'In_study_not_category': table[0][1],
            'Not_study_in_category': table[1][0],
            'Not_study_not_category': table[1][1]
        })

        # Print results
        significance = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else ""
        print(f"  {category:20s}: OR={or_val:6.2f}, p={p_val:.2e} {significance} (overlap: {overlap}/{len(study_genes)})")

print()

# Save enrichment results to file
enrichment_df = pd.DataFrame(enrichment_results)
enrichment_output = OUTPUT_DIR / 'enrichment_results.tsv'
enrichment_df.to_csv(enrichment_output, sep='\t', index=False)
print(f"Enrichment results saved to: {enrichment_output}")

# Create a summary of significant results (p < 0.05)
significant_results = enrichment_df[enrichment_df['P_value'] < 0.05].copy()
if len(significant_results) > 0:
    significant_results = significant_results.sort_values('P_value')
    sig_output = OUTPUT_DIR / 'significant_enrichments.tsv'
    significant_results.to_csv(sig_output, sep='\t', index=False)
    print(f"Significant enrichments saved to: {sig_output}")
    print(f"\nFound {len(significant_results)} significant enrichments (p < 0.05)")
else:
    print("\nNo significant enrichments found (p < 0.05)")

print()
print("="*70)
print("GO enrichment analysis with functional categories completed.")
print(f"\nAll results saved to: {OUTPUT_DIR}")