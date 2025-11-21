#!/usr/bin/env python3

"""
Overlap analysis with .fischer method on gene lists using focused gene sets from GO hierarchies
Scripts 
1. Loads focused gene sets by explore_GO_hierarchies.py
2. Performs significant overlap analysis for high confidence genes 
3. Uses X chromosome protein coding genes as background

Author: Mit Van Bruggen
Date: 2025-11-18
Updated: 2025-11-20

"""

import pandas as pd
from pathlib import Path
from geneinfo.genelist import GeneList as glist
import geneinfo.ontology as go
import shutil
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import numpy as np

go.email('au799024@uni.au.dk')

# Setup paths
PROJECT_DIR = Path('/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib')
RESULTS_DIR = PROJECT_DIR / 'results/analysis'
HIERARCHY_DIR= RESULTS_DIR / 'functional_enrichment/GO_hierarchies'
OUTPUT_DIR = RESULTS_DIR / 'functional_enrichment/focused_enrichment'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# My study genes (high-confidence X chromosome genes)
study_genes = ['DACH2', 'DIAPH2', 'FAAH2', 'GPC3', 'HTR2C',
               'IL1RAPL2', 'KDM6A', 'NAP1L2', 'NCBP2L',
               'TRPC5', 'ZMAT1', 'KLHL13', 'PCDH11X'] #, #'OPN1LW', 'OPN1MW', 'OPN1MW2', 'OPN1MW3'#

# Create a GeneList object
study_gl = glist(study_genes)

print("="*70) #nice output separator
print("GO Enrichment Analysis with Functional Categories")
print("="*70)
print(f"\nStudy genes: {', '.join(study_genes)}")
print(f"Total: {len(study_genes)} genes")
print()

#========================================================
# Step 1: Load focused gene sets from hierarchy analysis
#========================================================
#Load summary to see categories 
summary_file = HIERARCHY_DIR / 'subtree_summary.tsv'
summary_df=pd.read_csv(summary_file, sep='\t')
print("available focused gene sets:")
print(summary_df.to_string(index=False))
print()
#Load gene lists for each category 
category_genes={}
for _, row in summary_df.iterrows():
    category = row['category']
    genes_file = HIERARCHY_DIR / f'genes_unique_{category}.txt'
    if genes_file.exists():
        with open(genes_file, 'r') as f: 
            genes = [line.strip() for line in f if line.strip()]
        category_genes[category] = genes
        print(f" Loaded {category:30s}: {len(genes):4d} genes")
    else: 
        print(f" Warning: Could not find {genes_file}")
    print()

#========================================================
# Step 2: Define background gene set (X chromosome protein coding genes)
#========================================================
#Same as in GO_enrichment_with_terms.py
print("Step 2: Defining background gene set...")
print()

# Get all protein-coding genes
all_protein_coding = go.all_protein_coding()

# Get X chromosome gene symbols
gene_table = go.gene_annotation_table()
x_gene_symbols = set(gene_table[gene_table['chromosome'] == 'X']['Symbol'].tolist())

# Filter to X chromosome only
all_x_protein_coding = [gene for gene in all_protein_coding if gene in x_gene_symbols]

background_genes = all_x_protein_coding
print(f"  Background: All protein-coding X chromosome genes")
print(f"  Total background genes: {len(background_genes)}")
print()

#==============================================================================
# Step 3: Filter category genes to X chromosome only (step is unnecessary here because filtering already took place in explore_GO_hierarchies.py)
#==============================================================================

print("Step 3: Filtering category genes to X chromosome only...")
print()

background_set = set(background_genes)

category_genes_x = {}
for category, genes in category_genes.items():
    genes_x = [g for g in genes if g in background_set]
    category_genes_x[category] = genes_x
    print(f"  {category:30s}: {len(genes_x):4d} X genes (was {len(genes):4d} total)")

print()

#==============================================================================
# Step 4: Perform Fisher's exact tests (same as in GO_enrichment_with_terms.py)
#==============================================================================

print("Step 4: Performing Fisher's exact tests for enrichment...")
print("="*80)
print()

# Convert to GeneList objects
study_gl = glist(study_genes)
background_gl = glist(background_genes)

# Store results
enrichment_results = []

# Test each category
for category, category_genes_list in category_genes_x.items():
    if len(category_genes_list) == 0:
        print(f"{category:30s}: No X genes (skipped)")
        continue
    
    # Create GeneList
    category_gl = glist(category_genes_list)
    
    # Perform Fisher's exact test
    p_value = study_gl.fisher(category_gl, background=background_gl)
    
    # Calculate overlap
    study_set = set(study_genes)
    category_set = set(category_genes_list)
    overlap_genes = list(study_set & category_set)
    
    # Calculate odds ratio
    a = len(study_set & category_set)  # In study AND category
    b = len(study_set - category_set)  # In study NOT category
    c = len(category_set - study_set)  # In category NOT study
    d = len(background_set - study_set - category_set)  # In neither
    
    if b > 0 and c > 0:
        odds_ratio = (a * d) / (b * c)
    else:
        odds_ratio = float('inf') if a > 0 else 0.0
    
    # Store results
    enrichment_results.append({
        'Category': category,
        'Study_genes_total': len(study_genes),
        'Category_genes_X': len(category_genes_list),
        'Background_genes_total': len(background_genes),
        'Overlap': len(overlap_genes),
        'Overlapping_genes': ', '.join(sorted(overlap_genes)) if overlap_genes else 'None',
        'Odds_Ratio': odds_ratio,
        'P_value': p_value,
        'In_study_in_category': a,
        'In_study_not_category': b,
        'Not_study_in_category': c,
        'Not_study_not_category': d
    })
    
    # Print results
    sig = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else ""
    print(f"{category:30s}: OR={odds_ratio:6.2f}, p={p_value:.4e} {sig:3s}")
    print(f"  Overlap: {len(overlap_genes)}/{len(study_genes)} study genes")
    print(f"  Category size: {len(category_genes_list)} X genes")
    if overlap_genes:
        print(f"  Genes: {', '.join(sorted(overlap_genes))}")
    print()

print("="*80)

#==============================================================================
# Step 5: Apply FDR correction and save results
#==============================================================================

print("\nStep 5: Applying FDR correction...")
print()

# Create results dataframe
enrichment_df = pd.DataFrame(enrichment_results)

# Apply FDR correction
reject, pvals_corrected, _, _ = multipletests(
    enrichment_df['P_value'],
    alpha=0.05,
    method='fdr_bh'
)

enrichment_df['FDR_corrected'] = pvals_corrected
enrichment_df['Significant_FDR_0.05'] = reject

# Sort by p-value
enrichment_df = enrichment_df.sort_values('P_value')

# Save results
enrichment_output = OUTPUT_DIR / 'focused_enrichment_results.tsv'
enrichment_df.to_csv(enrichment_output, sep='\t', index=False)

print(f"Results saved to: {enrichment_output}")
print()

# Print summary with FDR
print("Summary with FDR correction:")
print("="*80)
for _, row in enrichment_df.iterrows():
    sig_fdr = "***" if row['FDR_corrected'] < 0.001 else "**" if row['FDR_corrected'] < 0.01 else "*" if row['FDR_corrected'] < 0.05 else ""
    print(f"{row['Category']:30s}: p={row['P_value']:.4e}, FDR={row['FDR_corrected']:.4e} {sig_fdr}")
print("="*80)



