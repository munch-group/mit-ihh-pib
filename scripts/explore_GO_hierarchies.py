#!/usr/bin/env python3
"""
GO enrichment analysis with functional categories and DAG visualization

This script:
1. Searches for GO terms matching functional categories
2. Performs enrichment analysis for high-confidence genes
3. Visualizes GO term relationships using DAGs

Author: MIT van Bruggen
Date: 2025-11-11
Updated: 2025-11-12
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
OUTPUT_DIR = RESULTS_DIR / 'functional_enrichment/GO_with_terms'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# My study genes (high-confidence X chromosome genes)
study_genes = ['DACH2', 'DIAPH2', 'FAAH2', 'GPC3', 'HTR2C',
               'IL1RAPL2', 'KDM6A', 'NAP1L2', 'NCBP2L',
               'TRPC5', 'ZMAT1', 'KLHL13', 'PCDH11X']

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

print("Step 1: Defining functional categories...") #had some help defining the keywords (Claude)
print()

KEYWORD_CATEGORIES = {
    'fertility': [
        'sperm', 'oocyte', 'spermatogenesis', 'oogenesis','ovary', 'testis', 'meiosis', 'flagellum motor', 'microtubules'
    ],

    'immunity': [
        'immune', 'inflammation','cytokine', 'antibody',
        'lymphocyte','T cell', 'B cell', 'macrophage', 'neutrophil',
    ],
    'neurodevelopment': [
        'neuron','brain', 'neural','neuronal', 'behavior', 'neurotransmitter',
        'dopamine', 'neurodegener', 'neurodevelopmental',
        'sensory', 'psychiatric'
    ],
    'development': [
        'development', 'differentiation', 'morphogenesis', 'organogenesis',
        'embryo',
        'embryonic', 'fetal','angiogenesis', 'vasculogenesis', 'cardiogenesis', 'myogenesis', 'osteogenesis',
        'chondrogenesis', 'adipogenesis', 'hematopoiesis'
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

