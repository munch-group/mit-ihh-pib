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
    'reproduction': [
        'fertility', 'sperm', 'oocyte', 'spermatogenesis', 'oogenesis', 'reproductive',
        'fertilization', 'ovary', 'testis', 'gamete', 'meiosis', 'gonad', 'gametogenesis',
        'ovulation', 'spermiogenesis', 'acrosome', 'zona pellucida', 'follicle',
        'seminiferous', 'leydig', 'sertoli', 'corpus luteum', 'placenta', 'pregnancy',
        'implantation', 'zygote', 'conception', 'sex determination', 'steroidogenesis',
        'testosterone', 'estrogen', 'progesterone', 'luteinizing', 'follicle-stimulating',
        'sexual maturation', 'puberty', 'reproductive behavior', 'mating'
    ],
    'immunity': [
        'immune', 'inflammation', 'defense', 'cytokine', 'antibody',
        'inflammatory', 'innate immune', 'adaptive immune', 'lymphocyte',
        'T cell', 'B cell', 'interferon', 'interleukin', 'antigen',
        'immunoglobulin', 'macrophage', 'dendritic cell', 'neutrophil', 'eosinophil',
        'basophil', 'monocyte', 'natural killer', 'complement', 'chemokine',
        'major histocompatibility', 'MHC', 'HLA', 'toll-like receptor', 'TLR',
        'inflammasome', 'pathogen', 'phagocytosis', 'opsonization', 'apoptosis',
        'immunological', 'leukocyte', 'cytotoxic', 'helper T cell', 'regulatory T cell',
        'plasma cell', 'germinal center', 'thymus', 'spleen', 'lymph node',
        'autoimmun', 'allergy', 'hypersensitivity', 'immunodeficiency'
    ],
    'neurodevelopment': [
        'neuron', 'synapse', 'cognition', 'brain', 'neural', 'axon',
        'dendrite', 'synaptic', 'neuronal', 'nervous system', 'learning',
        'memory', 'behavior', 'neurogenesis', 'neurotransmitter',
        'dopamine', 'serotonin', 'glutamate', 'GABA', 'acetylcholine', 'norepinephrine',
        'neuropeptide', 'receptor', 'ion channel', 'voltage-gated', 'ligand-gated',
        'plasticity', 'long-term potentiation', 'LTP', 'long-term depression', 'LTD',
        'cortex', 'hippocampus', 'amygdala', 'cerebellum', 'striatum', 'thalamus',
        'hypothalamus', 'brainstem', 'spinal cord', 'peripheral nerve',
        'myelination', 'oligodendrocyte', 'astrocyte', 'microglia', 'glia',
        'neuroplasticity', 'neuroprotection', 'neurodegener', 'neurodevelopmental',
        'neural tube', 'neural crest', 'migration', 'neurite', 'growth cone',
        'guidance', 'circuit', 'network', 'sensory', 'motor', 'psychiatric'
    ],
    'development': [
        'development', 'differentiation', 'morphogenesis', 'organogenesis',
        'embryo', 'pattern', 'specification', 'determination',
        'embryonic', 'fetal', 'gastrulation', 'somite', 'mesoderm', 'endoderm', 'ectoderm',
        'germ layer', 'axis', 'anterior-posterior', 'dorsal-ventral', 'left-right',
        'segmentation', 'limb', 'organ', 'tissue', 'cell fate',
        'lineage', 'progenitor', 'stem cell', 'pluripotent', 'multipotent',
        'proliferation', 'apoptosis', 'cell death', 'growth factor', 'signaling pathway',
        'Wnt', 'Notch', 'Hedgehog', 'TGF-beta', 'BMP', 'FGF', 'EGF',
        'transcription factor', 'HOX', 'homeobox', 'morphogen', 'gradient',
        'induction', 'competence', 'commitment', 'maturation', 'remodeling',
        'angiogenesis', 'vasculogenesis', 'cardiogenesis', 'myogenesis', 'osteogenesis',
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

#==============================================================================
# Step 5: Define background gene set - All protein-coding X chromosome genes
#==============================================================================

print("\n\nStep 5: Defining background gene set...")
print()

# Get all protein-coding genes, then filter for X chromosome
all_protein_coding = go.all_protein_coding()

# Get X chromosome gene symbols from gene annotation table
gene_table = go.gene_annotation_table()
x_gene_symbols = set(gene_table[gene_table['chromosome'] == 'X']['Symbol'].tolist())

# Filter protein-coding genes to only those on X chromosome
all_x_protein_coding = [gene for gene in all_protein_coding if gene in x_gene_symbols]

background_genes = all_x_protein_coding
print(f"  Background: All protein-coding X chromosome genes")
print(f"  Total background genes: {len(background_genes)}")
print()

#==============================================================================
# Step 6: Filter category genes to only include X chromosome genes
#==============================================================================

print("Step 6: Filtering category genes to X chromosome only...")
print()

# Only keep category genes that are on X chromosome (in background)
background_set = set(background_genes)

neurodevelopment_x = [g for g in neurodevelopment_genes if g in background_set]
reproduction_x = [g for g in reproduction_genes if g in background_set]
immunity_x = [g for g in immunity_genes if g in background_set]
development_x = [g for g in development_genes if g in background_set]

print(f"  Neurodevelopment (X only): {len(neurodevelopment_x)} genes")
print(f"  Reproduction (X only):     {len(reproduction_x)} genes")
print(f"  Immunity (X only):         {len(immunity_x)} genes")
print(f"  Development (X only):      {len(development_x)} genes")
print()

#==============================================================================
# Step 7: Perform Fisher's exact tests for enrichment
#==============================================================================

print("Step 7: Performing Fisher's exact tests for enrichment...")
print("Using GeneList.fisher() method as described by https://munch-group.org/geneinfo/pages/gene_lists.html")
print()

# Convert study genes to GeneList object
study_gl = glist(study_genes)

# Create GeneList objects for X-filtered categories
neurodevelopment_gl_x = glist(neurodevelopment_x)
reproduction_gl_x = glist(reproduction_x)
immunity_gl_x = glist(immunity_x)
development_gl_x = glist(development_x)

# Background: all protein-coding X chromosome genes
background_gl = glist(background_genes)

# Store results
enrichment_results = []

# Define functional categories
categories = {
    'neurodevelopment': (neurodevelopment_gl_x, neurodevelopment_x),
    'reproduction': (reproduction_gl_x, reproduction_x),
    'immunity': (immunity_gl_x, immunity_x),
    'development': (development_gl_x, development_x)
}

# Perform Fisher's exact test for each category
# Using: study_genes.fisher(category_genes, background=background_genes)
print("Results:")
print("-" * 80)
for category_name, (category_gl, category_genes_list) in categories.items():
    # Perform Fisher's exact test
    p_value = study_gl.fisher(category_gl, background=background_gl)

    # Calculate overlap
    study_set = set(study_genes)
    category_set = set(category_genes_list)
    overlap_genes = list(study_set & category_set)

    # Calculate odds ratio manually for reporting
    a = len(study_set & category_set)  # In study AND in category
    b = len(study_set - category_set)  # In study but NOT in category
    c = len(category_set - study_set)  # In category but NOT in study
    d = len(set(background_genes) - study_set - category_set)  # In neither

    # Odds ratio = (a/b) / (c/d) = (a*d) / (b*c)
    if b > 0 and c > 0:
        odds_ratio = (a * d) / (b * c)
    else:
        odds_ratio = float('inf') if a > 0 else 0.0

    # Store results
    enrichment_results.append({
        'Category': category_name,
        'Study_genes_total': len(study_genes),
        'Category_genes_total': len(category_genes_list),
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

    # Print results with significance markers
    sig = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else ""
    print(f"{category_name:20s}: OR={odds_ratio:6.2f}, p={p_value:.4e} {sig:3s}")
    print(f"  Overlap: {len(overlap_genes)}/{len(study_genes)} study genes")
    if overlap_genes:
        print(f"  Genes: {', '.join(sorted(overlap_genes))}")
    print()

print("-" * 80)

# Create results dataframe
enrichment_df = pd.DataFrame(enrichment_results)

# Apply FDR correction (Benjamini-Hochberg)
from statsmodels.stats.multitest import multipletests
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
enrichment_output = OUTPUT_DIR / 'fisher_exact_test_results.tsv'
enrichment_df.to_csv(enrichment_output, sep='\t', index=False)

print(f"\nResults saved to: {enrichment_output}")
print()

# Print summary table
print("Summary with FDR correction:")
print("=" * 80)
for _, row in enrichment_df.iterrows():
    sig_fdr = "***" if row['FDR_corrected'] < 0.001 else "**" if row['FDR_corrected'] < 0.01 else "*" if row['FDR_corrected'] < 0.05 else ""
    print(f"{row['Category']:20s}: p={row['P_value']:.4e}, FDR={row['FDR_corrected']:.4e} {sig_fdr}")

print("=" * 80)

#==============================================================================
# Step 8: Test enrichment for combined/intersecting categories
#==============================================================================

print("\n\nStep 8: Testing enrichment for combined functional categories...")
print("(genes that belong to multiple categories simultaneously)")
print()

# Create intersections using GeneList set operations (&)
combined_categories = {
    'reproduction_AND_development': (
        reproduction_gl_x & development_gl_x,
        'Reproduction AND Development'
    ),
    'reproduction_AND_neurodevelopment': (
        reproduction_gl_x & neurodevelopment_gl_x,
        'Reproduction AND Neurodevelopment'
    ),
    'reproduction_AND_immunity': (
        reproduction_gl_x & immunity_gl_x,
        'Reproduction AND Immunity'
    ),
    'neurodevelopment_AND_immunity': (
        neurodevelopment_gl_x & immunity_gl_x,
        'Neurodevelopment AND Immunity'
    ),
    'neurodevelopment_AND_development': (
        neurodevelopment_gl_x & development_gl_x,
        'Neurodevelopment AND Development'
    ),
    'immunity_AND_development': (
        immunity_gl_x & development_gl_x,
        'Immunity AND Development'
    ),
    # Triple combinations
    'reproduction_AND_neurodevelopment_AND_development': (
        reproduction_gl_x & neurodevelopment_gl_x & development_gl_x,
        'Reproduction AND Neurodevelopment AND Development'
    ),
    'reproduction_AND_immunity_AND_development': (
        reproduction_gl_x & immunity_gl_x & development_gl_x,
        'Reproduction AND Immunity AND Development'
    ),
    'neurodevelopment_AND_immunity_AND_development': (
        neurodevelopment_gl_x & immunity_gl_x & development_gl_x,
        'Neurodevelopment AND Immunity AND Development'
    ),
    'reproduction_AND_neurodevelopment_AND_immunity': (
        reproduction_gl_x & neurodevelopment_gl_x & immunity_gl_x,
        'Reproduction AND Neurodevelopment AND Immunity'
    ),
    # Quadruple combination - all four categories
    'reproduction_AND_neurodevelopment_AND_immunity_AND_development': (
        reproduction_gl_x & neurodevelopment_gl_x & immunity_gl_x & development_gl_x,
        'Reproduction AND Neurodevelopment AND Immunity AND Development (ALL FOUR)'
    ),
}

# Store combined results
combined_results = []

print("Combined category sizes (X chromosome only):")
for cat_name, (cat_gl, display_name) in combined_categories.items():
    print(f"  {display_name:55s}: {len(cat_gl):3d} genes")
print()

print("Results:")
print("-" * 80)

for cat_name, (category_gl, display_name) in combined_categories.items():
    # Skip if category is empty or too small
    if len(category_gl) == 0:
        print(f"{display_name:55s}: No genes in this combination (skipped)")
        print()
        continue

    # Perform Fisher's exact test using GeneList.fisher()
    p_value = study_gl.fisher(category_gl, background=background_gl)

    # Calculate overlap
    study_set = set(study_genes)
    cat_genes_set = set(list(category_gl))
    overlap_genes = list(study_set & cat_genes_set)

    # Calculate odds ratio
    a = len(study_set & cat_genes_set)
    b = len(study_set - cat_genes_set)
    c = len(cat_genes_set - study_set)
    d = len(set(background_genes) - study_set - cat_genes_set)

    if b > 0 and c > 0:
        odds_ratio = (a * d) / (b * c)
    else:
        odds_ratio = float('inf') if a > 0 else 0.0

    # Store results
    combined_results.append({
        'Category': display_name,
        'Study_genes_total': len(study_genes),
        'Category_genes_total': len(category_gl),
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
    print(f"{display_name:55s}: OR={odds_ratio:6.2f}, p={p_value:.4e} {sig:3s}")
    print(f"  Overlap: {len(overlap_genes)}/{len(study_genes)} study genes (category size: {len(category_gl)} genes)")
    if overlap_genes:
        print(f"  Genes: {', '.join(sorted(overlap_genes))}")
    print()

print("-" * 80)

# Create combined results dataframe
if len(combined_results) > 0:
    combined_df = pd.DataFrame(combined_results)

    # Apply FDR correction
    reject, pvals_corrected, _, _ = multipletests(
        combined_df['P_value'],
        alpha=0.05,
        method='fdr_bh'
    )

    combined_df['FDR_corrected'] = pvals_corrected
    combined_df['Significant_FDR_0.05'] = reject

    # Sort by p-value
    combined_df = combined_df.sort_values('P_value')

    # Save combined results
    combined_output = OUTPUT_DIR / 'fisher_exact_test_combined_categories.tsv'
    combined_df.to_csv(combined_output, sep='\t', index=False)

    print(f"\nCombined category results saved to: {combined_output}")
    print()

    # Print summary with FDR
    print("Summary with FDR correction (combined categories):")
    print("=" * 80)
    for _, row in combined_df.iterrows():
        sig_fdr = "***" if row['FDR_corrected'] < 0.001 else "**" if row['FDR_corrected'] < 0.01 else "*" if row['FDR_corrected'] < 0.05 else ""
        print(f"{row['Category']:55s}: p={row['P_value']:.4e}, FDR={row['FDR_corrected']:.4e} {sig_fdr}")
    print("=" * 80)

print("\nAnalysis complete!")

