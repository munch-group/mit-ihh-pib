#!/usr/bin/env python3
"""
Analyze gene overlaps between GO categories to identify redundancy
"""

import pandas as pd
from pathlib import Path
from itertools import combinations

# Setup paths
PROJECT_DIR = Path('/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib')
HIERARCHY_DIR = PROJECT_DIR / 'results/analysis/functional_enrichment/GO_hierarchies'

# Study genes
study_genes = {'DACH2', 'DIAPH2', 'FAAH2', 'GPC3', 'HTR2C',
               'IL1RAPL2', 'KDM6A', 'NAP1L2', 'NCBP2L',
               'TRPC5', 'ZMAT1', 'KLHL13', 'PCDH11X'}

# Load all categories
categories = {}
for genes_file in HIERARCHY_DIR.glob('genes_unique_*.txt'):
    category = genes_file.stem.replace('genes_unique_', '')
    with open(genes_file, 'r') as f:
        genes = set(line.strip() for line in f if line.strip())
    categories[category] = genes

print("="*80)
print("Category Gene Set Overlaps")
print("="*80)
print()

# Calculate pairwise overlaps
overlaps = []
for cat1, cat2 in combinations(sorted(categories.keys()), 2):
    genes1 = categories[cat1]
    genes2 = categories[cat2]

    overlap = genes1 & genes2
    union = genes1 | genes2

    # Jaccard similarity
    jaccard = len(overlap) / len(union) if union else 0

    # Overlap coefficients
    overlap_coef1 = len(overlap) / len(genes1) if genes1 else 0
    overlap_coef2 = len(overlap) / len(genes2) if genes2 else 0

    # Study gene overlap
    study_overlap1 = genes1 & study_genes
    study_overlap2 = genes2 & study_genes
    study_overlap_both = study_overlap1 & study_overlap2

    overlaps.append({
        'category1': cat1,
        'category2': cat2,
        'size1': len(genes1),
        'size2': len(genes2),
        'overlap_size': len(overlap),
        'jaccard': jaccard,
        'overlap_coef1': overlap_coef1,
        'overlap_coef2': overlap_coef2,
        'study_genes1': len(study_overlap1),
        'study_genes2': len(study_overlap2),
        'study_genes_both': len(study_overlap_both),
        'overlap_genes': ', '.join(sorted(list(overlap)[:10])) + ('...' if len(overlap) > 10 else '')
    })

df = pd.DataFrame(overlaps)
df = df.sort_values('jaccard', ascending=False)

print("\nTop 20 Most Overlapping Category Pairs (by Jaccard similarity):")
print("-"*80)
for _, row in df.head(20).iterrows():
    print(f"{row['category1']:35s} <-> {row['category2']:35s}")
    print(f"  Sizes: {row['size1']:4d} / {row['size2']:4d} genes")
    print(f"  Overlap: {row['overlap_size']:4d} genes (Jaccard: {row['jaccard']:.3f})")
    print(f"  Study genes: {row['study_genes1']}/{row['study_genes2']} (both: {row['study_genes_both']})")
    print()

# Check for highly redundant pairs (>50% overlap)
print("\n" + "="*80)
print("High Redundancy Pairs (Jaccard > 0.3 or one contains >50% of the other):")
print("="*80)
print()

high_redundancy = df[(df['jaccard'] > 0.3) | (df['overlap_coef1'] > 0.5) | (df['overlap_coef2'] > 0.5)]
if len(high_redundancy) > 0:
    for _, row in high_redundancy.iterrows():
        print(f"{row['category1']:35s} <-> {row['category2']:35s}")
        print(f"  Jaccard: {row['jaccard']:.3f}")
        print(f"  {row['category1']} contains {row['overlap_coef1']*100:.1f}% overlap")
        print(f"  {row['category2']} contains {row['overlap_coef2']*100:.1f}% overlap")
        print()
else:
    print("No high redundancy pairs found.")

# Individual category analysis
print("\n" + "="*80)
print("Individual Category Analysis:")
print("="*80)
print()

for category in sorted(categories.keys()):
    genes = categories[category]
    study_overlap = genes & study_genes

    print(f"{category:35s}: {len(genes):4d} genes, {len(study_overlap)} study genes")
    if study_overlap:
        print(f"  Study genes: {', '.join(sorted(study_overlap))}")

print("\n" + "="*80)
print("Study Gene Distribution Across Categories:")
print("="*80)
print()

for gene in sorted(study_genes):
    in_categories = [cat for cat, genes in categories.items() if gene in genes]
    print(f"{gene:15s}: {len(in_categories):2d} categories - {', '.join(in_categories)}")
