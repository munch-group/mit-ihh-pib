#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Retrieve detailed information for high-confidence candidate genes.

This script:
1. Reads the list of high-confidence genes
2. Uses geneinfo to fetch gene information for each
3. Extracts key details: function, location, aliases, type
4. Saves detailed gene information table
5. Highlights genes related to reproduction/immunity

Author: MIT van Bruggen
Date: October 29, 2025
"""

import pandas as pd
import numpy as np
from pathlib import Path
import geneinfo.information as gi
import sys
import re

# Input/output paths
PROJECT_DIR = Path('/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib')
GENE_LIST_FILE = PROJECT_DIR / 'results/analysis/high_significance_genes/high_confidence_X_genes.txt'
OUTPUT_DIR = PROJECT_DIR / 'results/analysis/high_significance_genes'
OUTPUT_FILE = OUTPUT_DIR / 'gene_details.tsv'

# Genome assembly
ASSEMBLY = 'hg38'

# Keywords for functional classification
REPRODUCTIVE_KEYWORDS = [
    'gonad', 'gamete', 'sperm', 'oocyte', 'testis', 'ovary', 'fertility',
    'meiosis', 'spermatogenesis', 'oogenesis', 'fertilization', 'reproduction',
    'sex', 'sexual', 'germ cell'
]

IMMUNE_KEYWORDS = [
    'immune', 'immunity', 'inflammation', 'inflammatory', 'cytokine',
    'interferon', 'antibody', 'antigen', 'T cell', 'B cell', 'lymphocyte',
    'defense', 'innate', 'adaptive', 'toll-like', 'MHC', 'HLA'
]

NEURAL_KEYWORDS = [
    'neuron', 'neural', 'brain', 'synapse', 'synaptic', 'cognitive',
    'neurodevelopment', 'neurological', 'receptor', 'serotonin', 'dopamine',
    'neurotransmitter', 'axon', 'dendrite'
]


def read_gene_list(file_path):
    """Read gene list from file."""
    with open(file_path, 'r') as f:
        genes = [line.strip() for line in f if line.strip()]
    return genes


def classify_function(summary, aliases=''):
    """
    Classify gene function based on keywords in summary.

    Returns list of categories.
    """
    categories = []
    text = (summary + ' ' + aliases).lower()

    if any(keyword in text for keyword in REPRODUCTIVE_KEYWORDS):
        categories.append('reproductive')

    if any(keyword in text for keyword in IMMUNE_KEYWORDS):
        categories.append('immune')

    if any(keyword in text for keyword in NEURAL_KEYWORDS):
        categories.append('neural')

    if not categories:
        categories.append('other')

    return ', '.join(categories)


def parse_gene_info(gene_name):
    """
    Fetch and parse gene information for a single gene.

    Returns dict with gene details.
    """
    try:
        # Get gene information
        info = gi.gene_info(gene_name)

        # The gi.gene_info returns a formatted string, we need to parse it
        # Let's capture it by redirecting output
        import io
        from contextlib import redirect_stdout

        f = io.StringIO()
        with redirect_stdout(f):
            gi.gene_info(gene_name)
        info_text = f.getvalue()

        # Parse the output
        lines = info_text.strip().split('\n')

        # Extract components
        symbol_line = lines[0] if len(lines) > 0 else ''
        full_name = lines[1] if len(lines) > 1 else ''
        summary = lines[2] if len(lines) > 2 else ''
        location = ''

        # Parse symbol line for type and aliases
        gene_type = 'unknown'
        aliases = ''
        if 'protein-coding' in symbol_line:
            gene_type = 'protein-coding'
        elif 'lncRNA' in symbol_line or 'long non-coding' in symbol_line.lower():
            gene_type = 'lncRNA'
        elif 'ncRNA' in symbol_line:
            gene_type = 'ncRNA'

        if 'Aliases:' in symbol_line:
            aliases = symbol_line.split('Aliases:')[1].strip()

        # Extract genomic position
        for line in lines:
            if 'Human genomic position:' in line:
                location = line.replace('Human genomic position:', '').strip()
                # Extract hg38 coordinates
                if 'hg38' in location:
                    location = location.split('(hg38)')[0].strip()
                break

        # Parse chromosome and coordinates from location
        chrom, start, end = '', '', ''
        if location and ':' in location:
            parts = location.split(':')
            chrom = parts[0]
            if len(parts) > 1 and '-' in parts[1]:
                coords = parts[1].split('-')
                start = coords[0]
                end = coords[1]

        # Get gene coordinates using geneinfo
        try:
            coords = gi.gene_coords([gene_name], assembly=ASSEMBLY)
            if coords:
                chrom = coords[0][0]
                start = coords[0][1]
                end = coords[0][2]
        except:
            pass

        # Classify function
        func_category = classify_function(summary + ' ' + full_name, aliases)

        return {
            'gene_symbol': gene_name,
            'gene_type': gene_type,
            'full_name': full_name[:200] if len(full_name) <= 200 else full_name[:197] + '...',
            'summary': summary[:500] if len(summary) <= 500 else summary[:497] + '...',
            'chromosome': chrom,
            'start': start,
            'end': end,
            'gene_length': int(end) - int(start) if start and end else np.nan,
            'aliases': aliases,
            'functional_category': func_category,
            'status': 'success'
        }

    except Exception as e:
        print(f"  WARNING: Could not fetch info for {gene_name}: {e}")
        return {
            'gene_symbol': gene_name,
            'gene_type': 'unknown',
            'full_name': 'N/A',
            'summary': f'Error: {str(e)}',
            'chromosome': '',
            'start': '',
            'end': '',
            'gene_length': np.nan,
            'aliases': '',
            'functional_category': 'unknown',
            'status': 'failed'
        }


def main():
    """Main analysis pipeline."""
    print(f"\n{'='*70}")
    print(f"Fetching Detailed Gene Information")
    print(f"{'='*70}\n")

    # Read gene list
    print(f"Reading gene list from: {GENE_LIST_FILE}")
    genes = read_gene_list(GENE_LIST_FILE)
    print(f"  - Found {len(genes)} genes\n")

    # Fetch information for each gene
    print("Fetching gene information from Ensembl/NCBI...")
    gene_details = []

    for idx, gene in enumerate(genes, 1):
        print(f"  Processing {idx}/{len(genes)}: {gene}")
        details = parse_gene_info(gene)
        gene_details.append(details)

    # Create dataframe
    df = pd.DataFrame(gene_details)

    # Sort by chromosome position
    df['chrom_sort'] = df['chromosome'].str.replace('chr', '').str.replace('X', '23').str.replace('Y', '24')
    df['chrom_sort'] = pd.to_numeric(df['chrom_sort'], errors='coerce')
    df['start_int'] = pd.to_numeric(df['start'], errors='coerce')
    df = df.sort_values(['chrom_sort', 'start_int']).reset_index(drop=True)
    df = df.drop(['chrom_sort', 'start_int'], axis=1)

    # Save results
    print(f"\nSaving results to: {OUTPUT_FILE}")
    df.to_csv(OUTPUT_FILE, sep='\t', index=False)

    # Generate summary
    print(f"\n{'='*70}")
    print("GENE INFORMATION SUMMARY")
    print(f"{'='*70}\n")

    # 1. Gene types
    print("1. Gene Types:")
    type_counts = df['gene_type'].value_counts()
    for gene_type, count in type_counts.items():
        print(f"   - {gene_type:20s}: {count}")

    # 2. Functional categories
    print("\n2. Functional Categories:")

    # Expand categories (some genes may have multiple)
    all_categories = []
    for cats in df['functional_category']:
        all_categories.extend([c.strip() for c in str(cats).split(',')])

    from collections import Counter
    cat_counts = Counter(all_categories)
    for category, count in cat_counts.most_common():
        print(f"   - {category:20s}: {count} genes")

    # 3. Reproductive genes
    print("\n3. Reproductive Function Genes:")
    repro_genes = df[df['functional_category'].str.contains('reproductive', na=False)]
    if len(repro_genes) > 0:
        for _, row in repro_genes.iterrows():
            print(f"   - {row['gene_symbol']:15s}: {row['full_name'][:60]}")
    else:
        print("   - None identified")

    # 4. Immune genes
    print("\n4. Immune Function Genes:")
    immune_genes = df[df['functional_category'].str.contains('immune', na=False)]
    if len(immune_genes) > 0:
        for _, row in immune_genes.iterrows():
            print(f"   - {row['gene_symbol']:15s}: {row['full_name'][:60]}")
    else:
        print("   - None identified")

    # 5. Neural genes
    print("\n5. Neural Function Genes:")
    neural_genes = df[df['functional_category'].str.contains('neural', na=False)]
    if len(neural_genes) > 0:
        for _, row in neural_genes.iterrows():
            print(f"   - {row['gene_symbol']:15s}: {row['full_name'][:60]}")
    else:
        print("   - None identified")

    # 6. Gene length distribution
    print("\n6. Gene Length Distribution:")
    if df['gene_length'].notna().sum() > 0:
        print(f"   - Mean length:   {df['gene_length'].mean():>10,.0f} bp")
        print(f"   - Median length: {df['gene_length'].median():>10,.0f} bp")
        print(f"   - Min length:    {df['gene_length'].min():>10,.0f} bp")
        print(f"   - Max length:    {df['gene_length'].max():>10,.0f} bp")

    # 7. Highlight key candidates
    print("\n7. Key Candidate Genes (Reproductive/Immune):")
    key_genes = df[df['functional_category'].str.contains('reproductive|immune', na=False)]
    if len(key_genes) > 0:
        print(f"\n   {len(key_genes)} genes with reproductive or immune function:\n")
        for _, row in key_genes.iterrows():
            print(f"   {row['gene_symbol']:15s} ({row['functional_category']})")
            print(f"   {row['chromosome']}:{row['start']}-{row['end']}")
            print(f"   {row['full_name']}")
            print()

    # 8. Success rate
    print(f"\n8. Data Retrieval Success:")
    success = len(df[df['status'] == 'success'])
    failed = len(df[df['status'] == 'failed'])
    print(f"   - Successfully retrieved: {success}/{len(df)} ({100*success/len(df):.1f}%)")
    if failed > 0:
        print(f"   - Failed: {failed}")
        print(f"   - Failed genes: {', '.join(df[df['status'] == 'failed']['gene_symbol'].tolist())}")

    print(f"\n{'='*70}")
    print("Analysis complete!")
    print(f"{'='*70}\n")

    return df


if __name__ == '__main__':
    try:
        df = main()
    except Exception as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)
