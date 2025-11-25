#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Retrieve detailed information for high-confidence candidate genes using MyGene API.

This script:
1. Reads the list of high-confidence genes
2. Uses geneinfo's mygene_get_gene_info to fetch structured gene data
3. Extracts key details: function, location, aliases, type
4. Classifies genes by function (reproductive, immune, neural, other)
5. Saves detailed gene information table

Author: MIT van Bruggen
Date: October 29, 2025
"""

import pandas as pd
import numpy as np
from pathlib import Path
# Import the internal function that returns structured data
from geneinfo.information import mygene_get_gene_info
from geneinfo.coords import gene_coords
import sys

# Input/output paths
PROJECT_DIR = Path('/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib')
GENE_LIST_FILE = PROJECT_DIR / 'results/ihs_fixed_par1/analysis/high_significance_genes/high_confidence_X_genes.txt'
OUTPUT_DIR = PROJECT_DIR / 'results/ihs_fixed_par1/analysis/high_significance_genes'
OUTPUT_FILE = OUTPUT_DIR / 'gene_details_complete.tsv'

# Genome assembly
ASSEMBLY = 'hg38'

# Keywords for functional classification
REPRODUCTIVE_KEYWORDS = [
    'gonad', 'gamete', 'sperm', 'oocyte', 'testis', 'ovary', 'fertility',
    'meiosis', 'spermatogenesis', 'oogenesis', 'fertilization', 'reproduction',
    'sex', 'sexual', 'germ cell', 'embryo'
]

IMMUNE_KEYWORDS = [
    'immune', 'immunity', 'inflammation', 'inflammatory', 'cytokine',
    'interferon', 'antibody', 'antigen', 'T cell', 'B cell', 'lymphocyte',
    'defense', 'innate', 'adaptive', 'toll-like', 'MHC', 'HLA', 'interleukin'
]

NEURAL_KEYWORDS = [
    'neuron', 'neural', 'brain', 'synapse', 'synaptic', 'cognitive',
    'neurodevelopment', 'neurological', 'receptor', 'serotonin', 'dopamine',
    'neurotransmitter', 'axon', 'dendrite', 'cerebral', 'hippocampus'
]

DEVELOPMENTAL_KEYWORDS = [
    'development', 'developmental', 'morphogenesis', 'differentiation',
    'organogenesis', 'pattern', 'axis', 'embryonic'
]


def read_gene_list(file_path):
    """Read gene list from file."""
    with open(file_path, 'r') as f:
        genes = [line.strip() for line in f if line.strip()]
    return genes


def classify_function(summary, name='', aliases=''):
    """
    Classify gene function based on keywords.
    Returns list of categories.
    """
    categories = []
    text = (str(summary) + ' ' + str(name) + ' ' + str(aliases)).lower()

    if any(keyword in text for keyword in REPRODUCTIVE_KEYWORDS):
        categories.append('reproductive')

    if any(keyword in text for keyword in IMMUNE_KEYWORDS):
        categories.append('immune')

    if any(keyword in text for keyword in NEURAL_KEYWORDS):
        categories.append('neural')

    if any(keyword in text for keyword in DEVELOPMENTAL_KEYWORDS):
        categories.append('developmental')

    if not categories:
        categories.append('other')

    return '; '.join(categories)


def parse_gene_info(gene_name):
    """
    Fetch and parse gene information for a single gene using MyGene API.
    Returns dict with gene details.
    """
    try:
        # Use the mygene API via geneinfo
        info = mygene_get_gene_info(
            gene_name,
            species='human',
            scopes='symbol',
            fields='symbol,alias,name,type_of_gene,summary,genomic_pos,genomic_pos_hg19'
        )

        if info is None:
            raise Exception(f"Gene {gene_name} not found in MyGene database")

        # Extract fields
        symbol = info.get('symbol', gene_name)
        full_name = info.get('name', '')
        gene_type = info.get('type_of_gene', 'unknown')
        summary = info.get('summary', 'No summary available')

        # Handle aliases (can be string or list)
        aliases = info.get('alias', [])
        if isinstance(aliases, str):
            aliases = [aliases]
        elif aliases is None:
            aliases = []
        aliases_str = ', '.join(aliases) if aliases else ''

        # Get genomic position (hg38)
        genomic_pos = info.get('genomic_pos', {})
        if isinstance(genomic_pos, list):
            genomic_pos = genomic_pos[0] if genomic_pos else {}

        chrom = genomic_pos.get('chr', '')
        start = genomic_pos.get('start', '')
        end = genomic_pos.get('end', '')
        strand = genomic_pos.get('strand', '')

        # Calculate gene length
        gene_length = np.nan
        if start and end:
            try:
                gene_length = int(end) - int(start)
            except:
                pass

        # Classify function
        func_category = classify_function(summary, full_name, aliases_str)

        # Truncate long fields
        summary_short = summary[:300] + '...' if len(str(summary)) > 300 else summary

        return {
            'gene_symbol': symbol,
            'gene_type': gene_type,
            'full_name': full_name,
            'summary': summary_short,
            'chromosome': chrom,
            'start': start,
            'end': end,
            'strand': strand,
            'gene_length': gene_length,
            'aliases': aliases_str,
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
            'strand': '',
            'gene_length': np.nan,
            'aliases': '',
            'functional_category': 'unknown',
            'status': 'failed'
        }


def main():
    """Main analysis pipeline."""
    print(f"\n{'='*70}")
    print(f"Fetching Detailed Gene Information (MyGene API)")
    print(f"{'='*70}\n")

    # Read gene list
    print(f"Reading gene list from: {GENE_LIST_FILE}")
    genes = read_gene_list(GENE_LIST_FILE)
    print(f"  - Found {len(genes)} genes\n")

    # Fetch information for each gene
    print("Fetching gene information from MyGene.info API...")
    gene_details = []

    for idx, gene in enumerate(genes, 1):
        print(f"  Processing {idx}/{len(genes)}: {gene:20s}", end='')
        details = parse_gene_info(gene)
        gene_details.append(details)
        print(f" [{details['status']}] - {details['gene_type']}")

    # Create dataframe
    df = pd.DataFrame(gene_details)

    # Sort by chromosome position
    df['chrom_sort'] = df['chromosome'].astype(str).str.replace('X', '23').str.replace('Y', '24')
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

    # Success rate
    success = len(df[df['status'] == 'success'])
    failed = len(df[df['status'] == 'failed'])
    print(f"Data Retrieval: {success}/{len(df)} successful ({100*success/len(df):.1f}%)")
    if failed > 0:
        failed_genes = ', '.join(df[df['status'] == 'failed']['gene_symbol'].tolist())
        print(f"Failed genes: {failed_genes}\n")

    # Filter to successful retrievals for analysis
    df_success = df[df['status'] == 'success']

    if len(df_success) == 0:
        print("No genes successfully retrieved!")
        return df

    # 1. Gene types
    print("\n1. Gene Types:")
    type_counts = df_success['gene_type'].value_counts()
    for gene_type, count in type_counts.items():
        pct = 100 * count / len(df_success)
        print(f"   - {gene_type:20s}: {count:2d} ({pct:5.1f}%)")

    # 2. Functional categories
    print("\n2. Functional Categories:")
    all_categories = []
    for cats in df_success['functional_category']:
        all_categories.extend([c.strip() for c in str(cats).split(';')])

    from collections import Counter
    cat_counts = Counter(all_categories)
    for category, count in cat_counts.most_common():
        pct = 100 * count / len(df_success)
        print(f"   - {category:20s}: {count:2d} ({pct:5.1f}%)")

    # 3. Gene length distribution
    print("\n3. Gene Length Distribution:")
    if df_success['gene_length'].notna().sum() > 0:
        lengths = df_success['gene_length'].dropna()
        print(f"   - Mean length:   {lengths.mean():>10,.0f} bp")
        print(f"   - Median length: {lengths.median():>10,.0f} bp")
        print(f"   - Min length:    {lengths.min():>10,.0f} bp ({df_success[df_success['gene_length'] == lengths.min()]['gene_symbol'].values[0]})")
        print(f"   - Max length:    {lengths.max():>10,.0f} bp ({df_success[df_success['gene_length'] == lengths.max()]['gene_symbol'].values[0]})")

    # 4. Reproductive genes
    print("\n4. Reproductive Function Genes:")
    repro_genes = df_success[df_success['functional_category'].str.contains('reproductive', na=False)]
    if len(repro_genes) > 0:
        for _, row in repro_genes.iterrows():
            print(f"   * {row['gene_symbol']:15s} - {row['full_name'][:50]}")
    else:
        print("   - None identified")

    # 5. Immune genes
    print("\n5. Immune Function Genes:")
    immune_genes = df_success[df_success['functional_category'].str.contains('immune', na=False)]
    if len(immune_genes) > 0:
        for _, row in immune_genes.iterrows():
            print(f"   * {row['gene_symbol']:15s} - {row['full_name'][:50]}")
    else:
        print("   - None identified")

    # 6. Neural genes
    print("\n6. Neural/Behavioral Genes:")
    neural_genes = df_success[df_success['functional_category'].str.contains('neural', na=False)]
    if len(neural_genes) > 0:
        for _, row in neural_genes.iterrows():
            print(f"   * {row['gene_symbol']:15s} - {row['full_name'][:50]}")
    else:
        print("   - None identified")

    # 7. Developmental genes
    print("\n7. Developmental Genes:")
    dev_genes = df_success[df_success['functional_category'].str.contains('developmental', na=False)]
    if len(dev_genes) > 0:
        for _, row in dev_genes.iterrows():
            print(f"   * {row['gene_symbol']:15s} - {row['full_name'][:50]}")
    else:
        print("   - None identified")

    # 8. Key highlights
    print("\n8. KEY CANDIDATE GENES TO INVESTIGATE:")
    key_genes = df_success[df_success['functional_category'].str.contains('reproductive|immune|neural', na=False)]
    if len(key_genes) > 0:
        print(f"\n   {len(key_genes)} genes with reproductive, immune, or neural function:\n")
        for _, row in key_genes.iterrows():
            print(f"   {row['gene_symbol']:15s} ({row['functional_category']})")
            print(f"      Location: chr{row['chromosome']}:{row['start']}-{row['end']}")
            print(f"      Function: {row['full_name']}")
            summary_preview = row['summary'][:150].replace('\n', ' ')
            print(f"      Summary:  {summary_preview}...")
            print()
    else:
        print("   - No obvious reproductive/immune/neural candidates")

    print(f"{'='*70}")
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
