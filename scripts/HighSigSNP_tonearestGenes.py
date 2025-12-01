#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Map high-significance SNPs to nearest genes using geneinfo module.

This script:
1. Reads SNP positions from chrX_high_significance.tsv
2. For each SNP, queries genes in surrounding region using geneinfo
3. Calculates distance to nearest gene(s)
4. Determines overlap type (inside gene, promoter, distal)
5. Saves results and generates summary statistics

Author: MIT van Bruggen
Date: October 29, 2025
"""

import pandas as pd
import numpy as np
from pathlib import Path
from geneinfo.coords import gene_coords_region, all_gene_coords
import sys

# Configuration
ASSEMBLY = 'hg38'
SEARCH_WINDOW = 500000  # Search �500kb around each SNP
PROMOTER_DISTANCE = 2000  # Promoter defined as <2kb from TSS
ENHANCER_DISTANCE = 50000  # Enhancer/regulatory region <50kb

# Input/output paths
PROJECT_DIR = Path('/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib')
INPUT_FILE = PROJECT_DIR / 'results/analysis/candidates/per_chromosome/chrX_candidates_stringent.tsv'
OUTPUT_DIR = PROJECT_DIR / 'results/ihs_revised_HCsnps_logic/analysis/high_significance_genes'
OUTPUT_FILE = OUTPUT_DIR / 'snp_to_gene_mapping.tsv'

# Create output directory
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def get_genes_in_region(chrom, pos, window=SEARCH_WINDOW, assembly=ASSEMBLY):
    """
    Get all genes within a window around a position.

    Parameters
    ----------
    chrom : str
        Chromosome (e.g., 'X' or 'chrX')
    pos : int
        Position on chromosome
    window : int
        Window size (�window around pos)
    assembly : str
        Genome assembly (default: 'hg38')

    Returns
    -------
    list
        List of gene information tuples
    """
    # Ensure chromosome format is 'chrX'
    if not chrom.startswith('chr'):
        chrom = f'chr{chrom}'

    start = max(0, pos - window)
    end = pos + window

    # Query genes using geneinfo
    genes = gene_coords_region(chrom=chrom, start=start, end=end, assembly=assembly)

    return genes


def calculate_distance_to_gene(snp_pos, gene_start, gene_end, gene_strand):
    """
    Calculate distance from SNP to gene.

    Returns distance and overlap type:
    - If SNP is inside gene: distance=0, type='inside_gene'
    - If SNP is upstream/downstream: calculate distance to TSS/TES

    Parameters
    ----------
    snp_pos : int
        SNP position
    gene_start : int
        Gene start position
    gene_end : int
        Gene end position
    gene_strand : str
        Gene strand ('+' or '-')

    Returns
    -------
    tuple
        (distance, overlap_type, position_relative_to_gene)
    """
    # Check if SNP is inside gene
    if gene_start <= snp_pos <= gene_end:
        return 0, 'inside_gene', 'intronic/exonic'

    # Determine TSS (transcription start site)
    if gene_strand == '+':
        tss = gene_start
        tes = gene_end
    else:  # strand == '-'
        tss = gene_end
        tes = gene_start

    # Calculate distance to TSS and TES
    dist_to_tss = abs(snp_pos - tss)
    dist_to_tes = abs(snp_pos - tes)

    # Determine position relative to gene
    if snp_pos < gene_start:
        if gene_strand == '+':
            position = 'upstream'
            distance = dist_to_tss
        else:
            position = 'downstream'
            distance = dist_to_tes
    else:  # snp_pos > gene_end
        if gene_strand == '+':
            position = 'downstream'
            distance = dist_to_tes
        else:
            position = 'upstream'
            distance = dist_to_tss

    # Determine overlap type based on distance
    if distance < PROMOTER_DISTANCE:
        overlap_type = 'promoter'
    elif distance < ENHANCER_DISTANCE:
        overlap_type = 'regulatory_region'
    else:
        overlap_type = 'distal'

    return distance, overlap_type, position


def map_snp_to_genes(snp_row):
    """
    Map a single SNP to nearest gene(s).

    Parameters
    ----------
    snp_row : pd.Series
        Row from SNP dataframe

    Returns
    -------
    list
        List of dicts with gene mapping information
    """
    chrom = snp_row['chr']
    pos = snp_row['pos']

    # Get genes in region
    genes = get_genes_in_region(chrom, pos)

    if not genes:
        return [{
            'snp_location': snp_row['Location'],
            'chr': chrom,
            'pos': pos,
            'nearest_gene': 'NO_GENE_FOUND',
            'gene_start': np.nan,
            'gene_end': np.nan,
            'gene_strand': np.nan,
            'distance': np.nan,
            'overlap_type': 'no_nearby_gene',
            'position_relative': 'N/A',
            'neg_log10_p': snp_row['neg_log10_p'],
            'std_ihs': snp_row['Std iHS']
        }]

    # Calculate distance to each gene
    results = []
    for gene_info in genes:
        # Debug: print first gene_info to see structure
        if len(results) == 0:
            print(f"\n  DEBUG: gene_info structure: {gene_info}", file=sys.stderr)

        # gene_coords_region returns: (gene_name, chrom, start, end, strand, ...)
        # But the actual structure needs to be determined
        try:
            if len(gene_info) >= 5 and isinstance(gene_info[2], (int, str)):
                # Structure: (name, chrom, start, end, strand, ...)
                gene_name = gene_info[0]
                gene_start = int(gene_info[2])
                gene_end = int(gene_info[3])
                gene_strand = gene_info[4] if len(gene_info) > 4 else '+'
            else:
                # Try alternative structure: (name, start, end, strand)
                gene_name = gene_info[0]
                gene_start = int(gene_info[1])
                gene_end = int(gene_info[2])
                gene_strand = gene_info[3] if len(gene_info) > 3 else '+'
        except (ValueError, IndexError) as e:
            print(f"\n  WARNING: Could not parse gene_info {gene_info}: {e}", file=sys.stderr)
            continue

        distance, overlap_type, position_relative = calculate_distance_to_gene(
            int(pos), gene_start, gene_end, gene_strand  # Ensure pos is also int
        )

        results.append({
            'snp_location': snp_row['Location'],
            'chr': chrom,
            'pos': pos,
            'nearest_gene': gene_name,
            'gene_start': gene_start,
            'gene_end': gene_end,
            'gene_strand': gene_strand,
            'distance': distance,
            'overlap_type': overlap_type,
            'position_relative': position_relative,
            'neg_log10_p': snp_row['neg_log10_p'],
            'std_ihs': snp_row['Std iHS']
        })

    # Sort by distance to find nearest
    results.sort(key=lambda x: x['distance'])

    return results


def main():
    """Main analysis pipeline."""
    print(f"\n{'='*70}")
    print(f"Mapping High-Significance SNPs to Nearest Genes")
    print(f"{'='*70}\n")

    # Read SNP data
    print(f"Reading SNP data from: {INPUT_FILE}")
    snps_df = pd.read_csv(INPUT_FILE, sep='\t')
    print(f"  - Found {len(snps_df)} high-significance SNPs\n")

    # Process each SNP
    print("Mapping SNPs to genes...")
    all_mappings = []

    for idx, snp_row in snps_df.iterrows():
        print(f"  Processing SNP {idx+1}/{len(snps_df)}: {snp_row['Location']}", end='\r')
        mappings = map_snp_to_genes(snp_row)
        all_mappings.extend(mappings)

    print(f"\n  - Found {len(all_mappings)} total SNP-gene associations")

    # Create results dataframe
    results_df = pd.DataFrame(all_mappings)

    # For each SNP, keep only the nearest gene(s)
    # Filter out entries with no genes found (where distance is NaN)
    results_with_genes = results_df[results_df['distance'].notna()].copy()

    if len(results_with_genes) > 0:
        nearest_df = results_with_genes.loc[results_with_genes.groupby('snp_location')['distance'].idxmin()]
    else:
        nearest_df = pd.DataFrame()

    # Add back SNPs with no genes found
    no_genes = results_df[results_df['distance'].isna()].copy()
    if len(no_genes) > 0:
        nearest_df = pd.concat([nearest_df, no_genes], ignore_index=True)

    print(f"  - {len(nearest_df)} unique SNP-to-nearest-gene mappings\n")

    # Save full results
    print(f"Saving results to: {OUTPUT_FILE}")
    results_df.to_csv(OUTPUT_FILE, sep='\t', index=False)

    # Save nearest genes only
    nearest_file = OUTPUT_DIR / 'snp_to_nearest_gene_only.tsv'
    nearest_df.to_csv(nearest_file, sep='\t', index=False)
    print(f"Saving nearest genes only to: {nearest_file}\n")

    # Generate summary statistics
    print(f"{'='*70}")
    print("SUMMARY STATISTICS")
    print(f"{'='*70}\n")

    # 1. Overlap type distribution
    print("1. SNP Location Categories:")
    overlap_counts = nearest_df['overlap_type'].value_counts()
    for overlap_type, count in overlap_counts.items():
        pct = 100 * count / len(nearest_df)
        print(f"   - {overlap_type:20s}: {count:3d} ({pct:5.1f}%)")

    # 2. Distance distribution
    print("\n2. Distance Distribution (for SNPs outside genes):")
    outside_genes = nearest_df[nearest_df['overlap_type'] != 'inside_gene']
    if len(outside_genes) > 0:
        distances = outside_genes['distance']
        print(f"   - Median distance: {distances.median():,.0f} bp")
        print(f"   - Mean distance:   {distances.mean():,.0f} bp")
        print(f"   - Min distance:    {distances.min():,.0f} bp")
        print(f"   - Max distance:    {distances.max():,.0f} bp")

    # 3. Genes with multiple SNPs nearby
    print("\n3. Genes with Multiple High-Significance SNPs:")
    gene_counts = nearest_df['nearest_gene'].value_counts()
    multi_snp_genes = gene_counts[gene_counts > 1]
    if len(multi_snp_genes) > 0:
        print(f"   - {len(multi_snp_genes)} genes have >1 nearby SNP:")
        for gene, count in multi_snp_genes.items():
            print(f"     * {gene}: {count} SNPs")
    else:
        print("   - No genes have multiple nearby SNPs")

    # 4. Unique genes
    print(f"\n4. Unique Genes Identified:")
    unique_genes = nearest_df['nearest_gene'].nunique()
    print(f"   - {unique_genes} unique genes near high-significance SNPs")

    # Save unique gene list
    gene_list_file = OUTPUT_DIR / 'high_confidence_X_genes.txt'
    unique_gene_list = sorted(nearest_df['nearest_gene'].unique())
    with open(gene_list_file, 'w') as f:
        for gene in unique_gene_list:
            f.write(f"{gene}\n")
    print(f"   - Saved to: {gene_list_file}")

    # 5. Top signals
    print(f"\n5. Top 10 Highest-Significance SNPs and Their Nearest Genes:")
    top_snps = nearest_df.nlargest(10, 'neg_log10_p')
    for idx, row in top_snps.iterrows():
        print(f"   {row['snp_location']:25s} -> {row['nearest_gene']:15s} "
              f"(-log10 p = {row['neg_log10_p']:.2f}, dist = {row['distance']:>7,.0f} bp, {row['overlap_type']})")

    print(f"\n{'='*70}")
    print("Analysis complete!")
    print(f"{'='*70}\n")

    return results_df, nearest_df


if __name__ == '__main__':
    try:
        results_df, nearest_df = main()
    except Exception as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)
