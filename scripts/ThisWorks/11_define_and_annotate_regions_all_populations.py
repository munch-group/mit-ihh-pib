#!/usr/bin/env python3
"""
Automated Selection Region Definition and Gene Annotation for All Populations

Purpose:
    Run the complete region definition and gene annotation pipeline for all 1000 Genomes populations.
    This script automates what was previously done manually for individual populations.

Key features:
    - Automatically calculates population-specific q99 threshold from raw iHS data
    - Defines selection regions using sliding window approach (20kb windows, 4kb steps)
    - Annotates regions with genes from GENCODE v44
    - Includes X-inactivation status 
    - Processes X chromosome only (as per original pipeline)

Based on:
    - 02_define_selection_regions.py (region definition logic)
    - 03_annotate_regions_with_genes.py (gene annotation logic)

Author: Mit Van Bruggen 
Date: 2025-12-18
"""

import pandas as pd
import numpy as np
import gzip
import re
import urllib.request
import sys
from pathlib import Path
from collections import defaultdict
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Superpopulation structure
SUPERPOPULATIONS = {
    'AFR': ['ACB', 'ASW', 'ESN', 'GWD', 'LWK', 'MSL', 'YRI'],
    'AMR': ['CLM', 'MXL', 'PEL', 'PUR'],
    'EAS': ['CDX', 'CHB', 'CHS', 'JPT', 'KHV'],
    'EUR': ['CEU', 'FIN', 'GBR', 'IBS', 'TSI'],
    'SAS': ['BEB', 'GIH', 'ITU', 'PJL', 'STU']
}

# Get all populations
ALL_POPULATIONS = [pop for pops in SUPERPOPULATIONS.values() for pop in pops]

# Paths
base_dir = Path("/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib")
data_dir = base_dir / "data/annotations"
data_dir.mkdir(parents=True, exist_ok=True)

# GENCODE configuration
GENCODE_VERSION = "44"
gencode_gtf = data_dir / f"gencode.v{GENCODE_VERSION}.annotation.gtf.gz"
GENCODE_URL = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{GENCODE_VERSION}/gencode.v{GENCODE_VERSION}.annotation.gtf.gz"

# X-inactivation escapee genes from Tukiainen et al. (2017) Nature
X_INACTIVATION_ESCAPEES = {
    'AKAP17A', 'ASMTL', 'CD99', 'CRLF2', 'CSF2RA', 'DHRSX', 'DDX3X', 'EIF1AX',
    'GYG2', 'HDHD1A', 'IL3RA', 'KDM5C', 'KDM6A', 'NLGN4X', 'P2RY8', 'PLCXD1',
    'PPP2R3B', 'PRKX', 'RPS4X', 'SHOX', 'SLC25A6', 'TBL1X', 'TXLNG', 'USP9X',
    'UTX', 'XIST', 'ZFX', 'DDX3X', 'KDM5C', 'PNPLA4', 'JPX', 'CHIC1', 'CD99',
    'XG', 'ZBED1', 'SSR4', 'PRKX', 'NLGN4X', 'VCX', 'ANOS1'
}

# Region definition configuration (from Script 02)
WINDOW_SIZE = 20000      # 20kb windows
STEP_SIZE = 4000         # 4kb steps (80% overlap)
MIN_SNPS = 20            # Minimum SNPs per region
PERCENTILE = 99          # 99th percentile for threshold

#nice output headers 
print("=" * 80)
print("Automated Selection Region Definition and Gene Annotation")
print("=" * 80)
print()
print(f"Configuration:")
print(f"  Window size: {WINDOW_SIZE:,} bp ({WINDOW_SIZE/1000:.0f} kb)")
print(f"  Step size: {STEP_SIZE:,} bp ({STEP_SIZE/1000:.0f} kb)")
print(f"  Window overlap: {100 * (1 - STEP_SIZE/WINDOW_SIZE):.0f}%")
print(f"  Minimum SNPs per region: {MIN_SNPS}")
print(f"  Threshold: {PERCENTILE}th percentile of |Std iHS| per population")
print(f"  Chromosome: X only")
print()
print(f"Populations to process: {len(ALL_POPULATIONS)}")
print()


# ============================================================================
# GENCODE Gene Loading Functions
# ============================================================================
# added 'if not already cached' for gencode since I re-ran this script multiple times
def download_gencode():
    """Download GENCODE GTF file if not already cached"""
    if gencode_gtf.exists():
        file_size = gencode_gtf.stat().st_size / (1024**2)  # MB
        print(f"GENCODE v{GENCODE_VERSION} already cached ({file_size:.1f} MB)")
        return

    try:
        def report_progress(block_num, block_size, total_size):
            downloaded = block_num * block_size
            percent = min(100, downloaded * 100 / total_size)
            mb_downloaded = downloaded / (1024**2)
            mb_total = total_size / (1024**2)
            print(f"\r  Progress: {percent:.1f}% ({mb_downloaded:.1f}/{mb_total:.1f} MB)",
                  end='', file=sys.stderr) #progress bar for downlad 

        urllib.request.urlretrieve(GENCODE_URL, gencode_gtf, reporthook=report_progress)
        print("\nDownload complete!")
        print()
    except Exception as e: #claude addition
        print(f"\n Error downloading GENCODE: {e}", file=sys.stderr)
        print(f"  Please download manually from: {GENCODE_URL}")
        raise 


def load_gencode_genes():
    """Load gene annotations from GENCODE GTF (cached)"""
    download_gencode()

    genes = []
    with gzip.open(gencode_gtf, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom = fields[0]
            feature_type = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8]

            # Only process genes on X chromosome
            if chrom != 'chrX' or feature_type != 'gene':
                continue

            # Parse attributes
            attr_dict = {}
            pattern = r'(\w+)\s+"([^"]+)"'
            for match in re.finditer(pattern, attributes):
                key, value = match.groups()
                attr_dict[key] = value

            gene_name = attr_dict.get('gene_name', '')
            gene_id = attr_dict.get('gene_id', '').split('.')[0]
            gene_type = attr_dict.get('gene_type', '')

            # X-inactivation status
            x_inactivation_status = 'escapee' if gene_name in X_INACTIVATION_ESCAPEES else 'subject'

            genes.append({
                'gene_name': gene_name,
                'gene_id': gene_id,
                'gene_type': gene_type,
                'chr': 'X',  # Standardize to 'X' (without 'chr' prefix)
                'start': start,
                'end': end,
                'strand': strand,
                'x_inactivation_status': x_inactivation_status
            })

    genes_df = pd.DataFrame(genes)
    print(f"  Loaded {len(genes_df):,} genes on X chromosome")
    print()

    return genes_df


# ============================================================================
# Region Definition Functions (from Script 02)
# ============================================================================
# get raw ihs data 
def load_ihs_data(pop_code: str):
    """Load raw iHS data for a population"""
    ihs_file = base_dir / f"results/ihs_{pop_code}/{pop_code}.chrX.ihs.tsv"

    if not ihs_file.exists(): #added this for reproducibility 
        print(f"  Warning: iHS file not found for {pop_code}")
        return None

    df = pd.read_csv(ihs_file, sep='\t')

    # Parse position from Location column
    df['pos'] = df['Location'].str.split(':').str[1].astype(int)
    df['chr'] = 'X'

    return df


def calculate_threshold(ihs_df: pd.DataFrame, percentile: float = 99) -> float:
    """Calculate percentile threshold from |Std iHS|"""
    # Remove NaN values (I had a couple in one population so just to be safe)
    clean_data = ihs_df[~ihs_df['Std iHS'].isna()].copy()

    if len(clean_data) == 0:
        return np.nan

    # Calculate absolute standardized iHS
    abs_std_ihs = np.abs(clean_data['Std iHS'])

    # Calculate threshold
    threshold = np.percentile(abs_std_ihs, percentile)

    return threshold


def cluster_candidates_to_regions(candidates: pd.DataFrame, window_size: int,
                                 step_size: int, min_snps: int):
    """
    Cluster candidates into selection regions using sliding windows
    (Villegas-Mirón et al. 2021 approach)
    """
    regions = []

    # Sort by position
    candidates = candidates.sort_values('pos').reset_index(drop=True)
    positions = candidates['pos'].values

    chr_start = int(positions.min())
    chr_end = int(positions.max())

    # Sliding window scan with fixed step size
    outlier_windows = []
    for win_start in range(chr_start, chr_end, step_size):
        win_end = win_start + window_size

        # Find all candidates in this window
        in_window = (positions >= win_start) & (positions < win_end)
        window_indices = np.where(in_window)[0]

        # Keep windows with sufficient candidates
        if len(window_indices) >= min_snps:
            outlier_windows.append({
                'start': win_start,
                'end': win_end,
                'indices': window_indices.tolist()
            })

    # Merge overlapping/consecutive windows into regions
    if len(outlier_windows) == 0:
        return pd.DataFrame()

    clusters = []
    current_cluster = {
        'start': outlier_windows[0]['start'],
        'end': outlier_windows[0]['end'],
        'indices': set(outlier_windows[0]['indices'])
    }

    for window in outlier_windows[1:]:
        # Check if window overlaps with current cluster
        if window['start'] <= current_cluster['end']:
            # Extend and merge
            current_cluster['end'] = max(current_cluster['end'], window['end'])
            current_cluster['indices'].update(window['indices'])
        else:
            # Save current cluster and start new one
            clusters.append(current_cluster)
            current_cluster = {
                'start': window['start'],
                'end': window['end'],
                'indices': set(window['indices'])
            }

    # Don't forget last cluster
    clusters.append(current_cluster)

    # Calculate statistics for each cluster (optional but I did it so I coule report on the max iHS scores in each defined region)
    for cluster in clusters:
        # Get candidates in this cluster
        indices_list = sorted(list(cluster['indices']))
        cluster_candidates = candidates.iloc[indices_list]

        # Use actual candidate positions for region boundaries
        cluster_start = cluster_candidates['pos'].min()
        cluster_end = cluster_candidates['pos'].max()
        cluster_length = cluster_end - cluster_start

        # Calculate statistics
        abs_std_ihs = np.abs(cluster_candidates['Std iHS'])
        mean_std_ihs = abs_std_ihs.mean()
        max_std_ihs = abs_std_ihs.max()

        # P-values (two-tailed) (didn't use this in the end)
        mean_pvalue = 2 * (1 - stats.norm.cdf(mean_std_ihs))
        min_pvalue = 2 * (1 - stats.norm.cdf(max_std_ihs))

        # Top SNP
        top_snp_idx = abs_std_ihs.idxmax()
        top_snp = cluster_candidates.loc[top_snp_idx]

        regions.append({
            'chr': 'X',
            'start': cluster_start,
            'end': cluster_end,
            'length': cluster_length,
            'n_snps': len(cluster_candidates),
            'mean_std_ihs': mean_std_ihs,
            'max_std_ihs': max_std_ihs,
            'mean_pvalue': mean_pvalue,
            'min_pvalue': min_pvalue,
            'top_snp_location': top_snp['Location'],
            'top_snp_std_ihs': abs_std_ihs.loc[top_snp_idx],
            'chr_type': 'X chromosome'
        })

    regions_df = pd.DataFrame(regions)

    return regions_df


def define_regions_for_population(pop_code: str):
    """
    Define selection regions for a single population
    Returns: (regions_df, threshold_used)
    """
    print(f"\n{pop_code}:")
    print("-" * 40)

    # Load iHS data
    print(f"  Loading iHS data...")
    ihs_df = load_ihs_data(pop_code)

    if ihs_df is None:
        print(f"  Skipping {pop_code}: No iHS data")
        return None, None

    print(f"    {len(ihs_df):,} total variants")

    # Calculate threshold
    threshold = calculate_threshold(ihs_df, PERCENTILE)
    print(f"    {PERCENTILE}th percentile |Std iHS|: {threshold:.4f}")

    # Filter candidates
    ihs_df['abs_std_ihs'] = np.abs(ihs_df['Std iHS'])
    candidates = ihs_df[ihs_df['abs_std_ihs'] >= threshold].copy()
    print(f"    {len(candidates):,} candidates above threshold")

    if len(candidates) == 0:
        print(f"  Skipping {pop_code}: No candidates above threshold")
        return None, threshold

    # Define regions
    print(f"  Clustering candidates into regions...")
    regions = cluster_candidates_to_regions(candidates, WINDOW_SIZE, STEP_SIZE, MIN_SNPS)

    if len(regions) == 0:
        print(f"    No regions with >= {MIN_SNPS} SNPs")
        return None, threshold

    print(f"    {len(regions)} regions defined")

    return regions, threshold


# ============================================================================
# Gene Annotation Functions (from Script 03)
# ============================================================================

def find_overlapping_genes(region_chr, region_start, region_end, genes_df, flank=50000):
    """Find genes overlapping a region (with optional flanking region)"""
    # Add flanking regions
    search_start = max(0, region_start - flank)
    search_end = region_end + flank

    # Filter to chromosome
    chr_genes = genes_df[genes_df['chr'] == str(region_chr)]

    # Find overlapping genes
    overlapping = chr_genes[
        (chr_genes['end'] >= search_start) &
        (chr_genes['start'] <= search_end)
    ].copy()

    # Calculate overlap type
    def get_overlap_type(row):
        if row['start'] >= region_start and row['end'] <= region_end:
            return 'contained'
        elif row['start'] < region_start and row['end'] > region_end:
            return 'contains_region'
        elif row['start'] < search_start or row['end'] > search_end:
            return 'flanking'
        else:
            return 'overlapping'

    if len(overlapping) > 0:
        overlapping['overlap_type'] = overlapping.apply(get_overlap_type, axis=1)

    return overlapping


def annotate_regions_with_genes(regions_df, genes_df, flank=50000):
    """
    Annotate all regions with overlapping genes
    Returns: (annotated_regions_df, genes_in_regions_df)
    """
    annotated_regions = []
    all_region_genes = []

    for idx, region in regions_df.iterrows():
        # Find overlapping genes
        overlapping_genes = find_overlapping_genes(
            region['chr'],
            region['start'],
            region['end'],
            genes_df,
            flank=flank
        )

        # Separate by overlap type
        if len(overlapping_genes) > 0 and 'overlap_type' in overlapping_genes.columns:
            contained = overlapping_genes[overlapping_genes['overlap_type'] == 'contained']
            overlapping_only = overlapping_genes[overlapping_genes['overlap_type'] == 'overlapping']
            flanking = overlapping_genes[overlapping_genes['overlap_type'] == 'flanking']
        else:
            contained = overlapping_genes[:]
            overlapping_only = overlapping_genes[:]
            flanking = overlapping_genes[:]

        # Create region annotation
        region_annotation = region.to_dict()
        region_annotation.update({
            'n_genes_total': len(overlapping_genes),
            'n_genes_contained': len(contained),
            'n_genes_overlapping': len(overlapping_only),
            'n_genes_flanking': len(flanking),
            'gene_names': ','.join(overlapping_genes['gene_name'].tolist()) if len(overlapping_genes) > 0 else 'None',
            'gene_ids': ','.join(overlapping_genes['gene_id'].tolist()) if len(overlapping_genes) > 0 else 'None',
            'protein_coding_genes': ','.join(
                overlapping_genes[overlapping_genes['gene_type'] == 'protein_coding']['gene_name'].tolist()
            ) if len(overlapping_genes) > 0 else 'None'
        })
        annotated_regions.append(region_annotation)

        # Add genes to master list with region info
        for _, gene in overlapping_genes.iterrows():
            gene_info = gene.to_dict()
            gene_info.update({
                'region_chr': region['chr'],
                'region_start': region['start'],
                'region_end': region['end'],
                'region_length': region['length'],
                'region_n_snps': region['n_snps'],
                'region_max_std_ihs': region['max_std_ihs']
            })
            all_region_genes.append(gene_info)

    annotated_regions_df = pd.DataFrame(annotated_regions)
    genes_in_regions_df = pd.DataFrame(all_region_genes) if all_region_genes else pd.DataFrame()

    return annotated_regions_df, genes_in_regions_df


# ============================================================================
# Main Processing Function (cleaned up by claude)
# ============================================================================

def process_population(pop_code: str, genes_df: pd.DataFrame):
    """
    Complete pipeline for one population: define regions + annotate with genes
    """
    # Create output directory
    output_dir = base_dir / f"results/ihs_{pop_code}/analysis/regions"
    output_dir.mkdir(parents=True, exist_ok=True)

    gene_output_dir = base_dir / f"results/ihs_{pop_code}/analysis/gene_annotation"
    gene_output_dir.mkdir(parents=True, exist_ok=True)

    # Define regions
    regions, threshold = define_regions_for_population(pop_code)

    if regions is None or len(regions) == 0:
        print(f"  No regions defined for {pop_code}")
        return {
            'population': pop_code,
            'threshold': threshold if threshold else np.nan,
            'n_regions': 0,
            'n_genes': 0,
            'status': 'no_regions'
        }

    # Save regions
    regions_file = output_dir / f"selection_regions_X_chromosome.tsv"
    regions.to_csv(regions_file, sep='\t', index=False)
    print(f"  Saved: {regions_file}")

    # Annotate with genes
    print(f"  Annotating regions with genes...")
    annotated_regions, genes_in_regions = annotate_regions_with_genes(regions, genes_df)

    # Save annotated regions
    annotated_file = gene_output_dir / f"annotated_regions.tsv"
    annotated_regions.to_csv(annotated_file, sep='\t', index=False)
    print(f"  Saved: {annotated_file}")

    # Save genes in regions
    if len(genes_in_regions) > 0:
        # Save all genes
        genes_file = gene_output_dir / f"genes_in_regions.tsv"
        genes_in_regions.to_csv(genes_file, sep='\t', index=False)
        print(f"  Saved: {genes_file}")

        # Filter and save protein-coding genes only
        protein_coding = genes_in_regions[genes_in_regions['gene_type'] == 'protein_coding']
        if len(protein_coding) > 0:
            protein_coding_file = gene_output_dir / f"genes_in_regions_protein_coding.tsv"
            protein_coding.to_csv(protein_coding_file, sep='\t', index=False)
            print(f"  Saved: {protein_coding_file}")
            print(f"    {len(protein_coding)} protein-coding genes ({protein_coding['gene_name'].nunique()} unique)")

        # Gene summary
        unique_genes = genes_in_regions['gene_name'].nunique()
        unique_protein_coding = protein_coding['gene_name'].nunique() if len(protein_coding) > 0 else 0
        print(f"    {unique_genes} unique genes total ({unique_protein_coding} protein-coding)")

        return {
            'population': pop_code,
            'threshold': threshold,
            'n_regions': len(regions),
            'n_genes': unique_genes,
            'n_protein_coding_genes': unique_protein_coding,
            'status': 'success'
        }
    else:
        print(f"    No genes found in regions")
        return {
            'population': pop_code,
            'threshold': threshold,
            'n_regions': len(regions),
            'n_genes': 0,
            'n_protein_coding_genes': 0,
            'status': 'success_no_genes'
        }


# ============================================================================
# Main Execution (added print statements to show progress / debug)
# ============================================================================

def main():
    # Load genes once (used for all populations)
    print("\n" + "=" * 80)
    print("Step 1: Loading GENCODE gene annotations")
    print("=" * 80 + "\n")

    genes_df = load_gencode_genes()

    # Process all populations
    print("\n" + "=" * 80)
    print("Step 2: Processing all populations")
    print("=" * 80)

    summary_stats = []

    for pop_code in ALL_POPULATIONS:
        result = process_population(pop_code, genes_df)
        summary_stats.append(result)

    # Save summary
    print("\n" + "=" * 80)
    print("Step 3: Saving summary")
    print("=" * 80 + "\n")

    summary_df = pd.DataFrame(summary_stats)
    summary_file = base_dir / "results/all_populations_region_summary.tsv"
    summary_df.to_csv(summary_file, sep='\t', index=False)
    print(f"Saved summary: {summary_file}")

    # Print summary statistics
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80 + "\n")

    print(summary_df.to_string(index=False))
    print()

    successful = summary_df[summary_df['status'].str.contains('success')]
    print(f"Successfully processed: {len(successful)}/{len(ALL_POPULATIONS)} populations")
    print(f"Total regions defined: {successful['n_regions'].sum():,}")
    print(f"Total unique genes across all populations: {summary_df['n_genes'].sum():,}")
    print()


if __name__ == "__main__":
    main()
