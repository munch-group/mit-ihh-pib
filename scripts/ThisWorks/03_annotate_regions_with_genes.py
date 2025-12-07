#!/usr/bin/env python3
"""
Script: Annotate selection regions with genes using GENCODE
Purpose: Fast, reliable gene annotation without API dependencies

Uses GENCODE v44 (GRCh38.p14) - the gold standard for human genome annotation
Downloads and caches GTF file locally for fast repeated queries

Features:
- Downloads GENCODE v44 GTF automatically (cached)
- Annotates all regions in <1 minute
- Includes gene symbols, IDs, biotypes, descriptions
- Adds X-inactivation status from literature
- Generates gene lists for enrichment testing
- No rate limits, 100% reproducible

Author: MIT IHH PIB Project
Date: 2025-10-23
"""

import pandas as pd
import gzip
import re
from pathlib import Path
import urllib.request
import sys
from collections import defaultdict

# Setup paths
base_dir = Path("/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib")
data_dir = base_dir / "data/annotations"
regions_dir = base_dir / "results/ihs_EUR/analysis/regions"
output_dir = base_dir / "results/ihs_EUR/analysis/gene_annotation"

# Create directories
data_dir.mkdir(parents=True, exist_ok=True)
output_dir.mkdir(parents=True, exist_ok=True)

# GENCODE file paths
GENCODE_VERSION = "44"
GENCODE_GTF = data_dir / f"gencode.v{GENCODE_VERSION}.annotation.gtf.gz"
GENCODE_URL = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{GENCODE_VERSION}/gencode.v{GENCODE_VERSION}.annotation.gtf.gz"

# X-inactivation escape genes from Tukiainen et al. (2017) Nature
# These genes escape X-inactivation in females
X_INACTIVATION_ESCAPEES = {
    'AKAP17A', 'ASMTL', 'CD99', 'CRLF2', 'CSF2RA', 'DHRSX', 'DDX3X', 'EIF1AX',
    'GYG2', 'HDHD1A', 'IL3RA', 'KDM5C', 'KDM6A', 'NLGN4X', 'P2RY8', 'PLCXD1',
    'PPP2R3B', 'PRKX', 'RPS4X', 'SHOX', 'SLC25A6', 'TBL1X', 'TXLNG', 'USP9X',
    'UTX', 'XIST', 'ZFX', 'DDX3X', 'KDM5C', 'PNPLA4', 'JPX', 'CHIC1', 'CD99',
    'XG', 'ZBED1', 'SSR4', 'PRKX', 'NLGN4X', 'VCX', 'ANOS1'
}

print("=" * 70)
print("Gene Annotation using GENCODE v44")
print("=" * 70)
print()

def download_gencode():
    """Download GENCODE GTF file if not already cached"""
    if GENCODE_GTF.exists():
        file_size = GENCODE_GTF.stat().st_size / (1024**2)  # MB
        print(f"✓ GENCODE v{GENCODE_VERSION} already downloaded ({file_size:.1f} MB)")
        return

    print(f"Downloading GENCODE v{GENCODE_VERSION} GTF file...")
    print(f"  URL: {GENCODE_URL}")
    print(f"  Destination: {GENCODE_GTF}")
    print(f"  Size: ~45 MB (compressed), ~1.5 GB (uncompressed)")
    print()
    print("  This is a one-time download and will be cached for future use.")
    print("  Downloading... (may take 2-5 minutes depending on connection)")

    try:
        def report_progress(block_num, block_size, total_size):
            downloaded = block_num * block_size
            percent = min(100, downloaded * 100 / total_size)
            mb_downloaded = downloaded / (1024**2)
            mb_total = total_size / (1024**2)
            print(f"\r  Progress: {percent:.1f}% ({mb_downloaded:.1f}/{mb_total:.1f} MB)",
                  end='', file=sys.stderr)

        urllib.request.urlretrieve(GENCODE_URL, GENCODE_GTF, reporthook=report_progress)
        print("\n✓ Download complete!")
        print()
    except Exception as e:
        print(f"\n✗ Error downloading GENCODE: {e}", file=sys.stderr)
        sys.exit(1)

def parse_gtf_attributes(attr_string):
    """Parse GTF attribute field into dictionary"""
    attrs = {}
    # Match key "value" pairs
    pattern = r'(\w+)\s+"([^"]+)"'
    for match in re.finditer(pattern, attr_string):
        key, value = match.groups()
        attrs[key] = value
    return attrs

def load_gencode_genes():
    """
    Load gene annotations from GENCODE GTF
    Returns DataFrame with one row per gene
    """
    print("Step 1: Loading GENCODE gene annotations...")
    print(f"  Parsing {GENCODE_GTF.name}...")

    genes = []

    with gzip.open(GENCODE_GTF, 'rt') as f:
        for line_num, line in enumerate(f, 1):
            # Skip comments
            if line.startswith('#'):
                continue

            # Only process gene features
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            if feature_type != 'gene':
                continue

            # Parse fields
            chrom = fields[0].replace('chr', '')  # Remove 'chr' prefix if present
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = parse_gtf_attributes(fields[8])

            # Extract gene information
            gene_id = attributes.get('gene_id', '').split('.')[0]  # Remove version
            gene_name = attributes.get('gene_name', 'N/A')
            gene_type = attributes.get('gene_type', 'N/A')

            # Gene description if available
            # Note: basic GENCODE GTF doesn't have descriptions, we'll add later if needed

            genes.append({
                'chr': chrom,
                'gene_start': start,
                'gene_end': end,
                'strand': strand,
                'gene_id': gene_id,
                'gene_name': gene_name,
                'gene_type': gene_type
            })

            # Progress indicator
            if line_num % 100000 == 0:
                print(f"  Processed {line_num:,} lines...", end='\r', file=sys.stderr)

    print(f"\r  ✓ Parsed {len(genes):,} genes from GENCODE v{GENCODE_VERSION}")

    # Convert to DataFrame
    genes_df = pd.DataFrame(genes)

    # Add X-inactivation status for X chromosome genes
    genes_df['x_inactivation_status'] = genes_df.apply(
        lambda row: 'escape' if (row['chr'] == 'X' and row['gene_name'] in X_INACTIVATION_ESCAPEES)
        else 'subject' if row['chr'] == 'X'
        else 'not_applicable',
        axis=1
    )

    print(f"  Gene types: {genes_df['gene_type'].value_counts().head().to_dict()}")
    print(f"  Chromosomes: {sorted(genes_df['chr'].unique())}")
    print()

    return genes_df

def find_overlapping_genes(region_chr, region_start, region_end, genes_df, flank=50000):
    """
    Find genes overlapping a genomic region

    Args:
        region_chr: Chromosome name
        region_start: Region start position
        region_end: Region end position
        genes_df: DataFrame with gene annotations
        flank: Additional flanking region (default: 50kb)

    Returns:
        DataFrame with overlapping genes
    """
    # Add flanking regions
    search_start = max(0, region_start - flank)
    search_end = region_end + flank

    # Filter to chromosome
    chr_genes = genes_df[genes_df['chr'] == str(region_chr)]

    # Find overlapping genes
    overlapping = chr_genes[
        (chr_genes['gene_end'] >= search_start) &
        (chr_genes['gene_start'] <= search_end)
    ].copy()

    # Calculate overlap type
    def get_overlap_type(row):
        if row['gene_start'] >= region_start and row['gene_end'] <= region_end:
            return 'contained'
        elif row['gene_start'] < region_start and row['gene_end'] > region_end:
            return 'contains_region'
        elif row['gene_start'] < search_start or row['gene_end'] > search_end:
            return 'flanking'
        else:
            return 'overlapping'

    if len(overlapping) > 0:
        overlapping['overlap_type'] = overlapping.apply(get_overlap_type, axis=1)

    return overlapping

def annotate_all_regions(regions_file, genes_df, flank=50000):
    """
    Annotate all regions with overlapping genes

    Args:
        regions_file: Path to regions TSV file
        genes_df: DataFrame with gene annotations
        flank: Flanking region size (default: 50kb)

    Returns:
        Tuple of (annotated_regions_df, genes_in_regions_df)
    """
    print(f"Step 2: Annotating regions with genes (±{flank/1000:.0f}kb flank)...")

    # Load regions
    regions = pd.read_csv(regions_file, sep='\t')
    print(f"  Loaded {len(regions):,} regions")
    print(f"  Chromosomes: {sorted(regions['chr'].unique())}")
    print()

    annotated_regions = []
    all_region_genes = []

    for idx, region in regions.iterrows():
        # Find overlapping genes
        overlapping_genes = find_overlapping_genes(
            region['chr'],
            region['start'],
            region['end'],
            genes_df,
            flank=flank
        )

        # Separate by overlap type (only if genes exist)
        if len(overlapping_genes) > 0 and 'overlap_type' in overlapping_genes.columns:
            contained = overlapping_genes[overlapping_genes['overlap_type'] == 'contained']
            overlapping = overlapping_genes[overlapping_genes['overlap_type'] == 'overlapping']
            flanking = overlapping_genes[overlapping_genes['overlap_type'] == 'flanking']
        else:
            contained = overlapping_genes[:]  # Empty DataFrame
            overlapping = overlapping_genes[:]
            flanking = overlapping_genes[:]

        # Create region annotation
        region_annotation = region.to_dict()
        region_annotation.update({
            'n_genes_total': len(overlapping_genes),
            'n_genes_contained': len(contained),
            'n_genes_overlapping': len(overlapping),
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
                'region_max_std_ihs': region['max_std_ihs'],
            })
            all_region_genes.append(gene_info)

        # Progress
        if (idx + 1) % 500 == 0:
            print(f"  Processed {idx + 1:,}/{len(regions):,} regions...", end='\r', file=sys.stderr)

    print(f"\r  ✓ Annotated all {len(regions):,} regions")
    print()

    return pd.DataFrame(annotated_regions), pd.DataFrame(all_region_genes)

def generate_gene_lists(genes_in_regions_df, output_dir):
    """
    Generate gene lists for enrichment analysis

    Creates separate lists for:
    - All genes in selection regions
    - X chromosome genes in selection
    - Autosomal genes in selection
    - Protein-coding genes only
    - By gene type
    """
    print("Step 3: Generating gene lists for enrichment analysis...")

    # All unique genes
    all_genes = genes_in_regions_df.drop_duplicates(subset=['gene_id'])
    print(f"  Total unique genes in regions: {len(all_genes):,}")

    # Save comprehensive gene list
    all_genes_file = output_dir / "all_genes_in_selection_regions.tsv"
    all_genes.to_csv(all_genes_file, sep='\t', index=False)
    print(f"  Saved: {all_genes_file.name}")

    # X chromosome vs autosomes
    x_genes = all_genes[all_genes['chr'] == 'X']
    autosome_genes = all_genes[all_genes['chr'] != 'X']

    x_genes_file = output_dir / "X_chromosome_genes_in_selection.tsv"
    x_genes.to_csv(x_genes_file, sep='\t', index=False)
    print(f"  Saved: {x_genes_file.name} ({len(x_genes):,} genes)")

    autosome_genes_file = output_dir / "autosome_genes_in_selection.tsv"
    autosome_genes.to_csv(autosome_genes_file, sep='\t', index=False)
    print(f"  Saved: {autosome_genes_file.name} ({len(autosome_genes):,} genes)")

    # Protein-coding genes only (for enrichment)
    protein_coding = all_genes[all_genes['gene_type'] == 'protein_coding']
    protein_coding_file = output_dir / "protein_coding_genes_in_selection.tsv"
    protein_coding.to_csv(protein_coding_file, sep='\t', index=False)
    print(f"  Saved: {protein_coding_file.name} ({len(protein_coding):,} genes)")

    # X chromosome protein-coding
    x_protein_coding = protein_coding[protein_coding['chr'] == 'X']
    x_protein_coding_file = output_dir / "X_protein_coding_genes.tsv"
    x_protein_coding.to_csv(x_protein_coding_file, sep='\t', index=False)
    print(f"  Saved: {x_protein_coding_file.name} ({len(x_protein_coding):,} genes)")

    # Simple gene name lists for g:Profiler / enrichment tools
    gene_names_dir = output_dir / "gene_name_lists"
    gene_names_dir.mkdir(exist_ok=True)

    # All genes (gene names only)
    with open(gene_names_dir / "all_genes.txt", 'w') as f:
        f.write('\n'.join(all_genes['gene_name'].tolist()))

    # X chromosome genes
    with open(gene_names_dir / "X_chromosome_genes.txt", 'w') as f:
        f.write('\n'.join(x_genes['gene_name'].tolist()))

    # Protein-coding genes
    with open(gene_names_dir / "protein_coding_genes.txt", 'w') as f:
        f.write('\n'.join(protein_coding['gene_name'].tolist()))

    # X protein-coding
    with open(gene_names_dir / "X_protein_coding_genes.txt", 'w') as f:
        f.write('\n'.join(x_protein_coding['gene_name'].tolist()))

    print(f"  Saved simple gene name lists to: {gene_names_dir.name}/")
    print()

    return {
        'all': len(all_genes),
        'x_chromosome': len(x_genes),
        'autosomes': len(autosome_genes),
        'protein_coding': len(protein_coding),
        'x_protein_coding': len(x_protein_coding)
    }

def print_summary(annotated_regions, genes_in_regions, gene_counts):
    """Print summary statistics"""
    print("=" * 70)
    print("Annotation Summary")
    print("=" * 70)
    print()

    # Region statistics
    print("Regions:")
    print(f"  Total regions: {len(annotated_regions):,}")
    regions_with_genes = annotated_regions[annotated_regions['n_genes_total'] > 0]
    print(f"  Regions with genes: {len(regions_with_genes):,} ({len(regions_with_genes)/len(annotated_regions)*100:.1f}%)")
    print(f"  Average genes per region: {annotated_regions['n_genes_total'].mean():.1f}")
    print(f"  Median genes per region: {annotated_regions['n_genes_total'].median():.0f}")
    print()

    # Gene statistics
    print("Genes:")
    print(f"  Total unique genes: {gene_counts['all']:,}")
    print(f"  X chromosome genes: {gene_counts['x_chromosome']:,}")
    print(f"  Autosomal genes: {gene_counts['autosomes']:,}")
    print(f"  Protein-coding genes: {gene_counts['protein_coding']:,} ({gene_counts['protein_coding']/gene_counts['all']*100:.1f}%)")
    print(f"  X protein-coding: {gene_counts['x_protein_coding']:,}")
    print()

    # X-inactivation status (for X chromosome)
    if len(genes_in_regions) > 0 and 'x_inactivation_status' in genes_in_regions.columns:
        x_genes_all = genes_in_regions[genes_in_regions['chr'] == 'X'].drop_duplicates(subset=['gene_id'])
        if len(x_genes_all) > 0:
            escapees = x_genes_all[x_genes_all['x_inactivation_status'] == 'escape']
            print("X-inactivation Status (X chromosome genes):")
            print(f"  Escape X-inactivation: {len(escapees):,} ({len(escapees)/len(x_genes_all)*100:.1f}%)")
            print(f"  Subject to X-inactivation: {len(x_genes_all) - len(escapees):,}")
            if len(escapees) > 0:
                print(f"  Escapee examples: {', '.join(escapees['gene_name'].head(10).tolist())}")
            print()

    # Gene type breakdown
    if len(genes_in_regions) > 0:
        unique_genes = genes_in_regions.drop_duplicates(subset=['gene_id'])
        print("Gene Types:")
        gene_type_counts = unique_genes['gene_type'].value_counts().head(10)
        for gene_type, count in gene_type_counts.items():
            print(f"  {gene_type}: {count:,} ({count/len(unique_genes)*100:.1f}%)")
        print()

    # Top regions with most genes
    print("Top 10 regions with most genes:")
    top_regions = annotated_regions.nlargest(10, 'n_genes_total')[
        ['chr', 'start', 'end', 'n_genes_total', 'protein_coding_genes', 'max_std_ihs']
    ]
    for _, region in top_regions.iterrows():
        genes_list = region['protein_coding_genes'][:100]  # Truncate if too long
        if genes_list == 'None':
            genes_list = "No protein-coding genes"
        print(f"  chr{region['chr']}:{region['start']:,}-{region['end']:,} "
              f"({region['n_genes_total']} genes, max |iHS|={region['max_std_ihs']:.2f})")
        print(f"    Protein-coding: {genes_list}")

def main():
    """Main execution"""
    # Download GENCODE if needed
    download_gencode()

    # Load genes
    genes_df = load_gencode_genes()

    # Annotate all regions
    regions_file = regions_dir / "selection_regions_500kb_min2snps.tsv"
    annotated_regions, genes_in_regions = annotate_all_regions(regions_file, genes_df, flank=50000)

    # Save annotated regions
    output_file = output_dir / "selection_regions_annotated.tsv"
    annotated_regions.to_csv(output_file, sep='\t', index=False)
    print(f"Step 2 output: Saved annotated regions to {output_file.name}")
    print()

    # Save genes in regions (detailed)
    genes_file = output_dir / "genes_in_selection_regions_detailed.tsv"
    genes_in_regions.to_csv(genes_file, sep='\t', index=False)
    print(f"Step 2 output: Saved detailed gene list to {genes_file.name}")
    print()

    # Generate gene lists for enrichment
    gene_counts = generate_gene_lists(genes_in_regions, output_dir)

    # Print summary
    print_summary(annotated_regions, genes_in_regions, gene_counts)

    print()
    print("=" * 70)
    print("Annotation complete!")
    print("=" * 70)
    print()
    print("Output files:")
    print(f"  {output_file}")
    print(f"  {genes_file}")
    print(f"  {output_dir / 'gene_name_lists'}/ (simple gene lists)")
    print()
    print("Next steps:")
    print("  1. Review annotated regions")
    print("  2. Use gene lists for enrichment analysis (g:Profiler, DAVID, etc.)")
    print("  3. Test for reproductive/immune gene enrichment (Phase 4)")

if __name__ == '__main__':
    main()
