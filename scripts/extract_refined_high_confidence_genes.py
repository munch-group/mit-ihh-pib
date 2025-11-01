#!/usr/bin/env python3
"""
Extract refined high-confidence gene sets from SNP-to-gene mappings.

This script creates multiple gene sets based on different criteria:
1. Stringent: genes overlapping high-sig SNPs (neg_log10_p >= 6.0)
2. Moderate: genes within 50kb of high-sig SNPs (neg_log10_p >= 6.0)
3. Liberal: genes with multiple high-sig SNPs or very high signals

This provides better statistical power for GO enrichment while maintaining
biological relevance.
"""

import pandas as pd
from pathlib import Path
from datetime import datetime

# Define paths
BASE_DIR = Path('/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib')
HIGH_SIG_DIR = BASE_DIR / 'results/analysis/high_significance_genes'
OUTPUT_DIR = HIGH_SIG_DIR
INPUT_FILE = HIGH_SIG_DIR / 'snp_to_gene_mapping.tsv'

print("=" * 80)
print("Extracting Refined High-Confidence Gene Sets")
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 80)
print()

# ============================================================================
# Load SNP-to-gene mappings
# ============================================================================

print("Step 1: Loading SNP-to-gene mappings...")
df = pd.read_csv(INPUT_FILE, sep='\t')
print(f"  Loaded {len(df)} SNP-to-gene mappings")
print(f"  Unique genes: {df['nearest_gene'].nunique()}")
print(f"  Unique SNPs: {df['snp_location'].nunique()}")
print()

# ============================================================================
# Analyze signal distribution
# ============================================================================

print("Step 2: Analyzing signal distribution...")
print(f"  neg_log10_p range: {df['neg_log10_p'].min():.2f} - {df['neg_log10_p'].max():.2f}")
print(f"  std_ihs range: {df['std_ihs'].min():.2f} - {df['std_ihs'].max():.2f}")
print()

print("  Distribution of overlap types:")
print(df['overlap_type'].value_counts())
print()

# ============================================================================
# Gene Set 1: STRINGENT - Overlapping genes with strong signal
# ============================================================================

print("Step 3: Extracting STRINGENT gene set...")
print("  Criteria:")
print("    - overlap_type = 'inside_gene'")
print("    - neg_log10_p >= 6.0")
print()

stringent = df[
    (df['overlap_type'] == 'inside_gene') &
    (df['neg_log10_p'] >= 6.0)
].copy()

stringent_genes = sorted(stringent['nearest_gene'].unique())
print(f"  Found {len(stringent_genes)} genes in stringent set")

# Save gene list
stringent_file = OUTPUT_DIR / 'high_confidence_genes_stringent.txt'
with open(stringent_file, 'w') as f:
    for gene in stringent_genes:
        f.write(f"{gene}\n")
print(f"  Saved: {stringent_file.name}")
print()

# ============================================================================
# Gene Set 2: MODERATE - Nearby genes with strong signal
# ============================================================================

print("Step 4: Extracting MODERATE gene set...")
print("  Criteria:")
print("    - overlap_type = 'inside_gene' OR distance < 50kb")
print("    - neg_log10_p >= 6.0")
print()

moderate = df[
    (
        (df['overlap_type'] == 'inside_gene') |
        (df['distance'].abs() < 50000)
    ) &
    (df['neg_log10_p'] >= 6.0)
].copy()

moderate_genes = sorted(moderate['nearest_gene'].unique())
print(f"  Found {len(moderate_genes)} genes in moderate set")

# Save gene list
moderate_file = OUTPUT_DIR / 'high_confidence_genes_moderate.txt'
with open(moderate_file, 'w') as f:
    for gene in moderate_genes:
        f.write(f"{gene}\n")
print(f"  Saved: {moderate_file.name}")
print()

# ============================================================================
# Gene Set 3: LIBERAL - Multiple SNPs or very high signal
# ============================================================================

print("Step 5: Extracting LIBERAL gene set...")
print("  Criteria:")
print("    - (overlap_type = 'inside_gene' OR distance < 100kb) AND")
print("    - (neg_log10_p >= 6.0 OR multiple high-sig SNPs)")
print()

# Count SNPs per gene
snps_per_gene = df[df['neg_log10_p'] >= 6.0].groupby('nearest_gene').size()
genes_with_multiple_snps = set(snps_per_gene[snps_per_gene >= 2].index)

liberal = df[
    (
        (df['overlap_type'] == 'inside_gene') |
        (df['distance'].abs() < 100000)
    ) &
    (
        (df['neg_log10_p'] >= 6.0) |
        (df['nearest_gene'].isin(genes_with_multiple_snps))
    )
].copy()

liberal_genes = sorted(liberal['nearest_gene'].unique())
print(f"  Found {len(liberal_genes)} genes in liberal set")
print(f"    Including {len(genes_with_multiple_snps)} genes with multiple high-sig SNPs")

# Save gene list
liberal_file = OUTPUT_DIR / 'high_confidence_genes_liberal.txt'
with open(liberal_file, 'w') as f:
    for gene in liberal_genes:
        f.write(f"{gene}\n")
print(f"  Saved: {liberal_file.name}")
print()

# ============================================================================
# Gene Set 4: PROTEIN-CODING ONLY (moderate criteria)
# ============================================================================

print("Step 6: Extracting PROTEIN-CODING subset (moderate criteria)...")

# We need to filter out LOC genes and non-coding genes
# Heuristic: exclude genes starting with LOC, MIR, LINC, or ending with -DT, -AS1, etc.
def is_likely_protein_coding(gene_name):
    """Simple heuristic to identify likely protein-coding genes."""
    if pd.isna(gene_name):
        return False

    # Exclude non-coding gene patterns
    exclude_patterns = ['LOC', 'MIR', 'LINC', 'SNORD', 'SNORA']
    exclude_suffixes = ['-DT', '-AS1', '-AS2', 'RNA', 'RNA-']

    gene_upper = str(gene_name).upper()

    # Check if starts with excluded pattern
    for pattern in exclude_patterns:
        if gene_upper.startswith(pattern):
            return False

    # Check if ends with excluded suffix
    for suffix in exclude_suffixes:
        if gene_upper.endswith(suffix):
            return False

    return True

moderate_protein_coding = [g for g in moderate_genes if is_likely_protein_coding(g)]
print(f"  Found {len(moderate_protein_coding)} protein-coding genes (filtered from {len(moderate_genes)})")

# Save gene list
protein_coding_file = OUTPUT_DIR / 'high_confidence_genes_protein_coding.txt'
with open(protein_coding_file, 'w') as f:
    for gene in moderate_protein_coding:
        f.write(f"{gene}\n")
print(f"  Saved: {protein_coding_file.name}")
print()

# ============================================================================
# Generate detailed statistics for each gene set
# ============================================================================

print("Step 7: Generating detailed statistics...")

for gene_set_name, genes in [
    ('Stringent', stringent_genes),
    ('Moderate', moderate_genes),
    ('Liberal', liberal_genes),
    ('Protein-coding', moderate_protein_coding)
]:
    # Get statistics for each gene
    gene_stats = []
    for gene in genes:
        gene_df = df[df['nearest_gene'] == gene]
        gene_stats.append({
            'Gene': gene,
            'N_SNPs': len(gene_df),
            'Min_log10p': gene_df['neg_log10_p'].min(),
            'Max_log10p': gene_df['neg_log10_p'].max(),
            'Mean_log10p': gene_df['neg_log10_p'].mean(),
            'Max_stdIHS': gene_df['std_ihs'].max(),
            'Has_Overlap': 'inside_gene' in gene_df['overlap_type'].values,
            'Min_Distance': gene_df['distance'].abs().min()
        })

    stats_df = pd.DataFrame(gene_stats)
    stats_df = stats_df.sort_values('Max_log10p', ascending=False)

    # Save statistics
    stats_file = OUTPUT_DIR / f'high_confidence_genes_{gene_set_name.lower()}_stats.tsv'
    stats_df.to_csv(stats_file, sep='\t', index=False)
    print(f"  Saved {gene_set_name} statistics: {stats_file.name}")

print()

# ============================================================================
# Summary report
# ============================================================================

print("=" * 80)
print("SUMMARY")
print("=" * 80)
print()
print(f"{'Gene Set':<25} {'N Genes':<10} {'File':<40}")
print("-" * 80)
print(f"{'Original (current)':<25} {13:<10} {'high_confidence_X_genes.txt':<40}")
print(f"{'Stringent':<25} {len(stringent_genes):<10} {'high_confidence_genes_stringent.txt':<40}")
print(f"{'Moderate':<25} {len(moderate_genes):<10} {'high_confidence_genes_moderate.txt':<40}")
print(f"{'Liberal':<25} {len(liberal_genes):<10} {'high_confidence_genes_liberal.txt':<40}")
print(f"{'Protein-coding (moderate)':<25} {len(moderate_protein_coding):<10} {'high_confidence_genes_protein_coding.txt':<40}")
print()

print("Recommendations:")
print("  - For GO enrichment: Use 'Moderate' or 'Protein-coding' (30-60 genes)")
print("  - For manuscript: Report 'Stringent' as high-confidence core set")
print("  - For exploration: Use 'Liberal' to capture broader patterns")
print()

# Save summary
summary_file = OUTPUT_DIR / 'high_confidence_gene_sets_summary.txt'
with open(summary_file, 'w') as f:
    f.write("High-Confidence Gene Set Extraction Summary\n")
    f.write("=" * 80 + "\n\n")
    f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    f.write(f"{'Gene Set':<25} {'N Genes':<10}\n")
    f.write("-" * 40 + "\n")
    f.write(f"{'Original':<25} {13:<10}\n")
    f.write(f"{'Stringent':<25} {len(stringent_genes):<10}\n")
    f.write(f"{'Moderate':<25} {len(moderate_genes):<10}\n")
    f.write(f"{'Liberal':<25} {len(liberal_genes):<10}\n")
    f.write(f"{'Protein-coding':<25} {len(moderate_protein_coding):<10}\n")

print(f"Summary saved: {summary_file}")
print()
print("=" * 80)
print("Gene set extraction complete!")
print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 80)
