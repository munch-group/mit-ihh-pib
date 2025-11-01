#!/usr/bin/env python3
"""
Comprehensive GO Enrichment Analysis for X Chromosome Selection Study

Performs three key analyses:
A. High-confidence genes vs all X genes
B. High-confidence genes vs X genes under selection
C. X genes under selection vs autosomal genes under selection

Uses pre-downloaded GO annotations from NCBI.
"""

import pandas as pd
from pathlib import Path
import geneinfo.ontology as go
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
import geneinfo
import sys
from datetime import datetime

# Set email for NCBI
go.email('au799024@uni.au.dk')

# Define paths
BASE_DIR = Path('/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib')
GENE_LISTS_DIR = BASE_DIR / 'results/analysis/gene_annotation/gene_name_lists'
HIGH_SIG_DIR = BASE_DIR / 'results/analysis/high_significance_genes'
ENRICHMENT_DIR = BASE_DIR / 'results/analysis/functional_enrichment'
OUTPUT_DIR = ENRICHMENT_DIR
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("Comprehensive GO Enrichment Analysis for X Chromosome Selection")
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 80)
print()

# ============================================================================
# Load GO annotations
# ============================================================================

print("Step 1: Loading GO annotations...")

# For memory efficiency, we'll load GO annotations on-demand for each analysis
# First, just load the X chromosome annotations file to verify it exists
x_go_file = ENRICHMENT_DIR / 'X_chromosome_GO_annotations_complete.tsv'
if not x_go_file.exists():
    print(f"  ERROR: X chromosome GO annotations not found: {x_go_file}")
    print(f"  Please run retrieve_GO_terms_X_genes.py first!")
    sys.exit(1)

print(f"  Found X chromosome GO annotations file: {x_go_file.name}")
print()

# ============================================================================
# Load GO DAG
# ============================================================================

print("Step 2: Loading GO database...")
geneinfo_dir = Path(geneinfo.__file__).parent
obo_file = geneinfo_dir / 'cache' / 'go-basic.obo'

if not obo_file.exists():
    print(f"  ERROR: GO database not found at {obo_file}")
    sys.exit(1)

godag = GODag(str(obo_file), optional_attrs='relationship')
print(f"  Loaded {len(godag)} GO terms")
print()

# ============================================================================
# Helper function to load gene-to-GO mapping
# ============================================================================

def load_gene2go_from_df(df, gene_col='Gene', go_col='GO_ID'):
    """Convert GO annotations dataframe to gene-to-GO mapping."""
    gene2go = {}
    for _, row in df.iterrows():
        gene = row[gene_col]
        go_id = row[go_col]
        if pd.notna(gene) and pd.notna(go_id):
            if gene not in gene2go:
                gene2go[gene] = set()
            gene2go[gene].add(go_id)
    return gene2go

# ============================================================================
# Load gene lists
# ============================================================================

print("Step 3: Loading gene lists...")

def load_gene_list(filepath):
    """Load gene list from file."""
    with open(filepath) as f:
        genes = [line.strip() for line in f if line.strip()]
    return genes

# Study and background gene lists
# Try to use the refined moderate gene set, fall back to original if not found
moderate_file = HIGH_SIG_DIR / 'high_confidence_genes_moderate.txt'
if moderate_file.exists():
    high_confidence_genes = load_gene_list(moderate_file)
    print(f"  Using REFINED high-confidence gene set (moderate criteria)")
else:
    high_confidence_genes = load_gene_list(HIGH_SIG_DIR / 'high_confidence_X_genes.txt')
    print(f"  Using ORIGINAL high-confidence gene set")

x_selection_genes = load_gene_list(GENE_LISTS_DIR / 'X_protein_coding_genes.txt')
all_x_genes = load_gene_list(GENE_LISTS_DIR / 'X_chromosome_genes.txt')
autosome_selection_genes = load_gene_list(GENE_LISTS_DIR / 'autosome_protein_coding_in_selection.txt')

print(f"  High-confidence X genes: {len(high_confidence_genes)}")
print(f"  X genes under selection: {len(x_selection_genes)}")
print(f"  All X chromosome genes: {len(all_x_genes)}")
print(f"  Autosomal genes under selection: {len(autosome_selection_genes)}")
print()

# ============================================================================
# Helper function for enrichment analysis
# ============================================================================

def run_enrichment(study_genes, population_genes, gene2go_mapping, analysis_name):
    """Run GO enrichment analysis and return results."""

    print(f"Running Analysis: {analysis_name}")
    print(f"  Study set: {len(study_genes)} genes")
    print(f"  Population: {len(population_genes)} genes")

    # Filter to genes with GO annotations
    study_with_go = [g for g in study_genes if g in gene2go_mapping]
    pop_with_go = [g for g in population_genes if g in gene2go_mapping]

    print(f"  Study genes with GO: {len(study_with_go)}/{len(study_genes)}")
    print(f"  Population genes with GO: {len(pop_with_go)}/{len(population_genes)}")

    if len(study_with_go) == 0:
        print(f"  WARNING: No study genes have GO annotations!")
        return None

    if len(pop_with_go) == 0:
        print(f"  WARNING: No population genes have GO annotations!")
        return None

    # Run enrichment
    goeaobj = GOEnrichmentStudy(
        pop_with_go,  # population (background)
        gene2go_mapping,  # gene-to-GO associations
        godag,  # GO DAG
        propagate_counts=True,
        alpha=0.05,  # significance threshold
        methods=['fdr_bh']  # Benjamini-Hochberg FDR correction
    )

    results = goeaobj.run_study(study_with_go)

    print(f"  Total GO terms tested: {len(results)}")
    sig_results = [r for r in results if r.p_fdr_bh < 0.05]
    print(f"  Significant GO terms (FDR < 0.05): {len(sig_results)}")
    print()

    return results

def save_results(results, output_prefix, analysis_name):
    """Save enrichment results to file."""

    if results is None:
        print(f"  Skipping save for {analysis_name} (no results)")
        return

    # Convert to DataFrame
    results_data = []
    for r in results:
        results_data.append({
            'GO_ID': r.GO,
            'GO_Name': r.name,
            'GO_Category': r.NS,
            'Enrichment': r.enrichment,
            'p_uncorrected': r.p_uncorrected,
            'p_fdr_bh': r.p_fdr_bh,
            'Study_Count': r.study_count,
            'Study_Total': r.study_n,
            'Pop_Count': r.pop_count,
            'Pop_Total': r.pop_n,
            'Study_Genes': ','.join(sorted(r.study_items)) if hasattr(r, 'study_items') else '',
            'Depth': r.depth,
        })

    df = pd.DataFrame(results_data)
    df = df.sort_values('p_fdr_bh')

    # Save all results
    all_file = OUTPUT_DIR / f'{output_prefix}_all.tsv'
    df.to_csv(all_file, sep='\t', index=False)
    print(f"  Saved all results: {all_file.name}")

    # Save significant only
    df_sig = df[df['p_fdr_bh'] < 0.05]
    if len(df_sig) > 0:
        sig_file = OUTPUT_DIR / f'{output_prefix}_significant.tsv'
        df_sig.to_csv(sig_file, sep='\t', index=False)
        print(f"  Saved significant results: {sig_file.name}")
    else:
        print(f"  No significant results to save")

    return df, df_sig

# ============================================================================
# Analysis A: High-confidence genes vs all X genes
# ============================================================================

print("=" * 80)
print("ANALYSIS A: High-confidence X genes vs All X chromosome genes")
print("=" * 80)
print()

# Load X chromosome GO annotations for Analyses A and B
print("Loading X chromosome GO annotations...")
x_go_annotations = pd.read_csv(x_go_file, sep='\t')
x_gene2go = load_gene2go_from_df(x_go_annotations)
print(f"  Created mapping for {len(x_gene2go)} X chromosome genes")
print()

results_a = run_enrichment(
    study_genes=high_confidence_genes,
    population_genes=all_x_genes,
    gene2go_mapping=x_gene2go,
    analysis_name="Analysis A"
)

df_a, df_a_sig = save_results(results_a, 'analysis_A_highconf_vs_all_X', 'Analysis A')

# ============================================================================
# Analysis B: High-confidence genes vs X genes under selection
# ============================================================================

print("=" * 80)
print("ANALYSIS B: High-confidence X genes vs X genes under selection")
print("=" * 80)
print()

results_b = run_enrichment(
    study_genes=high_confidence_genes,
    population_genes=x_selection_genes,
    gene2go_mapping=x_gene2go,
    analysis_name="Analysis B"
)

df_b, df_b_sig = save_results(results_b, 'analysis_B_highconf_vs_X_selection', 'Analysis B')

# ============================================================================
# Analysis C: X genes under selection vs autosomal genes under selection
# ============================================================================

print("=" * 80)
print("ANALYSIS C: X genes under selection vs Autosomal genes under selection")
print("=" * 80)
print()

# Free up memory from X-only annotations
del x_go_annotations
del x_gene2go
print("Freed X-specific GO annotations from memory")
print()

# Load all GO annotations for Analysis C
print("Loading all human GO annotations for Analysis C...")
print("  (This may take 1-2 minutes)")

try:
    all_go_annotations = go.go_annotation_table(taxid=9606)
    gene_info = go.gene_annotation_table(taxid=9606)
    gene_id_to_symbol = dict(zip(gene_info['GeneID'], gene_info['Symbol']))

    # Add gene symbols
    all_go_annotations['Symbol'] = all_go_annotations['GeneID'].map(gene_id_to_symbol)

    # Free gene_info to save memory
    del gene_info
    del gene_id_to_symbol

    print(f"  Retrieved {len(all_go_annotations)} total GO annotations")

    # Filter to genes in both study set and population (X genes + autosomal genes)
    # We need GO annotations for both groups to perform the enrichment test
    combined_genes = list(set(x_selection_genes + autosome_selection_genes))
    all_go_annotations = all_go_annotations[all_go_annotations['Symbol'].isin(combined_genes)]

    print(f"  Filtered to {len(all_go_annotations)} annotations for study and background genes")

    # Create mapping
    all_gene2go = load_gene2go_from_df(all_go_annotations, gene_col='Symbol')
    print(f"  Created mapping for {len(all_gene2go)} genes")
    print()

except Exception as e:
    print(f"  ERROR: Failed to load GO annotations for Analysis C: {e}")
    all_gene2go = None

if all_gene2go:
    # For Analysis C, the population must include BOTH study and background genes
    # This is required for goatools to work correctly
    # The enrichment test will compare X genes vs autosomal genes within this combined population
    combined_population = list(set(x_selection_genes + autosome_selection_genes))

    results_c = run_enrichment(
        study_genes=x_selection_genes,
        population_genes=combined_population,
        gene2go_mapping=all_gene2go,
        analysis_name="Analysis C"
    )

    df_c, df_c_sig = save_results(results_c, 'analysis_C_X_vs_autosome_selection', 'Analysis C')
else:
    results_c = None
    df_c = None
    df_c_sig = pd.DataFrame()

# ============================================================================
# Generate Summary Report
# ============================================================================

print("=" * 80)
print("SUMMARY REPORT")
print("=" * 80)
print()

summary_file = OUTPUT_DIR / 'GO_enrichment_summary_report.txt'
with open(summary_file, 'w') as f:
    f.write("=" * 80 + "\n")
    f.write("GO Enrichment Analysis Summary Report\n")
    f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write("=" * 80 + "\n\n")

    # Analysis A
    f.write("ANALYSIS A: High-confidence X genes vs All X chromosome genes\n")
    f.write("-" * 80 + "\n")
    f.write(f"Study genes: {len(high_confidence_genes)}\n")
    f.write(f"Background genes: {len(all_x_genes)}\n")
    if results_a:
        f.write(f"Total GO terms tested: {len(results_a)}\n")
        f.write(f"Significant GO terms (FDR < 0.05): {len(df_a_sig)}\n")
    f.write("\n")

    # Analysis B
    f.write("ANALYSIS B: High-confidence X genes vs X genes under selection\n")
    f.write("-" * 80 + "\n")
    f.write(f"Study genes: {len(high_confidence_genes)}\n")
    f.write(f"Background genes: {len(x_selection_genes)}\n")
    if results_b:
        f.write(f"Total GO terms tested: {len(results_b)}\n")
        f.write(f"Significant GO terms (FDR < 0.05): {len(df_b_sig)}\n")
    f.write("\n")

    # Analysis C
    f.write("ANALYSIS C: X genes under selection vs Autosomal genes under selection\n")
    f.write("-" * 80 + "\n")
    f.write(f"Study genes (X under selection): {len(x_selection_genes)}\n")
    f.write(f"Population (all genes under selection): X + autosomal = {len(x_selection_genes) + len(autosome_selection_genes)}\n")
    if results_c:
        f.write(f"Total GO terms tested: {len(results_c)}\n")
        f.write(f"Significant GO terms (FDR < 0.05): {len(df_c_sig)}\n")
    else:
        f.write(f"Analysis failed or was skipped\n")
    f.write("\n")

    f.write("=" * 80 + "\n")
    f.write("Output Files\n")
    f.write("=" * 80 + "\n")
    f.write("Analysis A: analysis_A_highconf_vs_all_X_*.tsv\n")
    f.write("Analysis B: analysis_B_highconf_vs_X_selection_*.tsv\n")
    f.write("Analysis C: analysis_C_X_vs_autosome_selection_*.tsv\n")

print(f"Summary report saved: {summary_file}")
print()

print("=" * 80)
print("GO Enrichment Analysis Complete!")
print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 80)
