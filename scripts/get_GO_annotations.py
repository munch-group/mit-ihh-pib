#!/usr/bin/env python3
"""
Retrieve and analyze GO (Gene Ontology) terms for candidate genes

This script:
1. Retrieves GO terms for the 11 overlap genes (high confidence)
2. Categorizes genes by functional categories (reproduction, immunity, neurodevelopment)
3. Creates summary statistics and visualizations

Author: MIT van Bruggen
Date: 2025-10-30
"""

import pandas as pd
import numpy as np
from pathlib import Path
import re
from collections import defaultdict

# Import geneinfo
import geneinfo.ontology as go

# Set email for GO database access
go.email('vanbruggenmit@birc.au.dk')

# Setup paths
PROJECT_DIR = Path('/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib')
RESULTS_DIR = PROJECT_DIR / 'results/analysis'
OUTPUT_DIR = RESULTS_DIR / 'functional_enrichment'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Overlap genes (highest confidence)
OVERLAP_GENES = ['DACH2', 'DIAPH2', 'FAAH2', 'GPC3', 'HTR2C',
                 'IL1RAPL2', 'KDM6A', 'NAP1L2', 'NCBP2L', 'TRPC5', 'ZMAT1']

print("="*70)
print("GO Term Analysis for High-Confidence Genes")
print("="*70)
print()

#==============================================================================
# 1. Retrieve GO terms for each gene
#==============================================================================

def get_go_terms_for_gene(gene_symbol):
    """
    Retrieve GO terms for a gene using geneinfo

    Returns dictionary with GO terms organized by category
    """
    try:
        # Get GO term IDs for the gene
        go_term_ids = go.get_go_terms_for_genes([gene_symbol])

        if go_term_ids is None or len(go_term_ids) == 0:
            print(f"  ⚠ No GO terms found for {gene_symbol}")
            return None

        # Access the GO DAG directly to get term names and namespaces
        from goatools.obo_parser import GODag
        from pathlib import Path
        import os

        # Get the GO DAG from geneinfo's cache
        import geneinfo
        import os
        geneinfo_dir = os.path.dirname(geneinfo.__file__)
        obo_file = Path(geneinfo_dir) / 'cache' / 'go-basic.obo'

        if not obo_file.exists():
            raise FileNotFoundError(f"GO database not found at {obo_file}")

        godag = GODag(str(obo_file), optional_attrs='relationship')

        # Organize by GO category
        result = {
            'gene': gene_symbol,
            'biological_process': [],
            'molecular_function': [],
            'cellular_component': [],
            'all_terms': []
        }

        # Get descriptions for each GO term
        for go_id in go_term_ids:
            try:
                if go_id in godag:
                    go_term = godag[go_id]
                    go_name = go_term.name
                    go_namespace = go_term.namespace  # biological_process, molecular_function, cellular_component
                    go_desc = go_term.defn if hasattr(go_term, 'defn') else ''

                    # Store term info
                    term_info = {
                        'id': go_id,
                        'name': go_name,
                        'category': go_namespace,
                        'description': go_desc
                    }

                    result['all_terms'].append(term_info)

                    # Categorize
                    if go_namespace == 'biological_process':
                        result['biological_process'].append(go_name)
                    elif go_namespace == 'molecular_function':
                        result['molecular_function'].append(go_name)
                    elif go_namespace == 'cellular_component':
                        result['cellular_component'].append(go_name)

            except Exception as e:
                # Skip terms that cause errors
                continue

        if len(result['all_terms']) > 0:
            return result
        else:
            print(f"  ⚠ No valid GO term descriptions for {gene_symbol}")
            return None

    except Exception as e:
        print(f"  ✗ Error retrieving GO terms for {gene_symbol}: {e}")
        import traceback
        traceback.print_exc()
        return None


print("Step 1: Retrieving GO terms for 11 overlap genes...")
print()

all_go_data = {}
for gene in OVERLAP_GENES:
    print(f"  Retrieving GO terms for {gene}...")
    go_data = get_go_terms_for_gene(gene)
    if go_data:
        all_go_data[gene] = go_data
        print(f"    ✓ Found {len(go_data['all_terms'])} GO terms")

print(f"\n  Retrieved GO terms for {len(all_go_data)}/{len(OVERLAP_GENES)} genes")

#==============================================================================
# 2. Categorize by functional keywords
#==============================================================================

print("\nStep 2: Categorizing genes by functional keywords...")
print()

# Define keyword categories
KEYWORD_CATEGORIES = {
    'reproduction': [
        'meiosis', 'gamete', 'fertility', 'gonad', 'sperm', 'oocyte',
        'spermatogenesis', 'oogenesis', 'germ cell', 'reproductive',
        'fertilization', 'ovary', 'testis'
    ],
    'immunity': [
        'immune', 'inflammation', 'defense', 'cytokine', 'antibody',
        'inflammatory', 'innate immune', 'adaptive immune', 'lymphocyte',
        'T cell', 'B cell', 'interferon', 'interleukin', 'antigen'
    ],
    'neurodevelopment': [
        'neuron', 'synapse', 'cognition', 'brain', 'neural', 'axon',
        'dendrite', 'synaptic', 'neuronal', 'nervous system', 'learning',
        'memory', 'behavior', 'neurogenesis', 'neurotransmitter'
    ],
    'development': [
        'development', 'differentiation', 'morphogenesis', 'organogenesis',
        'embryo', 'pattern', 'specification', 'determination'
    ]
}

def categorize_go_terms(go_data, categories):
    """
    Categorize GO terms based on keyword matching

    Returns dictionary of categories found
    """
    matches = defaultdict(list)

    for term_dict in go_data['all_terms']:
        term_name = term_dict['name'].lower()

        for category, keywords in categories.items():
            for keyword in keywords:
                if keyword.lower() in term_name:
                    matches[category].append(term_dict['name'])
                    break

    return dict(matches)


# Categorize each gene
gene_categories = {}
for gene, go_data in all_go_data.items():
    categories = categorize_go_terms(go_data, KEYWORD_CATEGORIES)
    gene_categories[gene] = categories

    print(f"  {gene}:")
    if categories:
        for cat, terms in categories.items():
            print(f"    - {cat}: {len(terms)} terms")
    else:
        print(f"    - No category keywords found")

#==============================================================================
# 3. Create comprehensive GO annotation table
#==============================================================================

print("\nStep 3: Creating comprehensive GO annotation table...")

# Flatten all GO terms into table format
go_records = []
for gene, go_data in all_go_data.items():
    categories = gene_categories.get(gene, {})

    for term_dict in go_data['all_terms']:
        # Determine which functional categories this term matches
        matched_categories = []
        term_name = term_dict['name'].lower()
        for cat, keywords in KEYWORD_CATEGORIES.items():
            for keyword in keywords:
                if keyword.lower() in term_name:
                    matched_categories.append(cat)
                    break

        go_records.append({
            'Gene': gene,
            'GO_ID': term_dict['id'],
            'GO_Name': term_dict['name'],
            'GO_Category': term_dict['category'],
            'Functional_Categories': ';'.join(matched_categories) if matched_categories else 'other'
        })

go_table = pd.DataFrame(go_records)
output_file = OUTPUT_DIR / 'gene_GO_annotations.tsv'
go_table.to_csv(output_file, sep='\t', index=False)
print(f"  Saved: {output_file}")
print(f"  Total GO term annotations: {len(go_table)}")

#==============================================================================
# 4. Create gene summary table
#==============================================================================

print("\nStep 4: Creating gene functional category summary...")

summary_records = []
for gene in OVERLAP_GENES:
    if gene in all_go_data:
        go_data = all_go_data[gene]
        categories = gene_categories.get(gene, {})

        summary_records.append({
            'Gene': gene,
            'Total_GO_Terms': len(go_data['all_terms']),
            'Biological_Process': len(go_data['biological_process']),
            'Molecular_Function': len(go_data['molecular_function']),
            'Cellular_Component': len(go_data['cellular_component']),
            'Reproduction_Terms': len(categories.get('reproduction', [])),
            'Immunity_Terms': len(categories.get('immunity', [])),
            'Neurodevelopment_Terms': len(categories.get('neurodevelopment', [])),
            'Development_Terms': len(categories.get('development', [])),
            'Primary_Function': max(categories.keys(), key=lambda k: len(categories[k])) if categories else 'other'
        })
    else:
        summary_records.append({
            'Gene': gene,
            'Total_GO_Terms': 0,
            'Biological_Process': 0,
            'Molecular_Function': 0,
            'Cellular_Component': 0,
            'Reproduction_Terms': 0,
            'Immunity_Terms': 0,
            'Neurodevelopment_Terms': 0,
            'Development_Terms': 0,
            'Primary_Function': 'no_GO_data'
        })

summary_table = pd.DataFrame(summary_records)
summary_file = OUTPUT_DIR / 'gene_functional_summary.tsv'
summary_table.to_csv(summary_file, sep='\t', index=False)
print(f"  Saved: {summary_file}")

#==============================================================================
# 5. Summary statistics
#==============================================================================

print("\n" + "="*70)
print("SUMMARY STATISTICS")
print("="*70)
print()

print("1. GO Term Coverage:")
genes_with_go = len(all_go_data)
genes_total = len(OVERLAP_GENES)
print(f"   Genes with GO terms: {genes_with_go}/{genes_total} ({genes_with_go/genes_total*100:.1f}%)")
print()

print("2. Functional Category Distribution:")
category_counts = {
    'Reproduction': 0,
    'Immunity': 0,
    'Neurodevelopment': 0,
    'Development': 0,
    'Other/Unknown': 0
}

for gene in OVERLAP_GENES:
    categories = gene_categories.get(gene, {})
    if not categories:
        category_counts['Other/Unknown'] += 1
    else:
        # Assign to primary category (most terms)
        primary = max(categories.keys(), key=lambda k: len(categories[k]))
        if primary == 'reproduction':
            category_counts['Reproduction'] += 1
        elif primary == 'immunity':
            category_counts['Immunity'] += 1
        elif primary == 'neurodevelopment':
            category_counts['Neurodevelopment'] += 1
        elif primary == 'development':
            category_counts['Development'] += 1
        else:
            category_counts['Other/Unknown'] += 1

for category, count in category_counts.items():
    pct = count / genes_total * 100
    print(f"   {category:20s}: {count:2d} genes ({pct:5.1f}%)")

print()
print("3. Genes by Primary Function:")
for idx, row in summary_table.sort_values('Total_GO_Terms', ascending=False).iterrows():
    print(f"   {row['Gene']:12s} - {row['Primary_Function']:20s} ({row['Total_GO_Terms']} GO terms)")

print()
print("="*70)
print("Analysis complete!")
print("="*70)
print()
print("Output files:")
print(f"  {output_file}")
print(f"  {summary_file}")
