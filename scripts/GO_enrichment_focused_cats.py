#!/usr/bin/env python3
""" 
GO enrichment analysis using geneinfo.go_enrichment()
Script
1. Takes study genes(13X genes under selection)
2. For each functional category that showed significant overlap with fischer test, gets relevant GO terms
3. Tests enrichment within those GO term subsets
4. Returns enriched GO terms with statistics 

Thus performs proper GO enrichment (using go.go_enrichment()) for: 
1. female fertility
2. calcium homeostasis 
3. neuroimmuna crosstalk 

Those were categories that showed enrichment in fischer's exact tests. 
Follows example from 
https://munch-group.org/geneinfo/pages/go_enrichment.html

Author: Mit Van Bruggen 
Updated: 2025-11-20
"""

import pandas as pd
from pathlib import Path
from geneinfo.genelist import GeneList as glist
import geneinfo.ontology as go
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import numpy as np

go.email('au799024@uni.au.dk')

# Setup paths
PROJECT_DIR = Path('/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib')
RESULTS_DIR = PROJECT_DIR / 'results/analysis'
OUTPUT_DIR = RESULTS_DIR / 'functional_enrichment/GO_proper_enrichment'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("="*80)
print("Proper GO Enrichment Analysis using geneinfo.go_enrichment()")
print("="*80)
print()

#==============================================================================
# Step 1: Define study genes
#==============================================================================

study_genes = ['DACH2', 'DIAPH2', 'FAAH2', 'GPC3', 'HTR2C',
               'IL1RAPL2', 'KDM6A', 'NAP1L2', 'NCBP2L',
               'TRPC5', 'ZMAT1', 'KLHL13', 'PCDH11X','ALG13','MID2','NCBP2L','TSC22D3','ATP6AP1','FAM50A','GDI1','RPL10']
print(f"Study genes ({len(study_genes)}):")
print(f"  {', '.join(study_genes)}")
print()

# Convert to GeneList
study_gl = glist(study_genes)

#==============================================================================
#Step 2: Define GO term search patterns for each category (inspired by explore_GO_hierarchies.py)
#==============================================================================

print("Step 1: Defining GO term patterns for significant categories...")
print()

# Focus on categories that showed enrichment
FOCUSED_GO_PATTERNS = {
    # PRIORITY 1: Significant after FDR (p < 0.05)
    'female_fertility': {
        'pattern': r'female.*gamete.*generation|oogenesis|female.*sexual.*characteristic',
        'description': 'Female fertility and reproduction ',
        'priority': 1,
        'expected_genes': ['DACH2', 'DIAPH2']
    },
    
    'calcium_homeostasis': {
        'pattern': r'calcium.*channel.*activity|calcium.*ion.*transmembrane.*transport|calcium.*channel.*activity',
        'description': 'Calcium signaling and homeostasis ',
        'priority': 1,
        'expected_genes': ['HTR2C', 'TRPC5']
    },
    
    # PRIORITY 2: Borderline significant (FDR < 0.10)
    'neuroimmune_crosstalk': {
        'pattern': r'cytokine.*signaling|immune.*response|innate.*immune|interleukin|cytokine-mediated.*signaling.*pathway',
        'description': 'Neuroimmune crosstalk (FDR=0.070, borderline)',
        'priority': 2,
        'expected_genes': ['GPC3', 'HTR2C', 'IL1RAPL2', 'OPN1LW', 'OPN1MW', 'OPN1MW2', 'OPN1MW3']
    },
    
    # PRIORITY 3: Suggestive (p < 0.15)
    'synapse_assembly': {
        'pattern': r'synapse|synaptic|glutamatergic|presynapse|postsynapse|neurotransmission',
        'description': 'Synaptic structure and function ',
        'priority': 3,
        'expected_genes': ['HTR2C', 'IL1RAPL2', 'TRPC5']
    },
    
    'axon_dendrite_morphogenesis': {
        'pattern': r'axon|dendrite|neurite|neuronal.*projection|growth.*cone|axonogenesis',
        'description': 'Neuron morphogenesis ',
        'priority': 3,
        'expected_genes': ['NAP1L2', 'TRPC5']
    },
    
    # PRIORITY 4: Exploratory (for completeness)
    'serotonin_signaling': {
        'pattern': r'serotonin|5.*HT|hydroxytryptamine',
        'description': 'Serotonin signaling (HTR2C)',
        'priority': 4,
        'expected_genes': ['HTR2C']
    },
    
    'chromatin_regulation': {
        'pattern': r'chromatin|nucleosome|histone|epigenetic',
        'description': 'Chromatin remodeling (NAP1L2)',
        'priority': 4,
        'expected_genes': ['NAP1L2']
    }
}

# Get X chromosome background
all_protein_coding = go.all_protein_coding()
gene_table = go.gene_annotation_table()
x_gene_symbols = set(gene_table[gene_table['chromosome'] == 'X']['Symbol'].tolist())
background_genes = [gene for gene in all_protein_coding if gene in x_gene_symbols]

print(f"Background: X chromosome protein-coding genes")
print(f"  Total: {len(background_genes)} genes")
print()

#==============================================================================
# Extract GO terms for each category (print statements added by claude)
#==============================================================================
category_terms = {}
category_info = {}

for category, config in FOCUSED_GO_PATTERNS.items():
    pattern = config['pattern']
    priority = config['priority']
    
    print(f"  [{priority}] {category}:")
    print(f"      {config['description']}")
    print(f"      Pattern: {pattern}")
    
    # Get GO terms
    terms = go.get_terms_for_go_regex(pattern)
    category_terms[category] = terms
    category_info[category] = config
    
    print(f"      Found {len(terms)} GO terms")
    print()

#==============================================================================
# Perform GO enrichment for each category (print statements by claude)
#==============================================================================

print("="*80)
print("Step 2: GO Enrichment Analysis (Priority Order)")
print("="*80)

all_results = []
enrichment_summary = []

# Process by priority
for priority in [1, 2, 3, 4]:
    categories_at_priority = {k: v for k, v in FOCUSED_GO_PATTERNS.items() if v['priority'] == priority}
    
    if not categories_at_priority:
        continue
    
    print(f"\n{'='*80}")
    print(f"PRIORITY {priority}")
    print('='*80)
    
    for category, config in categories_at_priority.items():
        terms = category_terms[category]
        
        if len(terms) == 0:
            print(f"\n{category}: No GO terms found (skipped)")
            continue
        
        print(f"\n{category.upper()}")
        print(f"  {config['description']}")
        print(f"  Testing {len(terms)} GO terms")
        print(f"  Expected genes: {', '.join(config['expected_genes'])}")
        print("-"*80)
        
        try:
            # Perform GO enrichment
            results = go.go_enrichment(
                gene_list=study_genes,
                terms=terms,
                background_genes=background_genes, #X chromosome background
                list_study_genes=True,
                alpha=0.05)
            
            if results.empty:
                print("  No enrichment found")
                enrichment_summary.append({
                    'category': category,
                    'priority': priority,
                    'n_go_terms': len(terms),
                    'n_enriched': 0,
                    'top_term': None,
                    'top_p_fdr': None,
                    'genes_found': None
                })
                continue
            
            # Add metadata
            results['category'] = category
            results['priority'] = priority
            
            # Count significant results
            n_sig_001 = len(results[results['p_fdr_bh'] < 0.001])
            n_sig_01 = len(results[results['p_fdr_bh'] < 0.01])
            n_sig_05 = len(results[results['p_fdr_bh'] < 0.05])
            
            print(f"\n  Significant terms:")
            print(f"    FDR < 0.001: {n_sig_001}")
            print(f"    FDR < 0.01:  {n_sig_01}")
            print(f"    FDR < 0.05:  {n_sig_05}")
            
            # Show top 5 terms
            if n_sig_05 > 0:
                print(f"\n  Top enriched GO terms:")
                top_results = results[results['p_fdr_bh'] < 0.05].head(5)
                
                for idx, row in top_results.iterrows():
                    term_id = row['term_id']
                    p_fdr = row['p_fdr_bh']
                    ratio = row.get('ratio', 'N/A')
                    genes = row.get('study_genes', '')
                    
                    # Get term name from GO DAG
                    try:
                        godag = go.GODag(go.cache_dir + '/go-basic.obo')
                        term_name = godag[term_id].name if term_id in godag else 'Unknown'
                    except:
                        term_name = 'Unknown'
                    
                    sig = "***" if p_fdr < 0.001 else "**" if p_fdr < 0.01 else "*"
                    
                    print(f"\n    {term_id} {sig}")
                    print(f"      {term_name}")
                    print(f"      FDR={p_fdr:.2e}, ratio={ratio:.3f}")
                    print(f"      Genes: {genes}")
            
            # Store summary
            if not results.empty:
                top_row = results.iloc[0]
                enrichment_summary.append({
                    'category': category,
                    'priority': priority,
                    'n_go_terms': len(terms),
                    'n_enriched': n_sig_05,
                    'top_term': top_row['term_id'],
                    'top_p_fdr': top_row['p_fdr_bh'],
                    'genes_found': top_row.get('study_genes', '')
                })
            
            # Store all results
            all_results.append(results)
            
        except Exception as e:
            print(f"  ✗ ERROR: {e}")
            continue
        
        print("-"*80)

#==============================================================================
# Save results
#==============================================================================

print("\n" + "="*80)
print("Step 3: Saving results...")
print("="*80)
print()

if len(all_results) > 0:
    # Combine all results
    combined_results = pd.concat(all_results, ignore_index=True)
    combined_results = combined_results.sort_values('p_fdr_bh')
    
    # Save full results
    full_file = OUTPUT_DIR / 'go_enrichment_all_results.tsv'
    combined_results.to_csv(full_file, sep='\t', index=False)
    print(f"Saved: {full_file}")
    
    # Save significant only
    sig_results = combined_results[combined_results['p_fdr_bh'] < 0.05]
    if not sig_results.empty:
        sig_file = OUTPUT_DIR / 'go_enrichment_significant_FDR05.tsv'
        sig_results.to_csv(sig_file, sep='\t', index=False)
        print(f"✓ Saved: {sig_file}")
        print(f"  {len(sig_results)} GO terms with FDR < 0.05")
    
    # Save per-category results (priority 1 & 2 only)
    for category in ['female_fertility', 'calcium_homeostasis', 'neuroimmune_crosstalk']:
        cat_results = combined_results[combined_results['category'] == category]
        if not cat_results.empty:
            cat_file = OUTPUT_DIR / f'go_enrichment_{category}.tsv'
            cat_results.to_csv(cat_file, sep='\t', index=False)
            print(f"✓ Saved: {cat_file}")

# Save summary
if enrichment_summary:
    summary_df = pd.DataFrame(enrichment_summary)
    summary_df = summary_df.sort_values('priority')
    summary_file = OUTPUT_DIR / 'enrichment_summary.tsv'
    summary_df.to_csv(summary_file, sep='\t', index=False)
    print(f" Saved: {summary_file}")

#==============================================================================
# Create visualization (claude)
#==============================================================================

print("\n" + "="*80)
print("Step 4: Creating visualizations...")
print("="*80)
print()

if len(all_results) > 0 and len(sig_results) > 0:
    # Plot: Top enriched terms
    fig, ax = plt.subplots(figsize=(14, 10))
    
    # Get top 15 terms
    top_terms = sig_results.head(15).copy()
    
    # Create labels
    labels = []
    for idx, row in top_terms.iterrows():
        try:
            godag = go.GODag(go.cache_dir + '/go-basic.obo')
            term_name = godag[row['term_id']].name if row['term_id'] in godag else row['term_id']
        except:
            term_name = row['term_id']
        
        # Truncate long names
        if len(term_name) > 50:
            term_name = term_name[:47] + "..."
        
        category = row['category']
        labels.append(f"{term_name}\n({category})")
    
    # Plot
    y_pos = np.arange(len(top_terms))
    
    # Color by category
    category_colors = {
        'female_fertility': '#e74c3c',
        'calcium_homeostasis': '#3498db',
        'neuroimmune_crosstalk': '#2ecc71',
        'synapse_assembly': '#f39c12',
        'axon_dendrite_morphogenesis': '#9b59b6',
        'serotonin_signaling': '#1abc9c',
        'chromatin_regulation': '#34495e'
    }
    
    colors = [category_colors.get(cat, 'gray') for cat in top_terms['category']]
    
    ax.barh(y_pos, -np.log10(top_terms['p_fdr_bh']), color=colors, alpha=0.8)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel('-log10(FDR-corrected p-value)', fontsize=12)
    ax.set_title('Top GO Terms Enriched in X Chromosome Genes Under Selection', 
                 fontsize=14, fontweight='bold')
    ax.axvline(-np.log10(0.05), color='red', linestyle='--', linewidth=1.5, label='FDR = 0.05')
    ax.axvline(-np.log10(0.01), color='darkred', linestyle='--', linewidth=1.5, label='FDR = 0.01')
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=color, label=cat.replace('_', ' ').title()) 
                      for cat, color in category_colors.items() 
                      if cat in top_terms['category'].values]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=9)
    
    plt.tight_layout()
    plot_file = OUTPUT_DIR / 'top_enriched_go_terms.png'
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {plot_file}")
    plt.close()

#==============================================================================
# Summary report (print statements by claude)
#==============================================================================

print("\n" + "="*80)
print("FINAL SUMMARY")
print("="*80)
print()

if enrichment_summary:
    print("Enrichment by category:")
    print()
    summary_df = pd.DataFrame(enrichment_summary).sort_values('priority')
    
    for _, row in summary_df.iterrows():
        priority_marker = "★" * row['priority']
        sig_marker = "✓" if row['n_enriched'] > 0 else ""
        
        print(f"  [{priority_marker}] {row['category']:30s}: {row['n_enriched']:3d} sig. terms {sig_marker}")
        if row['n_enriched'] > 0:
            print(f"      Top term: {row['top_term']} (FDR={row['top_p_fdr']:.2e})")
            print(f"      Genes: {row['genes_found']}")
    
    print()
    print("Priority legend: ★ = high priority (FDR < 0.05 in Fisher's test)")
    print("                 ★★ = moderate (FDR < 0.10)")
    print("                 ★★★ = exploratory")

print()
print("="*80)
print("Analysis Complete!")
print("="*80)
print()
print("Key files:")
print(f"  Main results: {OUTPUT_DIR / 'go_enrichment_all_results.tsv'}")
print(f"  Significant:  {OUTPUT_DIR / 'go_enrichment_significant_FDR05.tsv'}")
print(f"  Summary:      {OUTPUT_DIR / 'enrichment_summary.tsv'}")
print(f"  Plot:         {OUTPUT_DIR / 'top_enriched_go_terms.png'}")
print()

