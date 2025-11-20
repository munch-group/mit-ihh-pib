#!/usr/bin/env python3
"""
GO Hierarchy Explorer - Create focused gene sets from GO term subtrees

This script provides tools to:
1. Start with specific, biologically meaningful GO terms
2. Visualize GO hierarchies (DAGs)
3. Navigate and explore parent/child relationships
4. Extract coherent subtrees
5. Generate focused gene sets from selected subtrees

Author: MIT van Bruggen
Date: 2025-11-18
"""

import pandas as pd
from pathlib import Path
import geneinfo.ontology as go
from geneinfo.genelist import GeneList as glist
import matplotlib.pyplot as plt
#import networkx as nx
from typing import List, Dict, Set, Tuple
import json

# Setup
go.email('au799024@uni.au.dk')

# Load GO DAG (Gene Ontology Directed Acyclic Graph)
# This gives us access to the GO hierarchy
godag = go.GODag(go.cache_dir + '/go-basic.obo')

PROJECT_DIR = Path('/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib')
RESULTS_DIR = PROJECT_DIR / 'results/analysis'
OUTPUT_DIR = RESULTS_DIR / 'functional_enrichment/GO_hierarchies'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("GO Hierarchy Explorer")
print("=" * 80)
print()

#==============================================================================
# Core GO Terms - Specific, biologically meaningful starting points
#==============================================================================
#Used amigo.geneontology.org to find GO terms for each category (search terms next to GO terms)
#Used QuickGo to verify 

CORE_GO_TERMS: dict[str, dict[str, list[str] | str]] = {
    #Reproductive cluster
    'female_fertility': {
        'seed_terms': [
            'GO:0046545', #development of primary female sexual characteristics
            'GO:0007292', #female gamete generation (DIAPH2)
            'GO:0048477', #oogenesis (DIAPH2)
            'GO:0001556', #oocyte maturation
            'GO:0060065', #uterus development
            'GO:0001890', #placenta development
        ],
        'description': 'Female fertility, oocyte and reproductive organ development'
    },

    'ovarian_function': {
        'seed_terms': [
            'GO:0030716',  # oocyte fate determination
            'GO:0046620',  # regulation of organ growth (ovary context)
            'GO:0008585',  # female gonad development
        ],
        'description': 'Ovarian development and function, DIAPH2 linked to POF2A'
    },

    #Neuroendocrine function 
    'neuroendocrine_reproduction':{
        'seed_terms':[
            'GO:0007631', #feeding behavior (HTR2C)
            'GO:0032098', #regulation of appetite (HTR2C)
            'GO:0007610', #behavior(HTR2C, TRPC5) 
            'GO:0007200', #phospholipase C-activating GPCR signaling (HTR2C)
        ],
        'description': 'neuroendocrine control of reproduction and energy balance'
    },
    #Neurodevelopment
    'synapse_assembly':{
        'seed_terms':[
            'GO:0098978', #glutamatergic synapse (IL1RAPL2)
            'GO:1905606', #regulation of presynapse assembly (IL1RAPL2)
            'GO:0007399', #nervous system development (IL1RAPL2)
            'GO:0007417', #central nervous system development (IL1RAPL2)
            'GO:0007268', #chemical synaptic transmission
            'GO:0050808', #synapse organization
        ],
        'description':'Synaptic formation and function'
    }, 

    'axon_dendrite_morphogenesis': {
        'seed_terms': [
            'GO:0045773',  # positive regulation of axon extension (TRPC5)
            'GO:0050774',  # negative regulation of dendrite morphogenesis (TRPC5)
            'GO:0007409',  # axonogenesis
            'GO:0030182',  # neuron differentiation (NAP1L2, IL1RAPL2, TRPC5)
        ],
        'description': 'Neuron structure and differentiation',
    },

    #Calcium signaling 
    'calcium_homeostasis':{
        'seed_terms':[
            'GO:0005262',  # calcium channel activity (TRPC5)
            'GO:0051480',  # regulation of cytosolic calcium (TRPC5)
            'GO:0070588',  # calcium ion transmembrane transport (TRPC5)
            'GO:0015279',  # store-operated calcium channel (TRPC5)
        ],
        'description':'Calcium channel activity and signaling'
    }, 

    #Wnt/Developmental signaling (GPC3 cluster)
   'wnt_signaling': {
        'seed_terms': [
            'GO:0090263',  # positive regulation of canonical Wnt (GPC3)
            'GO:0060828',  # regulation of canonical Wnt pathway (GPC3)
            'GO:0001822',  # kidney development (GPC3)
            'GO:0072111',  # cell proliferation involved in kidney development (GPC3)
        ],
        'description': 'Wnt signaling and organogenesis'
    },
    #Immunity 
    'neuroimmune_crosstalk':{
        'seed_terms':[
            'GO:0019221', # cytokine-mediated signaling pathway (IL1RAPL2)
            'GO:0002758', #innate immune response-activating signaling 
            'GO:0007165', #signal transduction (IL1RAPL2)
        ], 
        'description':'Immune signaling (neuro-immune)'
    },
}   


#==============================================================================
# Function 1: Get GO term information
#==============================================================================

def get_go_term_info(go_id: str) -> Dict:
    """
    Get detailed information about a GO term
    Args:
        go_id: GO identifier (e.g., 'GO:0097722')
    Returns:
        Dictionary with GO term details
    """
    try:
        # Access term from GO DAG
        if go_id in godag:
            term = godag[go_id]
            info = {
                'go_id': go_id,
                'name': term.name,
                'exists': True
            }
        else:
            info = {
                'go_id': go_id,
                'name': 'Unknown',
                'exists': False
            }

        return info
    except Exception as e:
        # If there's an error accessing the term
        return {
            'go_id': go_id,
            'name': f'Error: {str(e)}',
            'exists': False
        }


#==============================================================================
# Function 2: Get all descendants of a GO term (subtree)
#==============================================================================

def get_go_descendants(go_id: str, include_self: bool = True) -> Set[str]:
    """
    Get all descendant GO terms (children, grandchildren, etc.)
    This will extract the entire subtree below a GO term, which will represent
    a coherent biological module.
    Args:
        go_id: GO identifier
        include_self: Whether to include the seed term itself
    Returns:
        Set of GO IDs in the subtree
    """
    try:
        descendants = set()

        # Check if term exists in GO DAG
        if go_id not in godag:
            print(f"Warning: {go_id} not found in GO DAG")
            return {go_id} if include_self else set()

        term = godag[go_id]

        # Get children (terms more specific than this one)
        children = term.children

        if children:
            # Add children IDs
            for child in children:
                descendants.add(child.id)
                # Recursively get descendants of children
                descendants.update(get_go_descendants(child.id, include_self=False))

        if include_self:
            descendants.add(go_id)

        return descendants

    except Exception as e:
        print(f"Warning: Could not get descendants for {go_id}: {e}")
        return {go_id} if include_self else set()


#==============================================================================
# Function 3: Get all ancestors of a GO term
#==============================================================================

def get_go_ancestors(go_id: str, include_self: bool = False) -> Set[str]:
    """
    Get all ancestor GO terms (parents, grandparents, etc.)
    Shows where a term fits in the broader hierarchy.
    Args: idem as above
    Returns:
        Set of GO IDs that are ancestors
    """
    try:
        ancestors = set()

        # Check if term exists in GO DAG
        if go_id not in godag:
            print(f"Warning: {go_id} not found in GO DAG")
            return {go_id} if include_self else set()

        term = godag[go_id]

        # Get parents (terms more general than this one)
        parents = term.parents

        if parents:
            # Add parent IDs
            for parent in parents:
                ancestors.add(parent.id)
                # Recursively get ancestors of parents
                ancestors.update(get_go_ancestors(parent.id, include_self=False))

        if include_self:
            ancestors.add(go_id)

        return ancestors

    except Exception as e:
        print(f"Warning: Could not get ancestors for {go_id}: {e}")
        return {go_id} if include_self else set()


#==============================================================================
# Function 4: Build subtree from multiple seed terms
#==============================================================================

def build_subtree_from_seeds(seed_terms: List[str],
                              max_depth: int = None) -> Dict[str, Set[str]]:
    """
    Build a coherent subtree from multiple seed terms
    This combines multiple related GO terms into a single focused gene set.
    Args:
        seed_terms: List of GO IDs to start from
        max_depth: Maximum depth to traverse (None = unlimited)
    Returns:
        Dict with:
        - 'seeds': the input seed terms
        - 'descendants': all descendant terms
        - 'all_terms': seeds + descendants
    """
    all_descendants = set()

    for seed in seed_terms:
        descendants = get_go_descendants(seed, include_self=False)
        all_descendants.update(descendants)

    return {
        'seeds': set(seed_terms),
        'descendants': all_descendants,
        'all_terms': set(seed_terms) | all_descendants
    }


#==============================================================================
# Function 5: Get genes for a subtree
#==============================================================================

def get_genes_for_subtree(go_terms: Set[str]) -> Tuple[pd.DataFrame, List[str]]:
    """
    Get all genes annotated to any GO term in a subtree, filtered for X chromosome only
    Args:
        go_terms: Set of GO term IDs
    Returns:
        Tuple of (DataFrame with all annotations for X chromosome genes, list of unique gene symbols)
    """
    if not go_terms:
        return pd.DataFrame(), []

    try:
        # Convert set to list for geneinfo
        terms_list = list(go_terms)

        # Get genes for GO terms
        genes_df = go.get_genes_for_go_terms(terms_list)

        # Get X chromosome gene symbols from gene annotation table
        gene_table = go.gene_annotation_table()
        x_gene_symbols = set(gene_table[gene_table['chromosome'] == 'X']['Symbol'].tolist())

        # Filter genes to only those on X chromosome
        if not genes_df.empty:
            genes_df = genes_df[genes_df['symbol'].isin(x_gene_symbols)]
            unique_genes = genes_df['symbol'].unique().tolist()
        else:
            unique_genes = []

        return genes_df, unique_genes

    except Exception as e: #thx claude
        print(f"Warning: Error getting genes: {e}")
        return pd.DataFrame(), []


#==============================================================================
# Function 6: Visualize GO hierarchy using geneinfo's DAG like in Kasper's guide at https://munch-group.org/geneinfo/pages/go_graphs.html
#==============================================================================

def visualize_go_dag(go_terms: List[str], output_file: Path = None,
                     title: str = "GO Term Hierarchy") -> None:
    """
    Visualize GO term hierarchy using geneinfo's built-in DAG plotting
    Args:
        go_terms: List of GO term IDs to visualize
        output_file: Path to save the visualization (PNG format)
        title: Title for the plot
    """
    try:
        # Use geneinfo's built-in DAG visualization
        print(f"  Creating GO DAG for {len(go_terms)} terms...")

        if output_file:
            # Save to file using plot_gos
            go.plot_gos(str(output_file), go_terms, godag)
            print(f"  Saved DAG to: {output_file}")
        else:
            # Display interactively (for in jupyter notebook where I test some things)
            go.show_go_dag_for_terms(go_terms)
            print(f"  GO DAG displayed")

    except Exception as e:
        print(f"  Warning: Could not create GO DAG: {e}")

#==============================================================================
# Function 7: Create summary report for a subtree
#==============================================================================

def create_subtree_report(category_name: str,
                         subtree_info: Dict,
                         genes_df: pd.DataFrame,
                         unique_genes: List[str]) -> pd.DataFrame:
    """
    Create a detailed report for a subtree-based gene set
    Returns:
        DataFrame with summary statistics
    """
    n_seeds = len(subtree_info['seeds'])
    n_descendants = len(subtree_info['descendants'])
    n_total_terms = len(subtree_info['all_terms'])
    n_genes = len(unique_genes)

    report = pd.DataFrame([{
        'category': category_name,
        'n_seed_terms': n_seeds,
        'n_descendant_terms': n_descendants,
        'n_total_terms': n_total_terms,
        'n_unique_genes': n_genes,
        'genes_per_term': n_genes / n_total_terms if n_total_terms > 0 else 0
    }])

    return report


#==============================================================================
# EXECUTION of the functions defined above (parts by claude code + print statements fully by claude code, as well as checks and warnings )
#==============================================================================

print("Step 1: Validating seed GO terms...")
print("-" * 80)

# Validate the seed terms
for category, config in CORE_GO_TERMS.items():
    print(f"\n{category}: {config['description']}")
    print(f"  Seed terms:")
    for term_id in config['seed_terms']:
        info = get_go_term_info(term_id)
        status = "✓" if info['exists'] else "✗"
        print(f"    {status} {term_id}: {info['name']}")


print("\n" + "=" * 80)
print("Step 2: Building subtrees from seed terms...")
print("-" * 80)

# Build subtrees for each category
category_subtrees = {}
category_genes = {}
all_reports = []

for category, config in CORE_GO_TERMS.items():
    print(f"\n{category}: {config['description']}")
    print(f"  Building subtree from {len(config['seed_terms'])} seed terms...")

    # Build subtree
    subtree = build_subtree_from_seeds(config['seed_terms'])
    category_subtrees[category] = subtree

    print(f"  → Found {len(subtree['descendants'])} descendant terms")
    print(f"  → Total: {len(subtree['all_terms'])} GO terms in subtree")

    # Get genes
    print(f"  Getting genes for all {len(subtree['all_terms'])} terms...")
    genes_df, unique_genes = get_genes_for_subtree(subtree['all_terms'])
    category_genes[category] = {
        'genes_df': genes_df,
        'unique_genes': unique_genes,
        'genelist': glist(unique_genes)
    }

    print(f"  → Found {len(unique_genes)} unique genes")

    # Create report
    report = create_subtree_report(category, subtree, genes_df, unique_genes)
    all_reports.append(report)


print("\n" + "=" * 80)
print("Step 3: Saving results...")
print("-" * 80)

# Save summary report
summary_df = pd.concat(all_reports, ignore_index=True)
summary_file = OUTPUT_DIR / 'subtree_summary.tsv'
summary_df.to_csv(summary_file, sep='\t', index=False)
print(f"\n✓ Saved summary: {summary_file}")

# Save subtree definitions (GO terms)
for category, subtree in category_subtrees.items():
    subtree_file = OUTPUT_DIR / f'subtree_terms_{category}.json'

    # Convert sets to lists for JSON
    subtree_json = {
        'category': category,
        'description': CORE_GO_TERMS[category]['description'],
        'seed_terms': sorted(list(subtree['seeds'])),
        'descendant_terms': sorted(list(subtree['descendants'])),
        'all_terms': sorted(list(subtree['all_terms']))
    }

    with open(subtree_file, 'w') as f:
        json.dump(subtree_json, f, indent=2)

    print(f"✓ Saved subtree: {subtree_file}")

# Save gene lists
for category, gene_data in category_genes.items():
    # Save full annotations
    annot_file = OUTPUT_DIR / f'genes_annotated_{category}.tsv'
    gene_data['genes_df'].to_csv(annot_file, sep='\t', index=False)

    # Save gene list
    genes_file = OUTPUT_DIR / f'genes_unique_{category}.txt'
    with open(genes_file, 'w') as f:
        f.write('\n'.join(sorted(gene_data['unique_genes'])))

    print(f"✓ Saved genes: {genes_file}")


print("\n" + "=" * 80)
print("Step 4: Comparison with keyword approach")
print("-" * 80)

# Compare sizes
print("\nGene set sizes (hierarchy-based approach):")
for category, gene_data in category_genes.items():
    print(f"  {category:30s}: {len(gene_data['unique_genes']):5d} genes")

print("\n" + "=" * 80)
print("Step 5: Visualizing GO hierarchies (DAG plots)")
print("-" * 80)

# Create DAG visualizations for each category
dag_dir = OUTPUT_DIR / 'dag_plots'
dag_dir.mkdir(exist_ok=True)

print("\nCreating GO DAG visualizations for each category...")
print("Note: These show the hierarchical relationships between GO terms")
print()

for category, subtree in category_subtrees.items():
    print(f"\n{category}:")

    # Get the seed terms for this category
    seed_terms = list(subtree['seeds'])

    print(f"  Visualizing {len(seed_terms)} seed terms...")

    try:
        # Create DAG for seed terms and save to file
        output_file = dag_dir / f'dag_{category}.png'
        visualize_go_dag(seed_terms, output_file=output_file, title=f"{category} - Seed Terms")
        print(f"  ✓ DAG visualization complete")
    except Exception as e:
        print(f"  ✗ Could not create DAG: {e}")


print("\n" + "=" * 80)
print("Step 6: Example - Show single term hierarchy")
print("-" * 80)

# Show hierarchy for sperm motility as example
print("\nExample: Detailed hierarchy for sperm motility (GO:0097722)")
print("This shows what genes are annotated to this specific term:\n")

try:
    go.show_go_dag_for_gene('DYNLT3')  # Example gene
except:
    print("  Could not show example gene DAG")

print("\n" + "=" * 80)
print("DONE!")
print("=" * 80)
print("\nNext steps:")
print("1. Review the subtree_summary.tsv to see gene set sizes")
print("2. Examine individual subtree JSON files to see which GO terms are included")
print("3. Look at the GO DAG visualizations to understand term relationships")
print("4. Adjust seed terms if needed to make sets more/less specific")
print("5. Use the gene lists for enrichment analysis with your study genes")
print()
print("Key files created:")
print(f"  - Summary: {OUTPUT_DIR / 'subtree_summary.tsv'}")
print(f"  - Gene lists: {OUTPUT_DIR / 'genes_unique_*.txt'}")
print(f"  - Subtree definitions: {OUTPUT_DIR / 'subtree_terms_*.json'}")
print(f"  - DAG plots: {dag_dir / 'dag_*.png'}")
print()
