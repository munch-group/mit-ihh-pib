#!/usr/bin/env python3

"""
Population Gene Overlap Analysis for iHS Selection Signals

Compares protein coding genes under selection (from iHS analysis) across three populations:
- EUR (European)
- EAS (East Asian)
- AFR (African)

Identifies:
1. Genes present in all 3 populations (shared selection)
2. Genes present in exactly 2 populations (pairwise shared)
3. Genes unique to each population (population-specific selection)

Author: Mit Van Bruggen
Date: 2025-12-10
"""

import pandas as pd
from pathlib import Path
from geneinfo.genelist import GeneList as glist
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib_venn import venn3, venn3_circles
import seaborn as sns

# Setup paths
PROJECT_DIR = Path('/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib')
RESULTS_BASE = PROJECT_DIR / 'results'

# Gene list paths for each population
GENE_LISTS = {
    'EUR': RESULTS_BASE / 'ihs_EUR/analysis/gene_annotation/gene_name_lists/X_protein_coding_genes.txt',
    'EAS': RESULTS_BASE / 'ihs_EAS/analysis/gene_annotation/gene_name_lists/X_protein_coding_genes.txt',
    'AFR': RESULTS_BASE / 'ihs_AFR/analysis/gene_annotation/gene_name_lists/X_protein_coding_genes.txt'
}

# Output directory
OUTPUT_DIR = PROJECT_DIR / 'results/iHSperpop_genes_overlap'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

#==============================================================================
# Step 1: Load gene lists for each population
#==============================================================================

gene_sets = {}
gene_lists_ginfo = {}

for pop, filepath in GENE_LISTS.items():
    if filepath.exists():
        with open(filepath, 'r') as f:
            genes = [line.strip() for line in f if line.strip()]
        gene_sets[pop] = set(genes)
        gene_lists_ginfo[pop] = glist(genes)
        print(f"  {pop}: {len(genes):3d} genes under selection")
    else:
        print(f"  WARNING: {pop} file not found at {filepath}")
        gene_sets[pop] = set()
        gene_lists_ginfo[pop] = glist([])

print()

#==============================================================================
# Step 2: Identify overlaps
#==============================================================================

# All three populations
all_three = gene_sets['EUR'] & gene_sets['EAS'] & gene_sets['AFR']

# Exactly two populations
eur_eas_only = (gene_sets['EUR'] & gene_sets['EAS']) - gene_sets['AFR']
eur_afr_only = (gene_sets['EUR'] & gene_sets['AFR']) - gene_sets['EAS']
eas_afr_only = (gene_sets['EAS'] & gene_sets['AFR']) - gene_sets['EUR']

# Unique to each population
eur_unique = gene_sets['EUR'] - gene_sets['EAS'] - gene_sets['AFR']
eas_unique = gene_sets['EAS'] - gene_sets['EUR'] - gene_sets['AFR']
afr_unique = gene_sets['AFR'] - gene_sets['EUR'] - gene_sets['EAS']

# Print summary
print("OVERLAP SUMMARY:")
print()
print(f"Genes in ALL 3 populations:         {len(all_three):3d}")
print(f"  {sorted(all_three)}")
print()

print(f"Genes in exactly 2 populations:     {len(eur_eas_only) + len(eur_afr_only) + len(eas_afr_only):3d}")
print(f"  EUR & EAS only:                   {len(eur_eas_only):3d}")
if eur_eas_only:
    print(f"    {sorted(eur_eas_only)}")
print(f"  EUR & AFR only:                   {len(eur_afr_only):3d}")
if eur_afr_only:
    print(f"    {sorted(eur_afr_only)}")
print(f"  EAS & AFR only:                   {len(eas_afr_only):3d}")
if eas_afr_only:
    print(f"    {sorted(eas_afr_only)}")
print()

print(f"Population-specific genes:          {len(eur_unique) + len(eas_unique) + len(afr_unique):3d}")
print(f"  EUR unique:                       {len(eur_unique):3d}")
if eur_unique:
    print(f"    {sorted(eur_unique)}")
print(f"  EAS unique:                       {len(eas_unique):3d}")
if eas_unique:
    print(f"    {sorted(eas_unique)}")
print(f"  AFR unique:                       {len(afr_unique):3d}")
if afr_unique:
    print(f"    {sorted(afr_unique)}")
print()

print("="*80)
print()

#==============================================================================
# Step 3: Create detailed results dataframe
#==============================================================================

# Collect all unique genes
all_genes = gene_sets['EUR'] | gene_sets['EAS'] | gene_sets['AFR']

results = []
for gene in sorted(all_genes):
    in_eur = gene in gene_sets['EUR']
    in_eas = gene in gene_sets['EAS']
    in_afr = gene in gene_sets['AFR']

    # Determine category
    n_pops = sum([in_eur, in_eas, in_afr])
    if n_pops == 3:
        category = 'All_3_populations'
        detail = 'EUR+EAS+AFR'
    elif n_pops == 2:
        if in_eur and in_eas:
            category = 'Two_populations'
            detail = 'EUR+EAS'
        elif in_eur and in_afr:
            category = 'Two_populations'
            detail = 'EUR+AFR'
        else:  # in_eas and in_afr
            category = 'Two_populations'
            detail = 'EAS+AFR'
    else:  # n_pops == 1
        category = 'Population_specific'
        if in_eur:
            detail = 'EUR_only'
        elif in_eas:
            detail = 'EAS_only'
        else:
            detail = 'AFR_only'

    results.append({
        'Gene': gene,
        'EUR': in_eur,
        'EAS': in_eas,
        'AFR': in_afr,
        'N_populations': n_pops,
        'Category': category,
        'Detail': detail
    })

results_df = pd.DataFrame(results)

# Save results
results_file = OUTPUT_DIR / 'population_overlap_results.tsv'
results_df.to_csv(results_file, sep='\t', index=False)

# Save gene lists by category
categories = {
    'all_three_populations': sorted(all_three),
    'EUR_EAS_only': sorted(eur_eas_only),
    'EUR_AFR_only': sorted(eur_afr_only),
    'EAS_AFR_only': sorted(eas_afr_only),
    'EUR_unique': sorted(eur_unique),
    'EAS_unique': sorted(eas_unique),
    'AFR_unique': sorted(afr_unique)
}

for category_name, genes in categories.items():
    if genes:
        output_file = OUTPUT_DIR / f'{category_name}.txt'
        with open(output_file, 'w') as f:
            f.write('\n'.join(genes))
        print(f"  {category_name:25s}: {len(genes):3d} genes -> {output_file.name}")

print()

#==============================================================================
# Step 4: Create visualization - Venn Diagram
#==============================================================================

# Set style
sns.set_style("whitegrid")

# Figure Venn diagram - larger to accommodate annotations
fig, ax = plt.subplots(figsize=(18, 14))

venn = venn3(
    [gene_sets['EUR'], gene_sets['EAS'], gene_sets['AFR']],
    set_labels=('EUR', 'EAS', 'AFR'),
    set_colors=('#3498db', '#e74c3c', '#2ecc71'),
    alpha=0.6,
    ax=ax
)
venn3_circles(
    [gene_sets['EUR'], gene_sets['EAS'], gene_sets['AFR']],
    linewidth=2,
    ax=ax
)

plt.title('Overlap of Genes Under Selection (iHS)\nAcross EUR, EAS, and AFR Populations',
          fontsize=16, fontweight='bold', pad=30)

# Add gene names with text boxes and arrows
# Each tuple: (gene_set, label, text_position, arrow_target, horizontal_align, vertical_align, arrow_curve)
annotations = [
    # EUR only - top left
    (eur_unique, 'EUR only', (-0.72, 0.30), (-0.28, 0.32), 'right', 'bottom', 'arc3,rad=0.2'),
    # EAS only - top right
    (eas_unique, 'EAS only', (0.72, 0.60), (0.28, 0.32), 'left', 'bottom', 'arc3,rad=-0.2'),
    # AFR only - bottom center
    (afr_unique, 'AFR only', (-0.20, -0.70), (0.0, -0.32), 'center', 'top', 'arc3,rad=0'),
    # EUR & EAS only - top center (skip if empty)
    (eur_eas_only, 'EUR & EAS', (0.0, 0.80), (0.0, 0.30), 'center', 'bottom', 'arc3,rad=0'),
    # EUR & AFR only - far left
    (eur_afr_only, 'EUR & AFR', (-0.82, -0.30), (-0.20, -0.10), 'right', 'center', 'arc3,rad=-0.2'),
    # EAS & AFR only - far right
    (eas_afr_only, 'EAS & AFR', (0.82, -0.15), (0.15, -0.05), 'left', 'center', 'arc3,rad=0.2'),
    # All 3 populations - TOP CENTER, arrow pointing straight down to center
    (all_three, 'All 3 pops', (-0.15, 0.75), (-0.05, 0.14), 'center', 'bottom', 'arc3,rad=0'),
]

for gene_set, label, text_pos, arrow_pos, ha, va, connection_style in annotations:
    if gene_set:  # Only annotate if there are genes in this region
        # Create gene list text
        gene_text = '\n'.join(sorted(gene_set))
        full_text = f"{label} (n={len(gene_set)}):\n{gene_text}"

        # Add text box with arrow
        ax.annotate(
            full_text,
            xy=arrow_pos,           # Arrow points to this location (inside Venn region)
            xytext=text_pos,        # Text box location
            fontsize=11,
            family='monospace',
            ha=ha,
            va=va,
            bbox=dict(
                boxstyle='round,pad=0.5',
                facecolor='wheat',
                alpha=0.9,
                edgecolor='#8B4513',
                linewidth=1.5
            ),
            arrowprops=dict(
                arrowstyle='-|>',
                connectionstyle=connection_style,
                color='#555555',
                lw=1.5,
                mutation_scale=12
            )
        )

# Adjust axis limits to accommodate all annotations
ax.set_xlim(-1.1, 1.1)
ax.set_ylim(-1.15, 0.95)
ax.axis('off')

venn_file = OUTPUT_DIR / 'venn_diagram_population_overlap.png'
plt.savefig(venn_file, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print(f"  Venn diagram saved to: {venn_file}")

#==============================================================================
# Step 5: Summary statistics
#==============================================================================

total_unique_genes = len(all_genes)
total_eur = len(gene_sets['EUR'])
total_eas = len(gene_sets['EAS'])
total_afr = len(gene_sets['AFR'])

print(f"Total unique genes under selection: {total_unique_genes}")
print()
print(f"Population totals:")
print(f"  EUR: {total_eur:3d} genes")
print(f"  EAS: {total_eas:3d} genes")
print(f"  AFR: {total_afr:3d} genes")
print()

print(f"Sharing statistics:")
print(f"  Genes shared by all 3 populations:   {len(all_three):3d} ({len(all_three)/total_unique_genes*100:5.1f}% of total unique genes)")
print(f"  Genes shared by any 2 populations:   {len(eur_eas_only) + len(eur_afr_only) + len(eas_afr_only):3d} ({(len(eur_eas_only) + len(eur_afr_only) + len(eas_afr_only))/total_unique_genes*100:5.1f}%)")
print(f"  Population-specific genes:           {len(eur_unique) + len(eas_unique) + len(afr_unique):3d} ({(len(eur_unique) + len(eas_unique) + len(afr_unique))/total_unique_genes*100:5.1f}%)")
print()

print(f"Pairwise overlaps:")
eur_eas_overlap = len(gene_sets['EUR'] & gene_sets['EAS'])
eur_afr_overlap = len(gene_sets['EUR'] & gene_sets['AFR'])
eas_afr_overlap = len(gene_sets['EAS'] & gene_sets['AFR'])

print(f"  EUR & EAS: {eur_eas_overlap:3d} genes ({eur_eas_overlap/min(total_eur, total_eas)*100:5.1f}% of smaller set)")
print(f"  EUR & AFR: {eur_afr_overlap:3d} genes ({eur_afr_overlap/min(total_eur, total_afr)*100:5.1f}% of smaller set)")
print(f"  EAS & AFR: {eas_afr_overlap:3d} genes ({eas_afr_overlap/min(total_eas, total_afr)*100:5.1f}% of smaller set)")
print()

