#!/usr/bin/env python3
"""
iHS Analysis Workflow Visualization

Visualizes the complete iHS analysis pipeline for all subpopulations:
1. SNPs with |Std.iHS| > q99 along the X chromosome
2. Selection regions defined by clustering significant SNPs
3. Genes annotated within selection regions

Creates multi-panel figures showing the analysis workflow for each population.

Author: Mit Van Bruggen
Date: 2025-12-22
"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import warnings
warnings.filterwarnings('ignore')

BASE_DIR = Path('/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib')
RESULTS_DIR = BASE_DIR / 'results'
OUTPUT_DIR = BASE_DIR / 'results/subpopulation_ihs_pipeline_viz'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# All subpopulations (excluding pooled AFR, EAS, EUR which were methodologically incorrect)
POPS = ['ACB', 'ASW', 'BEB', 'CDX', 'CEU', 'CHB', 'CHS', 'CLM', 'ESN', 'FIN',
        'GBR', 'GIH', 'GWD', 'IBS', 'ITU', 'JPT', 'KHV', 'LWK', 'MSL', 'MXL',
        'PEL', 'PJL', 'PUR', 'STU', 'TSI', 'YRI']

# Population colors by superpopulation affiliation
COLORS = {
    # African
    'ACB': '#2ecc71', 'ASW': '#27ae60', 'ESN': '#229954', 'GWD': '#1e8449',
    'LWK': '#196f3d', 'MSL': '#145a32', 'YRI': '#0e4b2b',
    # East Asian
    'CDX': '#e74c3c', 'CHB': '#cb4335', 'CHS': '#b03a2e', 'JPT': '#943126', 'KHV': '#78281f',
    # European
    'CEU': '#3498db', 'FIN': '#2e86c1', 'GBR': '#2874a6', 'IBS': '#21618c', 'TSI': '#1b4f72',
    # South Asian
    'BEB': '#9b59b6', 'GIH': '#884ea0', 'ITU': '#76448a', 'PJL': '#633974', 'STU': '#512e5f',
    # Admixed American
    'CLM': '#f39c12', 'MXL': '#e67e22', 'PEL': '#d68910', 'PUR': '#ca6f1e'
}

# Full population names
NAMES = {
    'ACB': 'African Caribbean (Barbados)', 'ASW': 'African American (SW USA)',
    'ESN': 'Esan (Nigeria)', 'GWD': 'Gambian (Mandinka)', 'LWK': 'Luhya (Kenya)',
    'MSL': 'Mende (Sierra Leone)', 'YRI': 'Yoruba (Nigeria)',
    'CDX': 'Chinese Dai (Xishuangbanna)', 'CHB': 'Han Chinese (Beijing)',
    'CHS': 'Han Chinese (South)', 'JPT': 'Japanese (Tokyo)', 'KHV': 'Kinh (Vietnam)',
    'CEU': 'Utah Residents (CEPH)', 'FIN': 'Finnish', 'GBR': 'British (England/Scotland)',
    'IBS': 'Iberian (Spain)', 'TSI': 'Toscani (Italy)',
    'BEB': 'Bengali (Bangladesh)', 'GIH': 'Gujarati Indian (Houston)',
    'ITU': 'Indian Telugu (UK)', 'PJL': 'Punjabi (Lahore)', 'STU': 'Sri Lankan Tamil (UK)',
    'CLM': 'Colombian (Medellin)', 'MXL': 'Mexican (Los Angeles)',
    'PEL': 'Peruvian (Lima)', 'PUR': 'Puerto Rican'
}

print("="*80)
print("iHS Workflow Visualization")
print("="*80)
print()

def load_data(pop):
    """Load population data"""
    print(f"Loading {pop}...")
    pdir = RESULTS_DIR / f'ihs_{pop}'

    # Load ALL X chromosome SNPs to calculate proper q99
    all_snps_file = pdir / f'{pop}.chrX.ihs.tsv'
    rfile = pdir / 'analysis/regions/selection_regions_X_chromosome.tsv'
    gfile = pdir / 'analysis/gene_annotation/genes_in_regions.tsv'

    # Load all SNPs
    if all_snps_file.exists():
        all_snps = pd.read_csv(all_snps_file, sep='\t')
        # Parse Location column to get position
        location_parts = all_snps['Location'].str.split(':', expand=True)
        all_snps['chr'] = location_parts[0]
        all_snps['pos'] = location_parts[1].astype(int)
        print(f"  {len(all_snps):,} total X chromosome SNPs")

        # Calculate q99 from ALL SNPs
        q99 = np.abs(all_snps['Std iHS']).quantile(0.99)
        print(f"  q99 threshold: {q99:.3f}")

        # Filter to SNPs above q99 for visualization
        snps_above_q99 = all_snps[np.abs(all_snps['Std iHS']) >= q99].copy()
        print(f"  {len(snps_above_q99):,} SNPs ≥ q99")
    else:
        print(f"  WARNING: {all_snps_file} not found")
        snps_above_q99 = None
        q99 = None

    # Load regions and genes
    reg = pd.read_csv(rfile, sep='\t') if rfile.exists() else None
    gen = pd.read_csv(gfile, sep='\t') if gfile.exists() else None

    if reg is not None:
        print(f"  {len(reg)} regions")
    if gen is not None:
        print(f"  {len(gen):,} gene-region associations")

    return {'cand': snps_above_q99, 'reg': reg, 'gen': gen, 'q99': q99}

def plot_workflow(pop, data):
    """Plot workflow for one population"""
    print(f"Plotting {pop}...")

    if data['cand'] is None or data['reg'] is None:
        print(f"  Skipping {pop} - missing data")
        return

    fig, axes = plt.subplots(3, 1, figsize=(20, 12),
                            gridspec_kw={'height_ratios': [2, 1, 1.5]}, sharex=True)

    # Panel 1: SNPs (only those >= q99)
    ax1 = axes[0]
    q99 = data['q99']
    snps_above_q99 = data['cand']

    ax1.scatter(snps_above_q99['pos'], np.abs(snps_above_q99['Std iHS']),
               c=COLORS[pop], alpha=0.6, s=20, edgecolors='none',
               label=f'{len(snps_above_q99):,} SNPs ≥ q99')
    ax1.axhline(q99, color='red', linestyle='--', linewidth=1.5,
                label=f'q99 threshold = {q99:.3f}', alpha=0.7)
    ax1.set_ylabel('|Std.iHS|', fontsize=12, fontweight='bold')
    ax1.set_title(f'{NAMES[pop]} - iHS Analysis on X Chromosome',
                  fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.grid(alpha=0.3)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Panel 2: Regions - use vertical bars at midpoint
    ax2 = axes[1]
    for _, r in data['reg'].iterrows():
        # Calculate intensity based on max |Std.iHS|
        intensity = min(r['max_std_ihs']/8, 1)
        color = plt.cm.Reds(0.5 + 0.5*intensity)  # Use darker reds

        # Draw vertical bar at region midpoint
        midpoint = (r['start'] + r['end']) / 2
        ax2.plot([midpoint, midpoint], [0, 1],
                color=color, linewidth=2, alpha=0.9, solid_capstyle='butt')

    ax2.set_xlim(0, 155_000_000)
    ax2.set_ylim(0, 1)
    ax2.set_ylabel('Selection\nRegions', fontsize=11, fontweight='bold')
    ax2.set_yticks([])
    ax2.grid(alpha=0.3, axis='x')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)

    # Add text showing number of regions
    ax2.text(0.98, 0.95, f'{len(data["reg"])} regions',
             transform=ax2.transAxes, ha='right', va='top',
             fontsize=10, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Panel 3: Genes
    ax3 = axes[2]
    if data['gen'] is not None and len(data['gen']) > 0:
        ug = data['gen'].drop_duplicates(subset=['gene_id'])
        tracks = {}
        my = 0

        # First pass: assign genes to horizontal tracks
        for _, g in ug.iterrows():
            y = 0
            while any((y == oy and not (g['end'] < os or g['start'] > oe))
                     for os, oe, oy in tracks.values()):
                y += 1
            tracks[g['gene_name']] = (g['start'], g['end'], y)
            my = max(my, y)

        # Draw gene rectangles
        for _, g in ug.iterrows():
            y = tracks[g['gene_name']][2]
            col = '#2c3e50' if g['gene_type']=='protein_coding' else '#95a5a6'
            alph = 0.9 if g['gene_type']=='protein_coding' else 0.6
            rect = Rectangle((g['start'], y), g['end']-g['start'], 0.8,
                           facecolor=col, edgecolor='black', linewidth=0.3, alpha=alph)
            ax3.add_patch(rect)

        # Second pass: add labels for protein-coding genes with smart positioning (had help from Claude because they were often overlapping, they still do)
        protein_coding = ug[ug['gene_type'] == 'protein_coding'].copy()
        if len(protein_coding) > 0:
            # Sort by position
            protein_coding = protein_coding.sort_values('start')

            # Track label positions to avoid overlap
            label_positions = []  # (center_x, label_height)
            min_x_spacing = 8_000_000  # Minimum 8Mb spacing to avoid overlap

            for _, g in protein_coding.iterrows():
                y = tracks[g['gene_name']][2]
                gene_center = (g['start'] + g['end']) / 2

                # Determine label height - stagger if too close to previous labels
                label_height = 1.0  # Base height above gene
                for prev_x, prev_height in label_positions:
                    if abs(gene_center - prev_x) < min_x_spacing:
                        # Too close - use alternating heights
                        if prev_height == 1.0:
                            label_height = 1.8
                        else:
                            label_height = 1.0
                        break

                label_positions.append((gene_center, label_height))

                # Place label above the gene
                ax3.annotate(g['gene_name'],
                           xy=(gene_center, y + 0.8),  # Point to top of gene rectangle
                           xytext=(gene_center, y + 0.8 + label_height),  # Staggered label position
                           ha='center', va='bottom',
                           fontsize=7, fontweight='bold',
                           bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                                   edgecolor='black', linewidth=0.5, alpha=0.9),
                           arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0',
                                         lw=0.5, color='black'))

        # Add extra space above for labels (labels can go up to y+0.8+1.8, add padding)
        ax3.set_ylim(0, my+3.5)
        ax3.legend(handles=[
            mpatches.Patch(color='#2c3e50', label='Protein-coding'),
            mpatches.Patch(color='#95a5a6', label='Other', alpha=0.6)
        ], loc='upper right')
    else:
        ax3.text(0.5, 0.5, 'No gene data available', ha='center', va='center',
                transform=ax3.transAxes, fontsize=12, style='italic')
        ax3.set_ylim(0, 1)

    ax3.set_xlim(0, 155_000_000)
    ax3.set_ylabel('Genes in\nRegions', fontsize=11, fontweight='bold')
    ax3.set_yticks([])
    ax3.grid(alpha=0.3, axis='x')
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['left'].set_visible(False)
    ax3.set_xticks(np.arange(0, 160_000_000, 20_000_000))
    ax3.set_xticklabels([f'{int(x/1e6)}' for x in np.arange(0, 160_000_000, 20_000_000)])
    ax3.set_xlabel('X Chromosome Position (Mb)', fontsize=12, fontweight='bold')

    plt.tight_layout()
    outfile = OUTPUT_DIR / f'{pop}_ihs_workflow.png'
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outfile.name}")


# Main execution
all_data = {}
print("\nLoading data for all populations...")
for p in POPS:
    all_data[p] = load_data(p)


successful_plots = []
for p in POPS:
    plot_workflow(p, all_data[p])
    if all_data[p]['cand'] is not None and all_data[p]['reg'] is not None:
        successful_plots.append(p)

print()
print("="*80)
print("Complete")
print("="*80)
print(f"\nAll outputs saved to: {OUTPUT_DIR}")
print(f"\nSuccessfully generated {len(successful_plots)} workflow visualizations:")
for pop in successful_plots:
    print(f"  - {pop}_ihs_workflow.png ({NAMES[pop]})")
if len(successful_plots) < len(POPS):
    skipped = set(POPS) - set(successful_plots)
    print(f"\nSkipped {len(skipped)} populations due to missing data: {', '.join(sorted(skipped))}")
