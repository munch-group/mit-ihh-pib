#!/usr/bin/env python3
"""
iHS Analysis Workflow Visualization

Visualizes the complete iHS analysis pipeline for EUR, EAS, and AFR populations:
1. SNPs with |Std.iHS| > q99 along the X chromosome
2. Selection regions defined by clustering significant SNPs
3. Genes annotated within selection regions

Creates multi-panel figures showing the analysis workflow for each population.

Author: Mit Van Bruggen
Date: 2025-12-10
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
OUTPUT_DIR = BASE_DIR / 'results/iHSperpop_genes_overlap/ihs_workflow_viz'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

POPS = ['EUR', 'EAS', 'AFR']
COLORS = {'EUR': '#3498db', 'EAS': '#e74c3c', 'AFR': '#2ecc71'}
NAMES = {'EUR': 'European', 'EAS': 'East Asian', 'AFR': 'African'}

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
    rfile = pdir / 'analysis/regions/selection_regions_20kb_min20snps.tsv'
    gfile = pdir / 'analysis/gene_annotation/genes_in_selection_regions_detailed.tsv'

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
            while any((y == oy and not (g['gene_end'] < os or g['gene_start'] > oe))
                     for os, oe, oy in tracks.values()):
                y += 1
            tracks[g['gene_name']] = (g['gene_start'], g['gene_end'], y)
            my = max(my, y)

        # Draw gene rectangles
        for _, g in ug.iterrows():
            y = tracks[g['gene_name']][2]
            col = '#2c3e50' if g['gene_type']=='protein_coding' else '#95a5a6'
            alph = 0.9 if g['gene_type']=='protein_coding' else 0.6
            rect = Rectangle((g['gene_start'], y), g['gene_end']-g['gene_start'], 0.8,
                           facecolor=col, edgecolor='black', linewidth=0.3, alpha=alph)
            ax3.add_patch(rect)

        # Second pass: add labels for protein-coding genes with smart positioning
        protein_coding = ug[ug['gene_type'] == 'protein_coding'].copy()
        if len(protein_coding) > 0:
            # Sort by position
            protein_coding = protein_coding.sort_values('gene_start')

            # Track label positions to avoid overlap
            label_positions = []  # (center_x, label_height)
            min_x_spacing = 8_000_000  # Minimum 8Mb spacing to avoid overlap

            for _, g in protein_coding.iterrows():
                y = tracks[g['gene_name']][2]
                gene_center = (g['gene_start'] + g['gene_end']) / 2

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

def plot_comparison(all_data):
    """Compare all populations"""
    print("Creating population comparison...")

    fig, axes = plt.subplots(3, 1, figsize=(20, 10), sharex=True)

    for i, pop in enumerate(POPS):
        ax = axes[i]
        reg = all_data[pop]['reg']
        if reg is not None:
            for _, r in reg.iterrows():
                intensity = min(r['max_std_ihs']/8, 1)
                # Use population-specific color with intensity-based darkness
                color_map = plt.cm.Reds if pop == 'EUR' else (plt.cm.Blues if pop == 'EAS' else plt.cm.Greens)
                color = color_map(0.5 + 0.5*intensity)

                # Draw vertical bar at region midpoint
                midpoint = (r['start'] + r['end']) / 2
                ax.plot([midpoint, midpoint], [0, 1],
                       color=color, linewidth=2, alpha=0.9, solid_capstyle='butt')
            ax.set_ylabel(f'{NAMES[pop]}\n({len(reg)} regions)',
                         fontsize=11, fontweight='bold')
        ax.set_xlim(0, 155_000_000)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.grid(alpha=0.3, axis='x')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)

    axes[2].set_xlabel('X Chromosome Position (Mb)', fontsize=12, fontweight='bold')
    axes[2].set_xticks(np.arange(0, 160_000_000, 20_000_000))
    axes[2].set_xticklabels([f'{int(x/1e6)}' for x in np.arange(0, 160_000_000, 20_000_000)])

    plt.suptitle('Selection Regions Comparison Across Populations',
                fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout()
    outfile = OUTPUT_DIR / 'population_comparison.png'
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outfile.name}")

# Main execution
all_data = {}
for p in POPS:
    all_data[p] = load_data(p)

print()
print("="*80)
print("Creating Visualizations")
print("="*80)
print()

for p in POPS:
    plot_workflow(p, all_data[p])

print()
plot_comparison(all_data)

print()
print("="*80)
print("Complete!")
print("="*80)
print(f"\nAll outputs saved to: {OUTPUT_DIR}")
print("\nGenerated files:")
print("  - EUR_ihs_workflow.png")
print("  - EAS_ihs_workflow.png")
print("  - AFR_ihs_workflow.png")
print("  - population_comparison.png")
