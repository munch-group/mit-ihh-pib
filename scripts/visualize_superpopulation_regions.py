#!/usr/bin/env python3
"""
Visualize selection regions across superpopulations and their constituent subpopulations.

Creates three figures (AFR, EAS, EUR), each showing selection regions for:
- The UNION of all constituent subpopulation regions (top panel)
- Each constituent subpopulation

Author: Mit Van Bruggen
Date: 2025-12-22
"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

BASE_DIR = Path('/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib')
RESULTS_DIR = BASE_DIR / 'results'
OUTPUT_DIR = BASE_DIR / 'results/superpopulation_regions_viz'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Define superpopulations and their constituent subpopulations
# NOTE: Excluding the superpopulation codes themselves (AFR, EAS, EUR) as those
# represent pooled analyses which are methodologically incorrect
SUPERPOP_STRUCTURE = {
    'AFR': {
        'name': 'African',
        'subpops': ['ACB', 'ASW', 'ESN', 'GWD', 'LWK', 'MSL', 'YRI'],
        'color': '#2ecc71'
    },
    'EAS': {
        'name': 'East Asian',
        'subpops': ['CDX', 'CHB', 'CHS', 'JPT', 'KHV'],
        'color': '#e74c3c'
    },
    'EUR': {
        'name': 'European',
        'subpops': ['CEU', 'FIN', 'GBR', 'IBS', 'TSI'],
        'color': '#3498db'
    }
}

# Population labels (for display)
POP_LABELS = {
    'ACB': 'African Caribbean (Barbados)',
    'ASW': 'African American (SW USA)',
    'ESN': 'Esan (Nigeria)',
    'GWD': 'Gambian (Mandinka)',
    'LWK': 'Luhya (Kenya)',
    'MSL': 'Mende (Sierra Leone)',
    'YRI': 'Yoruba (Nigeria)',

    'CDX': 'Chinese Dai (Xishuangbanna)',
    'CHB': 'Han Chinese (Beijing)',
    'CHS': 'Han Chinese (South)',
    'JPT': 'Japanese (Tokyo)',
    'KHV': 'Kinh (Vietnam)',

    'CEU': 'Utah Residents (CEPH)',
    'FIN': 'Finnish',
    'GBR': 'British (England/Scotland)',
    'IBS': 'Iberian (Spain)',
    'TSI': 'Toscani (Italy)'
}

print("="*80)
print("Superpopulation Selection Regions Visualization")
print("="*80)
print()

def load_regions(pop):
    """Load selection regions for a population"""
    regions_file = RESULTS_DIR / f'ihs_{pop}' / 'analysis' / 'regions' / 'selection_regions_X_chromosome.tsv'

    if regions_file.exists():
        regions = pd.read_csv(regions_file, sep='\t')
        print(f"  {pop}: {len(regions)} regions")
        return regions
    else:
        print(f"  {pop}: No regions file found")
        return None

def create_union_regions(subpop_regions):
    """
    Create union of selection regions from multiple populations.
    Merges overlapping regions from different populations.
    """
    all_regions = []

    # Collect all regions from all subpopulations
    for pop, regions in subpop_regions.items():
        if regions is not None and len(regions) > 0:
            for _, r in regions.iterrows():
                all_regions.append({
                    'chr': r['chr'],
                    'start': r['start'],
                    'end': r['end'],
                    'mean_std_ihs': abs(r['mean_std_ihs']),
                    'source_pop': pop
                })

    if not all_regions:
        return None

    # Convert to DataFrame and sort by position
    df = pd.DataFrame(all_regions)
    df = df.sort_values('start').reset_index(drop=True)

    # Merge overlapping regions
    merged = []
    current_region = df.iloc[0].to_dict()

    for i in range(1, len(df)):
        next_region = df.iloc[i]

        # Check if regions overlap or are adjacent (within 1kb)
        if next_region['start'] <= current_region['end'] + 1000:
            # Merge: extend end, take max mean_std_ihs
            current_region['end'] = max(current_region['end'], next_region['end'])
            current_region['mean_std_ihs'] = max(current_region['mean_std_ihs'],
                                                  next_region['mean_std_ihs'])
        else:
            # Save current region and start new one
            merged.append(current_region)
            current_region = next_region.to_dict()

    # Add the last region
    merged.append(current_region)

    return pd.DataFrame(merged)

def plot_superpopulation(superpop_code, superpop_info):
    """Create visualization for one superpopulation and its subpopulations"""
    print(f"\nPlotting {superpop_info['name']}...")

    subpops = superpop_info['subpops']
    n_pops = len(subpops)

    # Load all regions
    all_regions = {}
    for pop in subpops:
        all_regions[pop] = load_regions(pop)

    # Create union of all subpopulation regions
    print(f"  Creating union of {len(subpops)} subpopulations...")
    union_regions = create_union_regions(all_regions)
    if union_regions is not None:
        print(f"  Union: {len(union_regions)} merged regions")

    # Create figure with n_pops + 1 panels (union + each subpop)
    fig, axes = plt.subplots(n_pops + 1, 1, figsize=(24, (n_pops + 1) * 1.5), sharex=True)

    # Define color gradients
    color_maps = {
        'AFR': plt.cm.Greens,
        'EAS': plt.cm.Reds,
        'EUR': plt.cm.Blues
    }
    color_map = color_maps[superpop_code]

    # Plot UNION in the first panel
    ax = axes[0]
    if union_regions is not None and len(union_regions) > 0:
        for _, r in union_regions.iterrows():
            # Calculate intensity based on mean |Std.iHS|
            intensity = min((r['mean_std_ihs'] - 3) / 5, 1)
            intensity = max(intensity, 0.3)

            color = color_map(0.4 + 0.6 * intensity)

            midpoint = (r['start'] + r['end']) / 2
            ax.plot([midpoint, midpoint], [0, 1],
                   color=color, linewidth=4, alpha=0.9, solid_capstyle='butt')

        ax.set_ylabel(f'{superpop_info["name"]}\nUNION\n({len(union_regions)} regions)',
                     fontsize=11, fontweight='bold')
    else:
        ax.set_ylabel(f'{superpop_info["name"]}\nUNION\n(0 regions)',
                     fontsize=11, fontweight='bold')
        ax.text(0.5, 0.5, 'No selection regions',
               ha='center', va='center', transform=ax.transAxes,
               fontsize=9, style='italic', color='gray')

    ax.set_xlim(0, 155_000_000)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.grid(alpha=0.3, axis='x')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_facecolor('#f8f9fa')  # Highlight union panel

    # Plot each subpopulation
    for i, pop in enumerate(subpops):
        ax = axes[i + 1]  # +1 because union is at index 0
        regions = all_regions[pop]

        if regions is not None and len(regions) > 0:
            for _, r in regions.iterrows():
                intensity = min((abs(r['mean_std_ihs']) - 3) / 5, 1)
                intensity = max(intensity, 0.3)

                color = color_map(0.4 + 0.6 * intensity)

                midpoint = (r['start'] + r['end']) / 2
                ax.plot([midpoint, midpoint], [0, 1],
                       color=color, linewidth=3, alpha=0.85, solid_capstyle='butt')

            label = POP_LABELS.get(pop, pop)
            ax.set_ylabel(f'{label}\n({len(regions)} regions)',
                         fontsize=10, fontweight='bold')
        else:
            label = POP_LABELS.get(pop, pop)
            ax.set_ylabel(f'{label}\n(0 regions)',
                         fontsize=10, fontweight='bold')
            ax.text(0.5, 0.5, 'No selection regions identified',
                   ha='center', va='center', transform=ax.transAxes,
                   fontsize=9, style='italic', color='gray')

        ax.set_xlim(0, 155_000_000)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.grid(alpha=0.3, axis='x')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)

    # X-axis labels - every 10 Mb
    axes[-1].set_xlabel('X Chromosome Position (Mb)', fontsize=12, fontweight='bold')
    axes[-1].set_xticks(np.arange(0, 160_000_000, 10_000_000))
    axes[-1].set_xticklabels([f'{int(x/1e6)}' for x in np.arange(0, 160_000_000, 10_000_000)])

    # Title
    plt.suptitle(f'{superpop_info["name"]} - Selection Regions on X Chromosome',
                fontsize=16, fontweight='bold', y=0.998)

    plt.tight_layout()

    # Save figure
    outfile = OUTPUT_DIR / f'{superpop_code}_subpopulation_regions.png'
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outfile.name}")

# Main execution
print("Loading and plotting data...")

for superpop_code, superpop_info in SUPERPOP_STRUCTURE.items():
    plot_superpopulation(superpop_code, superpop_info)

print()
print("="*80)
print("Complete!")
print("="*80)
print(f"\nAll outputs saved to: {OUTPUT_DIR}")
print("\nGenerated files:")
print("  - AFR_subpopulation_regions.png")
print("  - EAS_subpopulation_regions.png")
print("  - EUR_subpopulation_regions.png")
print()
