#!/usr/bin/env python3
"""
Visualization using geneinfo package 
Shows high-significance SNPs with gene context
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
from adjustText import adjust_text
sns.set_style('white')

# High-quality figures
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

def load_snp_data():
    """Load SNP data"""
    snp_file = '/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/ihs_revised_HCsnps_logic/analysis/high_significance_genes/snp_to_nearest_gene_only.tsv'
    df = pd.read_csv(snp_file, sep='\t')

    # Use 'pos' column for position (snp_location is the variant ID)
    df['snp_pos'] = pd.to_numeric(df['pos'], errors='coerce')

    # Use existing neg_log10_p column
    df['neg_log_p'] = df['neg_log10_p']

    return df

def create_full_chromosome_view():
    """
    Create full X chromosome view with all SNPs
    """
    from geneinfo.plot import ChromIdeogram
    
    df = load_snp_data()
    
    # Create ideogram
    g = ChromIdeogram('chrX', assembly='hg38', font_size=8, figsize=(16, 10))
    
    # Draw chromosome at specific position
    g.draw_chromosomes(base=3, height=1, facecolor='lightgray')
    
    # Get top 10 genes for labeling
    top_snps = df.nlargest(10, 'neg_log_p')
    gene_coords = [(row['snp_pos'], row['nearest_gene']) 
                   for _, row in top_snps.iterrows()]
    
    # Prepare highlight dict (color genes by significance)
    highlight = defaultdict(dict)
    top_3 = top_snps.head(3)
    for _, row in top_3.iterrows():
        gene = row['nearest_gene']
        highlight[gene].update(dict(
            weight='bold',
            color='red',
            bbox=dict(edgecolor='red', facecolor='none', linewidth=1)
        ))
    
    # Add gene labels
    label_data = [(row['snp_pos'], row['nearest_gene'], 'red' if row['nearest_gene'] in top_3['nearest_gene'].values else 'black')
                  for _, row in top_snps.iterrows()]
    
    g.add_labels(label_data, base=g.ideogram_base, min_height=g.ideogram_height*2)
    
    # Create scatter plot axis above chromosome
    ax1, ax2 = g.add_axes(2, height_ratio=0.8, hspace=0.3)
    
    # Plot 1: All SNPs colored by significance
    colors = plt.cm.Reds(plt.Normalize(
        vmin=df['neg_log_p'].min(),
        vmax=df['neg_log_p'].max()
    )(df['neg_log_p']))
    
    ax1.scatter(
        df['snp_pos'],
        df['neg_log_p'],
        c=colors,
        s=df['neg_log_p'] * 8,
        alpha=0.7,
        edgecolors='black',
        linewidths=0.3
    )
    ax1.axhline(y=6, color='red', linestyle='--', linewidth=1, alpha=0.5, label='p = 10⁻⁶')
    ax1.set_ylabel('-log₁₀(p-value)', fontsize=10, weight='bold')
    ax1.set_title('High-Significance SNPs on X Chromosome', fontsize=12, weight='bold')
    ax1.legend(loc='upper right', frameon=True, fontsize=9)
    ax1.grid(True, alpha=0.3, linestyle=':')
    
    # Plot 2: Density/coverage
    bins = 100
    ax2.hist(df['snp_pos'], bins=bins, color='darkred', alpha=0.6, edgecolor='black')
    ax2.set_ylabel('SNP Count', fontsize=10, weight='bold')
    ax2.set_title('SNP Distribution Density', fontsize=11, weight='bold')
    ax2.grid(True, alpha=0.3, linestyle=':', axis='y')
    
    # Clean up
    sns.despine(ax=ax1)
    sns.despine(ax=ax2)
    
    plt.tight_layout()
    
    output_path = '/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/ihs_revised_HCsnps_logic/analysis/high_significance_genes/figures/geneinfo_full_chromosome.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"✓ Saved full chromosome view: {output_path}")
    plt.close()

def create_regional_views():
    """
    Create zoomed regional views for clusters
    """
    from geneinfo.plot import gene_plot
    
    df = load_snp_data()
    
    # Define interesting regions (manually identified clusters)
    regions = [
        {'chr': 'chrX', 'start': 1_000_000, 'end': 10_000_000, 'name': 'PAR1_region'},
        {'chr': 'chrX', 'start': 55_000_000, 'end': 65_000_000, 'name': 'Centromere_flanking'},
        {'chr': 'chrX', 'start': 80_000_000, 'end': 90_000_000, 'name': 'Mid_chromosome'},
        {'chr': 'chrX', 'start': 95_000_000, 'end': 105_000_000, 'name': 'DIAPH2_region'},
        {'chr': 'chrX', 'start': 115_000_000, 'end': 125_000_000, 'name': 'Telomeric_region'},
    ]
    
    for region in regions:
        # Get SNPs in this region
        region_snps = df[
            (df['snp_pos'] >= region['start']) &
            (df['snp_pos'] <= region['end'])
        ]
        
        if len(region_snps) == 0:
            continue
        
        print(f"\nCreating plot for {region['name']} ({len(region_snps)} SNPs)...")
        
        # Get genes to highlight
        highlight = defaultdict(dict)
        for _, row in region_snps.iterrows():
            gene = row['nearest_gene']
            sig = row['neg_log_p']
            
            if sig > 7.5:  # Very high significance
                highlight[gene].update(dict(
                    weight='bold',
                    color='red',
                    bbox=dict(edgecolor='red', facecolor='yellow', alpha=0.3, linewidth=1)
                ))
            elif sig > 7:  # High significance
                highlight[gene].update(dict(
                    weight='bold',
                    color='red'
                ))
            else:  # Moderate significance
                highlight[gene].update(dict(color='darkred'))
        
        # Create gene plot
        ax = gene_plot(
            region['chr'],
            region['start'],
            region['end'],
            assembly='hg38',
            highlight=highlight,
            figsize=(12, 5),
            exact_exons=False
        )
        
        # Overlay SNP data
        colors = plt.cm.Reds(plt.Normalize(
            vmin=region_snps['neg_log_p'].min(),
            vmax=region_snps['neg_log_p'].max()
        )(region_snps['neg_log_p']))
        
        ax.scatter(
            region_snps['snp_pos'],
            region_snps['neg_log_p'],
            c=colors,
            s=region_snps['neg_log_p'] * 15,
            alpha=0.7,
            edgecolors='black',
            linewidths=0.5,
            zorder=3
        )
        
        ax.axhline(y=6, color='red', linestyle='--', linewidth=1, alpha=0.5)
        ax.set_ylabel('-log₁₀(p-value)', fontsize=10, weight='bold')
        ax.set_title(
            f'{region["name"]} Region: {region["start"]/1e6:.1f}-{region["end"]/1e6:.1f} Mb '
            f'({len(region_snps)} SNPs)',
            fontsize=11,
            weight='bold'
        )
        ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)
        
        sns.despine(ax=ax)
        plt.tight_layout()
        
        output_path = f'/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/ihs_revised_HCsnps_logic/analysis/high_significance_genes/figures/geneinfo_{region["name"]}.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"  ✓ Saved: {output_path}")
        plt.close()

def create_simple_manhattan():
    """
    Simple Manhattan plot without chromosome bar - clean version with gene annotations
    """
    df = load_snp_data()

    # Create regular matplotlib figure (no ChromIdeogram)
    fig, ax = plt.subplots(figsize=(16, 8))

    # Color by significance
    colors = plt.cm.Reds(plt.Normalize(
        vmin=6,  # Start from threshold
        vmax=df['neg_log_p'].max()
    )(df['neg_log_p']))

    # Main scatter - plot all SNPs
    scatter = ax.scatter(
        df['snp_pos'] / 1e6,  # Convert to Mb
        df['neg_log_p'],
        c=colors,
        s=df['neg_log_p'] * 12,
        alpha=0.7,
        edgecolors='black',
        linewidths=0.5
    )

    # Threshold line
    ax.axhline(y=6, color='red', linestyle='--', linewidth=1.5, alpha=0.5,
               label='Threshold: p = 10⁻⁶')

    # Centromere region
    ax.axvspan(58.1, 63.8, alpha=0.1, color='gray', label='Centromere')

    # Annotate ALL SNPs with their gene names
    # Group SNPs by gene to avoid duplicate labels
    gene_groups = {}
    for _, row in df.iterrows():
        gene = row['nearest_gene']
        if gene not in gene_groups:
            gene_groups[gene] = []
        gene_groups[gene].append(row)

    # For each unique gene, create text labels (adjustText will position them)
    texts = []
    for gene, snps in gene_groups.items():
        # Use the most significant SNP's position for the label placement
        most_sig_snp = max(snps, key=lambda x: x['neg_log_p'])

        # Create text label with smaller font for clarity
        text = ax.text(
            most_sig_snp['snp_pos']/1e6,
            most_sig_snp['neg_log_p'],
            gene,
            fontsize=7,
            weight='bold' if most_sig_snp['neg_log_p'] > 7 else 'normal',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                     edgecolor='red' if most_sig_snp['neg_log_p'] > 7 else 'black',
                     linewidth=0.8 if most_sig_snp['neg_log_p'] > 7 else 0.4,
                     alpha=0.9),
            zorder=3
        )
        texts.append(text)

    # Use adjustText to prevent label overlap
    print(f"Adjusting {len(texts)} gene labels...")
    adjust_text(
        texts,
        arrowprops=dict(
            arrowstyle='-',
            color='black',
            lw=0.4,
            alpha=0.6,
            shrinkA=0,
            shrinkB=3
        ),
        expand_points=(1.5, 2.0),
        expand_text=(1.2, 1.5),
        force_points=(0.5, 0.6),
        force_text=(0.8, 0.9),
        min_arrow_len=10,
        only_move={'points': 'xy', 'text': 'xy'},
        avoid_self=True,
        ax=ax,
        time_lim=3.0
    )

    # Add additional arrows for genes with multiple SNPs
    for gene, snps in gene_groups.items():
        if len(snps) > 1:
            most_sig_snp = max(snps, key=lambda x: x['neg_log_p'])
            for snp in snps:
                if snp['snp_pos'] != most_sig_snp['snp_pos']:
                    ax.annotate(
                        '',
                        xy=(snp['snp_pos']/1e6, snp['neg_log_p']),
                        xytext=(most_sig_snp['snp_pos']/1e6, most_sig_snp['neg_log_p']),
                        textcoords='data',
                        arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.2',
                                      color='gray', lw=0.4, alpha=0.4)
                    )

    ax.set_xlabel('Position on X Chromosome (Mb)', fontsize=12, weight='bold')
    ax.set_ylabel('-log₁₀(p-value)', fontsize=12, weight='bold')
    ax.set_title(
        f'Distribution of {len(df)} High-Significance SNPs on X Chromosome\n'
        f'All {len(gene_groups)} unique genes labeled',
        fontsize=13,
        weight='bold',
        pad=15
    )

    ax.set_xlim(0, 156)
    ax.set_ylim(5.8, df['neg_log_p'].max() + 0.5)

    ax.legend(loc='upper right', frameon=True, fontsize=10)
    ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)

    # Add colorbar for significance
    sm = plt.cm.ScalarMappable(
        cmap=plt.cm.Reds,
        norm=plt.Normalize(vmin=6, vmax=df['neg_log_p'].max())
    )
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, pad=0.02, aspect=30, shrink=0.8)
    cbar.set_label('-log₁₀(p-value)', rotation=270, labelpad=20, weight='bold', fontsize=10)

    sns.despine(ax=ax)
    plt.tight_layout()

    output_path = '/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/ihs_revised_HCsnps_logic/analysis/high_significance_genes/figures/simple_manhattan_clean.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"✓ Saved simple Manhattan: {output_path}")
    plt.close()

def main():
    """Main execution"""
    print("\n" + "="*60)
    print("Creating Geneinfo-Style Visualizations")
    print("="*60 + "\n")
    
    # Create output directory
    import os
    output_dir = '/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/ihs_revised_HCsnps_logic/analysis/high_significance_genes/figures'
    os.makedirs(output_dir, exist_ok=True)
    
    print("1. Creating simple clean Manhattan plot...")
    create_simple_manhattan()
    
    print("\n2. Creating full chromosome view with gene context...")
    try:
        create_full_chromosome_view()
    except Exception as e:
        print(f"  Note: ChromIdeogram failed ({e}), skipping...")
    
    print("\n3. Creating regional zoom plots with gene annotations...")
    try:
        create_regional_views()
    except Exception as e:
        print(f"  Note: Gene plots failed ({e}), skipping...")
    
    print("\n" + "="*60)
    print("✓ Visualization complete!")
    print("="*60 + "\n")

if __name__ == '__main__':
    main()