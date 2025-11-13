#!/usr/bin/env python3
"""
Clean visualization of high-significance SNPs on X chromosome
Addresses label overlap and figure clarity issues
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import patches
from adjustText import adjust_text  # For smart label positioning
import warnings
from pathlib import Path
warnings.filterwarnings('ignore')

# Set style
sns.set_style('white')
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

def load_data():
    """Load high-significance SNPs data"""
    snp_file = '/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/analysis/high_significance_genes/snp_to_nearest_gene_only.tsv'
    df = pd.read_csv(snp_file, sep='\t')

    # Use 'pos' column for position (snp_location is the variant ID)
    # Keep snp_location for reference, but use pos for plotting
    df['snp_pos'] = pd.to_numeric(df['pos'], errors='coerce')

    # Calculate -log10(p-value)
    df['neg_log_p'] = df['neg_log10_p']

    # Sort by significance
    df = df.sort_values('neg_log_p', ascending=False)

    return df

def create_clean_manhattan(df, output_dir):
    """
    Create clean Manhattan plot with non-overlapping labels
    """
    # Create figure with more height for labels
    fig, ax = plt.subplots(figsize=(16, 8))
    
    # Get X chromosome length
    x_length = 156040895  # GRCh38 X chromosome length
    
    # Normalize colors by p-value (darker red = more significant)
    norm = plt.Normalize(vmin=df['neg_log_p'].min(), vmax=df['neg_log_p'].max())
    colors = plt.cm.Reds(norm(df['neg_log_p']))
    
    # Plot all SNPs
    scatter = ax.scatter(
        df['snp_pos'],
        df['neg_log_p'],
        c=colors,
        s=df['neg_log_p'] * 10,  # Size proportional to significance
        alpha=0.7,
        edgecolors='black',
        linewidths=0.5,
        zorder=2
    )
    
    # Add threshold line
    threshold = 6
    ax.axhline(
        y=threshold, 
        color='red', 
        linestyle='--', 
        linewidth=1.5,
        alpha=0.5,
        label=f'Threshold: p = 10⁻⁶',
        zorder=1
    )
    
    # Add centromere region (approximately)
    centromere_start = 58100000
    centromere_end = 63800000
    ax.axvspan(centromere_start, centromere_end, 
               alpha=0.1, color='gray', 
               label='Centromere region',
               zorder=0)
    
    # Label top 10 SNPs with smart positioning
    top_snps = df.nlargest(10, 'neg_log_p')
    
    texts = []
    for idx, row in top_snps.iterrows():
        # Create text label
        label = row['nearest_gene']
        x = row['snp_pos']
        y = row['neg_log_p']
        
        # Add connecting line from point to label
        text = ax.text(
            x, y, 
            label,
            fontsize=9,
            weight='bold',
            bbox=dict(
                boxstyle='round,pad=0.3',
                facecolor='white',
                edgecolor='black',
                linewidth=0.5,
                alpha=0.9
            ),
            zorder=3
        )
        texts.append(text)
    
    # Use adjustText to prevent overlap
    # This is the key to solving your overlap problem!
    from adjustText import adjust_text
    adjust_text(
        texts,
        arrowprops=dict(
            arrowstyle='-',
            color='black',
            lw=0.5,
            alpha=0.7
        ),
        expand_points=(1.5, 1.5),
        expand_text=(1.2, 1.2),
        force_points=(0.5, 0.5),
        force_text=(0.5, 0.5),
        ax=ax
    )
    
    # Styling
    ax.set_xlabel('Position on X Chromosome (Mb)', fontsize=12, weight='bold')
    ax.set_ylabel('-log₁₀(p-value)', fontsize=12, weight='bold')
    ax.set_title(
        'Distribution of 43 High-Significance SNPs on X Chromosome\n'
        f'(Total SNPs: {len(df)}, Threshold: -log₁₀(p) ≥ {threshold})',
        fontsize=14,
        weight='bold',
        pad=20
    )
    
    # Set x-axis to show full chromosome with proper scaling
    ax.set_xlim(0, x_length)
    ax.set_xticks(np.arange(0, x_length, 20_000_000))
    ax.set_xticklabels([f'{x/1e6:.0f}' for x in ax.get_xticks()])
    
    # Make y-axis start from threshold
    ax.set_ylim(threshold - 0.5, df['neg_log_p'].max() + 1)
    
    # Add legend
    ax.legend(
        loc='upper right',
        frameon=True,
        fancybox=True,
        shadow=True,
        fontsize=10
    )
    
    # Add colorbar for significance
    sm = plt.cm.ScalarMappable(cmap=plt.cm.Reds, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, pad=0.02, aspect=30)
    cbar.set_label('-log₁₀(p-value)', rotation=270, labelpad=20, weight='bold')
    
    # Grid for readability
    ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)
    
    # Remove top and right spines
    sns.despine(ax=ax)
    
    plt.tight_layout()
    
    # Save
    output_path = f'{output_dir}/highsig_snps_clean_manhattan.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"✓ Saved clean Manhattan plot: {output_path}")
    
    plt.close()

def create_density_plot(df, output_dir):
    """
    Create a density/coverage plot showing clustering of SNPs
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 8), 
                                     height_ratios=[1, 3],
                                     sharex=True)
    
    x_length = 156040895
    
    # Top panel: Density histogram
    bins = np.linspace(0, x_length, 100)
    ax1.hist(df['snp_pos'], bins=bins, color='darkred', alpha=0.7, edgecolor='black')
    ax1.set_ylabel('SNP Count', fontsize=10, weight='bold')
    ax1.set_title('SNP Density Distribution', fontsize=12, weight='bold')
    ax1.grid(True, alpha=0.3, linestyle=':', axis='y')
    sns.despine(ax=ax1)
    
    # Bottom panel: Scatter with gene context
    colors = plt.cm.Reds(plt.Normalize(
        vmin=df['neg_log_p'].min(), 
        vmax=df['neg_log_p'].max()
    )(df['neg_log_p']))
    
    scatter = ax2.scatter(
        df['snp_pos'],
        df['neg_log_p'],
        c=colors,
        s=df['neg_log_p'] * 15,
        alpha=0.7,
        edgecolors='black',
        linewidths=0.5
    )
    
    # Add threshold
    ax2.axhline(y=6, color='red', linestyle='--', linewidth=1.5, alpha=0.5)
    
    # Centromere
    centromere_start = 58100000
    centromere_end = 63800000
    ax2.axvspan(centromere_start, centromere_end, alpha=0.1, color='gray')
    
    ax2.set_xlabel('Position on X Chromosome (Mb)', fontsize=12, weight='bold')
    ax2.set_ylabel('-log₁₀(p-value)', fontsize=12, weight='bold')
    ax2.set_xlim(0, x_length)
    ax2.set_xticks(np.arange(0, x_length, 20_000_000))
    ax2.set_xticklabels([f'{x/1e6:.0f}' for x in ax2.get_xticks()])
    ax2.grid(True, alpha=0.3, linestyle=':')
    sns.despine(ax=ax2)
    
    plt.tight_layout()
    
    output_path = f'{output_dir}/highsig_snps_density_plot.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"✓ Saved density plot: {output_path}")
    
    plt.close()

def create_regional_plots(df, output_dir):
    """
    Create zoomed-in plots for regions with clustered SNPs
    """
    # Find regions with multiple SNPs (within 5 Mb windows)
    window_size = 5_000_000
    df_sorted = df.sort_values('snp_pos')

    regions = []
    i = 0
    while i < len(df_sorted):
        start_pos = df_sorted.iloc[i]['snp_pos']
        # Find SNPs within window
        window_snps = df_sorted[
            (df_sorted['snp_pos'] >= start_pos) &
            (df_sorted['snp_pos'] < start_pos + window_size)
        ]
        
        if len(window_snps) >= 3:  # At least 3 SNPs in region
            regions.append({
                'start': start_pos,
                'end': start_pos + window_size,
                'n_snps': len(window_snps),
                'snps': window_snps
            })
            i += len(window_snps)
        else:
            i += 1
    
    print(f"\n✓ Found {len(regions)} regions with ≥3 clustered SNPs")
    
    # Create plot for each region
    for idx, region in enumerate(regions[:5]):  # Top 5 regions
        fig, ax = plt.subplots(figsize=(12, 6))
        
        snps = region['snps']
        
        # Plot SNPs
        colors = plt.cm.Reds(plt.Normalize(
            vmin=snps['neg_log_p'].min(),
            vmax=snps['neg_log_p'].max()
        )(snps['neg_log_p']))
        
        scatter = ax.scatter(
            snps['snp_pos'],
            snps['neg_log_p'],
            c=colors,
            s=snps['neg_log_p'] * 30,
            alpha=0.7,
            edgecolors='black',
            linewidths=1
        )

        # Label all SNPs in this region
        texts = []
        for _, row in snps.iterrows():
            text = ax.text(
                row['snp_pos'],
                row['neg_log_p'],
                row['nearest_gene'],
                fontsize=9,
                weight='bold',
                bbox=dict(
                    boxstyle='round,pad=0.3',
                    facecolor='white',
                    edgecolor='black',
                    linewidth=0.5,
                    alpha=0.9
                )
            )
            texts.append(text)
        
        # Adjust labels
        from adjustText import adjust_text
        adjust_text(
            texts,
            arrowprops=dict(arrowstyle='-', color='black', lw=0.5),
            ax=ax
        )
        
        ax.axhline(y=6, color='red', linestyle='--', linewidth=1, alpha=0.5)
        ax.set_xlabel('Position (Mb)', fontsize=11, weight='bold')
        ax.set_ylabel('-log₁₀(p-value)', fontsize=11, weight='bold')
        ax.set_title(
            f'Regional View: {region["start"]/1e6:.1f}-{region["end"]/1e6:.1f} Mb '
            f'({region["n_snps"]} SNPs)',
            fontsize=13,
            weight='bold'
        )
        
        # Format x-axis
        ax.set_xlim(region['start'] - 500_000, region['end'] + 500_000)
        ax.ticklabel_format(style='plain', axis='x')
        ax.set_xticklabels([f'{x/1e6:.1f}' for x in ax.get_xticks()])
        
        ax.grid(True, alpha=0.3, linestyle=':')
        sns.despine(ax=ax)
        
        plt.tight_layout()
        
        output_path = f'{output_dir}/region_{idx+1}_zoom.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"  ✓ Saved regional plot {idx+1}: {output_path}")
        
        plt.close()

def main():
    """Main execution"""
    print("\n" + "="*60)
    print("Creating Clean High-Significance SNP Visualizations")
    print("="*60 + "\n")
    
    # Output directory
    output_dir = '/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/analysis/high_significance_genes/figures'

    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Load data
    print("Loading data...")
    df = load_data()
    print(f"✓ Loaded {len(df)} high-significance SNPs\n")
    
    # Create visualizations
    print("Creating visualizations...")
    print("\n1. Clean Manhattan plot with smart label positioning...")
    create_clean_manhattan(df, output_dir)
    
    print("\n2. Density plot showing SNP clustering...")
    create_density_plot(df, output_dir)
    
    print("\n3. Regional zoom plots for clustered SNPs...")
    create_regional_plots(df, output_dir)
    
    # Summary statistics
    print("\n" + "="*60)
    print("Summary Statistics")
    print("="*60)
    print(f"Total high-significance SNPs: {len(df)}")
    print(f"Range: {df['snp_pos'].min()/1e6:.1f} - {df['snp_pos'].max()/1e6:.1f} Mb")
    print(f"Max -log10(p): {df['neg_log_p'].max():.2f}")
    print(f"Mean -log10(p): {df['neg_log_p'].mean():.2f}")
    print(f"Unique genes: {df['nearest_gene'].nunique()}")
    print("\n" + "="*60)
    print("✓ All visualizations complete!")
    print("="*60 + "\n")

if __name__ == '__main__':
    main()