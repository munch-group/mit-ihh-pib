#!/usr/bin/env python3
"""
Phase 4b: Enrichment Results Visualization
==========================================

Creates publication-quality figures for enrichment analysis results.

Author: MIT van Bruggen
Date: October 24, 2025
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set style for publication-quality figures
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'sans-serif'


class EnrichmentVisualizer:
    """Creates visualizations for enrichment analysis results."""

    def __init__(self, results_dir: Path, output_dir: Path = None):
        """
        Initialize the visualizer.

        Args:
            results_dir: Directory containing enrichment results
            output_dir: Directory for output figures (default: results_dir/figures)
        """
        self.results_dir = Path(results_dir)
        self.output_dir = Path(output_dir) if output_dir else self.results_dir / "figures"
        self.output_dir.mkdir(parents=True, exist_ok=True)

        print(f"✓ Initialized enrichment visualizer")
        print(f"  Input directory: {self.results_dir}")
        print(f"  Output directory: {self.output_dir}")

    def load_results(self, filename: str) -> pd.DataFrame:
        """Load enrichment results file."""
        filepath = self.results_dir / filename
        if not filepath.exists():
            print(f"  Warning: {filename} not found")
            return pd.DataFrame()

        df = pd.read_csv(filepath, sep='\t')
        print(f"  Loaded {len(df)} terms from {filename}")
        return df

    def plot_top_terms_barplot(self,
                               results: pd.DataFrame,
                               top_n: int = 20,
                               title: str = "Top Enriched Terms",
                               filename: str = "top_terms_barplot.png",
                               color: str = "#3498db"):
        """
        Create horizontal barplot of top enriched terms.

        Args:
            results: Enrichment results DataFrame
            top_n: Number of top terms to show
            title: Plot title
            filename: Output filename
            color: Bar color
        """
        if results.empty:
            print(f"  Skipping {filename}: no results")
            return

        # Get top N terms by p-value
        plot_data = results.nsmallest(top_n, 'p_value').copy()

        # Calculate -log10(p-value) for visualization
        plot_data['neg_log10_pvalue'] = -np.log10(plot_data['p_value'])

        # Truncate long term names
        plot_data['short_name'] = plot_data['name'].apply(
            lambda x: x[:60] + '...' if len(x) > 60 else x
        )

        # Create figure
        fig, ax = plt.subplots(figsize=(10, max(6, top_n * 0.3)))

        # Create barplot
        y_pos = np.arange(len(plot_data))
        ax.barh(y_pos, plot_data['neg_log10_pvalue'], color=color, alpha=0.8)

        # Customize
        ax.set_yticks(y_pos)
        ax.set_yticklabels(plot_data['short_name'][::-1])
        ax.invert_yaxis()
        ax.set_xlabel('-log₁₀(p-value)', fontsize=12)
        ax.set_title(title, fontsize=14, fontweight='bold', pad=20)

        # Add significance threshold line
        ax.axvline(x=-np.log10(0.05), color='red', linestyle='--',
                  linewidth=1, alpha=0.5, label='p=0.05')
        ax.legend()

        # Add grid
        ax.grid(axis='x', alpha=0.3)

        plt.tight_layout()

        # Save
        output_path = self.output_dir / filename
        plt.savefig(output_path, bbox_inches='tight')
        print(f"✓ Saved {output_path}")
        plt.close()

    def plot_comparison_scatter(self,
                               comparison_df: pd.DataFrame,
                               title: str = "X Chromosome vs Autosomes Enrichment",
                               filename: str = "x_vs_autosome_scatter.png"):
        """
        Create scatter plot comparing X vs autosome enrichment.

        Args:
            comparison_df: Comparison results DataFrame
            title: Plot title
            filename: Output filename
        """
        if comparison_df.empty:
            print(f"  Skipping {filename}: no comparison data")
            return

        fig, ax = plt.subplots(figsize=(10, 10))

        # Calculate -log10 p-values
        x_scores = -np.log10(comparison_df['x_pvalue'])
        auto_scores = -np.log10(comparison_df['autosome_pvalue'])

        # Create scatter plot
        scatter = ax.scatter(auto_scores, x_scores,
                           alpha=0.6, s=50, c='#3498db', edgecolors='black', linewidths=0.5)

        # Add diagonal line (equal enrichment)
        max_val = max(x_scores.max(), auto_scores.max())
        ax.plot([0, max_val], [0, max_val], 'r--', alpha=0.5, linewidth=2, label='Equal enrichment')

        # Add significance thresholds
        sig_threshold = -np.log10(0.05)
        ax.axhline(y=sig_threshold, color='gray', linestyle=':', alpha=0.5)
        ax.axvline(x=sig_threshold, color='gray', linestyle=':', alpha=0.5)

        # Label quadrants
        ax.text(0.95, 0.05, 'Autosome-specific', transform=ax.transAxes,
               ha='right', va='bottom', fontsize=10, alpha=0.6)
        ax.text(0.05, 0.95, 'X-specific', transform=ax.transAxes,
               ha='left', va='top', fontsize=10, alpha=0.6)

        # Customize
        ax.set_xlabel('Autosome Enrichment (-log₁₀ p-value)', fontsize=12)
        ax.set_ylabel('X Chromosome Enrichment (-log₁₀ p-value)', fontsize=12)
        ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
        ax.legend()

        # Make square
        ax.set_aspect('equal', adjustable='box')

        plt.tight_layout()

        # Save
        output_path = self.output_dir / filename
        plt.savefig(output_path, bbox_inches='tight')
        print(f"✓ Saved {output_path}")
        plt.close()

    def plot_enrichment_heatmap(self,
                               results_dict: dict,
                               top_n: int = 30,
                               title: str = "Enrichment Heatmap",
                               filename: str = "enrichment_heatmap.png"):
        """
        Create heatmap showing enrichment across different gene sets.

        Args:
            results_dict: Dictionary of {analysis_name: results_df}
            top_n: Number of top terms to include
            title: Plot title
            filename: Output filename
        """
        # Collect all unique terms across analyses
        all_terms = {}

        for name, df in results_dict.items():
            if df.empty:
                continue
            for _, row in df.head(top_n).iterrows():
                term_id = row['native']
                term_name = row['name']
                if term_id not in all_terms:
                    all_terms[term_id] = {
                        'name': term_name,
                        'enrichment': {}
                    }
                all_terms[term_id]['enrichment'][name] = -np.log10(row['p_value'])

        if not all_terms:
            print(f"  Skipping {filename}: no terms to plot")
            return

        # Create matrix
        term_ids = list(all_terms.keys())
        analysis_names = list(results_dict.keys())

        matrix = np.zeros((len(term_ids), len(analysis_names)))

        for i, term_id in enumerate(term_ids):
            for j, analysis_name in enumerate(analysis_names):
                matrix[i, j] = all_terms[term_id]['enrichment'].get(analysis_name, 0)

        # Create figure
        fig, ax = plt.subplots(figsize=(12, max(8, len(term_ids) * 0.3)))

        # Create heatmap
        im = ax.imshow(matrix, cmap='YlOrRd', aspect='auto')

        # Set ticks
        ax.set_xticks(np.arange(len(analysis_names)))
        ax.set_yticks(np.arange(len(term_ids)))

        # Set labels
        ax.set_xticklabels(analysis_names, rotation=45, ha='right')
        term_labels = [all_terms[tid]['name'][:60] + '...' if len(all_terms[tid]['name']) > 60
                      else all_terms[tid]['name'] for tid in term_ids]
        ax.set_yticklabels(term_labels)

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('-log₁₀(p-value)', rotation=270, labelpad=20)

        # Title
        ax.set_title(title, fontsize=14, fontweight='bold', pad=20)

        plt.tight_layout()

        # Save
        output_path = self.output_dir / filename
        plt.savefig(output_path, bbox_inches='tight')
        print(f"✓ Saved {output_path}")
        plt.close()

    def plot_source_distribution(self,
                                results: pd.DataFrame,
                                title: str = "Enrichment by Database Source",
                                filename: str = "source_distribution.png"):
        """
        Create pie chart showing distribution of enriched terms by source.

        Args:
            results: Enrichment results DataFrame
            title: Plot title
            filename: Output filename
        """
        if results.empty or 'source' not in results.columns:
            print(f"  Skipping {filename}: no source data")
            return

        # Count by source
        source_counts = results['source'].value_counts()

        fig, ax = plt.subplots(figsize=(8, 8))

        # Create pie chart
        colors = sns.color_palette("Set2", len(source_counts))
        wedges, texts, autotexts = ax.pie(source_counts.values,
                                          labels=source_counts.index,
                                          autopct='%1.1f%%',
                                          colors=colors,
                                          startangle=90)

        # Customize
        for text in texts:
            text.set_fontsize(11)
        for autotext in autotexts:
            autotext.set_color('white')
            autotext.set_fontweight('bold')

        ax.set_title(title, fontsize=14, fontweight='bold', pad=20)

        plt.tight_layout()

        # Save
        output_path = self.output_dir / filename
        plt.savefig(output_path, bbox_inches='tight')
        print(f"✓ Saved {output_path}")
        plt.close()

    def plot_enrichment_summary(self,
                               summary_df: pd.DataFrame,
                               filename: str = "enrichment_summary.png"):
        """
        Create summary visualization showing key metrics across analyses.

        Args:
            summary_df: Summary DataFrame from enrichment analysis
            filename: Output filename
        """
        if summary_df.empty:
            print(f"  Skipping {filename}: no summary data")
            return

        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        # 1. Number of significant terms
        ax1 = axes[0, 0]
        ax1.bar(summary_df['analysis'], summary_df['total_significant_terms'],
               color='#3498db', alpha=0.8)
        ax1.set_ylabel('Number of Significant Terms', fontsize=11)
        ax1.set_title('Significant Terms by Analysis', fontweight='bold')
        ax1.tick_params(axis='x', rotation=45)
        ax1.grid(axis='y', alpha=0.3)

        # 2. Mean enrichment score
        ax2 = axes[0, 1]
        ax2.bar(summary_df['analysis'], summary_df['mean_enrichment_score'],
               color='#e74c3c', alpha=0.8)
        ax2.set_ylabel('Mean -log₁₀(p-value)', fontsize=11)
        ax2.set_title('Average Enrichment Strength', fontweight='bold')
        ax2.tick_params(axis='x', rotation=45)
        ax2.grid(axis='y', alpha=0.3)

        # 3. Minimum p-value
        ax3 = axes[1, 0]
        min_p_transformed = -np.log10(summary_df['min_pvalue'])
        ax3.bar(summary_df['analysis'], min_p_transformed,
               color='#2ecc71', alpha=0.8)
        ax3.set_ylabel('Max -log₁₀(p-value)', fontsize=11)
        ax3.set_title('Strongest Signal by Analysis', fontweight='bold')
        ax3.tick_params(axis='x', rotation=45)
        ax3.grid(axis='y', alpha=0.3)

        # 4. Top terms (text)
        ax4 = axes[1, 1]
        ax4.axis('off')
        top_terms_text = "Top Terms:\n\n"
        for _, row in summary_df.iterrows():
            term = row['top_term'][:40] + '...' if len(row['top_term']) > 40 else row['top_term']
            top_terms_text += f"{row['analysis']}:\n  {term}\n\n"
        ax4.text(0.1, 0.9, top_terms_text, transform=ax4.transAxes,
                fontsize=9, verticalalignment='top', family='monospace')

        plt.suptitle('Enrichment Analysis Summary', fontsize=16, fontweight='bold', y=0.995)
        plt.tight_layout()

        # Save
        output_path = self.output_dir / filename
        plt.savefig(output_path, bbox_inches='tight')
        print(f"✓ Saved {output_path}")
        plt.close()


def main():
    """Main visualization workflow."""

    print("\n" + "="*70)
    print("ENRICHMENT RESULTS VISUALIZATION")
    print("="*70 + "\n")

    # Paths
    base_dir = Path("/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib")
    results_dir = base_dir / "results" / "analysis" / "enrichment"

    if not results_dir.exists():
        print(f"Error: Results directory not found: {results_dir}")
        print("Please run enrichment analysis first (04_functional_enrichment_analysis.py)")
        return

    # Initialize visualizer
    viz = EnrichmentVisualizer(results_dir)

    # Load results
    print("\nLoading enrichment results...")
    x_results = viz.load_results("enrichment_X_chromosome_full.tsv")
    autosome_results = viz.load_results("enrichment_autosome_full.tsv")
    x_reproductive = viz.load_results("enrichment_x_reproductive.tsv")
    x_immune = viz.load_results("enrichment_x_immune.tsv")
    comparison_df = viz.load_results("x_vs_autosome_comparison.tsv")
    summary_df = viz.load_results("enrichment_summary.tsv")

    # Create visualizations
    print("\n" + "-"*70)
    print("Creating visualizations...")
    print("-"*70)

    # 1. Top X chromosome terms
    if not x_results.empty:
        viz.plot_top_terms_barplot(
            x_results, top_n=20,
            title="Top 20 Enriched Terms - X Chromosome Genes",
            filename="x_chromosome_top_terms.png",
            color="#9b59b6"
        )

    # 2. Top autosome terms
    if not autosome_results.empty:
        viz.plot_top_terms_barplot(
            autosome_results, top_n=20,
            title="Top 20 Enriched Terms - Autosomal Genes",
            filename="autosome_top_terms.png",
            color="#3498db"
        )

    # 3. Reproductive terms
    if not x_reproductive.empty:
        viz.plot_top_terms_barplot(
            x_reproductive, top_n=15,
            title="Reproductive Function Enrichment - X Chromosome",
            filename="x_reproductive_terms.png",
            color="#e74c3c"
        )

    # 4. Immune terms
    if not x_immune.empty:
        viz.plot_top_terms_barplot(
            x_immune, top_n=15,
            title="Immune Function Enrichment - X Chromosome",
            filename="x_immune_terms.png",
            color="#2ecc71"
        )

    # 5. X vs Autosome comparison scatter
    if not comparison_df.empty:
        viz.plot_comparison_scatter(
            comparison_df,
            title="X Chromosome vs Autosomes: Common Enriched Terms",
            filename="x_vs_autosome_scatter.png"
        )

    # 6. Source distribution
    if not x_results.empty:
        viz.plot_source_distribution(
            x_results,
            title="X Chromosome Enrichment by Database Source",
            filename="x_source_distribution.png"
        )

    # 7. Summary visualization
    if not summary_df.empty:
        viz.plot_enrichment_summary(
            summary_df,
            filename="enrichment_summary.png"
        )

    # 8. Heatmap (if we have multiple result sets)
    results_dict = {
        'X Chromosome': x_results,
        'Autosomes': autosome_results,
    }

    if not x_reproductive.empty:
        results_dict['X Reproductive'] = x_reproductive
    if not x_immune.empty:
        results_dict['X Immune'] = x_immune

    viz.plot_enrichment_heatmap(
        results_dict, top_n=25,
        title="Top Enriched Terms Across Analyses",
        filename="enrichment_heatmap.png"
    )

    print("\n" + "="*70)
    print("VISUALIZATION COMPLETE!")
    print("="*70)
    print(f"\nAll figures saved to: {viz.output_dir}")
    print("\nGenerated figures:")
    for fig_file in sorted(viz.output_dir.glob("*.png")):
        print(f"  - {fig_file.name}")
    print()


if __name__ == "__main__":
    main()
