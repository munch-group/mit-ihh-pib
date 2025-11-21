#!/usr/bin/env python3
"""
Compare iHS and RELATE selection signals on chromosome X

This script analyzes the overlap and differences between:
- iHS results: Strong selection signals (|Std iHS| > 2, p < 1e-6)
- RELATE results: Selection signals over variant lifetime (p < 1e-6)

Author: Mit van Bruggenmit
Date: 2025-11-21
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Tuple, Dict, List
import argparse
from scipy import stats

class SelectionComparator:
    """Compare iHS and RELATE selection signals"""

    def __init__(self, ihs_file: str, relate_file: str, output_dir: str):
        """
        Initialize the comparator

        Args:
            ihs_file: Path to iHS results TSV
            relate_file: Path to RELATE results CSV
            output_dir: Directory for output files
        """
        self.ihs_file = ihs_file
        self.relate_file = relate_file
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.ihs_data = None
        self.relate_data = None
        self.overlaps = None
        self.stats = {}

    def load_data(self) -> None:
        """Load iHS and RELATE data"""
        print("Loading data...")

        # Load iHS data
        self.ihs_data = pd.read_csv(self.ihs_file, sep='\t')
        print(f"  Loaded {len(self.ihs_data)} iHS candidates")

        # Load RELATE data
        self.relate_data = pd.read_csv(
            self.relate_file,
            names=['Population', 'Chromosome', 'Position', 'P_50pct', 'P_lifetime']
        )

        # Convert RELATE log10 p-values to -log10(p) for consistency
        self.relate_data['negLog10_P_lifetime'] = -self.relate_data['P_lifetime']
        self.relate_data['negLog10_P_50pct'] = -self.relate_data['P_50pct']

        print(f"  Loaded {len(self.relate_data)} RELATE signals across {self.relate_data['Population'].nunique()} populations")
        print(f"  Populations: {', '.join(sorted(self.relate_data['Population'].unique()))}")

    def find_overlaps(self, window: int = 50000, exact_match: bool = False) -> pd.DataFrame:
        """
        Find overlapping signals between iHS and RELATE

        Args:
            window: Window size in bp for considering overlap (default: 50kb)
            exact_match: If True, only exact position matches count

        Returns:
            DataFrame of overlapping signals
        """
        print(f"\nFinding overlaps (window={window} bp, exact_match={exact_match})...")

        overlaps = []

        for idx, ihs_row in self.ihs_data.iterrows():
            ihs_pos = ihs_row['pos']

            if exact_match:
                # Exact position match
                nearby = self.relate_data[self.relate_data['Position'] == ihs_pos]
            else:
                # Window-based match
                nearby = self.relate_data[
                    (self.relate_data['Position'] >= ihs_pos - window) &
                    (self.relate_data['Position'] <= ihs_pos + window)
                ]

            if len(nearby) > 0:
                for _, relate_row in nearby.iterrows():
                    distance = abs(ihs_pos - relate_row['Position'])

                    overlaps.append({
                        # Position info
                        'ihs_pos': ihs_pos,
                        'relate_pos': relate_row['Position'],
                        'distance_bp': distance,

                        # iHS info
                        'location': ihs_row['Location'],
                        'ref': ihs_row['ref'],
                        'alt': ihs_row['alt'],
                        'ihs': ihs_row['iHS'],
                        'std_ihs': ihs_row['Std iHS'],
                        'ihs_neg_log10_p': ihs_row['neg_log10_p'],

                        # RELATE info
                        'population': relate_row['Population'],
                        'relate_neg_log10_p_lifetime': relate_row['negLog10_P_lifetime'],
                        'relate_neg_log10_p_50pct': relate_row['negLog10_P_50pct'],

                        # iHH scores for haplotype analysis
                        'iHH_0': ihs_row['iHH_0'],
                        'iHH_1': ihs_row['iHH_1'],
                    })

        self.overlaps = pd.DataFrame(overlaps)

        if len(self.overlaps) > 0:
            # Classify signal types
            self.overlaps['signal_type'] = self.overlaps.apply(self._classify_signal, axis=1)

            # Add signal strength categories
            self.overlaps['ihs_strength'] = pd.cut(
                self.overlaps['std_ihs'],
                bins=[0, 2, 3, 4, np.inf],
                labels=['Moderate', 'Strong', 'Very Strong', 'Extreme']
            )

            self.overlaps['relate_strength'] = pd.cut(
                self.overlaps['relate_neg_log10_p_lifetime'],
                bins=[0, 6, 8, 10, np.inf],
                labels=['Moderate', 'Strong', 'Very Strong', 'Extreme']
            )

        print(f"  Found {len(self.overlaps)} overlapping signals")
        if len(self.overlaps) > 0:
            print(f"  Unique iHS positions with RELATE support: {self.overlaps['ihs_pos'].nunique()}")
            print(f"  Unique RELATE positions matching iHS: {self.overlaps['relate_pos'].nunique()}")

        return self.overlaps

    def _classify_signal(self, row: pd.Series) -> str:
        """Classify signal type based on iHS and RELATE strength"""
        ihs_strong = row['std_ihs'] > 3
        ihs_moderate = row['std_ihs'] > 2
        relate_strong = row['relate_neg_log10_p_lifetime'] > 8
        relate_moderate = row['relate_neg_log10_p_lifetime'] > 6

        if ihs_strong and relate_strong:
            return 'Both_Strong'
        elif ihs_strong and relate_moderate:
            return 'iHS_Strong_RELATE_Moderate'
        elif ihs_moderate and relate_strong:
            return 'iHS_Moderate_RELATE_Strong'
        elif ihs_moderate and relate_moderate:
            return 'Both_Moderate'
        else:
            return 'Weak'

    def find_ihs_only(self, window: int = 50000) -> pd.DataFrame:
        """Find iHS signals with NO RELATE support"""
        if self.overlaps is None:
            self.find_overlaps(window=window)

        # Get iHS positions that have RELATE support
        if len(self.overlaps) > 0:
            ihs_with_relate = set(self.overlaps['ihs_pos'].unique())
        else:
            ihs_with_relate = set()

        # Find iHS-only signals
        ihs_only = self.ihs_data[~self.ihs_data['pos'].isin(ihs_with_relate)].copy()

        print(f"\nFound {len(ihs_only)} iHS-only signals (no RELATE support within {window} bp)")

        return ihs_only

    def find_relate_only(self, window: int = 50000) -> pd.DataFrame:
        """Find RELATE signals with NO iHS support"""
        if self.overlaps is None:
            self.find_overlaps(window=window)

        # Get RELATE positions that have iHS support
        if len(self.overlaps) > 0:
            relate_with_ihs = set(self.overlaps['relate_pos'].unique())
        else:
            relate_with_ihs = set()

        # Find RELATE-only signals
        relate_only = self.relate_data[~self.relate_data['Position'].isin(relate_with_ihs)].copy()

        print(f"\nFound {len(relate_only)} RELATE-only signals (no iHS support within {window} bp)")
        print("  By population:")
        for pop in sorted(relate_only['Population'].unique()):
            count = len(relate_only[relate_only['Population'] == pop])
            print(f"    {pop}: {count}")

        return relate_only

    def compute_statistics(self, window: int = 50000) -> Dict:
        """Compute summary statistics"""
        print("\n" + "="*70)
        print("SUMMARY STATISTICS")
        print("="*70)

        if self.overlaps is None:
            self.find_overlaps(window=window)

        stats_dict = {
            'n_ihs_candidates': len(self.ihs_data),
            'n_relate_signals': len(self.relate_data),
            'n_relate_positions': self.relate_data['Position'].nunique(),
            'n_populations': self.relate_data['Population'].nunique(),
            'n_overlaps': len(self.overlaps),
            'n_unique_ihs_with_relate': self.overlaps['ihs_pos'].nunique() if len(self.overlaps) > 0 else 0,
            'n_unique_relate_with_ihs': self.overlaps['relate_pos'].nunique() if len(self.overlaps) > 0 else 0,
            'pct_ihs_with_relate': (self.overlaps['ihs_pos'].nunique() / len(self.ihs_data) * 100) if len(self.overlaps) > 0 else 0,
            'pct_relate_with_ihs': (self.overlaps['relate_pos'].nunique() / self.relate_data['Position'].nunique() * 100) if len(self.overlaps) > 0 else 0,
        }

        # Add iHS-only and RELATE-only counts
        ihs_only = self.find_ihs_only(window=window)
        relate_only = self.find_relate_only(window=window)

        stats_dict['n_ihs_only'] = len(ihs_only)
        stats_dict['n_relate_only'] = len(relate_only)

        # Print summary
        print(f"\nInput data:")
        print(f"  iHS candidates: {stats_dict['n_ihs_candidates']}")
        print(f"  RELATE signals: {stats_dict['n_relate_signals']} ({stats_dict['n_relate_positions']} unique positions)")
        print(f"  Populations: {stats_dict['n_populations']}")

        print(f"\nOverlap (within {window} bp window):")
        print(f"  Total overlaps: {stats_dict['n_overlaps']}")
        print(f"  iHS positions with RELATE support: {stats_dict['n_unique_ihs_with_relate']} ({stats_dict['pct_ihs_with_relate']:.1f}%)")
        print(f"  RELATE positions with iHS support: {stats_dict['n_unique_relate_with_ihs']} ({stats_dict['pct_relate_with_ihs']:.1f}%)")

        print(f"\nMethod-specific signals:")
        print(f"  iHS-only (no RELATE): {stats_dict['n_ihs_only']}")
        print(f"  RELATE-only (no iHS): {stats_dict['n_relate_only']}")

        if len(self.overlaps) > 0:
            print(f"\nSignal type distribution (overlaps):")
            signal_counts = self.overlaps['signal_type'].value_counts()
            for signal_type, count in signal_counts.items():
                print(f"  {signal_type}: {count}")

            # Correlation between iHS and RELATE
            corr = self.overlaps['std_ihs'].corr(self.overlaps['relate_neg_log10_p_lifetime'])
            stats_dict['correlation_ihs_relate'] = corr
            print(f"\nCorrelation (Std iHS vs RELATE -log10 p): {corr:.3f}")

            # Population-specific overlaps
            print(f"\nOverlaps by population:")
            pop_counts = self.overlaps.groupby('population').agg({
                'ihs_pos': 'nunique',
                'relate_pos': 'nunique'
            }).sort_values('ihs_pos', ascending=False)
            for pop, row in pop_counts.iterrows():
                print(f"  {pop}: {int(row['ihs_pos'])} unique iHS positions")

        self.stats = stats_dict
        return stats_dict

    def save_results(self, window: int = 50000) -> None:
        """Save all results to files"""
        print(f"\nSaving results to {self.output_dir}/...")

        # Save overlaps
        if self.overlaps is not None and len(self.overlaps) > 0:
            overlap_file = self.output_dir / 'ihs_relate_overlaps.tsv'
            self.overlaps.to_csv(overlap_file, sep='\t', index=False)
            print(f"  Saved overlaps: {overlap_file}")

        # Save iHS-only
        ihs_only = self.find_ihs_only(window=window)
        if len(ihs_only) > 0:
            ihs_only_file = self.output_dir / 'ihs_only_signals.tsv'
            ihs_only.to_csv(ihs_only_file, sep='\t', index=False)
            print(f"  Saved iHS-only signals: {ihs_only_file}")

        # Save RELATE-only
        relate_only = self.find_relate_only(window=window)
        if len(relate_only) > 0:
            relate_only_file = self.output_dir / 'relate_only_signals.tsv'
            relate_only.to_csv(relate_only_file, sep='\t', index=False)
            print(f"  Saved RELATE-only signals: {relate_only_file}")

        # Save summary statistics
        if self.stats:
            stats_file = self.output_dir / 'comparison_statistics.tsv'
            stats_df = pd.DataFrame([self.stats])
            stats_df.to_csv(stats_file, sep='\t', index=False)
            print(f"  Saved statistics: {stats_file}")

        # Save detailed summary by position
        if self.overlaps is not None and len(self.overlaps) > 0:
            summary_by_pos = self.overlaps.groupby('ihs_pos').agg({
                'std_ihs': 'first',
                'ihs_neg_log10_p': 'first',
                'location': 'first',
                'population': lambda x: ','.join(sorted(set(x))),
                'relate_neg_log10_p_lifetime': ['max', 'mean', 'count'],
                'signal_type': lambda x: x.value_counts().index[0]
            }).reset_index()

            summary_by_pos.columns = [
                'position', 'std_ihs', 'ihs_neg_log10_p', 'location',
                'populations', 'max_relate_p', 'mean_relate_p',
                'n_populations', 'dominant_signal_type'
            ]

            summary_file = self.output_dir / 'overlap_summary_by_position.tsv'
            summary_by_pos.to_csv(summary_file, sep='\t', index=False)
            print(f"  Saved position summary: {summary_file}")

    def create_visualizations(self) -> None: #Claude addition
        """Create comparison visualizations"""
        print(f"\nCreating visualizations...")

        if self.overlaps is None or len(self.overlaps) == 0:
            print("  No overlaps to visualize!")
            return

        # Set style
        sns.set_style("whitegrid")

        # 1. Scatter plot: iHS vs RELATE
        fig, ax = plt.subplots(figsize=(10, 8))

        scatter = ax.scatter(
            self.overlaps['std_ihs'],
            self.overlaps['relate_neg_log10_p_lifetime'],
            c=self.overlaps['distance_bp'],
            cmap='viridis',
            alpha=0.6,
            s=50
        )

        ax.axhline(y=8, color='red', linestyle='--', alpha=0.5, label='RELATE p=1e-8')
        ax.axvline(x=3, color='blue', linestyle='--', alpha=0.5, label='Std iHS=3')

        ax.set_xlabel('Standardized iHS', fontsize=12)
        ax.set_ylabel('RELATE -log10(p) [lifetime]', fontsize=12)
        ax.set_title('iHS vs RELATE Selection Signals (Chromosome X)', fontsize=14, fontweight='bold')
        ax.legend()

        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Distance (bp)', fontsize=10)

        plt.tight_layout()
        scatter_file = self.output_dir / 'ihs_vs_relate_scatter.png'
        plt.savefig(scatter_file, dpi=300, bbox_inches='tight')
        print(f"  Saved scatter plot: {scatter_file}")
        plt.close()

        # 2. Signal type distribution
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

        signal_counts = self.overlaps['signal_type'].value_counts()
        ax1.bar(range(len(signal_counts)), signal_counts.values)
        ax1.set_xticks(range(len(signal_counts)))
        ax1.set_xticklabels(signal_counts.index, rotation=45, ha='right')
        ax1.set_ylabel('Count', fontsize=12)
        ax1.set_title('Distribution of Signal Types', fontsize=12, fontweight='bold')
        ax1.grid(axis='y', alpha=0.3)

        # Population support
        pop_counts = self.overlaps.groupby('ihs_pos')['population'].nunique()
        ax2.hist(pop_counts, bins=range(1, pop_counts.max()+2), edgecolor='black', alpha=0.7)
        ax2.set_xlabel('Number of populations with RELATE signal', fontsize=12)
        ax2.set_ylabel('Number of iHS positions', fontsize=12)
        ax2.set_title('Population Support for iHS Signals', fontsize=12, fontweight='bold')
        ax2.grid(axis='y', alpha=0.3)

        plt.tight_layout()
        dist_file = self.output_dir / 'signal_distributions.png'
        plt.savefig(dist_file, dpi=300, bbox_inches='tight')
        print(f"  Saved distribution plots: {dist_file}")
        plt.close()

        # 3. Genomic distribution (Manhattan-style)
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

        # iHS track
        ax1.scatter(self.ihs_data['pos'], self.ihs_data['neg_log10_p'],
                   alpha=0.5, s=20, color='blue', label='All iHS')
        if len(self.overlaps) > 0:
            overlap_positions = self.overlaps.groupby('ihs_pos').agg({
                'ihs_neg_log10_p': 'first'
            }).reset_index()
            ax1.scatter(overlap_positions['ihs_pos'], overlap_positions['ihs_neg_log10_p'],
                       alpha=0.8, s=50, color='red', label='With RELATE support', zorder=5)
        ax1.axhline(y=6, color='gray', linestyle='--', alpha=0.5)
        ax1.set_ylabel('iHS -log10(p)', fontsize=11)
        ax1.set_title('Genomic Distribution of Selection Signals', fontsize=13, fontweight='bold')
        ax1.legend(loc='upper right')
        ax1.grid(axis='y', alpha=0.3)

        # RELATE track (by max p-value per position)
        relate_max = self.relate_data.groupby('Position')['negLog10_P_lifetime'].max().reset_index()
        ax2.scatter(relate_max['Position'], relate_max['negLog10_P_lifetime'],
                   alpha=0.5, s=20, color='green', label='All RELATE')
        if len(self.overlaps) > 0:
            overlap_relate = self.overlaps.groupby('relate_pos').agg({
                'relate_neg_log10_p_lifetime': 'max'
            }).reset_index()
            ax2.scatter(overlap_relate['relate_pos'], overlap_relate['relate_neg_log10_p_lifetime'],
                       alpha=0.8, s=50, color='red', label='With iHS support', zorder=5)
        ax2.axhline(y=6, color='gray', linestyle='--', alpha=0.5)
        ax2.set_xlabel('Position on Chromosome X (bp)', fontsize=11)
        ax2.set_ylabel('RELATE -log10(p)', fontsize=11)
        ax2.legend(loc='upper right')
        ax2.grid(axis='y', alpha=0.3)

        plt.tight_layout()
        manhattan_file = self.output_dir / 'genomic_distribution.png'
        plt.savefig(manhattan_file, dpi=300, bbox_inches='tight')
        print(f"  Saved genomic distribution: {manhattan_file}")
        plt.close()

        # 4. Venn diagram data (saved as text since we don't have venn library)
        venn_data = {
            'iHS_only': len(self.find_ihs_only()),
            'RELATE_only': len(self.find_relate_only()),
            'Both': self.overlaps['ihs_pos'].nunique()
        }

        venn_file = self.output_dir / 'venn_diagram_data.txt'
        with open(venn_file, 'w') as f:
            f.write("Venn Diagram Data\n")
            f.write("=" * 50 + "\n")
            for key, value in venn_data.items():
                f.write(f"{key}: {value}\n")
        print(f"  Saved Venn data: {venn_file}")

    def run_full_analysis(self, window: int = 50000) -> None:
        """Run complete analysis pipeline"""
        print("\n" + "="*70)
        print("iHS vs RELATE COMPARISON ANALYSIS")
        print("="*70)

        self.load_data()
        self.find_overlaps(window=window)
        self.compute_statistics(window=window)
        self.save_results(window=window)
        self.create_visualizations()

        print("\n" + "="*70)
        print("ANALYSIS COMPLETE!")
        print("="*70)
        print(f"\nResults saved to: {self.output_dir}")


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description='Compare iHS and RELATE selection signals',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python compare_ihs_relate.py

  # Custom window size (100kb)
  python compare_ihs_relate.py --window 100000

  # Exact position matching only
  python compare_ihs_relate.py --exact-match

  # Custom output directory
  python compare_ihs_relate.py --output results/comparison
        """
    )

    parser.add_argument(
        '--ihs-file',
        default='/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/analysis/candidates/chrX_high_significance.tsv',
        help='Path to iHS results file (default: chrX_high_significance.tsv)'
    )

    parser.add_argument(
        '--relate-file',
        default='/home/vanbruggenmit/mit-ihh-pib/data/relate/relateData_KM/all_chrX_snps_below_p_6.csv',
        help='Path to RELATE results file (default: all_chrX_snps_below_p_6.csv)'
    )

    parser.add_argument(
        '--output',
        default='/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/analysis/ihs_relate_comparison',
        help='Output directory for results (default: results/analysis/ihs_relate_comparison)'
    )

    parser.add_argument(
        '--window',
        type=int,
        default=50000,
        help='Window size in bp for overlap detection (default: 50000)'
    )

    parser.add_argument(
        '--exact-match',
        action='store_true',
        help='Only consider exact position matches (ignores --window)'
    )

    args = parser.parse_args()

    # Create comparator and run analysis
    comparator = SelectionComparator(
        ihs_file=args.ihs_file,
        relate_file=args.relate_file,
        output_dir=args.output
    )

    if args.exact_match:
        print("Using exact position matching")
        comparator.run_full_analysis(window=0)
    else:
        comparator.run_full_analysis(window=args.window)


if __name__ == '__main__':
    main()
