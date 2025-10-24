#!/usr/bin/env python3
"""
Phase 4: Functional Enrichment Analysis
========================================

Performs functional enrichment analysis on genes in selection regions to test
the hypothesis that X chromosome genes under selection are enriched for
reproductive and immune functions.

Uses g:Profiler API for enrichment testing with GO, KEGG, and Reactome databases.

Author: MIT van Bruggen
Date: October 24, 2025
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json
import sys
from typing import List, Dict, Tuple
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Try to import gprofiler - will provide installation instructions if missing
try:
    from gprofiler import GProfiler
except ImportError:
    print("ERROR: gprofiler-official package not installed.")
    print("Please install with: pip install gprofiler-official")
    sys.exit(1)


class EnrichmentAnalyzer:
    """Performs functional enrichment analysis on gene lists."""

    def __init__(self, base_dir: Path):
        """
        Initialize the enrichment analyzer.

        Args:
            base_dir: Base directory containing results/analysis/gene_annotation/
        """
        self.base_dir = Path(base_dir)
        self.gene_list_dir = self.base_dir / "results" / "analysis" / "gene_annotation" / "gene_name_lists"
        self.output_dir = self.base_dir / "results" / "analysis" / "enrichment"
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Initialize g:Profiler client
        self.gp = GProfiler(return_dataframe=True)

        print(f"✓ Initialized enrichment analyzer")
        print(f"  Input directory: {self.gene_list_dir}")
        print(f"  Output directory: {self.output_dir}")

    def load_gene_list(self, filename: str) -> List[str]:
        """Load a gene list from file."""
        filepath = self.gene_list_dir / filename
        with open(filepath, 'r') as f:
            genes = [line.strip() for line in f if line.strip()]
        print(f"  Loaded {len(genes)} genes from {filename}")
        return genes

    def run_enrichment(self,
                      query_genes: List[str],
                      background_genes: List[str] = None,
                      sources: List[str] = None,
                      name: str = "query") -> pd.DataFrame:
        """
        Run g:Profiler enrichment analysis.

        Args:
            query_genes: List of gene symbols to test
            background_genes: Optional background gene list
            sources: Data sources to query (default: GO:BP, KEGG, REAC)
            name: Name for this analysis (for output files)

        Returns:
            DataFrame with enrichment results
        """
        if sources is None:
            sources = ['GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC']

        print(f"\n{'='*70}")
        print(f"Running enrichment analysis: {name}")
        print(f"{'='*70}")
        print(f"Query genes: {len(query_genes)}")
        if background_genes:
            print(f"Background genes: {len(background_genes)}")
        print(f"Data sources: {', '.join(sources)}")

        try:
            result = self.gp.profile(
                organism='hsapiens',
                query=query_genes,
                background=background_genes,
                sources=sources,
                user_threshold=0.05,  # Significance threshold
                significance_threshold_method='fdr',  # FDR correction
                domain_scope='custom' if background_genes else 'annotated',
                ordered=False
            )

            if isinstance(result, pd.DataFrame) and not result.empty:
                print(f"✓ Found {len(result)} significant terms (FDR < 0.05)")

                # Add analysis name
                result['analysis'] = name

                # Sort by p-value
                result = result.sort_values('p_value')

                # Save results
                output_file = self.output_dir / f"enrichment_{name}.tsv"
                result.to_csv(output_file, sep='\t', index=False)
                print(f"✓ Saved results to {output_file}")

                return result
            else:
                print(f"✗ No significant enrichment found")
                return pd.DataFrame()

        except Exception as e:
            print(f"✗ Error during enrichment analysis: {e}")
            return pd.DataFrame()

    def filter_by_category(self, results: pd.DataFrame, keywords: List[str],
                          name: str = "") -> pd.DataFrame:
        """
        Filter enrichment results by keywords in term names.

        Args:
            results: Enrichment results DataFrame
            keywords: List of keywords to search for (case-insensitive)
            name: Category name for output

        Returns:
            Filtered DataFrame
        """
        if results.empty:
            return results

        # Create regex pattern for OR matching
        pattern = '|'.join(keywords)

        # Filter by name (case-insensitive)
        filtered = results[results['name'].str.contains(pattern, case=False, na=False)]

        if not filtered.empty:
            print(f"\n{name} terms: {len(filtered)} found")
            output_file = self.output_dir / f"enrichment_{name.lower().replace(' ', '_')}.tsv"
            filtered.to_csv(output_file, sep='\t', index=False)

        return filtered

    def create_summary_table(self, results_dict: Dict[str, pd.DataFrame]) -> pd.DataFrame:
        """
        Create summary table comparing enrichment across analyses.

        Args:
            results_dict: Dictionary of {analysis_name: results_df}

        Returns:
            Summary DataFrame
        """
        summaries = []

        for name, df in results_dict.items():
            if df.empty:
                continue

            summary = {
                'analysis': name,
                'total_significant_terms': len(df),
                'unique_sources': df['source'].nunique() if 'source' in df.columns else 0,
                'min_pvalue': df['p_value'].min() if 'p_value' in df.columns else np.nan,
                'mean_enrichment_score': -np.log10(df['p_value']).mean() if 'p_value' in df.columns else np.nan,
                'top_term': df.iloc[0]['name'] if len(df) > 0 and 'name' in df.columns else 'None',
                'top_pvalue': df.iloc[0]['p_value'] if len(df) > 0 and 'p_value' in df.columns else np.nan
            }
            summaries.append(summary)

        summary_df = pd.DataFrame(summaries)

        # Save summary
        output_file = self.output_dir / "enrichment_summary.tsv"
        summary_df.to_csv(output_file, sep='\t', index=False)
        print(f"\n✓ Created summary table: {output_file}")

        return summary_df

    def compare_x_vs_autosome(self,
                             x_results: pd.DataFrame,
                             autosome_results: pd.DataFrame,
                             top_n: int = 20) -> pd.DataFrame:
        """
        Compare enrichment between X chromosome and autosomes.

        Args:
            x_results: Enrichment results for X chromosome genes
            autosome_results: Enrichment results for autosomal genes
            top_n: Number of top terms to show

        Returns:
            Comparison DataFrame
        """
        if x_results.empty or autosome_results.empty:
            print("Cannot compare: one or both result sets are empty")
            return pd.DataFrame()

        print(f"\n{'='*70}")
        print("Comparing X chromosome vs Autosome enrichment")
        print(f"{'='*70}")

        # Get term IDs that appear in both
        x_terms = set(x_results['native'].tolist())
        autosome_terms = set(autosome_results['native'].tolist())

        common_terms = x_terms & autosome_terms
        x_specific = x_terms - autosome_terms
        autosome_specific = autosome_terms - x_terms

        print(f"Common terms: {len(common_terms)}")
        print(f"X-specific terms: {len(x_specific)}")
        print(f"Autosome-specific terms: {len(autosome_specific)}")

        # Create comparison for common terms
        comparisons = []

        for term_id in common_terms:
            x_row = x_results[x_results['native'] == term_id].iloc[0]
            auto_row = autosome_results[autosome_results['native'] == term_id].iloc[0]

            comparison = {
                'term_id': term_id,
                'term_name': x_row['name'],
                'source': x_row['source'],
                'x_pvalue': x_row['p_value'],
                'autosome_pvalue': auto_row['p_value'],
                'x_enrichment_score': -np.log10(x_row['p_value']),
                'autosome_enrichment_score': -np.log10(auto_row['p_value']),
                'enrichment_ratio': -np.log10(x_row['p_value']) / -np.log10(auto_row['p_value']),
                'x_intersection_size': x_row.get('intersection_size', 0),
                'autosome_intersection_size': auto_row.get('intersection_size', 0)
            }
            comparisons.append(comparison)

        comparison_df = pd.DataFrame(comparisons)

        if not comparison_df.empty:
            # Sort by enrichment ratio (X vs autosome)
            comparison_df = comparison_df.sort_values('enrichment_ratio', ascending=False)

            # Save full comparison
            output_file = self.output_dir / "x_vs_autosome_comparison.tsv"
            comparison_df.to_csv(output_file, sep='\t', index=False)
            print(f"✓ Saved comparison to {output_file}")

            # Show top terms more enriched in X
            print(f"\nTop {min(top_n, len(comparison_df))} terms MORE enriched in X chromosome:")
            print(comparison_df.head(top_n)[['term_name', 'enrichment_ratio', 'x_pvalue']].to_string())

        return comparison_df

    def print_top_results(self, results: pd.DataFrame, n: int = 10, title: str = ""):
        """Print top N enrichment results in a nice format."""
        if results.empty:
            print(f"\n{title}: No results to display")
            return

        print(f"\n{'='*70}")
        print(f"{title}")
        print(f"{'='*70}")

        for i, row in results.head(n).iterrows():
            print(f"\n{i+1}. {row['name']}")
            print(f"   Source: {row['source']}")
            print(f"   P-value: {row['p_value']:.2e}")
            if 'intersection_size' in row and 'term_size' in row:
                print(f"   Genes: {row['intersection_size']}/{row['term_size']}")
            if 'intersections' in row:
                genes = row['intersections'].split(',') if isinstance(row['intersections'], str) else []
                if len(genes) <= 5:
                    print(f"   Gene list: {', '.join(genes)}")
                else:
                    print(f"   Gene list: {', '.join(genes[:5])} ... (+{len(genes)-5} more)")


def main():
    """Main analysis workflow."""

    print("\n" + "="*70)
    print("PHASE 4: FUNCTIONAL ENRICHMENT ANALYSIS")
    print("="*70)
    print("\nTesting hypothesis: X chromosome genes under selection are enriched")
    print("for reproductive and immune functions\n")

    # Initialize analyzer
    base_dir = Path("/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib")
    analyzer = EnrichmentAnalyzer(base_dir)

    # Load gene lists
    print("\n" + "-"*70)
    print("Loading gene lists...")
    print("-"*70)

    x_protein_coding = analyzer.load_gene_list("X_protein_coding_genes.txt")
    all_protein_coding = analyzer.load_gene_list("protein_coding_genes.txt")

    # Create autosome gene list (all protein-coding minus X)
    x_set = set(x_protein_coding)
    autosome_protein_coding = [g for g in all_protein_coding if g not in x_set]
    print(f"  Created autosome list: {len(autosome_protein_coding)} genes (all - X)")

    # Also create a proper background (all genes in genome, not just in selection)
    # For now, we'll use all protein-coding genes as background
    background_genes = all_protein_coding

    # ==========================================================================
    # Analysis 1: X chromosome genes - full enrichment
    # ==========================================================================

    x_results = analyzer.run_enrichment(
        query_genes=x_protein_coding,
        background_genes=background_genes,
        name="X_chromosome_full"
    )

    if not x_results.empty:
        analyzer.print_top_results(x_results, n=15,
                                   title="Top 15 Enriched Terms - X Chromosome Genes")

    # ==========================================================================
    # Analysis 2: Autosome genes - full enrichment (for comparison)
    # ==========================================================================

    autosome_results = analyzer.run_enrichment(
        query_genes=autosome_protein_coding,
        background_genes=background_genes,
        name="autosome_full"
    )

    if not autosome_results.empty:
        analyzer.print_top_results(autosome_results, n=15,
                                   title="Top 15 Enriched Terms - Autosomal Genes")

    # ==========================================================================
    # Analysis 3: Filter for reproductive terms
    # ==========================================================================

    print("\n" + "="*70)
    print("FILTERING FOR REPRODUCTIVE TERMS")
    print("="*70)

    reproductive_keywords = [
        'reproduc', 'fertil', 'pregnan', 'ovulation', 'spermatogenesis',
        'oogenesis', 'meiosis', 'gamete', 'embryo', 'gonad', 'testis',
        'ovary', 'hormone', 'steroid'
    ]

    x_reproductive = analyzer.filter_by_category(
        x_results, reproductive_keywords, "X_reproductive"
    )
    autosome_reproductive = analyzer.filter_by_category(
        autosome_results, reproductive_keywords, "autosome_reproductive"
    )

    if not x_reproductive.empty:
        analyzer.print_top_results(x_reproductive, n=10,
                                   title="Reproductive Terms - X Chromosome")

    # ==========================================================================
    # Analysis 4: Filter for immune terms
    # ==========================================================================

    print("\n" + "="*70)
    print("FILTERING FOR IMMUNE TERMS")
    print("="*70)

    immune_keywords = [
        'immun', 'inflamm', 'cytokine', 'interferon', 'interleukin',
        'antibody', 'antigen', 'lymphocyte', 'T cell', 'B cell',
        'natural killer', 'complement', 'toll-like', 'TLR', 'MHC',
        'response to virus', 'response to bacteria', 'defense response'
    ]

    x_immune = analyzer.filter_by_category(
        x_results, immune_keywords, "X_immune"
    )
    autosome_immune = analyzer.filter_by_category(
        autosome_results, immune_keywords, "autosome_immune"
    )

    if not x_immune.empty:
        analyzer.print_top_results(x_immune, n=10,
                                   title="Immune Terms - X Chromosome")

    # ==========================================================================
    # Analysis 5: Compare X vs Autosomes
    # ==========================================================================

    if not x_results.empty and not autosome_results.empty:
        comparison_df = analyzer.compare_x_vs_autosome(x_results, autosome_results)

    # ==========================================================================
    # Create summary statistics
    # ==========================================================================

    print("\n" + "="*70)
    print("CREATING SUMMARY STATISTICS")
    print("="*70)

    results_dict = {
        'X_chromosome': x_results,
        'Autosomes': autosome_results,
        'X_reproductive': x_reproductive,
        'X_immune': x_immune,
        'Autosome_reproductive': autosome_reproductive,
        'Autosome_immune': autosome_immune
    }

    summary = analyzer.create_summary_table(results_dict)
    print("\n" + summary.to_string(index=False))

    # ==========================================================================
    # Statistical test: Are X genes MORE enriched for reproductive/immune?
    # ==========================================================================

    print("\n" + "="*70)
    print("STATISTICAL TEST: X vs AUTOSOME ENRICHMENT")
    print("="*70)

    # Count reproductive genes
    x_repro_genes = set(x_reproductive['intersections'].str.split(',').sum() if not x_reproductive.empty else [])
    auto_repro_genes = set(autosome_reproductive['intersections'].str.split(',').sum() if not autosome_reproductive.empty else [])

    print(f"\nReproductive genes:")
    print(f"  X chromosome: {len(x_repro_genes)} / {len(x_protein_coding)} ({100*len(x_repro_genes)/len(x_protein_coding):.1f}%)")
    print(f"  Autosomes: {len(auto_repro_genes)} / {len(autosome_protein_coding)} ({100*len(auto_repro_genes)/len(autosome_protein_coding):.1f}%)")

    # Fisher's exact test
    if len(x_repro_genes) > 0 and len(auto_repro_genes) > 0:
        contingency_table = [
            [len(x_repro_genes), len(x_protein_coding) - len(x_repro_genes)],
            [len(auto_repro_genes), len(autosome_protein_coding) - len(auto_repro_genes)]
        ]
        oddsratio, pvalue = stats.fisher_exact(contingency_table)
        print(f"\nFisher's exact test:")
        print(f"  Odds ratio: {oddsratio:.3f}")
        print(f"  P-value: {pvalue:.2e}")

        if oddsratio > 1:
            print(f"  → X chromosome has {oddsratio:.2f}x higher odds of reproductive enrichment")
        else:
            print(f"  → Autosomes have {1/oddsratio:.2f}x higher odds of reproductive enrichment")

    # Count immune genes
    x_immune_genes = set(x_immune['intersections'].str.split(',').sum() if not x_immune.empty else [])
    auto_immune_genes = set(autosome_immune['intersections'].str.split(',').sum() if not autosome_immune.empty else [])

    print(f"\nImmune genes:")
    print(f"  X chromosome: {len(x_immune_genes)} / {len(x_protein_coding)} ({100*len(x_immune_genes)/len(x_protein_coding):.1f}%)")
    print(f"  Autosomes: {len(auto_immune_genes)} / {len(autosome_protein_coding)} ({100*len(auto_immune_genes)/len(autosome_protein_coding):.1f}%)")

    # Fisher's exact test
    if len(x_immune_genes) > 0 and len(auto_immune_genes) > 0:
        contingency_table = [
            [len(x_immune_genes), len(x_protein_coding) - len(x_immune_genes)],
            [len(auto_immune_genes), len(autosome_protein_coding) - len(auto_immune_genes)]
        ]
        oddsratio, pvalue = stats.fisher_exact(contingency_table)
        print(f"\nFisher's exact test:")
        print(f"  Odds ratio: {oddsratio:.3f}")
        print(f"  P-value: {pvalue:.2e}")

        if oddsratio > 1:
            print(f"  → X chromosome has {oddsratio:.2f}x higher odds of immune enrichment")
        else:
            print(f"  → Autosomes have {1/oddsratio:.2f}x higher odds of immune enrichment")

    print("\n" + "="*70)
    print("ANALYSIS COMPLETE!")
    print("="*70)
    print(f"\nAll results saved to: {analyzer.output_dir}")
    print("\nKey output files:")
    print("  - enrichment_X_chromosome_full.tsv")
    print("  - enrichment_autosome_full.tsv")
    print("  - enrichment_X_reproductive.tsv")
    print("  - enrichment_X_immune.tsv")
    print("  - x_vs_autosome_comparison.tsv")
    print("  - enrichment_summary.tsv")
    print()


if __name__ == "__main__":
    main()
