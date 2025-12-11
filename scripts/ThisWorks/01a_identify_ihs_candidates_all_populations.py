#!/usr/bin/env python3
"""
Run iHS candidate identification for all populations.
This processes all 26 populations to identify significant iHS candidates
with p-values based on standard normal distribution.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Base directory
base_dir = Path("/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib")

# All populations
POPULATIONS = {
    'AFR': ['ACB', 'ASW', 'ESN', 'GWD', 'LWK', 'MSL', 'YRI'],
    'AMR': ['CLM', 'MXL', 'PEL', 'PUR'],
    'EAS': ['CDX', 'CHB', 'CHS', 'JPT', 'KHV'],
    'EUR': ['CEU', 'FIN', 'GBR', 'IBS', 'TSI'],
    'SAS': ['BEB', 'GIH', 'ITU', 'PJL', 'STU']
}

# Define thresholds
thresholds = {
    "liberal": 2.0,
    "moderate": 2.5,
    "stringent": 3.0
}

# Only process X chromosome
chroms = ['X']


def process_population(pop_code: str):
    """Process iHS data for a single population."""
    print(f"\n{'='*70}")
    print(f"Processing {pop_code}")
    print(f"{'='*70}")

    ihs_dir = base_dir / f"results/ihs_{pop_code}"
    output_dir = base_dir / f"results/ihs_{pop_code}/analysis/candidates"

    # Create output directories
    (output_dir / "candidates/per_chromosome").mkdir(parents=True, exist_ok=True)

    # Process X chromosome only
    chr_name = 'X'
    file = ihs_dir / f"{pop_code}.chr{chr_name}.ihs.tsv"

    if not file.exists():
        print(f"  WARNING: Missing file for {pop_code} chr{chr_name}")
        return None

    print(f"  Reading chr{chr_name}...", end=" ")

    # Read chromosome data
    chr_data = pd.read_csv(file, sep='\t')

    # Parse Location column (format: chr:pos:ref:alt or X:pos:ref:alt)
    location_parts = chr_data['Location'].str.split(':', expand=True)
    chr_data['chr'] = location_parts[0]
    chr_data['pos'] = location_parts[1].astype(int)
    chr_data['ref'] = location_parts[2]
    chr_data['alt'] = location_parts[3]

    print(f"{len(chr_data):,} variants")

    # Add p-values based on standard normal assumption
    # Two-tailed p-value: P(|Z| > |iHS|) where Z ~ N(0,1)
    chr_data['p_value'] = 2 * stats.norm.sf(np.abs(chr_data['Std iHS']))
    chr_data['neg_log10_p'] = -np.log10(chr_data['p_value'])

    # Extract candidates at each threshold
    print(f"  Extracting candidates:")
    for level, thresh in thresholds.items():
        candidates = chr_data[np.abs(chr_data['Std iHS']) >= thresh].copy()
        candidates = candidates.sort_values('Std iHS', key=abs, ascending=False)

        if len(candidates) > 0:
            # Save per-chromosome candidates
            outfile = output_dir / f"candidates/per_chromosome/chr{chr_name}_candidates_{level}.tsv"
            candidates.to_csv(outfile, sep='\t', index=False)
            print(f"    {level}: {len(candidates):,} variants (saved)")

    # Calculate statistics
    chr_stats = {
        'population': pop_code,
        'chr': chr_name,
        'n_variants': len(chr_data),
        'mean_std_iHS': chr_data['Std iHS'].mean(),
        'sd_std_iHS': chr_data['Std iHS'].std(),
        'max_abs_std_iHS': np.abs(chr_data['Std iHS']).max(),
        'q99': np.abs(chr_data['Std iHS']).quantile(0.99),
        'q995': np.abs(chr_data['Std iHS']).quantile(0.995),
        'q999': np.abs(chr_data['Std iHS']).quantile(0.999),
        'n_liberal': (np.abs(chr_data['Std iHS']) >= 2.0).sum(),
        'n_moderate': (np.abs(chr_data['Std iHS']) >= 2.5).sum(),
        'n_stringent': (np.abs(chr_data['Std iHS']) >= 3.0).sum()
    }

    return chr_stats


def main():
    """Process all populations."""
    print(f"{'='*70}")
    print("iHS Candidate Identification - All Populations")
    print(f"{'='*70}")
    print("\nProcessing chromosome X for all populations")

    # Get all populations
    all_populations = []
    for pops in POPULATIONS.values():
        all_populations.extend(pops)

    print(f"\nTotal populations: {len(all_populations)}")
    print(f"Populations: {', '.join(all_populations)}")

    # Process each population
    all_stats = []
    for pop_code in all_populations:
        stats = process_population(pop_code)
        if stats is not None:
            all_stats.append(stats)

    # Save summary statistics
    if all_stats:
        stats_df = pd.DataFrame(all_stats)
        output_file = base_dir / "results/ihs_candidates_summary_all_populations.tsv"
        stats_df.to_csv(output_file, sep='\t', index=False)

        print(f"\n{'='*70}")
        print("SUMMARY STATISTICS")
        print(f"{'='*70}\n")
        print(stats_df.to_string(index=False))

        print(f"\n\nSaved summary: {output_file}")

        print(f"\n{'='*70}")
        print("Processing Complete!")
        print(f"{'='*70}")
        print(f"\nProcessed {len(all_stats)} populations")
        print("\nOutput files created for each population:")
        print("  results/ihs_{POP}/analysis/candidates/candidates/per_chromosome/")
        print("    - chrX_candidates_liberal.tsv (|Std iHS| >= 2.0)")
        print("    - chrX_candidates_moderate.tsv (|Std iHS| >= 2.5)")
        print("    - chrX_candidates_stringent.tsv (|Std iHS| >= 3.0)")
        print("\nEach file includes:")
        print("  - Location, chr, pos, ref, alt")
        print("  - iHH_0, iHH_1, iHS, Std iHS")
        print("  - p_value, neg_log10_p (from standard normal)")


if __name__ == '__main__':
    main()
