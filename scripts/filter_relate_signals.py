#!/usr/bin/env python3
"""
Filter RELATE selection signals based on -log10(p-value) threshold.

This script filters variants from RELATE output based on the p-value for selection
in the entire lifetime of the variant (column 5).

Input: CSV file with columns: Population, Chromosome, Position, log10(P-value_50%), log10(P-value_lifetime)
Output: Filtered CSV with variants where -log10(p-value) >= 6
"""

import pandas as pd
import os
from pathlib import Path

# Input and output paths
input_file = "/home/vanbruggenmit/mit-ihh-pib/data/relate/relateData_KM/all_chrX_snps_below_p_6.csv"
output_dir = "/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/analysis/ihs_relate_comparison/relate_filtering"
output_file = os.path.join(output_dir, "chrX_signals_log10p_ge_6.csv")

# Threshold for filtering
LOG10_P_THRESHOLD = -6  # Since values are negative log10(p-values), we filter <= -6

def main():
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Read the CSV file using pandas
    df = pd.read_csv(input_file, header=None,
                     names=['Population', 'Chromosome', 'Position',
                            'log10_p_50pct', 'log10_p_lifetime'])
    # Filter based on the lifetime p-value (column 5)
    # I want -log10(p-value) >= 6, which means log10(p-value) <= -6
    filtered_df = df[df['log10_p_lifetime'] <= LOG10_P_THRESHOLD].copy() #make a copy of the filtered dataframe

    print(f"Variants after filtering (-log10(p) >= 6): {len(filtered_df)}")

    # Add a column with the actual -log10(p-value) for clarity in output file
    filtered_df['minus_log10_p'] = -filtered_df['log10_p_lifetime']

    # Save filtered results
    filtered_df.to_csv(output_file, index=False)

    # Print summary statistics
    print("\n=== Summary Statistics ===")
    print(f"Number of populations: {filtered_df['Population'].nunique()}")
    print(f"Populations: {', '.join(sorted(filtered_df['Population'].unique()))}")
    print(f"\nRange of -log10(p-values): {filtered_df['minus_log10_p'].min():.2f} to {filtered_df['minus_log10_p'].max():.2f}")
    print(f"Mean -log10(p-value): {filtered_df['minus_log10_p'].mean():.2f}")
    print(f"Median -log10(p-value): {filtered_df['minus_log10_p'].median():.2f}")

    # Show variants with strongest signals
    print("\n=== Top 10 strongest signals ===")
    top_signals = filtered_df.nlargest(10, 'minus_log10_p')[
        ['Population', 'Position', 'minus_log10_p']
    ]
    print(top_signals.to_string(index=False))

if __name__ == "__main__":
    main()
