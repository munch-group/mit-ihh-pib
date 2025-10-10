#!/usr/bin/env python3
"""
Create Manhattan plots for iHS analysis
"""
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import re

def parse_location(location):
    """Parse position from location string (handles both formats)"""
    # Format 1: chr21:14158687_T_C
    match = re.search(r':(\d+)_', location)
    if match:
        return int(match.group(1))
    # Format 2: 22:17071944:C:CT
    match = re.search(r'^\d+:(\d+):', location)
    if match:
        return int(match.group(1))
    return None

def load_ihs_data(filepath):
    """Load iHS TSV file"""
    df = pd.read_csv(filepath, sep='\t')
    df['Position'] = df['Location'].apply(parse_location)
    df = df.dropna(subset=['Position'])
    df['Position_Mb'] = df['Position'] / 1_000_000
    return df

def plot_manhattan(df, chr_name, output_file, threshold=2):
    """Create Manhattan plot"""
    fig, ax = plt.subplots(figsize=(14, 6))
    
    # Plot all points
    ax.scatter(df['Position_Mb'], df['Std iHS'], 
               alpha=0.6, s=10, c='#3b82f6', label='SNPs')
    
    # Highlight extreme values
    extreme = df[df['Std iHS'].abs() >= threshold]
    ax.scatter(extreme['Position_Mb'], extreme['Std iHS'], 
               alpha=0.8, s=20, c='red', label=f'|Std iHS| ≥ {threshold}')
    
    # Add threshold lines
    ax.axhline(y=threshold, color='red', linestyle='--', alpha=0.5, linewidth=1)
    ax.axhline(y=-threshold, color='red', linestyle='--', alpha=0.5, linewidth=1)
    ax.axhline(y=0, color='gray', linestyle='-', alpha=0.3, linewidth=0.5)
    
    ax.set_xlabel('Position (Mb)', fontsize=12)
    ax.set_ylabel('Standardized iHS', fontsize=12)
    ax.set_title(f'iHS Manhattan Plot - Chromosome {chr_name}', fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    print(f"Saved Manhattan plot to {output_file}")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Create Manhattan plots for iHS data')
    parser.add_argument('--input', required=True, help='Input iHS TSV file')
    parser.add_argument('--chr', required=True, help='Chromosome name (e.g., 21, 22)')
    parser.add_argument('--output', required=True, help='Output PNG file')
    parser.add_argument('--threshold', type=float, default=2.0, 
                        help='Threshold for extreme values (default: 2.0)')
    args = parser.parse_args()
    
    print(f"Loading data from {args.input}...")
    df = load_ihs_data(args.input)
    print(f"Loaded {len(df)} SNPs")
    
    extreme_count = len(df[df['Std iHS'].abs() >= args.threshold])
    print(f"Found {extreme_count} SNPs with |Std iHS| ≥ {args.threshold}")
    
    plot_manhattan(df, args.chr, args.output, args.threshold)

if __name__ == '__main__':
    main()