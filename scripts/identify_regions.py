#!/usr/bin/env python3
"""
Identify candidate regions under selection from iHS analysis
Groups nearby extreme SNPs into regions
"""
import pandas as pd
import argparse
import re

def parse_location(location):
    """Parse position from location string"""
    match = re.search(r':(\d+)[_:]', location)
    if match:
        return int(match.group(1))
    return None

def load_ihs_data(filepath):
    """Load iHS TSV file"""
    df = pd.read_csv(filepath, sep='\t')
    df['Position'] = df['Location'].apply(parse_location)
    df = df.dropna(subset=['Position'])
    return df

def identify_regions(df, chr_name, threshold=2.0, window=100000):
    """
    Identify candidate regions by grouping nearby extreme SNPs
    
    Args:
        df: DataFrame with iHS data
        chr_name: Chromosome name
        threshold: Minimum |Std iHS| value
        window: Maximum distance (bp) between SNPs in same region
    """
    # Filter for extreme values
    extreme = df[df['Std iHS'].abs() >= threshold].copy()
    extreme = extreme.sort_values('Position')
    
    if len(extreme) == 0:
        return []
    
    regions = []
    current_region = {
        'chr': chr_name,
        'start': extreme.iloc[0]['Position'],
        'end': extreme.iloc[0]['Position'],
        'snps': [extreme.iloc[0]],
        'max_abs_std_ihs': abs(extreme.iloc[0]['Std iHS'])
    }
    
    for idx in range(1, len(extreme)):
        snp = extreme.iloc[idx]
        
        # If within window, extend current region
        if snp['Position'] - current_region['end'] <= window:
            current_region['end'] = snp['Position']
            current_region['snps'].append(snp)
            if abs(snp['Std iHS']) > current_region['max_abs_std_ihs']:
                current_region['max_abs_std_ihs'] = abs(snp['Std iHS'])
        else:
            # Save current region and start new one
            regions.append(current_region)
            current_region = {
                'chr': chr_name,
                'start': snp['Position'],
                'end': snp['Position'],
                'snps': [snp],
                'max_abs_std_ihs': abs(snp['Std iHS'])
            }
    
    # Add last region
    regions.append(current_region)
    
    # Calculate additional stats
    for region in regions:
        region['size_kb'] = (region['end'] - region['start']) / 1000
        region['n_snps'] = len(region['snps'])
        
        # Determine selection direction (based on majority of SNPs)
        positive = sum(1 for snp in region['snps'] if snp['Std iHS'] > 0)
        negative = len(region['snps']) - positive
        region['direction'] = 'positive' if positive > negative else 'negative'
        
        # Get peak SNP
        peak_snp = max(region['snps'], key=lambda s: abs(s['Std iHS']))
        region['peak_location'] = peak_snp['Location']
        region['peak_std_ihs'] = peak_snp['Std iHS']
    
    # Sort by strongest signal
    regions.sort(key=lambda r: r['max_abs_std_ihs'], reverse=True)
    
    return regions

def main():
    parser = argparse.ArgumentParser(description='Identify candidate selection regions')
    parser.add_argument('--input', required=True, help='Input iHS TSV file')
    parser.add_argument('--chr', required=True, help='Chromosome name')
    parser.add_argument('--output', required=True, help='Output BED/TSV file')
    parser.add_argument('--threshold', type=float, default=2.0,
                        help='Minimum |Std iHS| threshold (default: 2.0)')
    parser.add_argument('--window', type=int, default=100000,
                        help='Maximum distance between SNPs in region (default: 100kb)')
    parser.add_argument('--format', choices=['bed', 'tsv'], default='tsv',
                        help='Output format')
    args = parser.parse_args()
    
    print(f"Loading data from {args.input}...")
    df = load_ihs_data(args.input)
    print(f"Loaded {len(df)} SNPs")
    
    extreme_count = len(df[df['Std iHS'].abs() >= args.threshold])
    print(f"Found {extreme_count} SNPs with |Std iHS| ≥ {args.threshold}")
    
    print(f"Identifying regions (window = {args.window} bp)...")
    regions = identify_regions(df, args.chr, args.threshold, args.window)
    print(f"Identified {len(regions)} candidate regions")
    
    # Write output
    if args.format == 'bed':
        # BED format: chr, start, end, name, score
        with open(args.output, 'w') as f:
            for i, region in enumerate(regions):
                name = f"region_{i+1}"
                score = int(min(region['max_abs_std_ihs'] * 100, 1000))
                f.write(f"chr{region['chr']}\t{region['start']}\t{region['end']}\t"
                       f"{name}\t{score}\n")
        print(f"Saved BED file to {args.output}")
    else:
        # TSV format with detailed information
        output_data = []
        for i, region in enumerate(regions):
            output_data.append({
                'region_id': i + 1,
                'chr': region['chr'],
                'start': region['start'],
                'end': region['end'],
                'size_kb': round(region['size_kb'], 2),
                'n_snps': region['n_snps'],
                'max_abs_std_ihs': round(region['max_abs_std_ihs'], 4),
                'direction': region['direction'],
                'peak_location': region['peak_location'],
                'peak_std_ihs': round(region['peak_std_ihs'], 4)
            })
        
        output_df = pd.DataFrame(output_data)
        output_df.to_csv(args.output, sep='\t', index=False)
        print(f"Saved TSV file to {args.output}")
        
        # Print summary
        print("\nTop 10 regions by signal strength:")
        print(output_df.head(10).to_string(index=False))

if __name__ == '__main__':
    main()