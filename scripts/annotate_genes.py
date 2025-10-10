#!/usr/bin/env python3
"""
Annotate candidate regions with genes using Ensembl REST API
Requires internet connection to query Ensembl
"""
import pandas as pd
import requests
import time
import argparse
import sys

def query_ensembl_genes(chr_name, start, end, species='human'):
    """
    Query Ensembl REST API for genes in a genomic region
    
    Args:
        chr_name: Chromosome (e.g., '21', '22')
        start: Start position (bp)
        end: End position (bp)
        species: Species name (default: 'human')
    
    Returns:
        List of gene dictionaries
    """
    server = "https://rest.ensembl.org"
    ext = f"/overlap/region/{species}/{chr_name}:{start}-{end}?feature=gene"
    
    try:
        response = requests.get(server + ext, headers={"Content-Type": "application/json"})
        
        if response.status_code == 429:
            # Rate limited - wait and retry
            retry_after = int(response.headers.get('Retry-After', 60))
            print(f"Rate limited. Waiting {retry_after} seconds...", file=sys.stderr)
            time.sleep(retry_after)
            return query_ensembl_genes(chr_name, start, end, species)
        
        if response.status_code != 200:
            print(f"Warning: Failed to query chr{chr_name}:{start}-{end} "
                  f"(status {response.status_code})", file=sys.stderr)
            return []
        
        data = response.json()
        
        # Extract relevant gene information
        genes = []
        for item in data:
            if item.get('feature_type') == 'gene':
                genes.append({
                    'gene_id': item.get('id'),
                    'gene_name': item.get('external_name', 'N/A'),
                    'gene_biotype': item.get('biotype'),
                    'gene_start': item.get('start'),
                    'gene_end': item.get('end'),
                    'strand': '+' if item.get('strand') == 1 else '-'
                })
        
        return genes
    
    except Exception as e:
        print(f"Error querying Ensembl: {e}", file=sys.stderr)
        return []

def annotate_regions(regions_file, output_file, delay=0.5):
    """
    Annotate regions with genes from Ensembl
    
    Args:
        regions_file: TSV file with candidate regions
        output_file: Output file with gene annotations
        delay: Delay between API calls (seconds)
    """
    # Load regions
    df = pd.read_csv(regions_file, sep='\t')
    print(f"Loaded {len(df)} regions")
    
    annotated_regions = []
    
    for idx, row in df.iterrows():
        print(f"Querying region {idx + 1}/{len(df)}: "
              f"chr{row['chr']}:{row['start']}-{row['end']}", file=sys.stderr)
        
        genes = query_ensembl_genes(row['chr'], row['start'], row['end'])
        
        if genes:
            gene_names = [g['gene_name'] for g in genes]
            gene_ids = [g['gene_id'] for g in genes]
            gene_biotypes = [g['gene_biotype'] for g in genes]
            
            annotated_regions.append({
                **row.to_dict(),
                'n_genes': len(genes),
                'gene_names': ','.join(gene_names),
                'gene_ids': ','.join(gene_ids),
                'gene_biotypes': ','.join(gene_biotypes)
            })
        else:
            annotated_regions.append({
                **row.to_dict(),
                'n_genes': 0,
                'gene_names': 'None',
                'gene_ids': 'None',
                'gene_biotypes': 'None'
            })
        
        # Be nice to Ensembl servers
        time.sleep(delay)
    
    # Save results
    output_df = pd.DataFrame(annotated_regions)
    output_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nSaved annotated regions to {output_file}")
    
    # Print summary
    regions_with_genes = output_df[output_df['n_genes'] > 0]
    print(f"\nSummary:")
    print(f"  Total regions: {len(output_df)}")
    print(f"  Regions with genes: {len(regions_with_genes)}")
    print(f"  Total genes found: {output_df['n_genes'].sum()}")
    
    if len(regions_with_genes) > 0:
        print(f"\nTop 10 regions with genes:")
        print(regions_with_genes[['region_id', 'chr', 'start', 'end', 
                                   'n_genes', 'gene_names', 'max_abs_std_ihs']]
              .head(10).to_string(index=False))

def main():
    parser = argparse.ArgumentParser(
        description='Annotate candidate regions with genes using Ensembl REST API'
    )
    parser.add_argument('--regions', required=True, 
                        help='Input TSV file with candidate regions')
    parser.add_argument('--output', required=True,
                        help='Output TSV file with gene annotations')
    parser.add_argument('--delay', type=float, default=0.5,
                        help='Delay between API calls in seconds (default: 0.5)')
    args = parser.parse_args()
    
    print("Querying Ensembl REST API for gene annotations...")
    print("This may take several minutes depending on the number of regions.")
    
    annotate_regions(args.regions, args.output, args.delay)

if __name__ == '__main__':
    main()