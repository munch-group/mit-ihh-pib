#!/usr/bin/env python3
"""
Split haplotype files by population for XP-EHH analysis.
This script takes a combined .hap file and splits it into population-specific files.
"""

import argparse
import sys
from pathlib import Path


def load_sample_list(sample_file):
    """Load sample IDs from a file (one per line)."""
    with open(sample_file) as f:
        return [line.strip() for line in f if line.strip()]


def get_sample_indices_from_vcf(vcf_file, target_samples):
    """
    Get column indices for target samples from VCF header.
    Returns dict mapping sample ID to (hap1_idx, hap2_idx) in the .hap file.
    """
    import gzip

    open_func = gzip.open if vcf_file.endswith('.gz') else open

    with open_func(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                # VCF header line with sample names
                fields = line.strip().split('\t')
                vcf_samples = fields[9:]  # Samples start at column 9

                # Map sample to haplotype column indices
                sample_to_hap_cols = {}
                for sample in target_samples:
                    if sample in vcf_samples:
                        vcf_idx = vcf_samples.index(sample)
                        # Each sample has 2 haplotypes in .hap file
                        hap1_idx = vcf_idx * 2
                        hap2_idx = vcf_idx * 2 + 1
                        sample_to_hap_cols[sample] = (hap1_idx, hap2_idx)

                return sample_to_hap_cols, len(vcf_samples) * 2  # Total haplotypes

    raise ValueError(f"Could not find #CHROM header in {vcf_file}")


def split_hap_file(hap_file, map_file, vcf_file, output_prefix, populations):
    """
    Split haplotype file by populations.

    Args:
        hap_file: Input .hap file (all populations)
        map_file: Input .map file (shared)
        vcf_file: Original VCF file (for sample order)
        output_prefix: Output file prefix
        populations: Dict of {pop_name: sample_list_file}
    """
    print(f"Reading haplotype file: {hap_file}")
    print(f"Reading map file: {map_file}")
    print(f"Reading VCF header from: {vcf_file}")

    # Load sample lists for each population
    pop_samples = {}
    pop_indices = {}

    for pop_name, sample_file in populations.items():
        print(f"\nLoading {pop_name} samples from {sample_file}")
        samples = load_sample_list(sample_file)
        print(f"  Found {len(samples)} samples")
        pop_samples[pop_name] = samples

    # Get column indices for each population's samples
    print(f"\nExtracting sample indices from VCF header...")

    for pop_name, samples in pop_samples.items():
        sample_to_hap, total_haps = get_sample_indices_from_vcf(vcf_file, samples)

        # Flatten to list of column indices to extract
        indices = []
        for sample in samples:
            if sample in sample_to_hap:
                hap1, hap2 = sample_to_hap[sample]
                indices.extend([hap1, hap2])

        pop_indices[pop_name] = sorted(indices)
        print(f"  {pop_name}: {len(indices)} haplotypes ({len(sample_to_hap)} samples found in VCF)")

    # Open output files
    out_files = {}
    for pop_name in populations.keys():
        out_hap = f"{output_prefix}.{pop_name}.hap"
        out_files[pop_name] = open(out_hap, 'w')
        print(f"\nCreating output file: {out_hap}")

    # Process haplotype file line by line
    print(f"\nProcessing haplotype file...")
    line_count = 0

    with open(hap_file) as f:
        for line in f:
            line_count += 1
            if line_count % 10000 == 0:
                print(f"  Processed {line_count} variants...", end='\r')

            # Split line into haplotypes (space-separated)
            haps = line.strip().split()

            # Extract haplotypes for each population
            for pop_name, indices in pop_indices.items():
                pop_haps = [haps[i] for i in indices if i < len(haps)]
                out_files[pop_name].write(' '.join(pop_haps) + '\n')

    print(f"\n  Total variants processed: {line_count}")

    # Close output files
    for f in out_files.values():
        f.close()

    # Copy map file for each population
    print(f"\nCopying map file for each population...")
    for pop_name in populations.keys():
        out_map = f"{output_prefix}.{pop_name}.map"
        import shutil
        shutil.copy(map_file, out_map)
        print(f"  Created: {out_map}")

    print("\nDone!")

    # Print summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    for pop_name, indices in pop_indices.items():
        print(f"{pop_name}:")
        print(f"  Haplotypes: {len(indices)}")
        print(f"  Output files: {output_prefix}.{pop_name}.{{hap,map}}")


def main():
    parser = argparse.ArgumentParser(
        description='Split haplotype files by population for XP-EHH analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  %(prog)s \\
    --hap /path/to/chrX.hap \\
    --map /path/to/chrX.map \\
    --vcf /path/to/chrX.vcf.gz \\
    --output /path/to/output_prefix \\
    --eur samples_EUR.txt \\
    --afr samples_AFR.txt \\
    --eas samples_EAS.txt
        """
    )

    parser.add_argument('--hap', required=True, help='Input .hap file (all populations)')
    parser.add_argument('--map', required=True, help='Input .map file')
    parser.add_argument('--vcf', required=True, help='Original VCF file (for sample order)')
    parser.add_argument('--output', required=True, help='Output file prefix')
    parser.add_argument('--eur', help='EUR sample list file')
    parser.add_argument('--afr', help='AFR sample list file')
    parser.add_argument('--eas', help='EAS sample list file')
    parser.add_argument('--amr', help='AMR sample list file')
    parser.add_argument('--sas', help='SAS sample list file')

    args = parser.parse_args()

    # Build populations dict
    populations = {}
    if args.eur:
        populations['EUR'] = args.eur
    if args.afr:
        populations['AFR'] = args.afr
    if args.eas:
        populations['EAS'] = args.eas
    if args.amr:
        populations['AMR'] = args.amr
    if args.sas:
        populations['SAS'] = args.sas

    if not populations:
        print("ERROR: At least one population must be specified", file=sys.stderr)
        sys.exit(1)

    split_hap_file(args.hap, args.map, args.vcf, args.output, populations)


if __name__ == '__main__':
    main()
