#!/usr/bin/env python3
"""
Fix chromosome X haplotype data for hapbin compatibility.

Males should have dummy '-' characters for their second haplotype in non-PAR regions.
This ensures all lines have the same number of fields, which hapbin requires.

Key difference from fix_chrX_haplotypes.py:
- This version outputs 2 haplotypes per individual for ALL regions
- Males in non-PAR get: real_haplotype + '-' (dummy)
- Females everywhere get: haplotype1 + haplotype2
- PAR regions: everyone gets 2 real haplotypes

This follows IMPUTE2 chrX conventions that hapbin expects.

GRCh38 PAR boundaries:
- PAR1: X:10001-2781479
- PAR2: X:155701383-156030895
"""

import sys
import gzip
import argparse
from collections import Counter

# GRCh38 PAR boundaries
PAR1_START = 10001
PAR1_END = 2781479
PAR2_START = 155701383
PAR2_END = 156030895

def auto_open(path, mode='rt'):
    """Open regular or gzipped files."""
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)

def is_par_region(pos):
    """Check if position is in a pseudoautosomal region."""
    return (PAR1_START <= pos <= PAR1_END) or (PAR2_START <= pos <= PAR2_END)

def infer_sex_from_vcf(vcf_path, sample_size=1000, verbose=True):
    """
    Infer sample sex by analyzing non-PAR X chromosome genotypes.
    Males should be mostly homozygous (0|0 or 1|1) in non-PAR regions.

    Returns: list of booleans, True = male, False = female
    """
    if verbose:
        sys.stderr.write("Inferring sample sex from VCF genotypes...\n")

    heterozygosity = []  # Track het rate per sample
    sample_ids = []
    n_variants_checked = 0

    with auto_open(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                # Extract sample IDs
                parts = line.rstrip('\n').split('\t')
                sample_ids = parts[9:]
                n_samples = len(sample_ids)
                heterozygosity = [[] for _ in range(n_samples)]
                if verbose:
                    sys.stderr.write(f"Found {n_samples} samples in VCF\n")
                continue

            # Process variant
            parts = line.rstrip('\n').split('\t')
            chrom = parts[0]
            pos = int(parts[1])

            # Skip PAR regions (diploid for everyone)
            if is_par_region(pos):
                continue

            # Get genotypes
            genotypes = parts[9:]
            for i, gt_field in enumerate(genotypes):
                gt = gt_field.split(':')[0]  # Get GT field
                if '|' in gt or '/' in gt:
                    alleles = gt.replace('|', '/').split('/')
                    if len(alleles) == 2 and alleles[0] != '.' and alleles[1] != '.':
                        is_het = (alleles[0] != alleles[1])
                        heterozygosity[i].append(is_het)

            n_variants_checked += 1
            if n_variants_checked >= sample_size:
                break

    if verbose:
        sys.stderr.write(f"Analyzed {n_variants_checked} non-PAR variants\n")

    # Calculate heterozygosity rate and infer sex
    is_male = []
    males = 0
    females = 0

    for i, het_list in enumerate(heterozygosity):
        if not het_list:
            # No data, assume female to be safe
            is_male.append(False)
            females += 1
            continue

        het_rate = sum(het_list) / len(het_list)
        # Males should have <5% heterozygosity in non-PAR X (mostly 0)
        # Females should have ~30-50% heterozygosity
        if het_rate < 0.05:
            is_male.append(True)
            males += 1
        else:
            is_male.append(False)
            females += 1

    if verbose:
        sys.stderr.write(f"Inferred: {males} males, {females} females\n")
        sys.stderr.write(f"Output format: 2 haplotypes per individual (6404 total)\n")
        sys.stderr.write(f"  Males in non-PAR: real + '-' (dummy)\n")
        sys.stderr.write(f"  Females: haplotype1 + haplotype2\n")
        sys.stderr.write(f"  PAR regions: everyone gets 2 real haplotypes\n")

    return is_male, sample_ids

def fix_hap_file_for_hapbin(hap_in_path, legend_in_path, hap_out_path, is_male, verbose=True):
    """
    Fix HAP file for hapbin compatibility.

    Output format: ALWAYS 2 columns per individual (6404 total)
    - Males in non-PAR: real haplotype + '-'
    - Females everywhere: haplotype1 + haplotype2
    - PAR regions: everyone gets 2 real haplotypes

    This ensures all lines have the same number of fields.
    """
    if verbose:
        sys.stderr.write(f"Processing HAP file: {hap_in_path}\n")
        sys.stderr.write(f"Output: {hap_out_path}\n")

    n_samples = len(is_male)
    n_males = sum(is_male)
    n_females = n_samples - n_males

    if verbose:
        sys.stderr.write(f"Samples: {n_samples} ({n_males} males, {n_females} females)\n")

    # Read legend to get positions
    if verbose:
        sys.stderr.write("Reading legend file for positions...\n")

    positions = []
    with auto_open(legend_in_path, 'rt') as f:
        header = f.readline()  # Skip header
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                pos = int(parts[1])  # position column
                positions.append(pos)

    if verbose:
        sys.stderr.write(f"Found {len(positions)} variants in legend\n")

    # Process HAP file
    with auto_open(hap_in_path, 'rt') as f_in, auto_open(hap_out_path, 'wt') as f_out:
        for line_no, line in enumerate(f_in):
            if line_no >= len(positions):
                sys.stderr.write(f"WARNING: More HAP lines than legend entries\n")
                break

            pos = positions[line_no]
            is_par = is_par_region(pos)

            haps = line.strip().split()

            # Verify we have expected number of haplotypes
            expected_haps = n_samples * 2
            if len(haps) != expected_haps:
                sys.stderr.write(f"WARNING: Line {line_no+1} has {len(haps)} haps, expected {expected_haps}\n")

            # Build output: ALWAYS 2 haplotypes per individual
            output_haps = []
            for sample_idx in range(n_samples):
                hap1_idx = sample_idx * 2
                hap2_idx = sample_idx * 2 + 1

                if hap1_idx >= len(haps):
                    # Missing data - fill with missing
                    output_haps.extend(['-', '-'])
                    continue

                hap1 = haps[hap1_idx] if hap1_idx < len(haps) else '-'
                hap2 = haps[hap2_idx] if hap2_idx < len(haps) else '-'

                if is_par:
                    # PAR region: both haplotypes are real for everyone
                    output_haps.extend([hap1, hap2])
                else:
                    # Non-PAR region
                    if is_male[sample_idx]:
                        # Male: first haplotype is real, second is dummy
                        output_haps.extend([hap1, '-'])
                    else:
                        # Female: both haplotypes are real
                        output_haps.extend([hap1, hap2])

            # Write line with consistent field count
            f_out.write(' '.join(output_haps) + '\n')

            if verbose and (line_no + 1) % 100000 == 0:
                par_status = "PAR" if is_par else "non-PAR"
                sys.stderr.write(f"Processed {line_no+1} variants ({par_status}, {len(output_haps)} fields)...\n")

    if verbose:
        sys.stderr.write(f"✓ Completed: {line_no+1} variants processed\n")
        sys.stderr.write(f"  All lines have {n_samples * 2} fields (2 per individual)\n")

def main():
    parser = argparse.ArgumentParser(
        description="Fix chrX HAP file for hapbin (constant field count with dummy '-' for males)",
        formatter_class=argparse.RawDescriptionHelpFormatter, 
        epilog=__doc__
    )
    parser.add_argument('--vcf', required=True, help='Input VCF file (to infer sex)')
    parser.add_argument('--hap-in', required=True, help='Input HAP file')
    parser.add_argument('--legend-in', required=True, help='Input LEGEND file')
    parser.add_argument('--hap-out', required=True, help='Output HAP file (fixed for hapbin)')
    parser.add_argument('--sample-size', type=int, default=1000,
                       help='Number of variants to use for sex inference (default: 1000)')
    parser.add_argument('--verbose', action='store_true', help='Verbose output')

    args = parser.parse_args()

    if args.verbose:
        sys.stderr.write("="*60 + "\n")
        sys.stderr.write("Fix chrX haplotypes for hapbin compatibility\n")
        sys.stderr.write("="*60 + "\n\n")

    # Step 1: Infer sex from VCF
    is_male, sample_ids = infer_sex_from_vcf(
        args.vcf,
        sample_size=args.sample_size,
        verbose=args.verbose
    )

    if args.verbose:
        sys.stderr.write("\n")

    # Step 2: Fix HAP file for hapbin
    fix_hap_file_for_hapbin(
        args.hap_in,
        args.legend_in,
        args.hap_out,
        is_male,
        verbose=args.verbose
    )

    if args.verbose:
        sys.stderr.write("\n" + "="*60 + "\n")
        sys.stderr.write("✓ chrX haplotype fixing complete!\n")
        sys.stderr.write("  Output has constant field count for hapbin compatibility\n")
        sys.stderr.write("="*60 + "\n")

if __name__ == "__main__":
    main()
