#!/usr/bin/env python3
"""
Subset chromosome X hapbin data to specific population for population-specific iHS analysis

Usage:
    python subset_population_hapbin.py <POP_CODE>

Example:
    python subset_population_hapbin.py YRI
    python subset_population_hapbin.py CEU
"""

from pathlib import Path
import sys

def subset_population(pop_code):
    """Subset hapbin data for a specific population."""

    # Paths
    base_dir = Path("/home/vanbruggenmit/mit-ihh-pib")
    data_dir = base_dir / "data/grch38/hapbin"
    output_dir = base_dir / f"data/grch38/hapbin_{pop_code}"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Sample files
    sample_pop_file = base_dir / "data/grch38/sample_population.txt"

    # Input files (October 2025 backup - what was actually used)
    # I tried to re-run iHS in multiple ways because at some point I thought the 1000G data wasn't suitable for the hapbin software
    # In the end it turned out my initial approach was fine, which is why I'm using these backup files here
    hap_file = data_dir / "chrX.hapbin.hap.backup_20251122_185513"
    map_file = data_dir / "chrX.hapbin.map.backup_20251122_185513"

    # Output files
    output_hap = output_dir / f"chrX_{pop_code}.hapbin.hap"
    output_map = output_dir / f"chrX_{pop_code}.hapbin.map"
    output_sample_list = output_dir / f"{pop_code}_samples.txt"

    print("=" * 80)
    print(f"Subsetting Chromosome X Hapbin Data to {pop_code} Samples")
    print("=" * 80)
    print()
    print("Using October 2025 backup files")
    print()

    # Load sample IDs for this population
    print(f"Loading {pop_code} sample IDs...")
    pop_samples = set()
    all_samples_ordered = []

    with open(sample_pop_file) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                sample_id, pop = parts[0], parts[1]
                all_samples_ordered.append(sample_id)
                if pop == pop_code:
                    pop_samples.add(sample_id)

    if len(pop_samples) == 0:
        print(f"ERROR: No samples found for population {pop_code}")
        print(f"Available populations in {sample_pop_file}:")
        pops = set()
        with open(sample_pop_file) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    pops.add(parts[1])
        for p in sorted(pops):
            print(f"  {p}")
        sys.exit(1)

    print(f"  Found {len(pop_samples)} {pop_code} samples")
    print()

    print(f"Total samples in order: {len(all_samples_ordered)}")

    # Determine which columns to keep
    print()
    print(f"Determining {pop_code} sample columns in HAP file...")

    # Read first line to check format
    with open(hap_file, 'r') as f:
        first_line = f.readline().strip().split()
        total_haplotypes = len(first_line)

    print(f"  Total haplotypes in file: {total_haplotypes}")

    # Build column indices for population samples
    pop_column_indices = []
    for idx, sample_id in enumerate(all_samples_ordered):
        if sample_id in pop_samples:
            # Assuming 2 haplotypes per sample
            col1 = idx * 2
            col2 = idx * 2 + 1
            if col1 < total_haplotypes:
                pop_column_indices.append(col1)
            if col2 < total_haplotypes:
                pop_column_indices.append(col2)

    print(f"  {pop_code} haplotype columns to extract: {len(pop_column_indices)}")
    print(f"  Expected for {len(pop_samples)} samples: ~{len(pop_samples)}-{len(pop_samples)*2} (depending on male/female for X)")
    print()

    # Save population sample list
    with open(output_sample_list, 'w') as f:
        for sample in sorted(pop_samples):
            f.write(f"{sample}\n")
    print(f"Saved {pop_code} sample list: {output_sample_list}")

    # Subset HAP file
    print()
    print("Subsetting HAP file (this may take several minutes)...")
    print(f"  Input: {hap_file}")
    print(f"  Input size: {hap_file.stat().st_size / (1024**3):.2f} GB")
    print(f"  Output: {output_hap}")

    lines_processed = 0
    with open(hap_file, 'r') as fin, open(output_hap, 'w') as fout:
        for line in fin:
            fields = line.strip().split()

            # Extract only population columns
            pop_fields = [fields[i] for i in pop_column_indices if i < len(fields)]

            # Write population-only line
            fout.write(' '.join(pop_fields) + '\n')

            lines_processed += 1
            if lines_processed % 10000 == 0:
                print(f"    Processed {lines_processed:,} variants...", end='\r')

    print(f"\n  Processed {lines_processed:,} variants")
    print(f"  Output: {output_hap}")
    print(f"  Output size: {output_hap.stat().st_size / (1024**3):.2f} GB")
    print()

    # Copy MAP file (same for all populations - just physical/genetic positions)
    print("Copying MAP file...")
    with open(map_file, 'r') as fin, open(output_map, 'w') as fout:
        map_content = fin.read()
        fout.write(map_content)
        map_lines = map_content.count('\n')

    print(f"  Copied: {output_map}")
    print(f"  MAP lines: {map_lines:,}")
    print()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python subset_population_hapbin.py <POP_CODE>")
        print()
        print("Example:")
        print("  python subset_population_hapbin.py YRI")
        print("  python subset_population_hapbin.py CEU")
        sys.exit(1)

    pop_code = sys.argv[1].upper()
    subset_population(pop_code)
