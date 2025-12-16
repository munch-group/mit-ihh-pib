#!/usr/bin/env python3
"""
Subset chromosomes 1-22 hapbin data to AFR, EUR, and EAS samples for population-specific iHS analysis

This script creates population-specific versions for all autosomes to identify
population-specific selection signals.

Populations:
- AFR: African ancestry
- EUR: European ancestry
- EAS: East Asian ancestry

Output directories:
- /home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin_AFR/
- /home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin_EUR/
- /home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin_EAS/
"""

from pathlib import Path
import sys

# Paths
base_dir = Path("/home/vanbruggenmit/mit-ihh-pib")
data_dir = base_dir / "data/grch38/hapbin"

# Populations to process
POPULATIONS = ["AFR", "EUR", "EAS"]

# Chromosomes to process
CHROMOSOMES = [f"chr{i}" for i in range(1, 23)]

# Sample files
sample_pop_file = base_dir / "data/grch38/sample_population.txt"

print("=" * 80)
print("Subsetting Chromosomes 1-22 Hapbin Data to AFR, EUR, and EAS Populations")
print("=" * 80)
print()
print(f"Populations: {', '.join(POPULATIONS)}")
print(f"Chromosomes: chr1-chr22 ({len(CHROMOSOMES)} chromosomes)")
print()

# Load sample-to-population mapping
print("Loading sample population mapping...")
all_samples_ordered = []
sample_to_pop = {}

with open(sample_pop_file) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            sample_id = parts[0]
            population = parts[1]
            all_samples_ordered.append(sample_id)
            sample_to_pop[sample_id] = population

print(f"  Total samples: {len(all_samples_ordered)}")
print()

# Load population-specific samples
pop_samples = {}
for pop in POPULATIONS:
    pop_file = base_dir / f"data/grch38/samples_{pop}.txt"
    samples = set()
    with open(pop_file) as f:
        for line in f:
            sample_id = line.strip()
            if sample_id:
                samples.add(sample_id)
    pop_samples[pop] = samples
    print(f"  {pop}: {len(samples)} samples")

print()

# Create output directories and sample lists
for pop in POPULATIONS:
    output_dir = base_dir / f"data/grch38/hapbin_{pop}"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save sample list
    output_sample_list = output_dir / f"{pop}_samples.txt"
    with open(output_sample_list, 'w') as f:
        for sample in sorted(pop_samples[pop]):
            f.write(f"{sample}\n")
    print(f"Created output directory and sample list for {pop}: {output_dir}")

print()

# Process each chromosome
for chrom in CHROMOSOMES:
    print("=" * 80)
    print(f"Processing {chrom}")
    print("=" * 80)
    print()

    # Input files
    hap_file = data_dir / f"{chrom}.hapbin.hap"
    map_file = data_dir / f"{chrom}.hapbin.map"

    # Check if input files exist
    if not hap_file.exists():
        print(f"  WARNING: HAP file not found: {hap_file}")
        print(f"  Skipping {chrom}")
        print()
        continue

    if not map_file.exists():
        print(f"  WARNING: MAP file not found: {map_file}")
        print(f"  Skipping {chrom}")
        print()
        continue

    # Read first line to get total haplotypes
    with open(hap_file, 'r') as f:
        first_line = f.readline().strip().split()
        total_haplotypes = len(first_line)

    print(f"  Input: {hap_file}")
    print(f"  Size: {hap_file.stat().st_size / (1024**3):.2f} GB")
    print(f"  Total haplotypes: {total_haplotypes}")
    print()

    # Process each population
    for pop in POPULATIONS:
        print(f"  Processing {pop}...")

        output_dir = base_dir / f"data/grch38/hapbin_{pop}"
        output_hap = output_dir / f"{chrom}_{pop}.hapbin.hap"
        output_map = output_dir / f"{chrom}_{pop}.hapbin.map"

        # Check if output already exists
        if output_hap.exists() and output_map.exists():
            print(f"    ✓ Output already exists, skipping")
            continue

        # Build column indices for this population
        pop_column_indices = []
        for idx, sample_id in enumerate(all_samples_ordered):
            if sample_id in pop_samples[pop]:
                # 2 haplotypes per sample
                col1 = idx * 2
                col2 = idx * 2 + 1
                if col1 < total_haplotypes:
                    pop_column_indices.append(col1)
                if col2 < total_haplotypes:
                    pop_column_indices.append(col2)

        print(f"    Haplotype columns: {len(pop_column_indices)}")

        # Subset HAP file
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
                    print(f"      Processed {lines_processed:,} variants...", end='\r', file=sys.stderr)

        print(f"      Processed {lines_processed:,} variants")

        # Copy MAP file (same for all populations)
        with open(map_file, 'r') as fin, open(output_map, 'w') as fout:
            fout.write(fin.read())

        print(f"    ✓ Created: {output_hap.name}")
        print(f"      Size: {output_hap.stat().st_size / (1024**3):.2f} GB")
        print()

    print()

print("=" * 80)
print("Subsetting Complete!")
print("=" * 80)
print()
print("Output directories:")
for pop in POPULATIONS:
    output_dir = base_dir / f"data/grch38/hapbin_{pop}"
    print(f"  {pop}: {output_dir}")
print()
print("Next step:")
print("  Run iHS analysis using:")
print("  sbatch scripts/ThisWorks/02_ihs_chr1-22_populations.slurm")
print()
