#!/usr/bin/env python3
"""
Subset chromosome X hapbin data to EAS samples only for population-specific iHS analysis

Input: chrX.hapbin.hap.backup_20251122_185513 (ALL populations, hapbin format)
Output: chrX_EAS.hapbin.hap (EAS samples only, hapbin format)

This script creates EAS-only versions to identify EAS-specific selection signals.
"""

from pathlib import Path
import sys

# Paths
base_dir = Path("/home/vanbruggenmit/mit-ihh-pib")
data_dir = base_dir / "data/grch38/hapbin"
output_dir = base_dir / "data/grch38/hapbin_EAS"
output_dir.mkdir(parents=True, exist_ok=True)

# Sample files
sample_pop_file = base_dir / "data/grch38/sample_population.txt"
eas_samples_file = base_dir / "data/grch38/samples_EAS.txt"

# Input files (October backup - what was actually used)
hap_file = data_dir / "chrX.hapbin.hap.backup_20251122_185513"
map_file = data_dir / "chrX.hapbin.map.backup_20251122_185513"

# Output files
output_hap = output_dir / "chrX_EAS.hapbin.hap"
output_map = output_dir / "chrX_EAS.hapbin.map"
output_sample_list = output_dir / "EAS_samples.txt"

print("=" * 80)
print("Subsetting Chromosome X Hapbin Data to EAS Samples")
print("=" * 80)
print()
print("Using October 2024 backup files (what was used for original iHS analysis)")
print()

# Load EAS sample IDs
print("Loading EAS sample IDs...")
eas_samples = set()
with open(eas_samples_file) as f:
    for line in f:
        sample_id = line.strip()
        if sample_id:
            eas_samples.add(sample_id)

print(f"  Found {len(eas_samples)} EAS samples")
print()

# Load sample-to-population mapping to get sample order
print("Loading sample population mapping...")
all_samples_ordered = []
with open(sample_pop_file) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            all_samples_ordered.append(parts[0])

print(f"  Total samples in order: {len(all_samples_ordered)}")

# Determine which columns to keep
print()
print("Determining EAS sample columns in HAP file...")

# Read first line to check format
with open(hap_file, 'r') as f:
    first_line = f.readline().strip().split()
    total_haplotypes = len(first_line)

print(f"  Total haplotypes in file: {total_haplotypes}")

# Build column indices for EAS samples
eas_column_indices = []
for idx, sample_id in enumerate(all_samples_ordered):
    if sample_id in eas_samples:
        # Assuming 2 haplotypes per sample
        col1 = idx * 2
        col2 = idx * 2 + 1
        if col1 < total_haplotypes:
            eas_column_indices.append(col1)
        if col2 < total_haplotypes:
            eas_column_indices.append(col2)

print(f"  EAS haplotype columns to extract: {len(eas_column_indices)}")
print(f"  Expected for 504 EAS samples: ~504-1008 (depending on male/female)")
print()

# Save EAS sample list
with open(output_sample_list, 'w') as f:
    for sample in sorted(eas_samples):
        f.write(f"{sample}\n")
print(f"Saved EAS sample list: {output_sample_list}")

# Subset HAP file
print()
print("Subsetting HAP file (this will take several minutes)...")
print(f"  Input: {hap_file}")
print(f"  Input size: {hap_file.stat().st_size / (1024**3):.2f} GB")
print(f"  Output: {output_hap}")

lines_processed = 0
with open(hap_file, 'r') as fin, open(output_hap, 'w') as fout:
    for line in fin:
        fields = line.strip().split()

        # Extract only EAS columns
        eas_fields = [fields[i] for i in eas_column_indices if i < len(fields)]

        # Write EAS-only line
        fout.write(' '.join(eas_fields) + '\n')

        lines_processed += 1
        if lines_processed % 10000 == 0:
            print(f"    Processed {lines_processed:,} variants...", end='\r', file=sys.stderr)

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

# Verify output
if output_hap.exists() and output_map.exists():
    with open(output_hap, 'r') as f:
        first_line_out = f.readline().strip().split()
        print(f"Verification:")
        print(f"  Output HAP haplotypes: {len(first_line_out)}")
        print(f"  Expected: ~{len(eas_column_indices)}")

    with open(output_map, 'r') as f:
        map_lines_out = sum(1 for _ in f)
        print(f"  Output MAP lines: {map_lines_out:,}")
        print(f"  Match with HAP: {map_lines_out == lines_processed}")

print()
print("=" * 80)
print("Subsetting Complete!")
print("=" * 80)
print()
print("Output files:")
print(f"  HAP: {output_hap}")
print(f"  MAP: {output_map}")
print(f"  Samples: {output_sample_list}")
print()
print("Next step:")
print("  Run iHS on EAS-only data using:")
print("  sbatch scripts/ThisWorks/02_ihs_chrX_EAS.slurm")
print()
print("Note:")
print("  This EAS-only dataset will allow you to identify EAS-specific selection")
print("  signals, which can then be validated with XP-EHH (EAS vs EUR/AFR)")
print("  Look for known signals like EDAR pathway (hair/tooth morphology)")
print()
