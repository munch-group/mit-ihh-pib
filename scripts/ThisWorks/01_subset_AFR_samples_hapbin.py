#!/usr/bin/env python3
"""
Subset chromosome X hapbin data to AFR samples only for population-specific iHS analysis

Input: chrX.hapbin.hap.backup_20251122_185513 (ALL populations, hapbin format)
Output: chrX_AFR.hapbin.hap (AFR samples only, hapbin format)

This script creates AFR-only versions to identify AFR-specific selection signals.
"""

from pathlib import Path
import sys

# Paths
base_dir = Path("/home/vanbruggenmit/mit-ihh-pib")
data_dir = base_dir / "data/grch38/hapbin"
output_dir = base_dir / "data/grch38/hapbin_AFR"
output_dir.mkdir(parents=True, exist_ok=True)

# Sample files
sample_pop_file = base_dir / "data/grch38/sample_population.txt"
afr_samples_file = base_dir / "data/grch38/samples_AFR.txt"

# Input files (October backup - what was actually used)
hap_file = data_dir / "chrX.hapbin.hap.backup_20251122_185513"
map_file = data_dir / "chrX.hapbin.map.backup_20251122_185513"

# Output files
output_hap = output_dir / "chrX_AFR.hapbin.hap"
output_map = output_dir / "chrX_AFR.hapbin.map"
output_sample_list = output_dir / "AFR_samples.txt"

print("=" * 80)
print("Subsetting Chromosome X Hapbin Data to AFR Samples")
print("=" * 80)
print()
print("Using October 2024 backup files (what was used for original iHS analysis)")
print()

# Load AFR sample IDs
print("Loading AFR sample IDs...")
afr_samples = set()
with open(afr_samples_file) as f:
    for line in f:
        sample_id = line.strip()
        if sample_id:
            afr_samples.add(sample_id)

print(f"  Found {len(afr_samples)} AFR samples")
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
print("Determining AFR sample columns in HAP file...")

# Read first line to check format
with open(hap_file, 'r') as f:
    first_line = f.readline().strip().split()
    total_haplotypes = len(first_line)

print(f"  Total haplotypes in file: {total_haplotypes}")

# Build column indices for AFR samples
afr_column_indices = []
for idx, sample_id in enumerate(all_samples_ordered):
    if sample_id in afr_samples:
        # Assuming 2 haplotypes per sample
        col1 = idx * 2
        col2 = idx * 2 + 1
        if col1 < total_haplotypes:
            afr_column_indices.append(col1)
        if col2 < total_haplotypes:
            afr_column_indices.append(col2)

print(f"  AFR haplotype columns to extract: {len(afr_column_indices)}")
print(f"  Expected for 504 AFR samples: ~504-1008 (depending on male/female)")
print()

# Save AFR sample list
with open(output_sample_list, 'w') as f:
    for sample in sorted(afr_samples):
        f.write(f"{sample}\n")
print(f"Saved AFR sample list: {output_sample_list}")

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

        # Extract only AFR columns
        afr_fields = [fields[i] for i in afr_column_indices if i < len(fields)]

        # Write AFR-only line
        fout.write(' '.join(afr_fields) + '\n')

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
        print(f"  Expected: ~{len(afr_column_indices)}")

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
print("  Run iHS on AFR-only data using:")
print("  sbatch scripts/ThisWorks/02_ihs_chrX_AFR.slurm")
print()
print("Note:")
print("  This AFR-only dataset will allow you to identify AFR-specific selection")
print("  signals, which can then be validated with XP-EHH (AFR vs EAS/EUR)")
print()
