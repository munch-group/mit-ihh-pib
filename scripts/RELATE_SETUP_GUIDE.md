# Relate Setup and Conversion Guide

This guide walks you through converting your VCF data to Relate format to verify your iHS results.

## Overview

Relate is a tool for inferring genome-wide genealogies and detecting selection. It requires:
- `.haps` files (haplotype data)
- `.sample` files (sample metadata)
- `.dist` files (genetic map/distances)

## Step 1: Install Relate

```bash
sbatch scripts/install_relate.sh
```

This will:
- Download Relate v1.2.1
- Install to: `software/relate/`
- Provide instructions to add to PATH

After installation, add to your PATH:
```bash
export PATH="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/software/relate/bin:$PATH"
```

Or add permanently to `~/.bashrc`:
```bash
echo 'export PATH="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/software/relate/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

## Step 2: Convert VCF to Relate Format

### For Autosomes (chromosomes 1-22):

```bash
# Convert all autosomes (runs as array job)
sbatch scripts/vcf_to_relate.sh

# Or convert a single chromosome for testing (e.g., chr22)
sbatch --array=22 scripts/vcf_to_relate.sh
```

### For Chromosome X:

```bash
sbatch scripts/vcf_to_relate_chrX.sh
```

**Important:** Chromosome X requires special handling:
- Males are hemizygous for non-PAR regions
- Use the 2/3 adjusted recombination map
- Relate handles hemizygosity if sex is specified in `.sample` file

## Step 3: Verify Conversion

Check that output files were created:

```bash
# For chr22 (test chromosome)
ls -lh /faststorage/project/mit-ihh-pib/data/grch38/relate/chr22.*

# Expected files:
# - chr22.haps (haplotype matrix)
# - chr22.sample (sample information)
```

## Step 4: Prepare Genetic Maps

Relate needs genetic maps in a specific format. You already have deCODE maps at:
- `/faststorage/project/mit-ihh-pib/data/grch38/maps/`

The maps need to be converted to Relate format (position in bp, genetic position in cM).

## Next Steps

After conversion, you can:

1. **Run Relate inference** to build genealogies
2. **Detect selection** using Relate's selection detection methods
3. **Compare with iHS results** to verify your findings

### Key Differences: Relate vs iHS

- **iHS**: Detects recent selection via extended haplotype homozygosity
- **Relate**: Infers full genealogies and detects selection via:
  - Branch length anomalies
  - Coalescence rate changes
  - Can detect older selection signals

Both methods should converge on strong selection signals, providing validation.

## Troubleshooting

### If RelateFileFormats not found:
```bash
# Check installation
ls -la /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/software/relate/bin/

# Verify PATH
echo $PATH

# Test manually
/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/software/relate/bin/RelateFileFormats --help
```

### If VCF not found:
Check that VCF files exist:
```bash
ls -la /faststorage/project/mit-ihh-pib/data/grch38/raw/chr*/chr*.vcf.gz
```

## File Locations

- **Input VCFs**: `/faststorage/project/mit-ihh-pib/data/grch38/raw/chr*/`
- **Output Relate files**: `/faststorage/project/mit-ihh-pib/data/grch38/relate/`
- **Genetic maps**: `/faststorage/project/mit-ihh-pib/data/grch38/maps/`
- **Scripts**: `/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/scripts/`
- **Logs**: `/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/`

## References

- Relate documentation: https://myersgroup.github.io/relate/
- Relate paper: Speidel et al. (2019) Nature Genetics
