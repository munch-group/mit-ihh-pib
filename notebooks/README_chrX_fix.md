# Chromosome X Haplotype Fix

## Problem

The original chrX conversion (`vcf2hap_chrX.sh` + `hapbin_chrX.sh`) treated all individuals as diploid, resulting in 6404 haplotypes (3202 samples × 2 haplotypes each).

However, **males are hemizygous for the X chromosome** in non-PAR regions:
- Males should contribute **1 haplotype** each
- Females should contribute **2 haplotypes** each

This caused `ihsbin` to report "Chromosomes per SNP: 6404" for chrX, when it should be ~3200.

## Solution

Two new scripts fix this issue:

### 1. `fix_chrX_haplotypes.py`

Python script that:
- Infers sample sex from VCF genotypes (males are homozygous in non-PAR X regions)
- Removes duplicate haplotypes for males in non-PAR regions
- Keeps both haplotypes for females everywhere
- Keeps both haplotypes for everyone in PAR regions

**PAR Boundaries (GRCh38):**
- PAR1: X:10,001-2,781,479
- PAR2: X:155,701,383-156,030,895

### 2. `reprocess_chrX_hapbin.sh`

SLURM script that:
1. Runs `fix_chrX_haplotypes.py` on the existing HAP file
2. Converts the fixed HAP to hapbin format
3. Backs up old hapbin files
4. Generates new corrected chrX.hapbin.hap and chrX.hapbin.map

## Usage

### Option 1: Run via SLURM (Recommended)

```bash
cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib
sbatch scripts/ThisWorks/reprocess_chrX_hapbin.sh
```

This will:
- Process the data (takes ~2-4 hours)
- Output fixed files to `data/grch38/hapbin/`
- Backup old files with `.backup` extension

### Option 2: Run manually

```bash
# Fix the haplotypes
python3 scripts/ThisWorks/fix_chrX_haplotypes.py \
    --vcf data/grch38/raw/chrX/chrX.vcf.gz \
    --hap-in data/grch38/work/chrX/chrX.impute.hap.gz \
    --legend-in data/grch38/work/chrX/chrX.impute.legend.gz \
    --hap-out data/grch38/work/chrX/chrX.impute.hap.fixed \
    --verbose

# Compress
gzip data/grch38/work/chrX/chrX.impute.hap.fixed

# Then run hapbin conversion (see reprocess_chrX_hapbin.sh for full command)
```

## Expected Results

**Before fix:**
- Chromosomes per SNP: 6404
- File size: ~24GB

**After fix:**
- Chromosomes per SNP: ~3200
- File size: ~12GB
- Males in non-PAR: 1 haplotype each
- Females everywhere: 2 haplotypes each
- PAR regions: 2 haplotypes for everyone

## Rerun iHS

After fixing, rerun the iHS calculation:

```bash
# Delete old output
rm results/ihs/ALL.chrX.ihs.tsv

# Rerun with job array (or manually for just chrX)
sbatch scripts/ThisWorks/ihs_all.slurm
```

Or manually for just chrX:

```bash
pixi run ihsbin \
    --hap data/grch38/hapbin/chrX.hapbin.hap \
    --map data/grch38/hapbin/chrX.hapbin.map \
    --out results/ihs/ALL.chrX.ihs.tsv
```

## Verification

Check that the fix worked:

```bash
# Should show ~3200 haplotypes per line (not 6404)
head -1 data/grch38/hapbin/chrX.hapbin.hap | wc -w

# Check iHS output
grep "Chromosomes per SNP" logs/ihs-*_22.out
```

## Technical Details

### Sex Inference

The script infers sex by analyzing heterozygosity in non-PAR X regions:
- **Males**: <5% heterozygosity (mostly homozygous 0|0 or 1|1)
- **Females**: 30-50% heterozygosity (normal diploid variation)

This is reliable because males only have one X chromosome copy in non-PAR regions.

### Files Modified

- **Input (unchanged):**
  - `data/grch38/raw/chrX/chrX.vcf.gz` - Original VCF
  - `data/grch38/work/chrX/chrX.impute.hap.gz` - Original HAP (6404 haplotypes)
  - `data/grch38/work/chrX/chrX.impute.legend.gz` - Legend file

- **New files created:**
  - `data/grch38/work/chrX/chrX.impute.hap.fixed.gz` - Fixed HAP (~3200 haplotypes)
  - `data/grch38/hapbin/chrX.hapbin.hap` - Final hapbin format (corrected)
  - `data/grch38/hapbin/chrX.hapbin.map` - Final hapbin map (unchanged)

- **Backups:**
  - `data/grch38/hapbin/chrX.hapbin.hap.backup` - Old (incorrect) hapbin file
  - `data/grch38/hapbin/chrX.hapbin.map.backup` - Old map file

## Author

Created to fix chrX hemizygosity issue in iHS analysis pipeline
Date: October 2025
