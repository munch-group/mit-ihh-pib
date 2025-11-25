# Comparison of Two chrX iHS Result Files

## Overview

Compared two iHS result files for chromosome X to understand why they show different values for the same SNPs and to count significant signals.

## Files Compared

### File 1: `ihs_fixed_par1_paraproc/ALL.chrX.ihs.tsv`
- **SNPs**: 196,860
- **Genomic range**: 31,516,278 - 154,690,904
- **Generated from**: Parallel computation (5 genomic regions)
- **Signals with |Std iHS| ≥ 2.0**: 10,973 (5.57%)
  - Positive (Std iHS ≥ 2.0): 4,487
  - Negative (Std iHS ≤ -2.0): 6,486

### File 2: `ihs/ALL.chrX.ihs.tsv`
- **SNPs**: 428,414
- **Genomic range**: 2,781,309 - 155,705,325
- **Generated from**: Single-pass computation
- **Signals with |Std iHS| ≥ 2.0**: 20,507 (4.79%)
  - Positive (Std iHS ≥ 2.0): 19,446
  - Negative (Std iHS ≤ -2.0): 1,061

## SNP Overlap Analysis

- **SNPs in both files**: 176,437 (89.63% of File 1)
- **SNPs only in File 1**: 20,423 (10.37%)
- **SNPs only in File 2**: 251,977 (58.82%)

File 2 starts ~29 million bp earlier, covering more of the PAR1 region.

## Comparison of Overlapping SNPs

For the 176,437 SNPs present in both files:
- **Correlation of Std iHS values**: r = 0.2888
- **Mean difference** (File 1 - File 2): -0.007
- **Standard deviation of difference**: 1.18
- **Maximum absolute difference**: 7.83

The low correlation indicates the two analyses are giving  different results for the same SNPs.

## Root Cause: Region 1 Failure

### What Happened

The parallel processing approach split chromosome X into 5 regions. Region 1 (positions 288,060 - 31,499,955) completely failed to produce valid iHS scores.

From the log file (`ihs_chrX_region1_13088628.out`):
```
Input: 581,406 SNPs
# valid loci: 0
# loci with MAF <= 0.05: 468,092 (80.5%)
# loci which reached the end of the chromosome: 99,974 (17.2%)
```

### Consequences

1. Region 1 output file contains only the header (33 bytes)
2. The merged file starts at region 2 (position 31,516,278)
3. **Missing ~29 million bp** including most/all of PAR1 region

### Why Region 1 Failed

- **80.5% of SNPs filtered** due to MAF ≤ 0.05
- **17.2% of SNPs** hit "end of chromosome" errors (haplotype extension failed)
- PAR1 recombines, which violates iHS assumptions and breaks LD structure
- Combined result: 0 valid iHS scores in this region

## Why the Std iHS Values Differ

### Different SNP Sets for Standardization

iHS standardization converts raw iHS to Std iHS via z-score transformation:
```
Std iHS = (iHS - mean) / SD
```

The two files computed mean and SD from **different sets of SNPs**:
- File 1: Based on regions 2-5 only (excludes region 1)
- File 2: Based on different genomic coverage starting at 2.7M

Same raw iHS value → Different Std iHS scores depending on the reference distribution.

### Distribution Differences

The signal distributions are very different:
- **File 1**: More negative signals (6,486 neg vs 4,487 pos)
- **File 2**: Heavily positive signals (19,446 pos vs 1,061 neg)

This indicates the analyses differ fundamentally, not just in standardization.

## Missing Coverage

### File 1 Missing Regions
- Positions 288,060 - 31,516,277 (≈29.2 Mb)
- Includes PAR1 region where RELATE detected strong selection signals

### File 2 Missing Regions
- Positions 288,060 - 2,781,308 (≈2.5 Mb)
- Partial PAR1 coverage, but still missing the earliest positions

## Findings

### 1. The Files Are Not Comparable
The correlation of 0.29 for overlapping SNPs means these are fundamentally different analyses, not just different standardizations of the same data.

### 2. Both Files Miss PAR1
Neither file includes the full PAR1 region:
- File 1: Completely missing (starts at 31.5M)
- File 2: Mostly missing (starts at 2.7M, PAR1 ends at 2.78M)

### 3. Different Input Processing
Based on the documentation in `notebooks/`:
- File 1: Generated from "fixed" hapbin files (after PAR1 format fix)
- File 2: Generated from earlier hapbin conversion
- Different hapbin conversions → Different input data → Different results

### 4. Parallel Processing Introduced Issues
The region splitting approach caused:
- Complete loss of region 1 data
- Different standardization (smaller SNP set)
- Cannot recover PAR1 analysis from this approach

## Answer to Original Question

**How many signals in `ihs_fixed_par1_paraproc/ALL.chrX.ihs.tsv` have Std iHS ≥ 2.0?**

**Total signals with |Std iHS| ≥ 2.0: 10,973**
- Positive signals (Std iHS ≥ 2.0): 4,487
- Negative signals (Std iHS ≤ -2.0): 6,486

Note: These signals are from positions 31,516,278 onward (non-PAR region only). PAR1 region is completely missing from this file.
