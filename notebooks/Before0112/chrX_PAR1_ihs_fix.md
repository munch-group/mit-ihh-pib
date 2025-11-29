# Fixing Missing PAR1 Region in chrX iHS Analysis

## Problem Summary

iHS analysis on chromosome X was missing all variants in the PAR1 (Pseudoautosomal Region 1) region (X:10,001-2,781,479), while RELATE analysis showed strong selection signals in this region. The iHS output file started at position 2,781,309, completely skipping PAR1.

## Root Cause Analysis

### Investigation Steps

1. **Compared RELATE vs iHS outputs**
   - RELATE showed selection signals primarily in PAR1 region (0-2.7 Mb)
   - iHS output started at position 2,781,309 (after PAR1 ends)

2. **Examined hapbin source code** (`/home/vanbruggenmit/mit-ihh-pib/hapbin/src/hapmap.cpp`)
   - Found critical limitation in lines 143-144:
   ```cpp
   std::getline(file, line);
   m_snpLength = (line.size()+1)/2;  // Assumes FIRST line determines haplotype count
   ```
   - hapbin assumes ALL lines in the HAP file have the same number of haplotypes as the first line

3. **Analyzed existing hapbin files**
   - Old chrX.hapbin.hap file had variable haplotype counts:
     - PAR1 variants (line 1): 6,404 haplotypes (diploid for all 3,202 individuals)
     - non-PAR variants (line 100,000): 4,468 haplotypes (haploid males + diploid females)
   - This variable count violated hapbin's assumption, causing it to fail on PAR1 variants

4. **Checked original IMPUTE2 files**
   - Original `chrX.impute.hap.gz` had **consistent 6,404 fields on every line**
   - This is the correct format: diploid representation for all individuals across entire chromosome
   - The previous hapbin conversion incorrectly created variable haplotype counts

## The Fix

### Solution Overview

Convert the original IMPUTE2 HAP file (which already has consistent field counts) directly to hapbin format, rather than trying to fix an already-broken hapbin file.

### Technical Challenges Resolved

1. **pysam dependency missing**
   - Problem: `convert_to_hapbin.py` required pysam, but it wasn't installed
   - Investigation: pysam was listed in `pyproject.toml` but not in the active `pixi.toml`
   - Solution: Added bioconda channel and pysam to `/faststorage/project/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/pixi.toml`
   ```toml
   channels = ["conda-forge", "munch-group", "bioconda"]
   dependencies:
     pysam = ">=0.23.3"
   ```

2. **Recombination map parsing error**
   - Problem: Map file loaded 0 rows
   - Investigation: Recombination map file is space-separated, not tab-separated
   - Solution: Added `--sep ' '` argument to conversion script

3. **SIGPIPE errors**
   - Problem: `zcat | head -1` with `set -o pipefail` caused exit code 13
   - Solution: Temporarily disable pipefail around the command:
   ```bash
   set +o pipefail
   FIELDS_LINE1=$(zcat "$HAP_IN" | head -1 | awk '{print NF}')
   set -o pipefail
   ```

### Implementation

Created conversion script: `/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/scripts/ThisWorks/convert_chrX_original_to_hapbin.sh`

Key features:
- Verifies input HAP file has 6,404 fields (diploid for all individuals)
- Uses original `chrX.impute.hap.gz` file
- Converts to hapbin format with consistent haplotype counts
- Backs up old hapbin files before overwriting

## Results

### Old hapbin file
- Started at position: 2,781,309 (after PAR1)
- Missing: All PAR1 variants

### New hapbin file
- Starts at position: **288,060** (within PAR1)
- Contains: **2,578,259 SNPs**
- Size: 24 GB
- Format: 6,404 fields per line (consistent throughout)
- Coverage: Includes full PAR1 region where RELATE detected selection signals

## Next Steps

1. Rerun iHS analysis using new hapbin files:
```bash
pixi run ihsbin --hap /home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin/chrX.hapbin.hap \
                --map /home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin/chrX.hapbin.map \
                --out results/ihs/ALL.chrX.ihs.tsv
```

2. Compare iHS and RELATE signals across entire chrX including PAR1

3. Investigate if iHS detects selection in PAR1 region (where RELATE shows strong signals)

## Key Lessons

1. **Source code analysis is critical** - Reading hapbin's source code revealed the fundamental constraint that wasn't documented elsewhere

2. **Verify assumptions at each processing step** - The variable haplotype count issue was introduced during an earlier conversion step

3. **Original data format matters** - The IMPUTE2 format with consistent diploid representation was actually correct for hapbin, despite seeming counterintuitive for male hemizygosity

4. **Different methods have different assumptions** - iHS assumes no recombination, which may affect its ability to detect selection in recombining regions like PAR1 compared to RELATE

## Files Modified

- `/faststorage/project/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/pixi.toml` - Added pysam dependency
- `/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/scripts/ThisWorks/convert_chrX_original_to_hapbin.sh` - Created conversion script
- `/home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin/chrX.hapbin.hap` - Regenerated hapbin file
- `/home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin/chrX.hapbin.map` - Regenerated map file
