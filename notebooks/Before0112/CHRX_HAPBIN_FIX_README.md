# Chromosome X hapbin Fix - Solution Documentation

## Problem Identified

**Root Cause:** hapbin cannot handle HAP files with variable field counts per line.

### Technical Details

1. **hapbin's assumption:** The `querySnpLengthAscii()` function reads the FIRST line of the HAP file and assumes ALL lines have the same number of fields:
   ```cpp
   std::getline(file, line);
   m_snpLength = (line.size()+1)/2;  // Uses FIRST line for entire file
   ```

2. **Your current format:**
   - PAR1 lines: 6,404 fields (2 per person × 3,202 people)
   - non-PAR lines: 4,468 fields (1 per male + 2 per female)
   - **This variable count breaks hapbin!**

3. **Result:** hapbin fails or produces invalid results for PAR1 variants, which is why iHS output starts at position 2,781,309 (right after PAR1 ends at 2,781,479).

## Solution

Create HAP files with **constant field count** (6,404 on every line) using IMPUTE2's chrX convention with dummy `-` characters for males in non-PAR regions.

### Output Format

**All lines have 6,404 fields (2 per individual):**

- **PAR regions:** Everyone gets 2 real haplotypes
  - Example: `0 1 0 0 1 1 0 0 ...` (all real alleles)

- **non-PAR regions:**
  - **Males:** real_haplotype + `-` (dummy)
  - **Females:** haplotype1 + haplotype2
  - Example: `0 - 0 1 1 - 0 0 ...` (males have `-` for second haplotype)

This follows IMPUTE2's chrX format specification that hapbin expects.

## Implementation

### Step 1: Create Fixed HAP File

Use the new script `fix_chrX_haplotypes_for_hapbin.py`:

```bash
python3 fix_chrX_haplotypes_for_hapbin.py \
    --vcf chr X.vcf.gz \
    --hap-in chrX.impute.hap.gz \
    --legend-in chrX.impute.legend.gz \
    --hap-out chrX.impute.hap.hapbin_format \
    --verbose
```

**What it does:**
1. Infers sex from VCF (males are homozygous in non-PAR)
2. Outputs HAP with constant 6,404 fields per line
3. Uses `-` dummy characters for male's second haplotype in non-PAR

### Step 2: Convert to hapbin Format

Use the existing `convert_to_hapbin.py` with the new HAP file:

```bash
python3 convert_to_hapbin.py \
    --hap chrX.impute.hap.hapbin_format.gz \
    --legend chrX.impute.legend.gz \
    --anc-fasta homo_sapiens_ancestor_X.fa \
    --chr X \
    --recomb-tsv chrX.decode.sexavg.cm.tsv \
    --no-header \
    --pos-col 1 \
    --cm-col 2 \
    --out-prefix chrX.hapbin \
    --verbose
```

### Step 3: Rerun iHS

```bash
pixi run ihsbin \
    --hap chrX.hapbin.hap \
    --map chrX.hapbin.map \
    --out ALL.chrX.ihs.tsv
```

**Expected result:** iHS output should now start from position ~288,060 (first variant in your data) instead of 2,781,309, **including PAR1 variants!**

## Quick Start

Run the all-in-one script:

```bash
sbatch /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/scripts/ThisWorks/reprocess_chrX_hapbin_v2.sh
```

This will:
1. Create hapbin-compatible HAP file
2. Convert to hapbin binary format
3. Verify consistent field counts
4. Backup old files automatically

## Verification

After running, verify the fix worked:

```bash
# Check field counts are consistent
head -1 chrX.hapbin.hap | awk '{print "PAR1:", NF}'
sed -n '100000p' chrX.hapbin.hap | awk '{print "non-PAR:", NF}'

# Both should show: 6404
```

Then check iHS output:

```bash
head -5 ALL.chrX.ihs.tsv
# Should start at position ~288,060 (in PAR1), not 2,781,309
```

## Key Files

### New Scripts
- `fix_chrX_haplotypes_for_hapbin.py` - Creates constant-field-count HAP
- `reprocess_chrX_hapbin_v2.sh` - All-in-one pipeline script

### Old Scripts (Don't Use)
- `fix_chrX_haplotypes.py` - Creates variable-field-count HAP (doesn't work with hapbin)
- `reprocess_chrX_hapbin.sh` - Uses old script

## Expected Outcomes

### Before Fix
- ❌ iHS output starts at position 2,781,309 (after PAR1)
- ❌ Missing ~86,882 PAR1 variants
- ❌ Cannot compare iHS vs RELATE in PAR1

### After Fix
- ✅ iHS output starts at position ~288,060 (includes PAR1)
- ✅ All ~86,882 PAR1 variants included
- ✅ Can compare iHS vs RELATE across entire X chromosome
- ✅ May reveal interesting PAR1 selection patterns!

## Important Notes

1. **PAR1 recombination:** Even with the fix, iHS results in PAR1 should be interpreted carefully because PAR1 recombines. The method assumes no recombination, so PAR1 results may be less reliable than non-PAR results.

2. **RELATE advantage:** RELATE models recombination explicitly, so its PAR1 signals are more reliable. The strong RELATE signals in PAR1 likely represent real selection that iHS struggles to detect.

3. **File sizes:** The new HAP file will be slightly larger because it includes dummy `-` characters, but the difference should be minimal when compressed.

## Questions?

If iHS still starts at 2,781,309 after the fix, check:
1. Are all lines in the new HAP file showing 6,404 fields?
2. Did the hapbin conversion succeed without errors?
3. Are you using the NEW hapbin files for iHS?

The hapbin source code review confirmed this is the issue - the fix should work!
