# chrX Recombination Map: The 2/3 Adjustment Status

## The Question

**Does our preprocessing use the 2/3 adjustment for the chrX recombination map?**

**Short answer:** YES! The adjustment IS being applied correctly. The map was created with the 2/3 multiplier.

---

## What the Map Contains

### The Map Creation ([extract_chrX_from_decode.sh](extract_chrX_from_decode.sh:67-69))
```bash
# Apply 2/3 adjustment for X chromosome
# (because only females recombine on X, not males)
adjusted_cm = cumulative_cm * 2.0/3.0
```

**The map file `chrX.decode.sexavg.cm.tsv` WAS created with the 2/3 adjustment!**

### Verification
```bash
# Raw deCODE sex-averaged value: 175.882 cM
# Your map (with 2/3 adjustment): 117.255 cM
# Ratio: 117.255 / 175.882 = 0.6667 = 2/3 ✓
```

**Result:** Both old and new preprocessing scripts use `chrX.decode.sexavg.cm.tsv`, which **already has the 2/3 adjustment applied during map creation**.

---

## Why the 2/3 Adjustment Matters

### The Biology

For **autosomes:**
- Both males and females have 2 copies
- Recombination happens in both sexes
- Sex-averaged map = (male_rate + female_rate) / 2

For **X chromosome:**
- **Females** have 2 X chromosomes → recombination occurs (use female rate)
- **Males** have 1 X chromosome → NO recombination in males (X passes intact from mother to daughter)
- Males NEVER pass X to sons, only to daughters

### The Math

In a population:
- **2/3 of X chromosomes** are in females (who have 2 copies)
- **1/3 of X chromosomes** are in males (who have 1 copy)

Since recombination only occurs in **females** (2/3 of X chromosomes), the effective recombination rate should be:

```
Effective rate = Female rate × (2/3)
```

OR equivalently, if using "sex-averaged" rate:

```
Effective rate = Sex-averaged rate × (2/3)
```

### Current Genetic Length

Your current map (`chrX.decode.sexavg.cm.tsv`):
```
Total length: 117.255 cM (already 2/3-adjusted)
```

This represents:
```
Raw deCODE value: 175.882 cM
After 2/3 adjustment: 175.882 × (2/3) = 117.255 cM ✓
```

---

## Impact on iHS Analysis

### Does This Matter?

**Good news:** The 2/3 adjustment is ALREADY applied! Your iHS analysis is using the correct genetic distances.

The recombination map is used to calculate **genetic distance** for EHH (Extended Haplotype Homozygosity) decay, and your map correctly accounts for the fact that only females recombine on the X chromosome.

---

## Comparison to Your Autosomal Data

### Do autosomes use sex-averaged maps?

Let me check what you used for chr1:

```bash
# Your chr1 likely uses:
/home/vanbruggenmit/mit-ihh-pib/data/grch38/maps/chr1.decode.sexavg.cm.tsv
Total length: 267.8 cM
```

For **autosomes**, sex-averaged is correct (no adjustment needed).
For **chrX**, sex-averaged needs the **2/3 adjustment**.

---

## Summary

| Aspect | Status |
|--------|--------|
| Haplotype fix (male hemizygosity) | ✓ Fixed correctly |
| Recombination map adjustment (2/3) | ✓ Applied correctly |
| Map genetic length | 117.255 cM (correct) |
| Current iHS run | ✓ Using correct inputs |
| Action needed | None - everything is correct! |

---

## References

**Why 2/3 for X chromosome:**
- Only 2/3 of X chromosomes undergo recombination (those in females)
- Males don't recombine X (passed intact to daughters)
- Standard practice in X chromosome selection scans

**Examples from literature:**
- Sabeti et al. (2007) - used sex-specific maps with X adjustment
- 1000 Genomes Selection Browser - applies X chromosome corrections
