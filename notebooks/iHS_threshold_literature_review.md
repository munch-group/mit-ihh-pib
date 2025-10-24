# iHS Threshold Selection: Literature Review

## Summary of Key Papers

### 1. Voight et al. (2006) - Foundational Study
**Finding:** Standardized iHS values are approximately standard normal (mean ~0, variance ~1) under neutral expectations when binned by allele frequency.

**Recommendation:**
- |iHS| > 2 represents ~2 standard deviations from mean
- However, emphasize **clustering of high values** rather than single SNPs

**Key Quote:** "standardized iHS will have mean 0 and variance 1. Therefore, |iHS|>2 will represent 2 variances away..."

**Citation:** Voight BF, et al. (2006). A map of recent positive selection in the human genome. PLoS Biol.

---

### 2. Paul et al. (2024) - Recent Application
**Threshold Used:** |iHS| ≥ 2.5

**Additional Criteria:**
- Combine with p-value threshold: -log10(p) > 4
- "arbitrarily considered that any extreme iHS absolute value (i.e. |iHS| ≥ 2.5) corresponds to a positive signature of selection"

**Method:** Define candidate regions, not just SNPs

**Citation:** Paul S, et al. (2024). Genome-wide scan for selection signatures reveals novel insights into the adaptive capacity in local Italian cattle. BMC Genomics.

---

### 3. Salazar-Tortosa et al. (2023) - Caution on Fixed Thresholds
**Key Finding:** Empirical distributions of iHS have **heavier tails** than simulated neutral distributions

**Implication:** Fixed normal cutoffs may be overly liberal/optimistic

**Recommendation:** Use empirical percentiles rather than assuming standard normal

**Citation:** Salazar-Tortosa D, et al. (2023). Selection scans revisited: An unbiased approach to detecting selection from genomic data. Mol Ecol.

---

### 4. Qanbari et al. (2011) - Cattle Study
**Approach:**
- No universal fixed threshold
- Compute windows and examine distribution of values
- Note: Once allele is fixed, iHS loses power

**Implication:** Very high |iHS| might reflect sweep in progress, not completed

**Citation:** Qanbari S, et al. (2011). A genome-wide scan for signatures of recent selection in Holstein cattle. Anim Genet.

---

### 5. Liu et al. (2013) - Emphasis on Clustering
**Definition:** Candidate regions require "uncharacteristic **clustering** of SNPs with high iHS statistics"

**Rationale:** Extended haplotype signals more credible when multiple contiguous SNPs show high |iHS|

**Citation:** Liu S, et al. (2013). Population genomics reveal recent speciation and rapid evolutionary adaptation in polar bears. Cell.

---

## Recommended Threshold Strategies

### Strategy 1: Fixed Thresholds (Assuming N(0,1))
| Threshold | Percentile | Interpretation | Source |
|-----------|------------|----------------|--------|
| \|iHS\| ≥ 2.0 | ~97.7% | Liberal, 2 SD | Voight et al. 2006 |
| \|iHS\| ≥ 2.5 | ~99.4% | Moderate | Paul et al. 2024 |
| \|iHS\| ≥ 3.0 | ~99.87% | Stringent, 3 SD | - |

### Strategy 2: P-value Based
| Threshold | Equivalent \|iHS\| | Use Case |
|-----------|-------------------|----------|
| p < 0.01 | \|iHS\| ≥ 2.58 | Standard significance |
| -log10(p) > 4 | \|iHS\| ≥ 3.71 | Paul et al. 2024 |
| p < 1e-6 | \|iHS\| ≥ 4.75 | Very stringent |

### Strategy 3: Empirical Percentiles (Recommended)
- Top 1% of genome-wide |iHS| values
- Top 0.5% for more stringent
- Top 0.1% for strongest signals

**Advantage:** Accounts for actual data distribution, doesn't assume normality

---

## Critical Recommendations from Literature

### 1. Require Clustering, Not Single SNPs
- Multiple contiguous SNPs with high |iHS| (Liu et al. 2013)
- Sliding windows (e.g., 500 kb) with minimum number of significant SNPs
- Reduces false positives from random outliers

### 2. Combine Multiple Criteria
**Example from Paul et al. (2024):**
- |iHS| ≥ 2.5 **AND**
- -log10(p) > 4 **AND**
- ≥30 SNPs per 500 kb window

### 3. Check Distribution Normality
- Verify mean ≈ 0, SD ≈ 1 for Std iHS
- If distributions deviate significantly, use empirical percentiles
- X chromosome may have different distribution than autosomes

### 4. X Chromosome Special Considerations
- Smaller effective population size
- Different recombination rate
- Consider computing X-specific empirical percentiles
- May need X-specific thresholds

---

## Practical Workflow for This Project

### Step 1: Check Normality
```r
# Verify Std iHS ~ N(0,1)
- Check mean ≈ 0
- Check SD ≈ 1
- Compare X vs autosomes
```

### Step 2: Choose Threshold Strategy

**Option A: Fixed Threshold**
- Use |Std iHS| ≥ 2.5 (following Paul et al. 2024)
- Good if distributions are approximately normal

**Option B: Empirical Percentile**
- Use top 0.5% or 1% of genome-wide |iHS|
- Better if distributions deviate from normality
- Safer approach (Salazar-Tortosa et al. 2023)

**Option C: Combined Criteria**
- |Std iHS| ≥ 2.5 **AND** -log10(p) > 4
- Most stringent, fewest false positives

### Step 3: Define Candidate Regions
Following Paul et al. (2024) approach:
- Sliding windows (e.g., 500 kb, overlap 10 kb)
- Require ≥2-3 SNPs per window above threshold
- This ensures clustering, not single outliers

### Step 4: Chromosome-Specific Analysis
- Compute separate thresholds for X if needed
- Compare X vs autosome distributions
- Consider X-specific evolutionary dynamics

---

## Output for This Analysis

The updated script (`01_identify_ihs_candidates.R`) provides:

1. **Normality diagnostics**
   - Mean and SD of Std iHS per chromosome
   - Check deviations from N(0,1)

2. **Empirical percentiles**
   - 99th, 99.5th, 99.9th percentiles by chromosome
   - Allows comparison to fixed thresholds

3. **P-values**
   - Two-tailed p-values for each SNP
   - -log10(p) for significance assessment

4. **Multiple threshold levels**
   - Liberal (|iHS| ≥ 2.0)
   - Moderate (|iHS| ≥ 2.5)
   - Stringent (|iHS| ≥ 3.0)

5. **X vs autosome comparison**
   - Separate statistics
   - Enrichment test

---

## References

1. Voight BF, Kudaravalli S, Wen X, Pritchard JK (2006). A map of recent positive selection in the human genome. PLoS Biol 4(3):e72.

---

## Decision Matrix for Your Project

| Condition | Recommended Approach |
|-----------|---------------------|
| Distributions approximately normal (mean≈0, SD≈1) | Fixed threshold: \|iHS\| ≥ 2.5 |
| Distributions deviate from normality | Empirical: top 0.5-1% |
| X chromosome differs from autosomes | X-specific empirical percentiles |
| Want to minimize false positives | Combined: \|iHS\| ≥ 2.5 AND -log10(p) > 4 |

**Next step:** Run analysis and examine normality diagnostics to decide!
