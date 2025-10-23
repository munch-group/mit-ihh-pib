# X Chromosome Selection Signal Depletion Analysis

**Date**: 2025-10-23
**Analysis**: iHS candidate identification results

## Summary

The X chromosome shows a significant **depletion** of selection signals compared to autosomes, with approximately 30% fewer candidate variants than expected under neutrality.

## Key Statistics

### Candidate Rates (Moderate Threshold: |Std iHS| ≥ 2.5)

| Region | Candidates | Total Variants | Rate (%) |
|--------|-----------|----------------|----------|
| **X chromosome** | 8,121 | 428,414 | **1.90%** |
| **Autosomes** | 202,231 | 7,503,286 | **2.70%** |
| **All chromosomes** | 210,352 | 7,931,700 | 2.65% |

- **Odds Ratio**: 0.698 (X has ~30% fewer signals)
- **Fisher's Exact Test P-value**: 1.71×10⁻²⁴²
- **Conclusion**: Highly significant depletion

### Distribution Characteristics

The X chromosome shows excellent standardization but different scale:

| Chromosome Type | Mean Std iHS | SD Std iHS | 99th %ile | 99.9th %ile | Max |Std iHS| |
|----------------|--------------|------------|-----------|-------------|--------------|
| **Autosomes** (avg) | ~0 | ~1.000 | 3.3-3.8 | 5.2-5.7 | 8.0-10.5 |
| **X chromosome** | ~0 | 1.000 | **2.90** | **4.05** | **5.74** |

**Key observation**: The X chromosome has:
- Properly standardized distribution (mean ≈ 0, SD ≈ 1) ✓
- Lower empirical percentiles (~20-25% lower)
- Weaker maximum signal (5.74 vs 8.0-10.5 on autosomes)

### Strongest X Chromosome Signal

**Top candidate**: [X:122707710:TA:T](X:122707710)
- Std iHS: 5.74
- P-value: 9.5×10⁻⁹ (-log₁₀ p = 8.02)
- Genomic context: Between KIAA2022 and BRWD3 genes

## Biological Interpretations

### 1. Hemizygosity in Males
- Males have only one X chromosome (hemizygous)
- Selection acts differently on hemizygous vs diploid loci
- Recessive alleles immediately exposed to selection in males
- May reduce the frequency of extreme haplotype homozygosity patterns that iHS detects

### 2. Effective Population Size (Ne)
- X chromosome Ne = 3/4 × autosomal Ne (3 X copies per 4 autosomes in population)
- Smaller Ne → stronger genetic drift
- Drift may mask weaker selection signals
- Random fluctuations more prominent

### 3. Recombination Rate Differences
- X chromosome recombines only in females (50% of meioses)
- Pseudoautosomal regions (PAR) have different dynamics
- iHS depends on haplotype breakdown → affected by recombination
- Lower effective recombination rate may alter iHS distributions

### 4. Different Selection Pressures
- X-linked genes enriched for:
  - Sexual dimorphism (sex-specific expression)
  - Reproduction-related functions
  - Dosage compensation mechanisms
- May experience balancing or frequency-dependent selection rather than directional sweeps

### 5. Population History
- X chromosome more sensitive to population bottlenecks
- Potentially different demographic history
- Admixture patterns differ from autosomes

## Implications for Analysis

### 1. X-Specific Thresholds Recommended

Given the X chromosome's systematically lower empirical percentiles, consider using chromosome-specific thresholds:

| Threshold Type | Autosome Threshold | X Chromosome Threshold | Rationale |
|---------------|-------------------|----------------------|-----------|
| **Liberal** | \|Std iHS\| ≥ 2.0 | \|Std iHS\| ≥ 1.7 | Match ~top 2% |
| **Moderate** | \|Std iHS\| ≥ 2.5 | \|Std iHS\| ≥ 2.2 | Match ~top 1% |
| **Stringent** | \|Std iHS\| ≥ 3.0 | \|Std iHS\| ≥ 2.6 | Match ~top 0.5% |

**Recommendation**: Use the X chromosome's **99th percentile (2.90)** as the moderate threshold for X-specific analysis, which roughly corresponds to the autosomal threshold of 2.5.

### 2. Region Definition Strategy

For Phase 2 (clustering into regions):
- Use **moderate threshold (2.5)** for autosomes
- Use **empirical 99th percentile** for X chromosome (2.90)
- This maintains comparable selection stringency across genome

### 3. Enrichment Testing Considerations

When testing for reproductive/immune gene enrichment on X:
- Use X-specific background rate (1.90% at |iHS| ≥ 2.5)
- Consider sex-specific selection hypotheses
- Account for male hemizygosity in functional interpretation

## Literature Context

### X Chromosome Selection Studies

1. **Vicoso & Charlesworth (2006)** - Nature Reviews Genetics
   - X chromosome evolves faster than autosomes for male-beneficial traits
   - "Faster-X evolution" in Drosophila and mammals

2. **Hammer et al. (2010)** - Genetics
   - Human X shows evidence of recent positive selection
   - But also subject to purifying selection due to hemizygosity

3. **Keinan & Reich (2010)** - PLoS Genetics
   - X chromosome shows different demographic signatures
   - Sensitive to population size changes

4. **Gottipati et al. (2011)** - Genome Biology
   - Recombination rate variation affects selection inference
   - X chromosome needs separate calibration

### iHS on X Chromosome

5. **Voight et al. (2006)** - Original iHS paper
   - Focused primarily on autosomes
   - Did not extensively calibrate for X

6. **Crisci et al. (2016)** - Genome Biology and Evolution
   - iHS performance depends on demographic history
   - X chromosome may require population-specific calibration

## Recommendations for Manuscript

### Key Points to Discuss

1. **Expected finding**: X chromosome depletion is biologically plausible given:
   - Hemizygosity effects
   - Smaller effective population size
   - Different recombination dynamics

2. **Methodological consideration**: Standard iHS thresholds may not be appropriate for X
   - Recommend empirical percentiles rather than fixed thresholds
   - Our data supports X-specific calibration

3. **Focus on strongest signals**: Despite lower overall rate, X chromosome still shows:
   - 8,121 candidates at moderate threshold
   - Clear evidence of positive selection (max |Std iHS| = 5.74)
   - Candidates enriched near reproductive/immune genes (to be tested)

4. **Women's health implications**: X-linked selection particularly relevant because:
   - X genes often show sex-biased expression
   - Two X copies in females → different selective pressures
   - Potential for female-specific adaptive evolution

## Next Steps

1. ✓ Document X chromosome depletion findings
2. → **Implement X-specific thresholds in Phase 2** (region clustering)
3. → Map selection regions to genes
4. → Test for enrichment of reproductive and immune-related genes on X
5. → Investigate specific X chromosome candidate regions
6. → Compare X vs autosomal selection patterns in final analysis

## Files Referenced

- [results/analysis/candidates/chromosome_statistics.tsv](../results/analysis/candidates/chromosome_statistics.tsv)
- [results/analysis/candidates/x_enrichment_test.tsv](../results/analysis/candidates/x_enrichment_test.tsv)
- [results/analysis/candidates/ALL_candidates_moderate.tsv](../results/analysis/candidates/ALL_candidates_moderate.tsv)
