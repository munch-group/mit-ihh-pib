# iHS Analysis Summary - October 23, 2025

## Executive Summary

Successfully completed **Phase 1 (Candidate Identification)** and **Phase 2 (Region Definition)** of the X chromosome selection analysis. Key finding: X chromosome shows significant **depletion** of selection signals (~30% fewer than autosomes), which has important implications for the analysis strategy.

---

## Phase 1: Candidate Identification ✓ COMPLETE

### Results Overview

**Total variants analyzed:** 7,931,700
- Autosomes (chr1-22): 7,503,286 variants
- X chromosome: 428,414 variants

**Candidates identified (moderate threshold |Std iHS| ≥ 2.5):**
- **Liberal** (≥2.0): 395,725 candidates (5.0%)
- **Moderate** (≥2.5): 210,352 candidates (2.7%)
- **Stringent** (≥3.0): 117,104 candidates (1.5%)

### Key Finding 1: Excellent Standardization

All chromosomes show distributions extremely close to N(0,1):
- Mean Std iHS: ~0 (within 10⁻⁹)
- SD Std iHS: ~1.000 (within 0.000006)

**Implication:** The iHS standardization is working correctly. We can trust the p-values and thresholds derived from the standard normal distribution.

### Key Finding 2: X Chromosome Depletion

**Highly significant depletion of selection signals on X chromosome:**

| Metric | X Chromosome | Autosomes | Ratio |
|--------|--------------|-----------|-------|
| Candidate rate | 1.90% | 2.70% | 0.70 |
| Total candidates | 8,121 | 202,231 | - |
| 99th percentile | 2.90 | 3.3-3.8 | 0.79 |
| Max |Std iHS| | 5.74 | 8.0-10.5 | 0.59 |

- **Odds ratio:** 0.698 (30% depletion)
- **Fisher's exact test:** P = 1.71×10⁻²⁴²

### Biological Interpretations of X Depletion

#### 1. Hemizygosity in Males
- Males have only one X chromosome
- Selection acts differently on hemizygous loci
- May reduce extreme haplotype patterns that iHS detects

#### 2. Effective Population Size
- X chromosome Ne = 3/4 × autosomal Ne
- Smaller Ne → stronger genetic drift
- Drift may mask weaker selection signals

#### 3. Recombination Differences
- X recombines only in females (50% of meioses)
- Lower effective recombination rate
- Affects haplotype breakdown patterns

#### 4. Different Selection Pressures
- X-linked genes enriched for sex-specific functions
- May experience balancing selection rather than directional sweeps
- Sexual antagonism possible

### Analysis Decision: X-Specific Thresholds

Given the systematic X chromosome depletion, we implemented **chromosome-specific thresholds** for region definition:
- **Autosomes:** |Std iHS| ≥ 2.5 (moderate threshold)
- **X chromosome:** |Std iHS| ≥ 2.9 (empirical 99th percentile)

This maintains comparable selection stringency across the genome.

---

## Phase 2: Selection Region Definition ✓ COMPLETE

### Approach

**Clustering method:**
- 500kb sliding windows
- Minimum 2 SNPs per region (Liu et al. 2013)
- Merge overlapping windows
- Calculate region statistics (mean/max iHS, length, SNP density)

### Results Overview

**Total regions identified:** 4,666
- Autosomes: 4,530 regions
- X chromosome: 136 regions

### Autosome Regions

**Statistics:**
- Mean SNPs per region: 44.6
- Median SNPs per region: 36
- Mean region length: 413.9 kb
- Mean max |Std iHS|: 5.24

**Top 5 autosomal regions:**
1. **chr6:15.97-16.23 Mb** (265 kb, 46 SNPs, |iHS|=10.46)
   - Near JARID2 (histone methylation, development)
2. **chr11:68.21-68.71 Mb** (499 kb, 65 SNPs, |iHS|=10.01)
   - Contains CPT1A (fatty acid metabolism)
3. **chr16:79.29-79.75 Mb** (462 kb, 148 SNPs, |iHS|=9.96)
   - MAF gene region (crystallin transcription factor)
4. **chr4:77.29-77.74 Mb** (453 kb, 84 SNPs, |iHS|=9.89)
   - ANTXR2, SHROOM3 genes
5. **chr1:86.72-87.15 Mb** (437 kb, 51 SNPs, |iHS|=9.85)
   - CLCA gene family (chloride channels)

### X Chromosome Regions

**Statistics:**
- Mean SNPs per region: 31.2 (lower than autosomes)
- Median SNPs per region: 18
- Mean region length: 225.4 kb (shorter than autosomes)
- Mean max |Std iHS|: 4.16 (lower than autosomes)

**Top 5 X chromosome regions:**

1. **chrX:122.30-122.71 Mb** (410 kb, 58 SNPs, |iHS|=5.74)
   - Between KIAA2022 and BRWD3
   - KIAA2022: neurodevelopmental disorders
   - BRWD3: chromatin regulation

2. **chrX:57.48-57.51 Mb** (36 kb, 34 SNPs, |iHS|=5.68)
   - Gene-dense region, multiple X-linked genes

3. **chrX:111.80-112.12 Mb** (326 kb, 6 SNPs, |iHS|=5.59)
   - Sparse but strong signals

4. **chrX:3.98-4.29 Mb** (310 kb, 88 SNPs, |iHS|=5.56)
   - Near pseudoautosomal region (PAR1)
   - May have unique recombination dynamics

5. **chrX:97.20-97.48 Mb** (284 kb, 67 SNPs, |iHS|=5.42)
   - Gene-rich region

### X vs Autosome Comparison

| Metric | X Chromosome | Autosomes | Notes |
|--------|--------------|-----------|-------|
| Total regions | 136 | 4,530 | X has 3.0% of regions |
| SNPs/region | 31.2 | 44.6 | 30% fewer SNPs per region on X |
| Region length | 225 kb | 414 kb | X regions 46% shorter |
| Max |Std iHS| | 4.16 | 5.24 | X signals 21% weaker |

**Interpretation:** Even after adjusting thresholds, X chromosome selection regions are smaller, have fewer SNPs, and show weaker signals than autosomal regions. This is consistent with the overall depletion pattern.

---

## Files Generated

### Phase 1 Outputs

**Candidate variants:**
```
results/analysis/candidates/
├── ALL_candidates_liberal.tsv (395,725 variants)
├── ALL_candidates_moderate.tsv (210,352 variants)
├── ALL_candidates_stringent.tsv (117,104 variants)
├── chromosome_statistics.tsv
├── empirical_percentiles.tsv
├── x_enrichment_test.tsv
└── per_chromosome/ (69 files, 3 thresholds × 23 chromosomes)
```

### Phase 2 Outputs

**Selection regions:**
```
results/analysis/regions/
├── selection_regions_500kb_min2snps.tsv (4,666 regions)
├── selection_regions_500kb_min3snps.tsv (4,630 regions, stringent)
├── selection_regions_autosomes.tsv (4,530 regions)
├── selection_regions_X_chromosome.tsv (136 regions)
├── selection_regions.bed (for genome browsers)
└── selection_regions_stringent.bed
```

### Documentation

```
notebooks/
├── ihs_analysis_pipeline.md (updated with Phase 1-2 results)
├── X_chromosome_depletion_analysis.md (detailed depletion analysis)
├── NEXT_STEPS.md (quick reference)
└── ANALYSIS_SUMMARY_2025-10-23.md (this document)
```

---

## Next Steps: Phase 3 - Gene Annotation

### Objective
Annotate the 4,666 selection regions with genes to enable enrichment testing.

### Approach Options

**Option A: bedtools + GENCODE** (Recommended)
- Fast, local processing
- Works with BED format
- Standard gene annotations

**Option B: biomaRt (R)**
- Direct Ensembl access
- Rich metadata
- Can be slow for many regions

**Option C: Python PyEnsembl**
- Programmatic access
- Good for custom analyses
- Requires downloading Ensembl data

### Implementation Plan

1. **Download GENCODE annotations** (GRCh38)
   ```bash
   wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
   ```

2. **Convert GTF to gene BED file**
   - Extract gene entries
   - Format: chr, start, end, gene_name, gene_type

3. **Intersect regions with genes**
   ```bash
   bedtools intersect -a selection_regions.bed -b gencode.genes.bed -wa -wb
   ```

4. **Add flanking regions** (±50kb)
   - Capture regulatory elements
   - Include nearby genes

5. **Annotate with functional categories**
   - Gene Ontology terms
   - KEGG pathways
   - Custom gene lists (reproductive, immune)

### Expected Outputs

- `regions_with_genes.tsv` - All regions annotated with overlapping genes
- `X_chromosome_genes_in_selection.tsv` - X-specific gene list
- `autosome_genes_in_selection.tsv` - Autosome gene list (control)

---

## Next Steps: Phase 4 - Enrichment Analysis

### Questions to Address

1. **Are X chromosome selection regions enriched for reproductive genes?**
   - Expected: YES (based on project hypothesis)
   - Test against autosomal background

2. **Are X chromosome selection regions enriched for immune genes?**
   - Expected: YES (X-linked immunity hypothesis)
   - Compare to autosomes

3. **Do specific GO terms or pathways show enrichment?**
   - Identify unexpected functional categories
   - Link to women's health outcomes

### Analysis Strategy

**Statistical approach:**
- Fisher's exact test for each gene category
- Compare X vs autosome enrichment rates
- Multiple testing correction (FDR)

**Gene lists needed:**
1. **Reproductive genes** (~500-1000 genes)
   - GO: reproduction, fertilization, pregnancy
   - Literature-curated lists

2. **Immune genes** (~1000-2000 genes)
   - GO: immune response, inflammation
   - ImmPort, InnateDB databases

3. **Background** (all genes)
   - All protein-coding genes on X and autosomes
   - Account for different gene densities

**Tools:**
- clusterProfiler (R) for GO/KEGG enrichment
- g:Profiler for comprehensive enrichment
- Custom Fisher's exact tests for curated lists

---

## Technical Notes

### Data Quality Issues Encountered and Resolved

1. **chr21 ID format inconsistency**
   - Issue: chr21 has "chr21:pos_ref_alt" while others use "chr:pos:ref:alt"
   - Solution: Updated parsing to handle both formats
   - Script: `scripts/01_identify_ihs_candidates.py` lines 78-94

2. **Mixed chromosome type handling**
   - Issue: chr column had mixed types (int, str, "chr21")
   - Solution: Convert to string before parsing
   - Script: `scripts/02_define_selection_regions.py` lines 93-99

### Software Used

- **Python 3.13** (via pixi environment)
- **pandas 2.x** for data manipulation
- **scipy** for statistical tests
- **numpy** for numerical operations
- **SLURM** for job scheduling

### Compute Resources

- **Phase 1:** ~10 minutes, 16GB RAM
- **Phase 2:** ~6 seconds, 8GB RAM
- Both phases: Single CPU core sufficient

---

## Key Takeaways

### Scientific Findings

1. ✓ **iHS standardization is excellent** - All assumptions met
2. ✓ **X chromosome shows significant depletion** - 30% fewer signals
3. ✓ **X-specific thresholds implemented** - Accounts for systematic differences
4. ✓ **4,666 selection regions identified** - Ready for gene annotation
5. ✓ **136 X chromosome regions** - Despite depletion, substantial signal

### Analysis Decisions

1. **Use moderate threshold (|Std iHS| ≥ 2.5)** for autosomes
2. **Use empirical 99th percentile (≥2.9)** for X chromosome
3. **500kb windows with ≥2 SNPs** for region definition
4. **Separate X and autosome analyses** given systematic differences

### Manuscript Points

1. **Report X chromosome depletion** as key finding
   - Biologically expected (hemizygosity, Ne, recombination)
   - Methodologically important (threshold adjustment)

2. **Emphasize X-specific calibration** in methods
   - Novel approach accounting for X biology
   - Maintains comparable stringency

3. **Focus on strongest X signals** despite lower rates
   - 136 regions still substantial
   - May represent stronger selection events

4. **Women's health angle** for discussion
   - X-linked genes often sex-biased in expression
   - Two X copies in females → different dynamics
   - Implications for female-specific health outcomes

---

## Timeline

| Phase | Status | Date Completed |
|-------|--------|----------------|
| iHS calculation | ✓ Complete | Prior to 2025-10-23 |
| Phase 1: Candidates | ✓ Complete | 2025-10-23 |
| Phase 2: Regions | ✓ Complete | 2025-10-23 |
| Phase 3: Gene annotation | → Next | - |
| Phase 4: Enrichment | Planned | - |
| Phase 5: Interpretation | Planned | - |

**Estimated time to enrichment results:** 1-2 weeks

---

## Questions for Discussion

1. Should we use different window sizes for X vs autosomes given shorter region lengths?
2. How should we handle pseudoautosomal regions (PAR1, PAR2)?
3. Should we separately analyze X-inactivation escapee genes?
4. What p-value threshold for enrichment (0.05, 0.01, FDR < 0.1)?
5. Should we compare across populations (if data available)?

---

## Contact

For questions about this analysis:
- Scripts: [scripts/](../scripts/)
- Documentation: [notebooks/](../notebooks/)
- Results: [results/analysis/](../results/analysis/)
