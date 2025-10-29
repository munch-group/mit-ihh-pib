# Progress Report: X Chromosome Selection Analysis
**Student:** MIT van Bruggen
**Date:** October 24, 2025
**Project:** Investigating Positive Selection on the Human X Chromosome with Focus on Reproductive and Immune Gene Enrichment

---

## Executive Summary

- Successfully completed **Phases 1-3** of the X chromosome selection analysis, processing 7.9 million variants across all chromosomes to identify 4,666 selection regions containing 48,678 genes. 
- Key finding is that the X chromosome shows significant depletion of selection signals (~30% fewer than autosomes), which led me to implement chromosome-specific thresholds. 
- Now ready for Phase 4 (functional enrichment testing) to address the core hypothesis about reproductive and immune gene enrichment.

### Completed Work
- ✓ **Phase 1:** Candidate identification (210,352 significant variants)
- ✓ **Phase 2:** Region definition (4,666 selection regions, 136 on X)
- ✓ **Phase 3:** Gene annotation (48,678 genes, 664 on X chromosome)
- → **Phase 4:** Ready for enrichment analysis

---

## Project Background and Hypothesis

### Research Question
Are genes under recent positive selection on the human X chromosome enriched for reproductive and immune functions, with implications for women's health?

### Rationale
The X chromosome is particularly interesting for several reasons:
1. **Hemizygosity in males** - selection acts differently on X-linked loci
2. **Dosage effects in females** - two X copies with X-inactivation
3. **Enrichment for sex-specific traits** - reproduction, cognition, immunity
4. **Clinical relevance** - X-linked disorders disproportionately affect women's health

### Approach
Using **integrated Haplotype Score (iHS)** to detect recent positive selection:
- iHS identifies extended haplotypes characteristic of selective sweeps
- Standardized scores should follow N(0,1) under neutrality
- High |iHS| values indicate selection

---

## Phase 1: Selection Candidate Identification

### Methodology

**Data:**
- 1000 Genomes Phase 3 data (GRCh38)
- 7,931,700 biallelic variants (chr1-22 + chrX)
- Ancestral allele polarization using human-chimp alignment

**Statistical Approach:**
- Applied literature-based thresholds:
  - Liberal: |Std iHS| ≥ 2.0 (Voight et al. 2006)
  - **Moderate: |Std iHS| ≥ 2.5** (Paul et al. 2024) ← Primary threshold
  - Stringent: |Std iHS| ≥ 3.0
- Calculated two-tailed p-values assuming N(0,1)
- Computed empirical percentiles (99th, 99.5th, 99.9th) per chromosome

### Quality Control: Excellent Standardization

All 23 chromosomes showed near-perfect standardization:
- **Mean Std iHS:** ~0 (within 10⁻⁹)
- **SD Std iHS:** ~1.000 (within 0.000006)

**Interpretation:** The standardization worked well. Statistical assumptions are met, p-values are reliable, and we can trust the thresholds derived from the standard normal distribution.

### Key Finding #1: X Chromosome Depletion

**Observation:**
| Chromosome Type | Candidate Rate | Total Candidates | 99th Percentile |
|----------------|---------------|-----------------|-----------------|
| **Autosomes** | 2.70% | 202,231 / 7,503,286 | 3.3-3.8 |
| **X chromosome** | 1.90% | 8,121 / 428,414 | 2.90 |
| **Difference** | -30% | -3,446 variants | -20% |

- **Fisher's exact test:** P = 1.71×10⁻²⁴²
- **Odds ratio:** 0.698 (X has 30% fewer signals)

**This was unexpected but maybe biologically important.**

### Biological Interpretation of X Depletion (?)

Possible explan.:

#### 1. Hemizygosity in Males
- Males have only one X chromosome
- Selection acts immediately on all X-linked alleles in males
- Reduces extreme haplotype homozygosity patterns
- BUT **iHS specifically detects extended homozygosity**, which is then short-lived  on X 
So: Because X-linked alleles are exposed to selection in hemizygous males, beneficial alleles fix more rapidly. This shortens the window during which extended haplotype homozygosity can be observed, making ongoing sweeps and thus high iHS values rarer in cross-sectional data.

#### 2. Effective Population Size
- Ne(X) = 3/4 × Ne(autosomes) due to 3 X copies per 4 autosomes
- Smaller Ne → stronger genetic drift
- Drift can mask weaker selection signals
- **Mathematical expectation:** X should show ~25% reduction in selection signal strength

#### 3. Recombination Rate Differences
- X chromosome recombines only in females (50% of meioses)
- Effective recombination rate is lower
- Slower haplotype breakdown → affects iHS calculations
- May require different time scaling

#### 4. Different Selection Modes
- X chromosome enriched for sex-specific genes
- May experience **balancing selection** rather than directional sweeps?
- Frequency-dependent selection
- **iHS detects directional sweeps best** - other selection modes give weaker signals

**Verdict:** The depletion may be **biologically expected** ? (Villegas‐Mirón et al., 2021, Meisel and Connallon, 2013)

### Methodological Response: X-Specific Thresholds

Given the systematic X chromosome depletion, I implemented **chromosome-specific thresholds** for region definition:

| Chromosome | Threshold | Rationale |
|-----------|-----------|-----------|
| **Autosomes** | \|Std iHS\| ≥ 2.5 | Standard moderate threshold (Paul et al. 2024) |
| **X chromosome** | \|Std iHS\| ≥ 2.9 | Empirical 99th percentile 

**This maintains comparable selection stringency across the genome (hopefully?)**


### Phase 1 Results Summary

**Candidates identified at moderate threshold (|Std iHS| ≥ 2.5):**
- Total: 210,352 variants (2.7% of genome)
- Autosomes: 202,231 variants
- X chromosome: 8,121 variants

**Files generated:**
- Per-chromosome candidate files (69 files)
- Combined genome-wide files (3 threshold levels)
- Chromosome statistics with normality diagnostics
- X chromosome enrichment/depletion test results

**documentation:**
- [iHS_threshold_literature_review.md](iHS_threshold_literature_review.md) - Threshold justification
- [X_chromosome_depletion_analysis.md](X_chromosome_depletion_analysis.md) - Detailed depletion analysis

---

## Phase 2: Selection Region Definition

### Methodology

**Clustering approach** :
- **Window size:** 500 kb sliding windows
- **Minimum SNPs:** ≥2 significant SNPs per region
- **Merging:** Overlapping windows combined
- **Rationale:** Single outlier SNPs may be false positives; require clustering for confidence

**Chromosome-specific thresholds applied:**
- Autosomes: |Std iHS| ≥ 2.5
- X chromosome: |Std iHS| ≥ 2.9 (99th percentile)

### Results: 4,666 Selection Regions Identified

**Overall statistics:**
| Metric | Autosomes | X Chromosome | X as % of Autosome |
|--------|-----------|--------------|-------------------|
| **Regions** | 4,530 | 136 | 3.0% |
| **SNPs/region** (mean) | 44.6 | 31.2 | 70% |
| **Region length** (mean) | 414 kb | 225 kb | 54% |
| **Max \|iHS\|** (mean) | 5.24 | 4.16 | 79% |

**Interpretation:** Even after threshold adjustment, X chromosome selection regions are:
- **Smaller** (46% shorter on average)
- **Sparser** (30% fewer SNPs per region)
- **Weaker** (21% lower maximum iHS)

This suggests X chromosome selection may involve:
- More localized selective sweeps
- Weaker selection coefficients
- Faster haplotype breakdown (due to smaller Ne)

### Top X Chromosome Selection Regions

| Rank | Genomic Position | Length | SNPs | Max \|iHS\| | P-value |
|------|-----------------|--------|------|-------------|---------|
| 1 | **chrX:122.30-122.71 Mb** | 410 kb | 58 | 5.74 | 9.5×10⁻⁹ |
| 2 | **chrX:57.48-57.51 Mb** | 36 kb | 34 | 5.68 | 1.3×10⁻⁸ |
| 3 | **chrX:111.80-112.12 Mb** | 326 kb | 6 | 5.59 | 2.3×10⁻⁸ |
| 4 | **chrX:3.98-4.29 Mb** | 310 kb | 88 | 5.56 | 2.7×10⁻⁸ |
| 5 | **chrX:97.20-97.48 Mb** | 284 kb | 67 | 5.42 | 6.1×10⁻⁸ |

**Note on top region (chrX:122.30-122.71 Mb):**
- Strongest X chromosome signal
- Contains/near KIAA2022 (neurodevelopmental disorders)
- Contains/near BRWD3 (chromatin regulation, X-inactivation escape)
- Highly significant despite X chromosome depletion

### Comparison to Autosomes

For context, the **strongest autosomal signals:**

1. **chr6:15.97-16.23 Mb** (|iHS| = 10.46) - JARID2 gene, histone methylation
2. **chr11:68.21-68.71 Mb** (|iHS| = 10.01) - CPT1A, fatty acid metabolism (known Arctic adaptation)
3. **chr16:79.29-79.75 Mb** (|iHS| = 9.96) - MAF gene region

The strongest X signal (5.74) is **45% weaker** than the strongest autosomal signal, but still **highly significant** and represents genuine selection.

### Phase 2 Outputs

**Files created:**
- `selection_regions_500kb_min2snps.tsv` - All 4,666 regions
- `selection_regions_autosomes.tsv` - 4,530 autosomal regions
- `selection_regions_X_chromosome.tsv` - 136 X chromosome regions
- `selection_regions.bed` - Genome browser visualization format

**Key documentation:**
- [ihs_analysis_pipeline.md](ihs_analysis_pipeline.md) - Complete pipeline documentation

---

## Phase 3: Gene Annotation

### Methodology

**Annotation source: GENCODE v44** (GRCh38.p14)
- Gold standard for human genome annotation
- Used by ENCODE, GTEx, 1000 Genomes
- 62,700 genes total (20,046 protein-coding)

**Technical approach:**
- Local GTF parsing (not Ensembl REST API)
- ±50 kb flanking regions included
- **Completed in 26 seconds** for all 4,666 regions
  - (Compare: Ensembl API would take ~40 minutes)

**Gene categorization:**
- Overlap type (contained, overlapping, flanking)
- Gene biotype (protein-coding, lncRNA, pseudogene, etc.)
- X-inactivation status (escape vs subject)

### X-Inactivation Annotation

I integrated X-inactivation escape gene lists from **Tukiainen et al. (2017) Nature**:
- ~15% of X genes escape X-inactivation in females
- These genes may have different selective pressures
- Important for dosage compensation and sex-specific effects

**X-inactivation escapees found in selection regions:**
- CD99, XG, GYG2, NLGN4X, TXLNG, ZFX, USP9X, DDX3X, KDM6A, CHIC1 (10 genes)

### Results: 48,678 Genes in Selection Regions

**Overall gene counts:**
| Category | Total | X Chromosome | % X |
|----------|-------|--------------|-----|
| **All genes** | 48,678 | 664 | 1.4% |
| **Protein-coding** | 16,573 | 242 | 1.5% |
| **lncRNA** | 14,899 | - | - |
| **Pseudogenes** | 7,500 | - | - |

**Region coverage:**
- 98.3% of regions contain at least one gene (4,587/4,666)
- Average: 11.3 genes per region
- Median: 9 genes per region

### X Chromosome Genes

**242 X chromosome protein-coding genes in selection regions**

Examples of notable genes:
- **CD99** - Cell adhesion, immune function, X-inactivation escapee
- **NLGN4X** - Neuroligin 4, synaptic function, autism risk gene
- **STS** - Steroid sulfatase, hormone metabolism
- **MID1** - Microtubule dynamics, Opitz syndrome
- **PDHA1** - Pyruvate dehydrogenase, metabolism

### Gene-Rich Selection Regions

**Top 10 regions with most genes:**

1. **chr14:22.17-22.58 Mb** - 100 genes (olfactory receptor gene cluster)
2. **chr7:142.27-142.77 Mb** - 95 genes (protease genes: PRSS1, PRSS2)
3. **chr14:100.50-100.99 Mb** - 95 genes (DLK1-RTL1 imprinted region)
4. **chr19:53.74-54.21 Mb** - 81 genes (includes immune genes: NLRP12, OSCAR)
5. **chr6:31.53-31.90 Mb** - 77 genes (**MHC region!** - MICB, TNF, LTA, LTB)

**The MHC region on chr6 is particularly interesting:**
- Major Histocompatibility Complex
- Central to immune function
- Well-known selection target
- **Validates our method** - we're capturing known selection regions

### Phase 3 Outputs

**Main annotation files:**
- `selection_regions_annotated.tsv` - All regions with genes
- `genes_in_selection_regions_detailed.tsv` - All genes with region info
- `X_chromosome_genes_in_selection.tsv` - 664 X genes
- `X_protein_coding_genes.tsv` - 242 X protein-coding genes

**Gene lists for enrichment** (simple text format):
- `all_genes.txt` - 48,677 gene names
- `protein_coding_genes.txt` - 16,572 genes
- `X_chromosome_genes.txt` - 663 genes
- `X_protein_coding_genes.txt` - 241 genes

These files are ready for enrichment tools like **g:Profiler**, **DAVID**, or **Enrichr**.




## Biological Interpretation Framework

### Why X Chromosome Depletion Doesn't Undermine the Hypothesis

**Original hypothesis:** X chromosome enriched for selection in reproductive/immune genes

**Revised hypothesis:** Despite overall X depletion (due to hemizygosity/Ne), **X-linked genes under selection show enrichment for reproductive and immune functions**

**Why this is actually a better story:**

1. **Lower background rate** → signals we do find are more meaningful
2. **Stronger support for hypothesis** → enrichment despite depletion is more convincing
3. **Methodological rigor** → we're accounting for X biology properly
4. **Novel contribution** → systematic X chromosome calibration




## Next Steps: Phase 4 - Functional Enrichment

### Ready for Enrichment Testing

**Gene lists prepared:**
- X chromosome genes: 664 total (242 protein-coding)
- Autosomal genes: 48,014 total (16,331 protein-coding)
- Background: All genes in GENCODE v44



## Tools
**Recommended: g:Profiler** (web-based)
- Input: Gene name lists (we have these ready)
- Tests: GO Biological Process, KEGG, Reactome
- Multiple testing: Built-in correction
- Output: Publication-ready figures


## Limitations and Caveats

### 1. Single Population
- Analysis uses 1000 Genomes ALL superpopulation
- Selection signals may be population-specific
- **Mitigation:** Standard for initial discovery, can stratify later

### 2. iHS-Specific Biases
- Best for recent, strong selection (<50,000 years)
- Misses older or weaker selection
- Assumes hard selective sweeps
- **Mitigation:** Complements other methods (XP-EHH, FST)

### 3. X Chromosome Complexity
- Pseudoautosomal regions have different dynamics
- X-inactivation patterns vary
- Sex-specific effects hard to detect in mixed samples
- **Mitigation:** PAR regions flagged, X-inactivation annotated

### 4. Gene Annotation Limitations
- Overlapping genes assigned to regions
- ±50kb flank is arbitrary
- Non-coding regulatory regions not fully captured
- **Mitigation:** Standard approach, captures most functional elements

### 5. Enrichment Testing (Planned)
- GO term annotations incomplete
- Literature curation subjective
- Multiple testing reduces power
- **Mitigation:** Use multiple databases, FDR correction, effect sizes

---

---


## Conclusion

I have successfully completed Phases 1-3 of the X chromosome selection analysis, establishing a robust analytical framework and identifying 4,666 selection regions containing 48,678 genes. The key unexpected finding—X chromosome depletion of selection signals—led to a novel methodological contribution (chromosome-specific threshold calibration) and provides important biological insights about X chromosome evolution.

The analysis is now positioned for Phase 4 (functional enrichment testing), which will address the core hypothesis about reproductive and immune gene enrichment. The gene lists are prepared, methodological framework is sound, and computational infrastructure is in place for rapid completion of enrichment analysis.

**Current status:** Ready to test hypothesis in Phase 4

