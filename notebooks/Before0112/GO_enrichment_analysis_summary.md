# GO Enrichment Analysis Results Summary 13112025 

**Script:** `scripts/GO_enrichment_with_terms.py`
**Date:** 13 November 2025
**Author:** MIT Van Bruggen

## Overview

This analysis tests whether a set of 13 high-confidence X chromosome genes show functional enrichment in specific Gene Ontology (GO) categories using Fisher's exact test.

## Study Genes (n=13)

High-confidence X chromosome genes identified from previous analysis:
- DACH2, DIAPH2, FAAH2, GPC3, HTR2C
- IL1RAPL2, KDM6A, NAP1L2, NCBP2L
- TRPC5, ZMAT1, KLHL13, PCDH11X

## Methodology

### 1. Functional Categories Defined

Four functional categories were defined using extensive keyword searches:
- **Reproduction** (108 X chr genes): fertility, sperm, oocyte, spermatogenesis, gonad, etc.
- **Immunity** (486 X chr genes): immune, inflammation, cytokine, T cell, antibody, etc.
- **Neurodevelopment** (419 X chr genes): neuron, synapse, brain, neural, learning, memory, etc.
- **Development** (664 X chr genes): embryo, differentiation, morphogenesis, growth factor, etc.

### 2. GO Term Search

- Used `geneinfo.ontology` to search GO database with keyword patterns
- Found thousands of GO terms matching each category:
  - Development: ~4,800 terms
  - Neurodevelopment: ~3,900 terms
  - Immunity: ~1,700 terms
  - Reproduction: ~425 terms

### 3. Background Gene Set

- **Background:** All protein-coding genes on X chromosome (n=852)
- This chromosome-specific background is appropriate according to https://munch-group.org/geneinfo/pages/gene_lists.html

### 4. Statistical Testing

- **Method:** Fisher's exact test using `GeneList.fisher()` method
- **Multiple testing correction:** Benjamini-Hochberg FDR correction (but is this needed??? does GeneList.fischer maybe already do this?)
- **Significance threshold:** FDR < 0.05

## Results

### A. Individual Functional Categories

| Category | Category Size | Overlap | Odds Ratio | P-value | FDR | Significant? |
|----------|--------------|---------|------------|---------|-----|--------------|
| **Immunity** | 486 | 10/13 | 2.54 | 0.118 | 0.472 | No |
| **Neurodevelopment** | 419 | 8/13 | 1.67 | 0.269 | 0.506 | No |
| **Development** | 664 | 11/13 | 1.57 | 0.426 | 0.506 | No |
| **Reproduction** | 108 | 2/13 | 1.26 | 0.506 | 0.506 | No |

**Key Findings:**
- **No significant enrichment** after FDR correction
- **Immunity** showed the strongest trend (p=0.118, OR=2.54) with 10/13 genes overlapping
- Most genes (11/13) are associated with general development processes
- Limited overlap with reproduction-specific functions (only 2/13 genes)

**Overlapping Genes by Category:**
- **Immunity (10):** DACH2, GPC3, HTR2C, IL1RAPL2, KDM6A, KLHL13, NAP1L2, PCDH11X, TRPC5, ZMAT1
- **Neurodevelopment (8):** DACH2, DIAPH2, GPC3, HTR2C, IL1RAPL2, KDM6A, NAP1L2, TRPC5
- **Development (11):** DACH2, DIAPH2, FAAH2, GPC3, HTR2C, IL1RAPL2, KDM6A, KLHL13, NAP1L2, TRPC5, ZMAT1
- **Reproduction (2):** DIAPH2, TRPC5

### B. Combined/Intersecting Categories

Testing genes that belong to multiple functional categories simultaneously:

| Category Combination | Size | Overlap | OR | P-value | FDR |
|---------------------|------|---------|----|---------|----|
| Neurodevelopment ∩ Immunity ∩ Development | 293 | 7/13 | 2.26 | 0.118 | 0.575 |
| Neurodevelopment ∩ Immunity | 294 | 7/13 | 2.24 | 0.120 | 0.575 |
| Neurodevelopment ∩ Development | 408 | 8/13 | 1.76 | 0.238 | 0.575 |
| Immunity ∩ Development | 483 | 9/13 | 1.73 | 0.265 | 0.575 |
| Reproduction ∩ Neurodevelopment ∩ Development | 73 | 2/13 | 1.97 | 0.308 | 0.575 |
| Reproduction ∩ Neurodevelopment | 74 | 2/13 | 1.94 | 0.314 | 0.575 |
| All other combinations | - | 1-2 | 0.85-1.97 | 0.48-0.71 | 0.70 |

**Key Findings:**
- No significant enrichment after FDR correction
- Strongest trend: genes involved in **all three** (Neurodevelopment ∩ Immunity ∩ Development) with 7/13 genes overlapping
- These 7 genes: DACH2, GPC3, HTR2C, IL1RAPL2, KDM6A, NAP1L2, TRPC5
- TRPC5 is the only gene appearing in all four categories

## Interpretation

### Main Conclusions

1. **No statistically significant enrichment** was detected after multiple testing correction
2. However, **biological trends** suggest these genes are involved in:
   - **Immunity** (strongest signal, 77% of genes)
   - **Neurodevelopment** (62% of genes)
   - **General development** (85% of genes)
   - Minimal representation in **reproduction** (15% of genes)

3. **Multi-functional genes:** Over half (7/13) of the study genes participate in neurodevelopment, immunity, AND development processes simultaneously

### Limitations --> maybe we can adjust

1. **Small sample size** (n=13 genes) limits statistical power (maybe test broader X gene set under selection)
2. **Multiple testing burden** (4 + 11 tests) reduces ability to detect modest effects (Can I make the background sets smaller and how)
3. **Broad categories** may dilute specific functional signals (how do I make these smaller, what key words should I use? )


## Output Files

All results saved to: `results/analysis/functional_enrichment/GO_with_terms/`

- `fisher_exact_test_results.tsv` - Individual category enrichment results
- `fisher_exact_test_combined_categories.tsv` - Combined category results
- `GO_terms_*.tsv` - GO term lists for each category (4 files)
- `GO_DAG_*.png` - GO term relationship visualizations (4 files)


