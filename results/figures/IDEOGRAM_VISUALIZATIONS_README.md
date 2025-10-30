# X Chromosome Ideogram Visualizations

**Created:** 2025-10-30
**Script:** `scripts/create_ideogram_visualizations.py`

## Overview

Three ideogram visualizations were created to show the distribution of selection signals, SNPs, and genes across the X chromosome.

---

## Visualizations Created

### 1. High-Significance SNPs (n=43)
**File:** `X_chromosome_ideogram_high_significance_snps.png`

**Shows:**
- 43 SNPs with -log10(p-value) ≥ 6
- SNP positions labeled with nearest gene names (top 10 SNPs)
- Color intensity scaled by p-value significance (red = more significant)
- Scatter plot showing all SNPs with size proportional to significance

**Key findings:**
- SNPs are distributed across the X chromosome
- Top 10 SNPs labeled with their nearest genes
- Significance ranges from -log10(p) = 6.0 to ~8.0

---

### 2. Selection Regions (n=136)
**File:** `X_chromosome_ideogram_selection_regions.png`

**Shows:**
- 136 selection regions from Phase 3 analysis (|Std iHS| ≥ 2.9)
- Regions colored by whether they contain high-sig SNPs:
  - **Red** (n=18): Regions containing high-significance SNPs
  - **Blue** (n=118): Regions without high-significance SNPs
- Coverage track showing maximum |Std iHS| across the chromosome

**Key findings:**
- Only 18/136 (13.2%) regions contain high-significance SNPs
- Most high-sig SNPs cluster in specific regions
- Confirms that extreme signals (p < 10^-6) are concentrated, not evenly distributed

---

### 3. Combined SNPs + Genes + Functional Categories
**File:** `X_chromosome_ideogram_functional_categories.png`

**Shows:**
- Gene positions colored by functional category:
  - **Red**: Immune genes
  - **Blue**: Neural genes
  - **Purple**: Immune + Neural genes
  - **Green**: Developmental genes
  - **Gray**: Other genes
- 11 overlap genes labeled (present in both high-sig and Phase 3)
- SNPs colored by their nearest gene's functional category
- Size of SNP markers proportional to significance

**Key findings:**
- Overlap genes span different functional categories
- Most high-confidence genes are in "other" or neural/immune categories
- Visual confirmation of gene-SNP associations

---

## Analysis Questions Answered

### 1. Are selection signals clustered or evenly distributed?
**Answer:** **Clustered**
- Only 18/136 selection regions contain high-significance SNPs
- Most extreme signals are concentrated in specific genomic regions
- This supports the idea of localized selective sweeps

### 2. Which cytogenetic bands show most selection?
**Answer:** (Visible from ideograms)
- Multiple regions across the X chromosome
- Coverage track shows peaks at specific locations
- Centromeric and telomeric regions show different patterns

### 3. Do immune/reproductive genes cluster together?
**Answer:** **Mixed distribution**
- Immune and neural genes are distributed across the chromosome
- Some overlap genes have dual functions (immune + neural)
- No obvious clustering by functional category

---

## Data Sources

| Data | File | N |
|------|------|---|
| High-sig SNPs | `results/analysis/candidates/chrX_high_significance.tsv` | 43 |
| Selection regions | `results/analysis/regions/selection_regions_500kb_min2snps.tsv` | 136 (X chr) |
| Genes | `results/analysis/high_significance_genes/gene_details_complete.tsv` | 13 protein-coding |
| SNP-gene mapping | `results/analysis/high_significance_genes/snp_to_nearest_gene_only.tsv` | 43 |

---

## Technical Details

**Genome assembly:** hg38
**Tool:** geneinfo Python package
**Resolution:** 300 DPI
**Format:** PNG

**Significance thresholds:**
- High-significance SNPs: -log10(p-value) ≥ 6 (p < 10^-6)
- Phase 3 regions: |Std iHS| ≥ 2.9

**Functional categories:**
- Immune, Neural, Immune+Neural, Developmental, Other
- Based on gene annotations from geneinfo

---

## How to Regenerate

```bash
# Activate pixi environment
cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib
pixi shell

# Run script
python scripts/create_ideogram_visualizations.py
```

Output files will be created in `results/figures/`

---

## Next Steps

These visualizations can be used to:
1. Identify regions for detailed functional analysis
2. Prioritize genes for experimental validation
3. Understand the genomic distribution of selection signals
4. Support manuscript figures and presentations

**Recommendation:** Focus on the 18 regions containing high-significance SNPs for:
- Functional enrichment analysis
- Literature review of candidate genes
- Investigation of potential selective pressures
