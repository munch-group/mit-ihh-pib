# Phase 4: Functional Enrichment Analysis - Quick Start Guide

**Date:** October 24, 2025
**Status:** Ready to run
**Estimated Runtime:** 10-30 minutes

---

## Overview

> **Hypothesis:** X chromosome genes under selection are enriched for reproductive and immune functions

The analysis uses **g:Profiler** (via Python API) to test enrichment in:
- GO Biological Process
- GO Molecular Function
- GO Cellular Component
- KEGG Pathways
- Reactome Pathways

---

## Quick Start (3 Steps)

### Step 1: Installing Required Package

The enrichment analysis uses the `gprofiler-official` Python package. Install it:

```bash
pip install --user gprofiler-official
```

Or if using conda -- i'm using pixi:
```bash
conda install -c conda-forge gprofiler-official
```

**Verification:**
```bash
python3 -c "import gprofiler; print('✓ g:Profiler installed successfully')"
```

---

### Step 2: Run Enrichment Analysis

#### Option A: Submit SLURM Job (Recommended)

```bash
cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/scripts

sbatch run_enrichment_analysis.slurm
```

Check job status:
```bash
squeue -u $USER
```

View output:
```bash
tail -f /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/enrichment_analysis_*.out
```

#### Option B: Run Directly (Interactive)

```bash
cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/scripts

python3 04_functional_enrichment_analysis.py
```

---

### Step 3: Generate Visualizations

After the enrichment analysis completes, create publication-quality figures:

```bash
python3 05_visualize_enrichment_results.py
```

---

## What the Analysis Does

### Main Analyses

1. **X Chromosome Full Enrichment**
   - Tests all 241 X chromosome protein-coding genes
   - Identifies all enriched functional categories
   - Uses all protein-coding genes as background

2. **Autosome Full Enrichment**
   - Tests autosomal genes in selection regions
   - Provides comparison baseline
   - Same background as X analysis

3. **Reproductive Terms Filter**
   - Filters results for reproduction-related keywords:
     - reproduction, fertility, pregnancy, ovulation
     - spermatogenesis, oogenesis, meiosis, gamete
     - embryo, gonad, hormone, steroid
   - Compares X vs autosomes

4. **Immune Terms Filter**
   - Filters results for immune-related keywords:
     - immune, inflammation, cytokine, interferon
     - antibody, antigen, lymphocyte, T cell, B cell
     - complement, toll-like, defense response
   - Compares X vs autosomes

5. **X vs Autosome Comparison**
   - Identifies common enriched terms
   - Calculates enrichment ratios
   - Highlights X-specific and autosome-specific terms

6. **Statistical Testing**
   - Fisher's exact test for reproductive gene enrichment
   - Fisher's exact test for immune gene enrichment
   - Calculates odds ratios with significance

### Key Statistics Calculated

- **P-values:** FDR-corrected significance
- **Enrichment scores:** -log₁₀(p-value)
- **Odds ratios:** X vs autosome enrichment
- **Gene overlap:** Number of genes in each category

---

## Expected Outputs

### Results Files

All outputs saved to: `/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/analysis/enrichment/`

**Main result files:**
- `enrichment_X_chromosome_full.tsv` - All X enrichment results
- `enrichment_autosome_full.tsv` - All autosome enrichment results
- `enrichment_x_reproductive.tsv` - X reproductive terms only
- `enrichment_x_immune.tsv` - X immune terms only
- `enrichment_autosome_reproductive.tsv` - Autosome reproductive terms
- `enrichment_autosome_immune.tsv` - Autosome immune terms
- `x_vs_autosome_comparison.tsv` - Direct comparison
- `enrichment_summary.tsv` - Summary statistics

### Visualization Files

All figures saved to: `results/analysis/enrichment/figures/`

**Figures generated:**
- `x_chromosome_top_terms.png` - Top 20 X chromosome terms (barplot)
- `autosome_top_terms.png` - Top 20 autosome terms (barplot)
- `x_reproductive_terms.png` - Reproductive enrichment (barplot)
- `x_immune_terms.png` - Immune enrichment (barplot)
- `x_vs_autosome_scatter.png` - Comparison scatter plot
- `x_source_distribution.png` - Enrichment by database (pie chart)
- `enrichment_heatmap.png` - Cross-analysis heatmap
- `enrichment_summary.png` - Multi-panel summary figure

---

## Interpreting Results

### Positive Results (Supporting Hypothesis)

If we see:
- ✓ **Reproductive terms significantly enriched in X genes**
  - Odds ratio > 1.0, p-value < 0.05
  - Examples: oogenesis, meiosis, fertility, hormone signaling

- ✓ **Immune terms significantly enriched in X genes**
  - Odds ratio > 1.0, p-value < 0.05
  - Examples: immune response, cytokine signaling, interferon

**Interpretation:** Confirms hypothesis that X chromosome selection targets reproductive and immune functions

### Alternative Results

If reproductive/immune enrichment is weak or absent:
- Check what IS enriched (may find other interesting patterns)
- X chromosome may be under selection for other functions:
  - Neurodevelopment (known X enrichment)
  - Chromatin regulation
  - Metabolism
  - Structural/housekeeping genes


### Statistical Thresholds

- **FDR < 0.05** - Significant enrichment (corrected for multiple testing)
- **FDR < 0.01** - Highly significant
- **FDR < 0.001** - Very strong signal
- **Odds ratio > 1.5** - Meaningful biological enrichment
- **Odds ratio > 2.0** - Strong enrichment

---

## Next Steps After Analysis

### 1. Review Results 

```bash
# Check summary
cat results/analysis/enrichment/enrichment_summary.tsv

# View top X chromosome terms
head -20 results/analysis/enrichment/enrichment_X_chromosome_full.tsv

# Check reproductive enrichment
cat results/analysis/enrichment/enrichment_x_reproductive.tsv

# Check immune enrichment
cat results/analysis/enrichment/enrichment_x_immune.tsv
```

### 2. Literature Review (1-2 Days)

For top enriched genes:
- Search PubMed for known functions
- Check OMIM for disease associations
- Look up GWAS hits (GWAS Catalog)
- Review GTEx for sex-biased expression

### 3. Create Manuscript Figures 

Priority figures:
1. **Figure 1:** Top X chromosome enrichment (barplot) - use `x_chromosome_top_terms.png`
2. **Figure 2:** Reproductive vs Immune enrichment (combined barplot)
3. **Figure 3:** X vs Autosome comparison (scatter plot)
4. **Table 1:** Summary statistics and Fisher's test results

### 4. Draft Results Section 

Key points to address:
- Overall enrichment patterns (what's enriched?)
- Reproductive gene enrichment (hypothesis test #1)
- Immune gene enrichment (hypothesis test #2)
- X vs autosome comparison (differential enrichment)
- Notable genes and pathways
- Statistical evidence (p-values, odds ratios)


## Alternative Enrichment Tools

If g:Profiler doesn't work :

### 1. Web-Based Tools (No Installation)

**g:Profiler Web** (https://biit.cs.ut.ee/gprofiler/gost)
- Upload: `X_protein_coding_genes.txt`
- Background: `protein_coding_genes.txt`
- Organism: Homo sapiens
- Click "Run query"

**DAVID** (https://david.ncifcrf.gov/)
- Upload gene list
- Select "OFFICIAL_GENE_SYMBOL"
- Run functional annotation

**Enrichr** (https://maayanlab.cloud/Enrichr/)
- Paste gene list
- Browse enrichment results
- Download figures

### 2. R-Based Tools

**clusterProfiler:**
```R
library(clusterProfiler)
library(org.Hs.eg.db)

# Read genes
x_genes <- read.table("X_protein_coding_genes.txt")[,1]

# Run enrichment
ego <- enrichGO(
  gene = x_genes,
  OrgDb = org.Hs.eg.db,
  keyType = 'SYMBOL',
  ont = 'BP',
  pAdjustMethod = 'fdr',
  pvalueCutoff = 0.05
)

# Plot
dotplot(ego, showCategory=20)
```



**Next command:**

```bash
# Install package
pip install --user gprofiler-official

# Run analysis
cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/scripts
sbatch run_enrichment_analysis.slurm

# Or run directly
python3 04_functional_enrichment_analysis.py
```

---

*Last updated: October 24, 2025*
*Author: MIT van Bruggen*
*Project: X Chromosome Selection Analysis*
