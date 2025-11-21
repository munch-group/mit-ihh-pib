# Progress Report: Functional Enrichment Analysis of X Chromosome Selection Candidates
Mit Van Bruggen 21/11/25
## Overview
Three complementary scripts to perform functional enrichment analysis on 13 high-confidence X chromosome genes showing signatures of positive selection (DACH2, DIAPH2, FAAH2, GPC3, HTR2C, IL1RAPL2, KDM6A, NAP1L2, NCBP2L, TRPC5, ZMAT1, KLHL13, PCDH11X).

---

## 1. explore_GO_hierarchies.py

**Script Location:** [scripts/explore_GO_hierarchies.py](../scripts/explore_GO_hierarchies.py)

### Purpose
Creates focused gene sets by extracting hierarchical GO term subtrees from biologically meaningful seed terms.

### Method
- Defined 8 functional categories with 3-6 GO seed terms each based on known gene functions
- Extracted complete GO subtrees (including all descendant terms) for each category using recursive traversal
- Filtered to X chromosome protein-coding genes only
- Generated DAG visualizations showing hierarchical term relationships

### Results

| Category | GO Terms | X Genes | Genes/Term | Description |
|----------|----------|---------|------------|-------------|
| Female fertility | 5 | 4 | 0.80 | Highly specific reproductive module |
| Ovarian function | 19 | 4 | 0.21 | Ovarian development |
| Neuroendocrine reproduction | 157 | 24 | 0.15 | Appetite/feeding control |
| Synapse assembly | 52 | 61 | 1.17 | Well-defined synaptic module |
| Axon/dendrite morphogenesis | 110 | 24 | 0.22 | Neuron structure development |
| Calcium homeostasis | 58 | 8 | 0.14 | Focused signaling module |
| Wnt signaling | 29 | 12 | 0.41 | Developmental signaling |
| Neuroimmune crosstalk | 585 | 160 | 0.27 | Broad immune signaling |

### Outputs
- `subtree_summary.tsv` - Summary statistics for each category
- `genes_unique_*.txt` - Gene lists for each category
- `subtree_terms_*.json` - Complete GO term hierarchies
- `dag_plots/dag_*.png` - Visual representations of GO hierarchies

---

## 2. Fischer_Overlap_Analysis.py

**Script Location:** [scripts/Fischer_Overlap_Analysis.py](../scripts/Fischer_Overlap_Analysis.py)

### Purpose
Tests whether study genes significantly overlap with each functional category using Fisher's exact test.

### Method
- Used focused gene sets from script #1 as functional categories
- Performed Fisher's exact test for each category
- Background: X chromosome protein-coding genes (n=852)
- Applied Benjamini-Hochberg FDR correction for multiple testing

### Key Results

#### Significant (FDR < 0.05)

1. **Female fertility**
   - Odds Ratio: 76.1
   - p-value: 0.0013
   - FDR: 0.010 
   - Overlapping genes: DACH2, DIAPH2
   - Overlap: 2/13 study genes

2. **Calcium homeostasis**
   - Odds Ratio: 25.2
   - p-value: 0.0057
   - FDR: 0.023 
   - Overlapping genes: HTR2C, TRPC5
   - Overlap: 2/13 study genes

#### Borderline Significant (FDR < 0.12)

3. **Axon/dendrite morphogenesis**
   - Odds Ratio: 6.8
   - p-value: 0.049
   - FDR: 0.119
   - Overlapping genes: NAP1L2, TRPC5
   - Overlap: 2/13 study genes

4. **Synapse assembly**
   - Odds Ratio: 4.0
   - p-value: 0.060
   - FDR: 0.119
   - Overlapping genes: HTR2C, IL1RAPL2, TRPC5
   - Overlap: 3/13 study genes

#### Not Significant - Maybe join them? 

- Wnt signaling: FDR = 0.271 (GPC3)
- Neuroendocrine reproduction: FDR = 0.416 (HTR2C)
- Neuroimmune crosstalk: FDR = 0.518 (GPC3, HTR2C, IL1RAPL2)
- Ovarian function: No overlap

###  Output
- `focused_enrichment_results.tsv` - Complete Fisher's test results with FDR correction

---

## 3. GO_enrichment_focused_cats.py

**Script Location:** [scripts/GO_enrichment_focused_cats.py](../scripts/GO_enrichment_focused_cats.py)

### Purpose
Performs proper GO enrichment analysis using the geneinfo package for categories showing Fisher's test enrichment.

### Method
- Tested 7 categories prioritized by Fisher's test results
- Used GO regex patterns to identify relevant terms for each category
- Applied standard GO enrichment with geneinfo.go_enrichment()
- Background: X chromosome protein-coding genes (n=852)
- Applied FDR correction
- Expanded gene list to 21 genes including additional candidates (OPN1LW, OPN1MW family, etc.)

### Results

#### Categories Tested

| Category | Priority | GO Terms | Enriched Terms (FDR<0.05) |
|----------|----------|----------|---------------------------|
| Female fertility | 1 | 14 | 0 |
| Calcium homeostasis | 1 | 48 | 0 |
| Neuroimmune crosstalk | 2 | 373 | 0 |
| Synapse assembly | 3 | 469 | 0 |
| Axon/dendrite morphogenesis | 3 | 352 | 0 |
| Serotonin signaling | 4 | 65 | 0 |
| Chromatin regulation | 4 | 353 | 0 |


### Key Outputs
- `go_enrichment_all_results.tsv` - All GO enrichment results
- `go_enrichment_significant_FDR05.tsv` - Significant terms only
- `enrichment_summary.tsv` - Summary by category
- `go_enrichment_*.tsv` - Results per category
- `top_enriched_go_terms.png` - Visualization of top terms

---

## Summary & Interpretation

### Strong Evidence For:

1. **Female fertility enrichment** (DACH2, DIAPH2)
   - Highly significant in Fisher's overlap analysis (FDR = 0.010)
   - Small, focused gene set (4 X genes total) - maybe too small 
   - High odds ratio (76.1) indicates strong association

2. **Calcium signaling enrichment** (HTR2C, TRPC5)
   - Significant in Fisher's overlap analysis (FDR = 0.023)
   - Moderate-sized gene set (8 X genes total)
   - Both genes have known calcium channel/signaling roles

## To do 
- Refine go.go_enrichment analysis 
- Refine study gene set 
   - p-value ≤ 10^-6 AND within +/- 100kb from detected SNP 
- Make ‘workflow’: all scripts that need to be run for my analysis
   - GWF
- Doing now: Got started on a script to see overlap/difference between relate signals and iHS signals 


