# Phase 4: Next Steps Using geneinfo Module
**Date:** October 29, 2025
**Student:** MIT van Bruggen
**Status:** Ready to implement functional enrichment analysis

---

## Overview

After completing Phase 3 (gene annotation) and extracting 43 high-significance SNPs (-log10 p ≥ 6), I'm now ready to use the **geneinfo module** for comprehensive functional analysis. The geneinfo module provides integrated tools for gene annotation, visualization, and functional enrichment - exactly what I need for Phase 4.

---

## What is geneinfo?

Based on the guide, geneinfo is a comprehensive Python package that provides:

1. **ID conversion** - Convert between gene IDs, symbols, etc.
2. **Getting gene information** - Retrieve gene details, coordinates, functions
3. **Gene context plotting** - Plot data over gene annotations
4. **Ideograms** - Chromosome visualization showing gene/SNP positions
5. **STRING networks** - Protein-protein interaction networks
6. **GO graphs** - Gene Ontology visualization
7. **GO information** - Retrieve GO terms for genes
8. **GO enrichment** - Statistical testing for functional enrichment
9. **Gene lists** - Manage and manipulate gene sets
10. **Caching** - Efficient data storage/retrieval

This replaces the need for multiple separate tools (g:Profiler, DAVID, manual STRING queries, etc.)!

---

## Current Status: What I Have

### Input Data Ready:

1. **High-significance SNPs** (43 SNPs, -log10 p ≥ 6)
   - File: `results/analysis/candidates/chrX_high_significance.tsv`
   - Columns: Location, iHH_0, iHH_1, iHS, Std iHS, chr, pos, ref, alt, p_value, neg_log10_p

2. **Phase 3 Gene Lists** (|Std iHS| ≥ 2.9, broader threshold)
   - X chromosome: 664 genes (242 protein-coding)
   - File: `results/analysis/gene_annotation/X_chromosome_genes_in_selection.tsv`
   - File: `results/analysis/gene_annotation/unique_Xgenes2410.txt`

3. **Autosomal genes under selection** (for background)
   - 48,014 autosomal genes
   - File: `results/analysis/gene_annotation/autosome_genes_in_selection.tsv`

---

##  Analysis Workflow Using geneinfo

### **Step 1: Map High-Significance SNPs to Nearest Genes**
*Priority: HIGH | Estimated time: 1 hour*

**Goal:** Find which genes are closest to the 43 strongest selection signals

**Input:**
- 43 SNP positions from `chrX_high_significance.tsv`
- Columns needed: chr (X), pos (genomic position)

**Using geneinfo:**
- Use **"Getting gene information"** feature to query genes by position
- Calculate distance from each SNP to nearest gene(s)
- Determine if SNP is inside gene, in promoter, or in distal regulatory region

**Output:**
- Table: [SNP_position, nearest_gene, distance, gene_type, overlap_type]
- Save as: `results/analysis/high_significance_genes/snp_to_gene_mapping.tsv`

**Expected result:**
- ~30-43 unique genes (some SNPs may be near the same gene)
- This is my **highest-confidence candidate gene list**

**Analysis questions:**
- How many SNPs are inside genes vs regulatory regions?
- What's the distribution of distances?
- Which genes appear multiple times (multiple SNPs nearby)?

---

### **Step 2: Extract and Compare Gene Lists**
*Priority: HIGH | Estimated time: 30 min*

**Goal:** Create focused gene list and compare with Phase 3 results

**Using geneinfo:**
- Use **"Gene lists"** feature to manage gene sets
- Use **"ID conversion"** if needed to standardize gene names

**Create three gene lists:**

1. **High-confidence genes** (from Step 1)
   - Genes near 43 high-significance SNPs
   - Save as: `high_confidence_X_genes.txt`

2. **Phase 3 X genes** (already have)
   - 242 protein-coding genes from |Std iHS| ≥ 2.9
   - File: `results/analysis/gene_annotation/gene_name_lists/X_protein_coding_genes.txt`

3. **Overlap genes**
   - Genes appearing in BOTH lists
   - These have highest confidence (strong signal + in selection region)

**Comparison analysis:**
```bash
# Extract unique genes from high-significance SNP mapping
cut -f2 snp_to_gene_mapping.tsv | sort | uniq > high_confidence_X_genes.txt

# Find overlap with Phase 3 genes
comm -12 <(sort high_confidence_X_genes.txt) <(sort X_protein_coding_genes.txt) > overlap_genes.txt
```

**Questions to answer:**
- How many high-confidence genes overlap with Phase 3? (Expected: most/all)
- Are there any high-confidence genes NOT in Phase 3? (Would be surprising)
- Which genes have the strongest combined evidence?

---

### **Step 3: Retrieve Gene Information**
*Priority: MEDIUM | Estimated time: 1 hour*

**Goal:** Get detailed information about candidate genes

**Using geneinfo:**
- Use **"Getting gene information"** feature for each candidate gene
- Retrieve: function, location, aliases, disease associations

**Information to collect:**
- Gene symbol and full name
- Chromosome location (start, end, strand)
- Gene type (protein-coding, lncRNA, etc.)
- Function description
- Known disease associations
- X-inactivation status (escape vs subject)

**Output:**
- Detailed gene table: `results/analysis/high_significance_genes/gene_details.tsv`

**Focus on:**
- Genes related to **reproduction** (fertility, meiosis, gamete function, sex determination)
- Genes related to **immunity** (immune response, inflammation, pathogen defense)
- X-inactivation escapees (enriched for dosage-sensitive functions)


---

### **Step 4: Create Ideogram Visualizations**
*Priority: HIGH | Estimated time: 2 hours*

**Goal:** Show where selection signals and genes are located on X chromosome

**Using geneinfo:**
- Use **"Ideograms"** feature to create chromosome visualization
- Plot SNP positions with marker sizes proportional to significance

**Visualizations to create:**

#### A. X Chromosome Ideogram with 43 High-Significance SNPs
- X-axis: Chromosome position (0-155 Mb)
- Markers at each SNP position
- Size/color scaled by -log10 p-value
- Gene names labeled for top signals

#### B. X Chromosome Ideogram with Selection Regions
- Show all 136 selection regions from Phase 3 as bands
- Highlight regions containing high-significance SNPs
- Color by maximum |iHS| in region

#### C. Combined View: SNPs + Genes + Functional Categories
- SNP positions as points
- Gene positions as bars
- Color genes by function (immune = red, reproductive = blue, other = gray)
- Show cytogenetic bands and centromere

**Output files:**
- `results/figures/X_chromosome_ideogram_high_significance_snps.png`
- `results/figures/X_chromosome_ideogram_selection_regions.png`
- `results/figures/X_chromosome_ideogram_functional_categories.png`

**Analysis questions:**
- Are selection signals clustered or evenly distributed?
- Which cytogenetic bands show most selection?
- Do immune/reproductive genes cluster together?

---

### **Step 5: Gene Context Plots**
*Priority: MEDIUM | Estimated time: 2 hours*

**Goal:** Show detailed gene annotation context around top SNPs

**Using geneinfo:**
- Use **"Gene context: plotting data over gene annotation"** feature
- Create detailed plots for top 5-10 signals

**For each top signal, create plot showing:**
- ±200kb window around SNP
- All genes in region (with exons/introns)
- SNP position marked
- iHS values for all variants in region
- Regulatory elements if available

**Example: Top signal X:122707710**
- Plot region X:122,507,710 - 122,907,710 (±200kb)
- Show KIAA2022 gene structure
- Mark the high-significance SNP
- Show iHS track for all variants
- Label nearby genes

**Output:**
- One detailed plot per top signal
- Save as: `results/figures/gene_context/chrX_122707710_context.png`

**This helps answer:**
- Is the SNP in an exon, intron, promoter?
- Are there multiple genes nearby?
- Do nearby genes also show selection?

---

### **Step 6: GO Information - Initial Exploration**
*Priority: HIGH | Estimated time: 1 hour*

**Goal:** Understand what GO terms are associated with candidate genes

**Using geneinfo:**
- Use **"GO information"** feature to retrieve GO terms for each gene
- Use **"GO graphs"** to visualize GO term relationships

**For high-confidence gene list:**
- Retrieve all GO terms (Biological Process, Molecular Function, Cellular Component)
- Save complete GO annotation table

**Analysis:**
- Count how many genes have GO terms related to:
  - Reproduction (keywords: meiosis, gamete, fertility, gonad, sperm, oocyte)
  - Immunity (keywords: immune, inflammation, defense, cytokine, antibody)
  - Neurodevelopment (keywords: neuron, synapse, cognition, brain)
  - Metabolism (other functions)

**Output:**
- `results/analysis/functional_enrichment/gene_GO_annotations.tsv`
- Initial counts of functional categories

**Hypothesis check:**
- Are more genes immune/reproductive than expected?
- What fraction have clear functional annotations?

---

### **Step 7: GO Enrichment Analysis**
*Priority: HIGH | Estimated time: 2 hours*

**Goal:** Statistical test for enrichment of reproductive and immune functions

**Using geneinfo:**
- Use **"GO enrichment"** feature for statistical testing

**Analysis A: High-confidence genes vs all X genes**
- **Test set:** 30-43 genes from high-significance SNPs
- **Background:** All X chromosome genes (~800 protein-coding genes)
- **Question:** Are high-significance genes enriched for specific functions?

**Analysis B: High-confidence genes vs X genes under selection**
- **Test set:** 30-43 genes from high-significance SNPs
- **Background:** 242 X protein-coding genes in selection regions (Phase 3)
- **Question:** Do strongest signals show functional enrichment vs moderate signals?

**Analysis C: X genes vs autosomal genes (both under selection)**
- **Test set:** 242 X protein-coding genes in selection
- **Background:** 16,331 autosomal protein-coding genes in selection
- **Question:** Is X chromosome enriched for reproductive/immune compared to autosomes?

**For each analysis:**
- Test GO Biological Process terms
- Apply multiple testing correction (FDR < 0.05)
- Report effect sizes (fold enrichment)
- Identify which genes drive each enrichment

**Output:**
- `results/analysis/functional_enrichment/GO_enrichment_high_confidence_vs_X.tsv`
- `results/analysis/functional_enrichment/GO_enrichment_X_vs_autosomes.tsv`
- Summary table of significantly enriched terms

**Focus on these GO categories:**
- Reproductive processes (GO:0000003)
- Immune system process (GO:0002376)
- Inflammatory response (GO:0006954)
- Meiotic processes
- Sexual reproduction
- Innate/adaptive immunity

---

### **Step 8: STRING Protein Interaction Networks**
*Priority: MEDIUM | Estimated time: 2 hours*

**Goal:** Identify protein interaction networks among candidate genes

**Using geneinfo:**
- Use **"STRING networks"** feature to query protein interactions
- Create network for high-confidence genes

**Network analysis:**
- Input: Gene list from high-significance SNPs
- Retrieve interactions from STRING database
- Filter by confidence score (e.g., > 0.7 for high confidence)

**Visualizations:**
- Network graph with genes as nodes, interactions as edges
- Color nodes by function (immune, reproductive, other)
- Size nodes by -log10 p-value
- Show interaction confidence as edge thickness

**Analysis questions:**
- Do candidate genes form connected networks?
- Are immune genes connected to each other?
- Are reproductive genes connected to each other?
- Any unexpected connections between functional categories?

**Known genes to check:**
- AXTRT1, AFF2 (professor mentioned)
- KIAA2022, BRWD3 (top region)
- Any X-inactivation escapees

**Output:**
- `results/figures/STRING_network_high_confidence_genes.png`
- Network statistics (clustering coefficient, modules, hubs)
- List of highly connected genes (potential key players)

---

### **Step 9: Focused Analysis on Specific Gene Categories**
*Priority: MEDIUM | Estimated time: 3 hours*

**Goal:** Deep dive into reproductive and immune genes specifically

#### A. Reproductive Gene Analysis

**Identify reproductive genes:**
- Use GO enrichment results
- Search for keywords in gene descriptions: fertility, gonad, gamete, meiosis, sperm, oocyte, testis, ovary
- Check overlap with known reproductive gene databases

**Questions:**
- How many reproductive genes in high-confidence list?
- Where are they located on X? (clustered or dispersed?)
- What specific reproductive processes? (spermatogenesis, oogenesis, fertilization?)
- Any known fertility disorders associated?

**Create:**
- Reproductive gene ideogram
- STRING network for just reproductive genes
- Distance to gene analysis for reproductive vs non-reproductive

#### B. Immune Gene Analysis

**Identify immune genes:**
- Use GO enrichment results
- Search for keywords: immune, inflammation, defense, cytokine, antibody, T cell, B cell
- Check overlap with ImmunoDB or similar

**Questions:**
- How many immune genes in high-confidence list?
- Which immune pathways? (innate, adaptive, inflammatory?)
- Any autoimmune disease associations?
- Compare X immune genes vs autosomal immune genes under selection

**Create:**
- Immune gene ideogram
- STRING network for just immune genes
- Functional subcategory breakdown

---

### **Step 10: Distance-to-Gene Analysis and Odds Ratios**
*Priority: MEDIUM | Estimated time: 2 hours*

**Goal:** Analyze relationship between SNP distance and gene function

**Now that we have gene mappings (from Step 1):**

**Analysis A: Distance distribution**
- Plot histogram: distance from SNP to nearest gene
- Categories: inside gene, <5kb (promoter), 5-50kb (enhancer), >50kb
- Compare to genome-wide expectations

**Analysis B: Odds ratios for functional categories**
- Calculate OR: immune genes vs non-immune for being close to high-significance SNPs
- Calculate OR: reproductive genes vs non-reproductive
- Test statistical significance

**Example calculation:**
```
Contingency table:
                    Near high-sig SNP | Not near
Immune genes        |       a          |    b
Non-immune genes    |       c          |    d

OR = (a/b) / (c/d)
```

**Questions:**
- Are immune/reproductive genes MORE likely to be near high-significance SNPs?
- Does distance from SNP relate to gene function?
- Are regulatory variants more common than coding variants?

**Output:**
- `results/analysis/distance_analysis/snp_distance_distribution.png`
- `results/analysis/distance_analysis/odds_ratios_by_function.tsv`
- Statistical tests for enrichment

---

### **Step 11: X-Inactivation Escapee Analysis**
*Priority: LOW-MEDIUM | Estimated time: 1 hour*

**Goal:** Test if X-inactivation escapees are enriched in selection signals

**Background:** ~15% of X genes escape X-inactivation (from Tukiainen et al. 2017)

**Analysis:**
- From Phase 3, I annotated X-inactivation status
- Test: Are escapee genes enriched in high-confidence list?

**Contingency table:**
```
                           In high-conf list | Not in list
Escapee genes              |        a        |      b
Subject-to-inactivation    |        c        |      d
```

**Hypothesis:**
- Escapees are dosage-sensitive (2 copies in females)
- May be under different selective pressures
- Could be enriched or depleted in selection signals

**Output:**
- Odds ratio and p-value for escapee enrichment
- List of escapee genes under selection
- Ideogram showing escapee gene positions

---

### **Step 12: Create Comprehensive Summary Figures**
*Priority: HIGH | Estimated time: 3 hours*

**Goal:** Publication-quality figures summarizing all results

#### Figure 1: X Chromosome Selection Landscape
- Multi-panel figure:
  - Panel A: Full X ideogram with all 136 selection regions
  - Panel B: 43 high-significance SNPs highlighted
  - Panel C: Gene density along X chromosome
  - Panel D: Functional category distribution

#### Figure 2: Functional Enrichment Results
- Multi-panel:
  - Panel A: GO enrichment bar plot (top 20 terms)
  - Panel B: Immune vs reproductive gene counts
  - Panel C: Comparison X vs autosomes
  - Panel D: Effect sizes (fold enrichment)

#### Figure 3: STRING Network
- Main network figure showing:
  - All high-confidence genes with interactions
  - Functional categories color-coded
  - Network modules/clusters highlighted
  - Key hub genes labeled

#### Figure 4: Gene Context Examples
- Top 3-4 selection signals with detailed context:
  - Gene structure
  - SNP position
  - iHS values
  - Functional annotation

**Output:**
- `results/figures/manuscript/Figure1_X_chromosome_landscape.pdf`
- `results/figures/manuscript/Figure2_functional_enrichment.pdf`
- `results/figures/manuscript/Figure3_STRING_network.pdf`
- `results/figures/manuscript/Figure4_gene_context.pdf`

---

### **Step 13: Generate Summary Tables**
*Priority: HIGH | Estimated time: 2 hours*

**Goal:** Create comprehensive supplementary tables

#### Table 1: High-Significance SNPs
- All 43 SNPs with: position, alleles, iHS, p-value, nearest gene, distance

#### Table 2: High-Confidence Candidate Genes
- All genes near high-sig SNPs with: name, position, function, GO terms, X-inactivation status

#### Table 3: GO Enrichment Results
- All significant GO terms with: GO ID, description, p-value, FDR, fold enrichment, genes

#### Table 4: STRING Network Statistics
- Network metrics, modules, hub genes

#### Table 5: Comparison Phase 3 vs High-Confidence
- Overlap analysis, functional comparison

**Output format:**
- Excel workbook: `results/analysis/supplementary_tables.xlsx`
- Individual TSV files for each table

---

## Implementation Plan: Order of Execution

### Week 1: Core Analysis (Priority: HIGH)
- **Day 1:** Step 1 - Map SNPs to genes
- **Day 2:** Step 2 - Gene list comparison + Step 3 - Gene information
- **Day 3:** Step 4 - Ideograms
- **Day 4:** Step 6 - GO information + Step 7 - GO enrichment
- **Day 5:** Step 8 - STRING networks

### Week 2: Detailed Analysis (Priority: MEDIUM)
- **Day 6:** Step 5 - Gene context plots
- **Day 7:** Step 9 - Reproductive/immune deep dive
- **Day 8:** Step 10 - Distance and odds ratios
- **Day 9:** Step 11 - X-inactivation analysis

### Week 3: Synthesis (Priority: HIGH)
- **Day 10-11:** Step 12 - Comprehensive figures
- **Day 12:** Step 13 - Summary tables
- **Day 13-14:** Write up results section
- **Day 15:** Final review and revisions

---

## Expected Outcomes

### Key Questions to Answer:

1. **Which genes are under the strongest selection on X?**
   - Expected: 30-43 high-confidence genes from high-sig SNPs

2. **Are these genes enriched for reproductive/immune functions?**
   - Hypothesis: YES, especially compared to autosomes
   - Will quantify with GO enrichment p-values

3. **What specific reproductive/immune processes?**
   - Spermatogenesis? Oogenesis? Innate immunity? Adaptive immunity?

4. **Do selection signals cluster on X chromosome?**
   - Ideograms will show spatial distribution

5. **Do candidate genes form functional networks?**
   - STRING analysis will reveal interactions

6. **How do high-significance genes compare to Phase 3 genes?**
   - Overlap analysis will show consistency

7. **Are X-inactivation escapees enriched?**
   - Statistical test will determine

### Deliverables:

1. **Manuscript results section** with:
   - Gene identification results
   - Functional enrichment statistics
   - Network analysis
   - Key candidate genes

2. **Publication-quality figures:**
   - 4 main figures
   - Multiple supplementary figures

3. **Supplementary tables:**
   - 5 comprehensive data tables

4. **Updated progress report** documenting Phase 4 completion

---

## Technical Notes

### Data Files Organization:

```
results/analysis/
├── high_significance_genes/
│   ├── snp_to_gene_mapping.tsv
│   ├── high_confidence_X_genes.txt
│   ├── overlap_genes.txt
│   └── gene_details.tsv
├── functional_enrichment/
│   ├── gene_GO_annotations.tsv
│   ├── GO_enrichment_high_confidence_vs_X.tsv
│   └── GO_enrichment_X_vs_autosomes.tsv
├── distance_analysis/
│   ├── snp_distance_distribution.png
│   └── odds_ratios_by_function.tsv
└── supplementary_tables.xlsx

results/figures/
├── X_chromosome_ideogram_high_significance_snps.png
├── X_chromosome_ideogram_selection_regions.png
├── X_chromosome_ideogram_functional_categories.png
├── STRING_network_high_confidence_genes.png
├── gene_context/
│   ├── chrX_122707710_context.png
│   └── [other top SNPs]
└── manuscript/
    ├── Figure1_X_chromosome_landscape.pdf
    ├── Figure2_functional_enrichment.pdf
    ├── Figure3_STRING_network.pdf
    └── Figure4_gene_context.pdf
```

### geneinfo Usage Tips:

1. **Use caching:** First time queries may be slow, but geneinfo caches results for faster subsequent access

2. **Batch queries:** Where possible, query multiple genes at once rather than one-by-one

3. **ID consistency:** Standardize on one gene ID system (gene symbols recommended) early

4. **Documentation:** Keep track of which geneinfo functions I use and their parameters

---

## Questions for Professor

Before starting:

1. **geneinfo installation/access:**
   - How do I install/load the geneinfo module?
   - Any specific configuration needed?
   - Documentation location?

2. **HDF5 file:**
   - Should I create an HDF5 file for my data?
   - Or does geneinfo handle data storage internally?

3. **Background set for GO enrichment:**
   - Does geneinfo allow custom backgrounds?
   - Recommend using all X genes or X genes in selection?

4. **STRING parameters:**
   - Recommended confidence threshold?
   - Include physical interactions only or also functional?

5. **Priority:**
   - Which analyses are most important?
   - Any specific genes/pathways to focus on?

---

## Success Metrics

This Phase 4 analysis will be successful if:

1. ✅ Identify 30-50 high-confidence candidate genes
2. ✅ Show statistically significant enrichment for reproductive/immune GO terms (p < 0.05 after correction)
3. ✅ Create clear visualizations showing selection signal distribution
4. ✅ Identify functional networks among candidate genes
5. ✅ Compare X vs autosomal selection patterns
6. ✅ Produce publication-quality figures and tables

**Timeline:** 3 weeks from geneinfo access to complete analysis

**Next immediate step:** Learn how to use geneinfo module and run Step 1 (SNP to gene mapping)

---

## References to Consult

- geneinfo module documentation (from professor)
- STRING database documentation
- Gene Ontology documentation (geneontology.org)
- Tukiainen et al. 2017 (X-inactivation escapees)
- My Phase 1-3 analysis results

---

## Notes to Self

- This is the core hypothesis-testing phase - everything leads to answering whether X genes are enriched for reproductive/immune functions
- geneinfo module integrates many analyses I was planning to do separately - huge time saver!
- Need to be careful about multiple testing correction across many GO terms
- Distance analysis will be key for understanding regulatory vs coding selection
- STRING networks may reveal unexpected connections between genes
- Keep autosomes as comparison/background throughout
