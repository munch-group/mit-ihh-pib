# Analysis Notes: High-Significance SNPs and geneinfo Workflow
**Date:** October 28, 2025
**Student:** MIT van Bruggen
**Context:** Meeting notes + planning next analysis steps

---

## Background

After completing Phase 3 (gene annotation using GENCODE with |Std iHS| ≥ 2.9), I had a meeting with K to discuss next steps. He gave me some notes that I needed to decode and understand how they fit into my overall analysis plan.

---

## Meeting notes

After meeting:
> "Just had a look at your file. I would start with the SNP with a -log10 pvalue above 6"
>
> ```bash
> cat ../vanbruggenmit/mit-ihh-pib/results/analysis/candidates/per_chromosome/chrX_candidates_liberal.tsv | awk '{ if ($11 >= 6 || $1 == "Location" ) print }'
> ```
>
> "And then use the geneinfo package to find the genes nearby"

Notes from the meeting:
- **Plots of odds ratios** - for closeness to gene
- **cut -f 6 ... | sort | uniq** - extract unique gene names
- **hdf5 file** - "p-values for X chromosome will be in data (last column, left one is the gene coordinates)"
- **geneinfo** - check for gene overlaps
- **Plot all significant genes** → ideograms → sizes relating to odds ratios
- **GO information and GO enrichment** - define background set (fertility or immunity)
- **STRING interaction network**
- **Gene examples:** AXTRT1, AFF2, AGG

---

## What I Did: Extracting High-Significance SNPs

### Command executed:
```bash
cat /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/analysis/candidates/per_chromosome/chrX_candidates_liberal.tsv | awk '{ if ($11 >= 6 || $1 == "Location" ) print }' > /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/analysis/candidates/chrX_high_significance.tsv
```

### Results:
- **43 SNPs** on X chromosome with -log10 pvalue ≥ 6
- Saved to: `results/analysis/candidates/chrX_high_significance.tsv`
- These SNPs have |Std iHS| values around **5.6-5.7** (much higher than the 2.9 threshold used in Phase 3!)

### Top signals:
| Location | Std iHS | -log10 p-value | Position |
|----------|---------|----------------|----------|
| X:122707710:TA:T | 5.74 | 8.02 | 122,707,710 |
| X:57480530:G:A | 5.68 | 7.87 | 57,480,530 |
| X:57482305:A:G | 5.68 | 7.87 | 57,482,305 |
| X:111803779:A:G | 5.59 | 7.63 | 111,803,779 |

The strongest signal (X:122707710) is in the region I identified in Phase 3 as containing KIAA2022 (neurodevelopmental disorders gene).

---

## Understanding the Workflow: Phase 3 vs geneinfo 
Just recapping so I don't get confused

### Two different approaches to gene mapping:

| Aspect | **Phase 3 (GENCODE)** | **Current (geneinfo)** |
|--------|----------------------|----------------------|
| **Input** | Selection regions (500kb windows) | Individual high-significance SNPs |
| **Threshold** | \|Std iHS\| ≥ 2.9 | -log10 p ≥ 6 (\|iHS\| ≥ ~5.6) |
| **Method** | Find ALL genes in/near regions (±50kb) | Find NEAREST gene(s) to each SNP |
| **X Results** | 136 regions → 664 genes (242 protein-coding) | 43 SNPs → ? genes (much fewer) |
| **Purpose** | Broad functional enrichment | Focused candidate gene identification |

### Why both approaches?

Think of it as **zooming in**?:

1. **Phase 3 (broad)**: "This 410kb region shows selection, here are all 10 genes in it"
   - Good for: Genome-wide GO enrichment, pathway analysis
   - Question: "Overall, what functions are enriched in selection regions?"

2. **geneinfo (precise)**: "This specific SNP at position 122707710 is near which gene?"
   - Good for: Focused functional analysis, STRING networks, candidate validation
   - Question: "What specific genes are driving the strongest selection signals?"

The geneinfo approach will give:
- **Distance** from SNP to gene (important for understanding mechanism!)
- **Causative gene candidates** - which gene is most likely driving selection
- **Regulatory context** - is SNP in coding sequence, promoter, enhancer?
- **Higher confidence** - focusing on strongest signals only

---

## Why Finding Genes Nearby Matters

### 1. SNPs don't have biological meaning by themselves
- A position like `X:122707710` doesn't tell us *what* is under selection
- Need to map to genes to understand biology and function

### 2. Distance tells something about the mechanism

| Distance from SNP to Gene | Likely Mechanism |
|---------------------------|------------------|
| **0 bp (inside gene)** | Coding variant, splice site, or intronic regulatory |
| **< 5 kb** | Promoter or proximal regulatory element |
| **5-50 kb** | Distal enhancer/silencer |
| **> 50 kb** | May not be the actual target gene |

This is why we analyze "closeness to gene"!

### 3. Focusing on strongest signals = highest confidence
- p-value < 10^-6 represents my **most confident** selection signals
- Will create a smaller, focused gene list for detailed follow-up
- More defensible results than using broader thresholds

### 4. Tests my core hypothesis
My research question:
> "Are genes under recent positive selection on the human X chromosome enriched for **reproductive and immune functions**?"

The workflow is:
1. Find strongest selection signals (43 SNPs, -log10 p ≥ 6)
2. → Map to nearby genes using geneinfo
3. → Check if those genes are enriched for reproductive/immune functions (GO enrichment)
4. → Validate with protein interaction networks (STRING)
5. → Create figures (ideograms showing where they are on X chromosome)

---

## Why I Analyzed All Chromosomes (Not Just X)

Initially I wondered today : if I'm only using X chromosome, why did I analyze chr1-22?

### Reasons autosomes are essential:

#### A. Quality Control & Validation
I need to verify iHS standardization worked correctly:
- All 23 chromosomes showed Mean ≈ 0, SD ≈ 1
- Proves my statistical assumptions are met

#### B. Discovered X Chromosome Depletion
Key finding: X has 30% fewer selection signals than autosomes (OR = 0.698)
- **Can't know X is unusual without comparing to autosomes!**

#### C. Determined X-Specific Threshold
- X chromosome: 99th percentile = 2.90
- Autosomes: 99th percentile = 3.3-3.8
- Used this to justify chromosome-specific thresholds

#### D. Background for GO Enrichment (future analysis)
For GO enrichment, I'll compare:
- **Test set**: X chromosome selected genes
- **Background set**: All genes under selection (autosomes + X)

This controls for selection-specific biases better than using "all genes" as background.

#### E. the Hypothesis Implicitly Requires Comparison
My question is really: "Are X genes more enriched for reproductive/immune functions **than autosomal genes under selection**?"

Need autosomes for proper statistical testing!

---

## Next Steps: Waiting for geneinfo Module

### Current Status:
- Extracted 43 high-significance SNPs (-log10 p ≥ 6)
- Data ready in `chrX_high_significance.tsv`

### Once I have geneinfo:
1. **Map SNPs to genes**
   - Input: 43 SNP positions
   - Output: Nearest gene(s) + distances

2. **Analyze distance distributions**
   - How many SNPs are in genes vs regulatory regions?
   - Plotting distance from SNP to nearest gene

3. **Extract unique genes**
   ```bash
   cut -f 6 geneinfo_output.tsv | sort | uniq
   ```

4. **Compare with Phase 3 results**
   - Which genes appear in both lists? (highest confidence!)
   - Are geneinfo genes a subset of Phase 3 genes?

5. **GO enrichment analysis**
   - Test for reproductive/immune function enrichment
   - Use autosomal selected genes as background

6. **Create visualizations**
   - Ideograms showing SNP positions on X chromosome
   - Marker sizes proportional to odds ratios or p-values
   - STRING protein interaction networks

7. **Focus on specific genes**
   - AXTRT1, AFF2, AGG (genes professor mentioned)
   - Check if these appear in my results

---

## Understanding "hdf5 file" Note

hdf5 file mentioned: "p-values for X chromosome will be in data (last column, left one is the gene coordinates)"

---


