# High-Confidence Candidate Gene Analysis Summary

**Date:** October 29, 2025
**Analysis:** Step 1 Complete - SNP to Gene Mapping
**Student:** MIT van Bruggen

---

## Overview

This analysis mapped 43 high-significance SNPs (-log10 p ≥ 6) from X chromosome iHS selection scans to their nearest genes using the geneinfo module. We identified **19 unique high-confidence candidate genes** and retrieved detailed functional information for each.

---

## Key Findings

### 1. SNP Location Distribution

- **Inside genes**: 14 SNPs (32.6%)
- **Regulatory regions** (<50kb): 12 SNPs (27.9%)
- **Distal** (>50kb): 17 SNPs (39.5%)

**Interpretation:** ~60% of high-significance SNPs are in or near genes (within 50kb), suggesting functional relevance rather than random drift.

### 2. Distance to Nearest Gene

For SNPs outside genes:
- **Median distance**: 132,715 bp (~133 kb)
- **Mean distance**: 192,858 bp (~193 kb)
- **Range**: 2,159 bp to 472,956 bp

The lower median vs mean indicates most SNPs are relatively close to genes, with a few distant outliers.

### 3. Genes with Multiple High-Significance SNPs

**10 genes have multiple nearby high-significance SNPs:**

| Gene | # SNPs | Region | Notes |
|------|--------|--------|-------|
| LOC101928201 | 7 | chrX:4.6M | Non-coding RNA |
| LOC124900496 | 6 | chrX:82.5M | Small nucleolar RNA |
| **DACH2** | 4 | chrX:86.1-86.8M | Developmental transcription factor |
| **HTR2C** | 4 | chrX:114.5-114.9M | Serotonin receptor |
| **FAAH2** | 3 | chrX:57.2-57.4M | Fatty acid metabolism |
| ZMAT1 | 2 | chrX:101.8-101.9M | Zinc finger protein |
| NAP1L2 | 2 | chrX:73.2M | Neuronal chromatin regulator |
| **DIAPH2** | 2 | chrX:96.6-97.6M | Ovarian development |
| TRPC5 | 2 | chrX:111.7-112.0M | Calcium channel |
| GPC3 | 2 | chrX:133.5-133.9M | Cell surface proteoglycan |

---

## Functional Classification of 19 Candidate Genes

### Gene Type Distribution

- **Protein-coding**: 10 genes (53%)
- **Non-coding RNA**: 5 genes (26%)
  - 3 lncRNAs/ncRNAs
  - 2 snoRNAs
- **Unknown**: 4 genes (21%)

### Functional Categories

Based on gene summaries and descriptions:

| Category | Count | Genes |
|----------|-------|-------|
| **Immune** | 3 | IL1RAPL2, TRPC5, KLHL13 |
| **Neural** | 4 | NAP1L2, IL1RAPL2, TRPC5, HTR2C |
| **Developmental** | 3 | DACH2, PCDH11X, DIAPH2 |
| **Other/Unknown** | 11 | (includes LOC genes, metabolic, structural) |

**Note:** Some genes have multiple functional roles (e.g., IL1RAPL2 is both immune and neural).

---

## Top 5 Priority Candidate Genes

### 1. **DIAPH2** - Diaphanous Related Formin 2
- **Location**: chrX:96,684,712-97,604,997 (920 kb)
- **Function**: Formin family protein involved in cytoskeleton organization
- **Key Evidence**:
  - **Explicitly involved in ovarian development and function**
  - **Linked to premature ovarian failure 2 (POF2A)**
  - Has **2 high-significance SNPs** nearby
- **Relevance**: **STRONGEST REPRODUCTIVE CANDIDATE** - direct link to female fertility
- **Hypothesis**: Selection on this gene could reflect adaptation in reproductive timing or fertility

---

### 2. **IL1RAPL2** - Interleukin 1 Receptor Accessory Protein Like 2
- **Location**: chrX:104,566,199-105,767,829 (1.2 Mb)
- **Function**: Interleukin 1 receptor family member
- **Key Evidence**:
  - Member of **immune signaling** pathway (interleukin receptors)
  - Similar to IL1RAPL1, which causes **intellectual disability** when mutated
  - Expressed in brain and immune tissues
- **Relevance**: Links **immunity and neural development** on X chromosome
- **Hypothesis**: Could be under selection for immune response or cognitive traits

---

### 3. **HTR2C** - 5-Hydroxytryptamine Receptor 2C (Serotonin Receptor)
- **Location**: chrX:114,584,078-114,910,061 (326 kb)
- **Function**: G-protein coupled receptor for serotonin (5-HT)
- **Key Evidence**:
  - Major **neurotransmitter receptor** affecting mood, appetite, behavior
  - Subject to **RNA editing** creating multiple functional isoforms
  - Has **4 high-significance SNPs** clustered around this gene
- **Relevance**: Selection could affect **behavior, mood regulation, or feeding**
- **Hypothesis**: Behavioral adaptation or sexual selection on temperament

---

### 4. **TRPC5** - Transient Receptor Potential Cation Channel Subfamily C Member 5
- **Location**: chrX:111,768,011-112,082,776 (315 kb)
- **Function**: Calcium-permeable cation channel
- **Key Evidence**:
  - Receptor-activated calcium channel
  - Involved in **immune cell signaling** and **neuronal function**
  - Has **2 high-significance SNPs** (one inside gene)
- **Relevance**: Dual role in immunity and neural signaling
- **Hypothesis**: Selection on calcium signaling pathways affecting multiple systems

---

### 5. **DACH2** - Dachshund Family Transcription Factor 2
- **Location**: chrX:86,148,451-86,832,604 (684 kb)
- **Function**: Transcription factor involved in developmental cell fate determination
- **Key Evidence**:
  - Homolog of Drosophila dachshund (controls eye, limb, **genital disc** development)
  - Developmental regulator
  - Has **4 high-significance SNPs** nearby
- **Relevance**: Developmental control, potentially including **reproductive organ development**
- **Hypothesis**: Selection on developmental timing or organ morphology

---

## Additional Notable Candidates

### **NAP1L2** - Nucleosome Assembly Protein 1 Like 2
- **Function**: Chromatin regulator specific to **neuronal cell proliferation**
- **Relevance**: Neural development, could affect brain size/function
- Small gene (2.9 kb), 2 SNPs nearby

### **FAAH2** - Fatty Acid Amide Hydrolase 2
- **Function**: Metabolizes endocannabinoids and other bioactive lipids
- **Relevance**: Has **3 high-significance SNPs** (2 inside gene), strong selection signal
- Lipid metabolism role unclear in selection context

### **PCDH11X** - Protocadherin 11 X-linked
- **Function**: Cell adhesion molecule, developmental
- **Relevance**: Large gene (844 kb), important in brain development
- Human-specific protocadherin with role in neural circuits

### **KDM6A** - Lysine Demethylase 6A
- **Function**: Histone demethylase (epigenetic regulator)
- **Relevance**: X-inactivation escapee, involved in chromatin remodeling
- Could affect gene expression genome-wide

---

## Genes with Uncertain/Unknown Function

Several LOC genes show strong selection signals but lack functional annotation:
- **LOC101928201** (7 SNPs) - non-coding RNA
- **LOC124900496** (6 SNPs) - snoRNA
- **LOC124905292** (top SNP) - snoRNA

These could be regulatory RNAs with important but uncharacterized functions.

---

## Evidence for Reproductive vs Immune Enrichment

### Reproductive Function: **MODERATE-STRONG Evidence**
- **DIAPH2**: Direct evidence (ovarian failure)
- **DACH2**: Indirect evidence (genital disc development in flies)
- 2 of 19 genes (11%) have clear reproductive links

### Immune Function: **MODERATE Evidence**
- **IL1RAPL2**: Interleukin receptor family
- **TRPC5**: Immune cell calcium signaling
- **KLHL13**: Ubiquitin ligase (protein regulation)
- 3 of 19 genes (16%) have immune-related functions

### Neural/Behavioral Function: **STRONG Evidence**
- **HTR2C**: Serotonin receptor (4 SNPs)
- **NAP1L2**: Neuronal chromatin regulation
- **IL1RAPL2**: Neural development (intellectual disability)
- **TRPC5**: Neuronal calcium signaling
- 4 of 19 genes (21%) have neural/behavioral roles

### Interpretation:
X chromosome selection signals show enrichment for:
1. **Neural/behavioral genes** (strongest signal)
2. **Immune function** (moderate signal)
3. **Reproductive function** (moderate signal, but with high-confidence candidate DIAPH2)

This could reflect:
- Sexual selection on behavior/cognition
- Adaptation to pathogens (immunity)
- Reproductive timing or fertility optimization

---

## Next Steps (Step 2-3 from Analysis Plan)

### Immediate Priorities:

1. **Compare with Phase 3 gene list**
   - How many of these 19 genes overlap with the 242 protein-coding genes from |Std iHS| ≥ 2.9?
   - Are high-significance genes consistently in broader selection regions?

2. **Gene Ontology (GO) enrichment analysis**
   - Test statistical enrichment for:
     - Reproductive processes (GO:0000003)
     - Immune system process (GO:0002376)
     - Nervous system development
   - Use all X genes as background

3. **Literature search for key candidates**
   - DIAPH2 + premature ovarian failure
   - HTR2C + behavior/sexual selection
   - IL1RAPL2 + X-linked cognitive disorders

4. **X-inactivation status**
   - Check if any candidates escape X-inactivation
   - KDM6A is known escapee - are others?

5. **STRING protein-protein interaction network**
   - Do these 19 genes interact?
   - Are they part of common pathways?

---

## Files Generated

| File | Description |
|------|-------------|
| `snp_to_gene_mapping.tsv` | All 357 SNP-gene associations |
| `snp_to_nearest_gene_only.tsv` | 43 nearest gene mappings |
| `high_confidence_X_genes.txt` | 19 unique gene names |
| `gene_details_complete.tsv` | Detailed gene information with functions |

---

## Summary Statement

**We identified 19 high-confidence candidate genes near the strongest X chromosome selection signals in our population. The top candidates include DIAPH2 (ovarian function), HTR2C (serotonin receptor), and IL1RAPL2 (immunity + neural development). These genes suggest selection on reproductive function, behavioral traits, and immune response on the X chromosome.**

**Priority #1: DIAPH2 - strongest evidence for reproductive selection (premature ovarian failure)**

---

## Questions for Discussion

1. Is the enrichment of neural/behavioral genes surprising? Could this reflect sexual selection?

2. DIAPH2 (ovarian failure) - could selection be acting on:
   - Age at menopause?
   - Fertility rate?
   - Egg quality?

3. Why are 4 SNPs clustered around HTR2C? Is this one selected variant or multiple?

4. Should we focus on protein-coding genes or also investigate the LOC non-coding RNAs?

5. How do these X chromosome candidates compare to autosomal selection signals?

---

**Analysis by:** MIT van Bruggen
**Date:** October 29, 2025
**Phase:** 4, Step 1 Complete
**Next:** Step 2 - Gene list comparison and GO enrichment
