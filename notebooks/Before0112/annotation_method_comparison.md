# Comparison: Phase 3 GENCODE vs High-Sig SNP Geneinfo Annotation

## Summary

**You are correct!** The geneinfo annotation (used for high-significance SNPs) is **much less strict distance-wise** than the GENCODE annotation used in Phase 3. This explains why KLHL13 and PCDH11X were identified by the high-sig SNP analysis but not by Phase 3.

---

## Key Differences

| Feature | Phase 3 GENCODE | High-Sig SNP Geneinfo |
|---------|-----------------|----------------------|
| **Search Window** | ±50 kb flank | ±500 kb window |
| **Approach** | Region-based overlap | Nearest gene to SNP |
| **Distance Limit** | 50 kb from region boundaries | 500 kb from SNP position |
| **Annotation Source** | GENCODE v44 GTF | geneinfo package |
| **Purpose** | Find genes in/near selection regions | Map individual SNPs to genes |

---

## Detailed Comparison

### Phase 3 GENCODE Annotation
**Script:** `03_annotate_regions_with_genes.py`

```python
def find_overlapping_genes(region_chr, region_start, region_end, genes_df, flank=50000):
    """
    Find genes overlapping a genomic region

    Args:
        flank: Additional flanking region (default: 50kb)
    """
    # Add flanking regions
    search_start = max(0, region_start - flank)
    search_end = region_end + flank

    # Find overlapping genes
    overlapping = chr_genes[
        (chr_genes['gene_end'] >= search_start) &
        (chr_genes['gene_start'] <= search_end)
    ]
```

**Key characteristics:**
- **±50 kb flank** around each selection region
- Only genes that **overlap** the extended region are included
- Conservative approach: genes must be relatively close to the selection signal
- Uses GENCODE v44 (gold standard annotation)

---

### High-Significance SNP Geneinfo Annotation
**Script:** `HighSigSNP_tonearestGenes.py`

```python
# Configuration
SEARCH_WINDOW = 500000  # Search ±500kb around each SNP
PROMOTER_DISTANCE = 2000  # Promoter defined as <2kb from TSS
ENHANCER_DISTANCE = 50000  # Enhancer/regulatory region <50kb

def get_genes_in_region(chrom, pos, window=SEARCH_WINDOW, assembly=ASSEMBLY):
    """
    Get all genes within a window around a position.

    Parameters:
        window : int
            Window size (±window around pos)
    """
    start = max(0, pos - window)
    end = pos + window

    # Query genes using geneinfo
    genes = gene_coords_region(chrom=chrom, start=start, end=end, assembly=assembly)
```

**Key characteristics:**
- **±500 kb search window** around each SNP (10x larger!)
- Finds the **nearest gene** even if very far away
- Liberal approach: ensures every SNP gets mapped to a gene
- Uses geneinfo package (likely based on RefSeq or similar)

---

## Why KLHL13 and PCDH11X Were Missed by Phase 3

### KLHL13
- **Gene position:** X:117,897,812-118,117,340
- **High-sig SNP:** X:117,424,856 (-log10 p = 6.18)
- **Distance:** 472,956 bp = **473 kb downstream**

**Phase 3 perspective:**
- Selection regions are defined by **sustained high iHS** across multiple SNPs
- KLHL13 was NOT in a Phase 3 selection region because:
  - The gene itself doesn't have multiple SNPs with high iHS
  - There's no 500kb selection region that extends 50kb to reach KLHL13
  - The single high-sig SNP is too far away to create a selection region that includes the gene

**High-sig SNP perspective:**
- The 500kb search window captures KLHL13 as the nearest gene
- Even though it's 473 kb away!

---

### PCDH11X
- **Gene position:** X:91,779,374-92,623,230
- **High-sig SNP:** X:93,062,356 (-log10 p = 7.18)
- **Distance:** 439,126 bp = **439 kb upstream**

**Phase 3 perspective:**
- Same reasoning as KLHL13
- The SNP is too far from the gene to be in the same selection region

**High-sig SNP perspective:**
- The 500kb search window captures PCDH11X as the nearest gene
- Even though it's 439 kb away!

---

## Implications

### The 10x Difference Matters

The **500 kb search window** (high-sig SNP) vs **50 kb flank** (Phase 3) is a **10-fold difference** in search distance.

**Visual representation:**

```
Phase 3 (±50 kb):
    ----[50kb]--[Selection Region]--[50kb]----
                     (genes here)

High-sig SNP (±500 kb):
[500kb]--------[SNP]--------[500kb]
                (any gene in this 1 Mb window)
```

### Why the Difference Exists

**Phase 3 logic:**
- Goal: Find genes likely under selection
- Method: Identify regions with **sustained** high iHS signals
- Conservative: Only include genes close to selection regions
- Reasoning: True selection targets should have multiple nearby SNPs with high iHS

**High-sig SNP logic:**
- Goal: Map every significant SNP to a gene for interpretation
- Method: Find **nearest gene** to each SNP
- Liberal: Use large window to ensure every SNP gets a gene assignment
- Reasoning: Even distal SNPs might affect gene regulation (long-range enhancers)

---

## Which Approach is "Better"?

Both approaches have merit for different purposes:

### Phase 3 (Conservative, Region-Based)
✅ **Better for:**
- Identifying high-confidence selection targets
- Genes with sustained selection signals
- Minimizing false positives

❌ **Limitation:**
- May miss genes affected by distal regulatory SNPs

### High-Sig SNP (Liberal, SNP-Based)
✅ **Better for:**
- Ensuring every significant SNP is interpretable
- Capturing potential long-range regulatory effects

❌ **Limitation:**
- Many "nearest gene" assignments may be biologically irrelevant
- High false positive rate for distal assignments (>400 kb)

---

## Biological Plausibility of Distal Associations

**Question:** Can a SNP 400-500 kb away truly affect a gene?

**Short answer:** Possible but unlikely for most cases.

**Long answer:**
- **Enhancers** can act at long distances (50-200 kb typical, up to 1 Mb in some cases)
- **Topologically Associating Domains (TADs)** can bring distant regions together
- **BUT** most functional variants are within 50 kb of their target genes
- Distances >400 kb are **statistically unlikely** to be causal without additional evidence

**For KLHL13 and PCDH11X:**
- Without chromatin conformation (Hi-C) or eQTL data, we cannot know if the distant SNPs truly affect these genes
- More likely: The high-sig SNPs are affecting **other genes** or **regulatory elements** closer to them
- The "nearest gene" may not be the **target gene**

---

## Conclusion

**Yes, the geneinfo annotation is much less strict distance-wise (10x larger search window).**

This explains why:
1. KLHL13 and PCDH11X appear in high-sig SNP results but not Phase 3
2. Both are >400 kb from their associated SNPs
3. Phase 3's 50 kb flank wouldn't reach them
4. High-sig SNP's 500 kb window captures them as "nearest gene"

**Recommendation:**
Focus on the **11 overlap genes** for downstream analysis, as they have:
- Close proximity to high-significance SNPs
- Sustained selection signals (Phase 3 regions)
- Much higher biological plausibility

The 2 non-overlapping genes (KLHL13, PCDH11X) likely represent **false positives** from the liberal "nearest gene" approach.

---

## Files for Reference

**Phase 3 annotation:**
- Script: `/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/scripts/03_annotate_regions_with_genes.py`
- Output: `/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/analysis/gene_annotation/X_chromosome_genes_in_selection.tsv`
- Flank: **50 kb**

**High-sig SNP annotation:**
- Script: `/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/scripts/HighSigSNP_tonearestGenes.py`
- Output: `/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/analysis/high_significance_genes/snp_to_gene_mapping.tsv`
- Window: **500 kb (±500 kb = 1 Mb total)**
