# iHS Analysis Pipeline for X Chromosome Selection Project

## Current Status (Updated 2025-10-23)
✓ iHS calculated for all chromosomes (chr1-22, chrX)
✓ Results formatted consistently across chromosomes
✓ Data quality validated
✓ **Phase 1 COMPLETE:** Candidates identified at 3 threshold levels
✓ **Phase 2 COMPLETE:** Selection regions defined with 500kb clustering
✓ **X chromosome depletion documented:** ~30% fewer signals than autosomes

## Completed Work

### Phase 1: Identify Selection Candidates ✓ COMPLETE

#### 1.1 Candidate Identification Results

**Script:** [scripts/01_identify_ihs_candidates.py](../scripts/01_identify_ihs_candidates.py)

```bash
sbatch scripts/run_01_identify_candidates.slurm
```

**Approach:**
- Processes each chromosome separately (memory efficient)
- Calculates p-values from standardized iHS (assuming N(0,1))
- Computes empirical percentiles (99th, 99.5th, 99.9th)
- Applies literature-based thresholds
- Tests for X chromosome enrichment/depletion

#### 1.2 Key Findings

**Normality Check:** ✓ EXCELLENT
- All chromosomes show mean Std iHS ≈ 0 (within 10⁻⁹)
- All chromosomes show SD Std iHS ≈ 1.000 (within 0.000006)
- Standardization is working correctly!

**X Chromosome Depletion:** Highly Significant
- X chromosome: 1.90% candidates (8,121/428,414)
- Autosomes: 2.70% candidates (202,231/7,503,286)
- Odds ratio: 0.698 (30% depletion)
- P-value: 1.71×10⁻²⁴²
- See detailed analysis: [notebooks/X_chromosome_depletion_analysis.md](X_chromosome_depletion_analysis.md)

**Empirical Percentiles:**
- Autosomes: 99th %ile ≈ 3.3-3.8, 99.9th ≈ 5.2-5.7
- X chromosome: 99th %ile = 2.90, 99.9th = 4.05

#### 1.3 Thresholds Used (Literature-Based)

Based on:
- Voight et al. (2006): Standardized iHS ~ N(0,1) under neutrality
- Paul et al. (2024): |Std iHS| ≥ 2.5 recommended
- Salazar-Tortosa et al. (2023): Check empirical distributions

**Applied thresholds:**
- **Liberal:** |Std iHS| ≥ 2.0 → 395,725 candidates (5.0% of genome)
- **Moderate:** |Std iHS| ≥ 2.5 → 210,352 candidates (2.7% of genome)
- **Stringent:** |Std iHS| ≥ 3.0 → 117,104 candidates (1.5% of genome)

#### 1.4 Output Files Created ✓

**Per-chromosome candidates:**
- [results/analysis/candidates/per_chromosome/](../results/analysis/candidates/per_chromosome/) (69 files)

**Combined genome-wide candidates:**
- [results/analysis/candidates/ALL_candidates_liberal.tsv](../results/analysis/candidates/ALL_candidates_liberal.tsv) (395,725 variants)
- [results/analysis/candidates/ALL_candidates_moderate.tsv](../results/analysis/candidates/ALL_candidates_moderate.tsv) (210,352 variants)
- [results/analysis/candidates/ALL_candidates_stringent.tsv](../results/analysis/candidates/ALL_candidates_stringent.tsv) (117,104 variants)

**Statistics and diagnostics:**
- [results/analysis/candidates/chromosome_statistics.tsv](../results/analysis/candidates/chromosome_statistics.tsv)
- [results/analysis/candidates/empirical_percentiles.tsv](../results/analysis/candidates/empirical_percentiles.tsv)
- [results/analysis/candidates/x_enrichment_test.tsv](../results/analysis/candidates/x_enrichment_test.tsv)

### Phase 2: Define Selection Regions ✓ COMPLETE

#### 2.1 Region Clustering Results

**Script:** [scripts/02_define_selection_regions.py](../scripts/02_define_selection_regions.py)

```bash
sbatch scripts/run_02_define_regions.slurm
```

**Approach:**
- Uses 500kb sliding windows
- Minimum 2 SNPs per region (Liu et al. 2013 recommendation)
- Merges overlapping windows
- **X-specific threshold:** Uses empirical 99th percentile (2.90) instead of fixed 2.5
- Accounts for X chromosome depletion

#### 2.2 Key Results

**Total Regions Identified:** 4,666
- **Autosomes:** 4,530 regions
  - Mean SNPs per region: 44.6
  - Median SNPs per region: 36
  - Mean region length: 413.9 kb
  - Mean max |Std iHS|: 5.24

- **X chromosome:** 136 regions
  - Mean SNPs per region: 31.2
  - Median SNPs per region: 18
  - Mean region length: 225.4 kb
  - Mean max |Std iHS|: 4.16

**Top 10 Selection Regions Genome-Wide:**
1. chr6:15,966,140-16,230,744 (265 kb, 46 SNPs, |Std iHS|=10.46)
2. chr11:68,207,213-68,706,406 (499 kb, 65 SNPs, |Std iHS|=10.01)
3. chr16:79,287,118-79,749,337 (462 kb, 148 SNPs, |Std iHS|=9.96)
4. chr4:77,291,111-77,744,353 (453 kb, 84 SNPs, |Std iHS|=9.89)
5. chr1:86,715,206-87,152,524 (437 kb, 51 SNPs, |Std iHS|=9.85)

**Top 5 X Chromosome Regions:**
1. chrX:122,303,906-122,714,169 (410 kb, 58 SNPs, |Std iHS|=5.74)
2. chrX:57,475,958-57,511,540 (36 kb, 34 SNPs, |Std iHS|=5.68)
3. chrX:111,798,727-112,124,702 (326 kb, 6 SNPs, |Std iHS|=5.59)
4. chrX:3,983,501-4,293,152 (310 kb, 88 SNPs, |Std iHS|=5.56)
5. chrX:97,195,776-97,479,943 (284 kb, 67 SNPs, |Std iHS|=5.42)

#### 2.3 Output Files Created ✓

**Region files:**
- [results/analysis/regions/selection_regions_500kb_min2snps.tsv](../results/analysis/regions/selection_regions_500kb_min2snps.tsv) (4,666 regions)
- [results/analysis/regions/selection_regions_500kb_min3snps.tsv](../results/analysis/regions/selection_regions_500kb_min3snps.tsv) (4,630 regions, stringent)
- [results/analysis/regions/selection_regions_autosomes.tsv](../results/analysis/regions/selection_regions_autosomes.tsv) (4,530 regions)
- [results/analysis/regions/selection_regions_X_chromosome.tsv](../results/analysis/regions/selection_regions_X_chromosome.tsv) (136 regions)

**Genome browser files:**
- [results/analysis/regions/selection_regions.bed](../results/analysis/regions/selection_regions.bed)
- [results/analysis/regions/selection_regions_stringent.bed](../results/analysis/regions/selection_regions_stringent.bed)

### Phase 3: Annotate Regions with Genes (NEXT)

#### 3.1 Gene Annotation Strategy
Extract genes overlapping or near selection regions:
- Use Ensembl BioMart or GENCODE annotations
- Include genes within ±50kb of region boundaries
- Annotate with gene biotype, description

#### 3.2 Functional Categorization
Categorize genes by function:
- **Reproductive:** GO terms for reproduction, fertility, pregnancy
- **Immune:** GO terms for immune response, inflammation
- **Other functional categories**

#### 3.3 Implementation Options

**Option A: Use biomaRt (R)**
```R
library(biomaRt)
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

genes <- getBM(
  attributes = c('chromosome_name', 'start_position', 'end_position',
                 'hgnc_symbol', 'gene_biotype', 'description'),
  filters = c('chromosome_name', 'start', 'end'),
  values = list(chr, region_start, region_end),
  mart = ensembl
)
```

**Option B: Use bedtools intersect**
```bash
# Download GENCODE annotations
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz

# Intersect regions with genes
bedtools intersect -a selection_regions.bed -b gencode.v44.genes.bed -wa -wb > regions_with_genes.txt
```

**Option C: Use Python with PyEnsembl**
```python
from pyensembl import EnsemblRelease
ensembl = EnsemblRelease(release=110, species='human')

for region in regions:
    genes = ensembl.genes_at_locus(
        contig=region.chr,
        position=region.start,
        end=region.end
    )
```

### Phase 4: Gene Enrichment Analysis (Week 3-4)

#### 4.1 Prepare Gene Lists
Extract genes from:
1. X chromosome selection peaks
2. Autosomal selection peaks (control)
3. Background (all genes)

#### 4.2 Test for Functional Enrichment

**Reproductive genes:**
- GO terms: reproduction, fertilization, gametogenesis
- KEGG pathways: oocyte meiosis, progesterone signaling
- HPO terms: infertility, pregnancy complications

**Immune genes:**
- GO terms: immune response, inflammation
- KEGG pathways: cytokine signaling, T/B cell receptor
- HPO terms: autoimmune disorders

**Tools to use:**
- g:Profiler (https://biit.cs.ut.ee/gprofiler/)
- DAVID (https://david.ncifcrf.gov/)
- Enrichr (https://maayanlab.cloud/Enrichr/)
- clusterProfiler (R package)

#### 4.3 Statistical Testing
- Fisher's exact test for enrichment
- Multiple testing correction (Bonferroni/FDR)
- Compare X vs autosome enrichment patterns

### Phase 5: Biological Interpretation (Week 4-5)

#### 5.1 Literature Review
For top candidate genes:
- Known functions
- Previous selection studies
- Disease associations
- Sex-specific effects

#### 5.2 X-Specific Considerations
- X-inactivation status
- Dosage compensation
- Sex-biased expression
- Male hemizygosity effects

#### 5.3 Health Implications
Link findings to:
- Pregnancy outcomes
- Fertility differences
- Sex-specific disease risk
- Immune response differences

### Phase 6: Extended Analyses (Optional, Week 5-6)

#### 6.1 XP-EHH Analysis
Cross-population extended haplotype homozygosity:
- Compare populations (e.g., African vs European)
- Identify population-specific selection
- Complement iHS findings

#### 6.2 Functional Validation
- eQTL analysis (GTEx data)
- GWAS overlap
- Expression in relevant tissues (ovary, testis, immune)

#### 6.3 Haplotype Analysis
- Reconstruct extended haplotypes
- Estimate selection coefficients
- Dating selection events

## Key Questions to Address

1. **Are there more/stronger selection signals on X vs autosomes?**
2. **Are X chromosome selection peaks enriched for reproductive genes?**
3. **Are X chromosome selection peaks enriched for immune genes?**
4. **Do selection signals differ between populations?**
5. **What are the health implications for women?**

## Tools & Resources Needed

### Software
- **R packages:** tidyverse, GenomicRanges, biomaRt, clusterProfiler
- **Python:** pandas, numpy, matplotlib, seaborn
- **Visualization:** ggplot2, manhattanly
- **Enrichment:** g:Profiler web interface or API

### Data Resources
- **Gene annotations:** Ensembl, UCSC Genome Browser
- **Functional databases:** Gene Ontology, KEGG, Reactome
- **Expression data:** GTEx, Human Protein Atlas
- **Previous studies:** HGDP selection browser, GWAS Catalog

### Reference Gene Lists
Create curated lists of:
- Reproductive genes (literature + GO terms)
- Immune genes (ImmPort, InnateDB)
- Sex-biased expression genes (GTEx)
- X-inactivation escape genes

## Timeline Summary

| Week | Phase | Deliverable |
|------|-------|-------------|
| 1-2 | Candidate identification | List of selection peaks |
| 2 | Visualization | Manhattan plots, distributions |
| 3 | Region definition | Annotated selection regions |
| 3-4 | Enrichment testing | Functional enrichment results |
| 4-5 | Interpretation | Biological interpretation document |
| 5-6 | Extended analyses | XP-EHH, validation (optional) |

## Immediate Action Items

1. **Create analysis directory structure:**
   ```bash
   mkdir -p results/analysis/{candidates,plots,enrichment,regions}
   ```

2. **Merge iHS results into genome-wide file**

3. **Create initial visualization:**
   - Histogram of Std iHS (X vs autosomes)
   - Basic Manhattan plot

4. **Define significance threshold:**
   - Calculate empirical distribution
   - Decide on conservative vs liberal approach

5. **Extract top 100-200 candidate regions** for initial exploration

## Expected Outputs

1. **Figures:**
   - Manhattan plots (genome-wide, X-specific)
   - Distribution comparisons
   - Enrichment bar plots

2. **Tables:**
   - Top selection candidates
   - Enriched functional categories
   - Gene lists with annotations

3. **Report/Manuscript:**
   - Methods section
   - Results with interpretations
   - Discussion of health implications

## Notes

- Focus on **X chromosome** as primary chromosome of interest
- Use **autosomes as controls** for enrichment comparisons
- Consider **sex-specific selection** mechanisms
- Link findings to **women's health outcomes**
- Be aware of **X chromosome-specific biases** in selection detection
