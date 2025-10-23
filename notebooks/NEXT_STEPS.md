# Next Steps: X Chromosome Selection Analysis

## Current Progress ✓

You have successfully calculated iHS scores for all chromosomes (chr1-22 and chrX) with proper variant IDs and validated data quality.

## Immediate Next Steps

### Step 1: Run Initial Candidate Identification (TODAY)

Run the analysis script I just created:

```bash
cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib
Rscript scripts/01_identify_ihs_candidates.R
```

**This will:**
- Load all chromosome iHS data (~11 million variants total)
- Calculate summary statistics
- Identify selection candidates at 3 thresholds
- Test if X chromosome is enriched for selection signals
- Create initial visualization plots

**Expected runtime:** 5-10 minutes

**Outputs:**
- List of candidate variants at different thresholds
- Statistical test for X enrichment
- Distribution plots (X vs autosomes)

### Step 2: Review Initial Results (TODAY/TOMORROW)

**Key questions to answer:**
1. How many selection candidates are on the X chromosome?
2. Is the X chromosome enriched compared to autosomes?
3. What threshold should you use (stringent/moderate/liberal)?
4. Do you see any obvious patterns in the plots?

### Step 3: Define Selection Regions (THIS WEEK)

Cluster nearby significant SNPs into contiguous regions:

**Tasks:**
- Merge SNPs within 100kb windows
- Count how many genes fall in each region
- Create a table of "selection peaks"

**I can help you create the next script for this**

### Step 4: Gene Annotation (THIS WEEK)

Map selection regions to genes:

**Tasks:**
- Download gene annotations (Ensembl/UCSC)
- Overlap selection peaks with gene coordinates
- Extract gene lists for X chromosome peaks

### Step 5: Enrichment Analysis (NEXT WEEK)

Test functional enrichment:

**For X chromosome genes in selection peaks, test enrichment for:**

**Reproductive functions:**
- GO:0000003 - reproduction
- GO:0007283 - spermatogenesis
- GO:0048477 - oogenesis
- GO:0007338 - single fertilization
- GO:0022414 - reproductive process

**Immune functions:**
- GO:0006955 - immune response
- GO:0002376 - immune system process
- GO:0006954 - inflammatory response
- GO:0045087 - innate immune response
- GO:0002250 - adaptive immune response

**Tools:**
- g:Profiler web interface (easiest)
- DAVID
- clusterProfiler R package

### Step 6: Literature Review & Interpretation (NEXT WEEK)

For top candidate genes:
- PubMed search for each gene
- Check GTEx for tissue-specific expression
- Look for GWAS associations
- Review previous selection studies

## Key Research Questions

Your project should address:

1. **Primary question:** Are X chromosome regions under selection enriched for reproductive and immune genes?

2. **Secondary questions:**
   - How does selection on X compare to autosomes?
   - Are there population-specific patterns?
   - Do selected genes show sex-biased expression?
   - What are the health implications for women?

## Timeline

| Week | Tasks | Deliverables |
|------|-------|-------------|
| **This week** | Run scripts 1-3, identify candidates | Candidate lists, initial plots |
| **Next week** | Gene annotation, enrichment analysis | Gene lists, enrichment results |
| **Week 3** | Literature review, interpretation | Summary of top genes |
| **Week 4** | Extended analysis (XP-EHH optional) | Additional evidence |
| **Week 5-6** | Write-up, final figures | Draft report/manuscript |

## Resources You'll Need

### Gene Databases
- **Ensembl:** https://www.ensembl.org/biomart/martview
- **UCSC Table Browser:** https://genome.ucsc.edu/cgi-bin/hgTables

### Enrichment Tools
- **g:Profiler:** https://biit.cs.ut.ee/gprofiler/gost
- **Enrichr:** https://maayanlab.cloud/Enrichr/
- **DAVID:** https://david.ncifcrf.gov/

### Expression Data
- **GTEx Portal:** https://gtexportal.org/
- **Human Protein Atlas:** https://www.proteinatlas.org/

### Previous Studies
- **HGDP Selection Browser:** http://hgdp.uchicago.edu/
- Search: "X chromosome selection humans"
- Search: "sex-biased genes immune reproduction"

## What I Can Help With Next

Let me know when you're ready for:

1. **Creating script 02** - Define selection regions/peaks
2. **Creating script 03** - Annotate regions with genes
3. **Help with enrichment analysis** - Prepare gene lists
4. **Create Manhattan plot** - Genome-wide visualization
5. **Statistical comparisons** - X vs autosomes

## Quick Start Commands

```bash
# Navigate to project directory
cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib

# Run initial analysis
Rscript scripts/01_identify_ihs_candidates.R

# Check outputs
ls -lh results/analysis/candidates/
ls -lh results/analysis/plots/

# View candidate counts
head results/analysis/candidates/ihs_candidates_moderate.tsv

# View enrichment test
cat results/analysis/candidates/x_enrichment_test.tsv
```

## Expected Key Findings

Based on the literature, you might expect to find:

1. **X chromosome shows stronger selection** than autosomes (due to hemizygosity in males)
2. **Enrichment for immune genes** (especially if using recent selection timeframes)
3. **Some reproductive genes** under selection
4. **Population differences** in selection signals
5. **Sex-biased expression** in selected genes

## Notes

- The X chromosome has unique evolutionary dynamics
- Consider X-inactivation in your interpretation
- Be aware of demographic effects (smaller Ne)
- Think about male vs female selection pressures
- Link findings to pregnancy/fertility outcomes

Good luck! Let me know when you're ready to proceed with the next steps.
