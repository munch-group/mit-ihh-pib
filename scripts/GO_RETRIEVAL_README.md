# GO Term Retrieval for Enrichment Analysis

## Overview

This directory contains scripts to retrieve GO (Gene Ontology) annotations for all X chromosome genes to enable proper GO enrichment analysis.

## Why This Is Needed

GO enrichment analysis requires GO annotations for both:
- **Test set**: 11 high-confidence genes (✓ already have)
- **Background set**: All 664 X chromosome genes (✗ need to retrieve)

Without background GO annotations, we cannot perform statistical enrichment testing.

## Files

### Python Scripts
- `retrieve_GO_terms_X_genes.py` - Retrieves GO terms for all X genes
- `GO_enrichment_analysis_fast.py` - Fast enrichment using pre-computed annotations

### Cluster Submission
- `submit_GO_retrieval.sh` - SLURM sbatch script for cluster submission

## How to Run

### Option 1: Submit to Cluster (Recommended)

```bash
cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib

# Submit the job
sbatch scripts/submit_GO_retrieval.sh

# Check job status
squeue -u vanbruggenmit

# Monitor progress (after job starts)
tail -f logs/GO_retrieval_X_*.out
```

### Option 2: Run Interactively (Not Recommended - takes 6-8 hours)

```bash
cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib
.pixi/envs/default/bin/python scripts/retrieve_GO_terms_X_genes.py
```

## Expected Runtime

- **Number of genes**: ~664 X chromosome genes
- **Estimated time**: 6-8 hours (approximately 3 seconds per gene)
- **Memory usage**: ~8 GB
- **Progress updates**: Every 10 genes

## Output Files

All files will be saved to: `results/analysis/functional_enrichment/`

1. **X_chromosome_GO_annotations_complete.tsv**
   - Main output: GO annotations for all X genes
   - Columns: Gene, GO_ID, GO_Name, GO_Category, In_X_Selection, In_All_X

2. **X_chromosome_GO_retrieval_summary.txt**
   - Summary statistics of retrieval
   - Gene counts, success rates, etc.

3. **X_chromosome_GO_retrieval_errors.tsv**
   - Log of any genes that failed (if any)

## After Completion

Once the GO term retrieval is complete, run the fast enrichment analysis:

```bash
cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib
.pixi/envs/default/bin/python scripts/GO_enrichment_analysis_fast.py
```

This will perform:
- **Analysis A**: High-confidence genes vs all X genes
- **Analysis B**: High-confidence genes vs X genes under selection
- **Analysis C**: X genes vs autosomal genes (requires additional retrieval)

## Troubleshooting

### Job not starting
```bash
# Check job queue
squeue -u vanbruggenmit

# Check job details
scontrol show job <job_id>
```

### Job failed
```bash
# Check error log
cat logs/GO_retrieval_X_*.err

# Check output log
cat logs/GO_retrieval_X_*.out
```

### Common issues
- **Network issues**: NCBI queries may timeout occasionally (script continues)
- **Memory**: Job uses ~8GB RAM for GO database loading
- **Time limit**: 12 hour limit should be sufficient for 664 genes

## Notes

- The script queries NCBI one gene at a time (required by geneinfo package)
- Progress is logged every 10 genes
- Errors are logged but don't stop the process
- The GO database (~44K terms) is loaded once at the start
