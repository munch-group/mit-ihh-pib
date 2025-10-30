# GO (Gene Ontology) Term Analysis

**Script:** `scripts/get_GO_annotations.py`
**Status:** Running (first run downloads GO database - may take 5-10 minutes)

## Purpose

Retrieve and analyze Gene Ontology (GO) terms for the 11 high-confidence overlap genes to understand their biological functions.

## Overlap Genes Analyzed

1. DACH2
2. DIAPH2
3. FAAH2
4. GPC3
5. HTR2C
6. IL1RAPL2
7. KDM6A
8. NAP1L2
9. NCBP2L
10. TRPC5
11. ZMAT1

## Analysis Steps

### 1. Retrieve GO Terms
For each gene, retrieve all associated GO terms using `geneinfo.ontology.get_go_terms_for_genes()`

### 2. Categorize by Function
Classify genes based on GO term keywords:

**Reproduction keywords:**
- meiosis, gamete, fertility, gonad, sperm, oocyte
- spermatogenesis, oogenesis, germ cell, reproductive
- fertilization, ovary, testis

**Immunity keywords:**
- immune, inflammation, defense, cytokine, antibody
- inflammatory, innate immune, adaptive immune, lymphocyte
- T cell, B cell, interferon, interleukin, antigen

**Neurodevelopment keywords:**
- neuron, synapse, cognition, brain, neural, axon
- dendrite, synaptic, neuronal, nervous system, learning
- memory, behavior, neurogenesis, neurotransmitter

**Development keywords:**
- development, differentiation, morphogenesis, organogenesis
- embryo, pattern, specification, determination

### 3. Generate Summary Statistics

**Questions answered:**
- How many genes have GO annotations?
- What fraction fall into each functional category?
- Which genes have the most GO terms?
- Are immune/reproductive functions enriched?

## Output Files

### 1. gene_GO_annotations.tsv
Complete table of all GO terms for all genes

**Columns:**
- `Gene`: Gene symbol
- `GO_ID`: GO term identifier (e.g., GO:0000776)
- `GO_Name`: GO term name
- `GO_Category`: biological_process, molecular_function, or cellular_component
- `Functional_Categories`: Matched functional categories (reproduction, immunity, etc.)

### 2. gene_functional_summary.tsv
Summary statistics per gene

**Columns:**
- `Gene`: Gene symbol
- `Total_GO_Terms`: Total number of GO terms
- `Biological_Process`: Number of BP terms
- `Molecular_Function`: Number of MF terms
- `Cellular_Component`: Number of CC terms
- `Reproduction_Terms`: Terms matching reproduction keywords
- `Immunity_Terms`: Terms matching immunity keywords
- `Neurodevelopment_Terms`: Terms matching neurodevelopment keywords
- `Development_Terms`: Terms matching development keywords
- `Primary_Function`: Most represented functional category

## Expected Insights

1. **Functional diversity**: Do the genes have diverse or similar functions?
2. **Immune enrichment**: Are immune-related genes over-represented?
3. **Reproductive function**: Any genes involved in reproduction/fertility?
4. **Neural function**: Connection to neurodevelopment (X-linked intellectual disability)?
5. **Annotation quality**: What fraction of genes have good GO coverage?

## Usage

```bash
# Run the analysis (first time will download GO database)
python scripts/get_GO_annotations.py

# View results
column -t -s $'\t' results/analysis/functional_enrichment/gene_GO_annotations.tsv | less -S
column -t -s $'\t' results/analysis/functional_enrichment/gene_functional_summary.tsv
```

## Next Steps After Completion

1. Review functional categories - do they support selection hypotheses?
2. Perform formal GO enrichment analysis (if patterns emerge)
3. Compare with background X chromosome genes
4. Integrate with literature review
5. Generate GO term visualizations using `geneinfo.ontology.show_go_dag_for_gene()`

## Notes

- First run downloads GO database (~44,000 terms) - this is normal and only happens once
- GO terms are standardized across all species
- Multiple GO terms per gene is expected (genes have multiple functions)
- Some genes may have limited GO annotations (especially newer/less-studied genes)
