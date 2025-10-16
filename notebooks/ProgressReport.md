# 120925 Meeting

Goal: detect regions of recent of positive selection on the human X chromosome using population-genomic data, and test whether they’re enriched for reproduction and immune genes, to see what that might mean for sex-specific health.   
Enriched: do genes in my set show up in certain categories more often than I’d expect by chance?  
**Step 1: You’ll have a list of genes that fall inside regions of the X chromosome with strong selection signals.**

**Step 2: You test whether those genes are enriched for particular functional categories, e.g.:**

* **Reproduction-related genes,**

* **Immune-related genes,**

* **Pathways from databases like GO (Gene Ontology), KEGG, or Reactome.**

**Step 3: If reproductive or immune genes are significantly more common in your list than in the whole genome → that’s enrichment.**

iHS: integrated haplotype score  
XP EHH: cross population extended haplotype homozygosity  
How to apply? Which tools? Which datasets? 

1000 Genomes project suitable? or other datasets? ancient DNA or Simons Genome Diversity Project

- publicly available   
- already phased   
- has global populations to compare

Paper 2024 Nam, Munch, Dutheil: Strong selective sweeps in ape X chromosomes: potential in focusing on reproduction and immunity 

Better to focus on one or two populations? Or compare across many?   
Are there specific genes/pathways on the X that would be especially interesting to test for enrichment?   
What are biggest pitfalls when interpreting sweep signals on X? 

abstract \+ intro \+ context   
iHS XP-EHH (multiple populations, cross) 

- african pop, asian, central america/PR 

What kind of input from 1000 genomes, what format does it need ← look into this   
Work on the cluster ,-- get account   
trees along genome,   
github thing \- Mit IHH PiB

Next week: git dingen, project installation, … 

**Bron: "Evolutionary\_and\_biomedical\_insights\_from\_a\_marmos.pdf"** Deze studie richtte zich op het genereren van een diploïde genoomassemblage voor de marmoset, met gebruikmaking van verschillende geavanceerde sequencingtechnologieën. Het onderzoek beoogde genetische varianten te identificeren, pseudoautosomale regio's op te lossen en nieuwe geslachts-differentiatiegebieden te ontdekken, terwijl ook potentiële pathogene humane plaatsen in het marmosetgenoom werden opgemerkt.  
• **Bron: "Fueling-the-microbiome\_2.pdf"** Deze thesis had als doel nieuwe putative pullulanase-coderende genen uit *Lactobacillus crispatus*\-stammen te identificeren en hun enzymatische activiteit functioneel te onderzoeken. Hoewel voorlopige resultaten een zekere mate van glycogeendegradatie suggereerden, kon definitieve enzymatische activiteit onder de geteste omstandigheden niet statistisch significant worden bevestigd, wat vragen oproept over stam-specificiteit of assaygevoeligheid.  
• **Bron: "Munch.pdf"** Deze studie onderzocht de divergentiepatronen van het X-chromosoom tussen mensen en chimpansees, met een focus op regio's met lage onvolledige lijn-splitsing (ILS). De conclusie was dat sterke selectieve *sweeps* in de mens-chimpansee voorouder, in plaats van achtergrondselectie of lagere mutatiesnelheden, de waargenomen lage divergentie en diversiteitswoestijnen op het X-chromosoom verklaren.  
• **Bron: "fgene-12-714491.pdf"** Deze studie voerde een uitgebreide analyse uit van recente positieve selectie op het X-chromosoom in diverse menselijke populaties, waarbij kenmerken van zowel 'harde' als 'zachte' selectieve *sweeps* in coderende en niet-coderende gebieden werden geïdentificeerd. Het concludeerde dat Sub-Sahara Afrikaanse populaties meer selectieve *sweeps* vertonen, er een verrijking is van neuraal- en vruchtbaarheidsgerelateerde genen, en dat regulerende veranderingen (zoals in enhancers) waarschijnlijk een belangrijke rol spelen in recente menselijke adaptatie.  
• **Bron: "nihms-111210.pdf"** Deze review stelt dat seks-specifieke genetische architectuur en genotype-geslacht-interacties wijdverspreid zijn bij mensen en de prevalentie, het verloop en de ernst van veelvoorkomende ziekten en kwantitatieve fenotypes beïnvloeden. Het benadrukt dat differentiële genregulatie tussen mannen en vrouwen, met name in sekssteroïd-responsieve genen en weefselspecifieke eQTLs, waarschijnlijk ten grondslag ligt aan de meeste fenotypische seksuele dimorfie en cruciaal is om te overwegen in genetische studies.  
• **Bron: "nihms905235.pdf"** Deze studie voerde een systematisch onderzoek uit naar het landschap van X-chromosoominactivatie (XCI) in verschillende menselijke weefsels, met behulp van drie complementaire RNA-sequencing benaderingen. Het toonde aan dat XCI onvolledig en variabel is over genen en individuen heen, waarbij 'ontsnapping' aan XCI vaak leidt tot seks-bevooroordeelde genexpressie die mogelijk bijdraagt aan seksverschillen in gezondheid en ziekte.  
• **Bron: "pbio.0040072.pdf"** Deze paper presenteert een genoomwijde kaart van onvolledige selectieve *sweeps* bij mensen, waarbij wijdverspreide signalen van recente positieve selectie in drie continentale populaties (Oost-Aziatisch, Europees en Yoruba) werden geïdentificeerd. Het concludeerde dat veel selectiesignalen regio-specifiek zijn, recent hebben plaatsgevonden (meestal in het Holoceen-tijdperk) en suggereren dat deze loci bijdragen aan aanzienlijke fenotypische variatie die relevant is voor complexe eigenschappen en lokale aanpassingen.

# 190925 Meeting → 220925

* Add me on github itself? → I can’t push documents yet  
* Slurm-jupyter   
  * popgen environment \= birc-project environment  
  * slurm-jupyter \-u vanbruggenmit \-A mit-ihh-pib \-e birc-project \--chrome  
    * birc-project environment wasn’t created yet   
    * So I did this and I installed jupyterlab on the cluster, because it wasn’t there yet.   
    * I also had to accept the conda terms of service for this  
  * I got it to work eventually so CHECK  
* Visual Studio code  
  * Done  
* Scheduling jobs  
  * mj command not working \- is it because job was really short?  
  * Job worked, succes.txt file created in scripts folder  
  * haven’t tried scancel yet  
  * cat success.txt  
  * firstjob.stdout and firstjob.stderr not working?  
* Copying files from and to cluster  
* Reproducible research   
  * ![][image1]  
* Quarto reporting   
  * executed all the commands  
  * R installation? (when checking quarto installation:  Unable to locate an installed version of R. Install R from [https://cloud.r-project.org/](https://cloud.r-project.org/))   
    * How to install R on cluster? Necessary?  
  * quarto render \--execute  
    * not working  
      * rm \-rf .quarto/quarto-session-temp\*  
      * conda install \-y numpy pandas matplotlib seaborn ipython  
      * cat global\_params.py  
        * error because: params \= load\_params("global\_params.yml") params.questions  
        * in global\_params.yml, there is no entry for questions. So questions never becomes an attribute of params.  
    * had to change some things in global\_params.yml → there were no questions in there, and those were needed.   
    * .html file created, opened by copying it to my laptop and opening from there   
      * easier way?   
      * scp vanbruggenmit@login.genome.au.dk:/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/\_book/thesis/index.html \~/Downloads/  
  * quarto publish gh-pages  
    * I can’t push things to github yet – no access. 

### Format 1000 Genomes → IHS format needed

* Input data from 1000 Genomes  
  * IHS requires phased haplotype data  
    * statistic is based on extended haplotype homozygosity  
  * From 1000 Genomes, you would typically use VCF files containing phased SNP genotypes  
    * Ideally, population-specific subsets (EUR, AFR, EAS)  
  * Also: need recombination map data for the same build of the genome  
    * HapMap/1000 Genomes recombination maps   
* Needed iHS format  
  * depends on software implementation   
    * rehh in R   
      * https://cran.r-project.org/web/packages/rehh/index.html  
    * selscan  
      * https://hpc.nih.gov/apps/selscan.html  
    * hapbin   
      * [https://pmc.ncbi.nlm.nih.gov/articles/PMC4651233/](https://pmc.ncbi.nlm.nih.gov/articles/PMC4651233/)  
      * much faster than the ones above  
  * In general:   
    * selscan: takes phased haplotypes in VCF or its own hap/legend/sample format   
    * rehh expects haplotype and SNP position files  
      * first convert VCF → haplotype matrix \+ map   
  * All require SNP positions (chromosome, bp coordinates), phased haplotypes for each individual, recombination map to convert physical distance into genetic distance. 

### 

### 

### Hapbin

* C++ toolkit for fast calculation of haplotype-based statistics (iHS, XP-EHH), designed for large-scale datasets like 1000 genomes → much faster  
  * [https://github.com/evotools/hapbin](https://github.com/evotools/hapbin)  
* Input requirements  
  * does not work directly on VCFs → needs converting to hapbin BINARY format   
1. Phased haplotypes: 1000 Genomes provides phased VCFs but Hapbin needs these in its own binary format .hap and .map files   
2. Genetic map: required for iHS to scale haplotype decay in terms of recombination distance, not just base pairs   
- 1000 genomes provides recombination maps in HapMap style .txt files: chrom, pos, rate   
- Hapbin has a tool convert\_map to turn these into its binary map format  
3. Population subset  
- Hapbin works per-population: need to filer 1000 Genomes individuals by population group   
- Usually done by making sample list and filtering VCFs before conversion   
* Workflow  
  * Download 1000 Genomes phased VCFs for populations of interest   
    * phased VCF   
    * make sure we’ve subset the individuals to population you want   
      * bcftools view \-S sample\_list.txt input.vcf.gz \-Oz \-o subset.vcf.gz  
  * Download recombination map for the same genome build   
  * Convert VCF → hapbin format   
    * convert\_vcf \--vcf input.vcf.gz \--map recomb.map \--hap output.hap \--legend output.legend \--sample output.sample  
    * Hapbin comes with a tool called vcf2hapmap or convert\_vcf  
    * check outputs: .hap, .map, .legend , .sample  
  * Run Hapbin iHS   
    * ihs \--hap output.hap \--map output.map \--out ihs\_results  
  * Normalize iHS: hapbin has a separate normalization utility, since raw iHS values depend on allele frequency → needs to happen before interpretation 

# 260925 Meeting

* Project structure setup: data in home/vanbruggenmit/mit-ihh-pib/data and scripts is home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/scripts
* Installed necessary tools: bcftools, htslib/tabix, plink, vcftools
    * bcftools to work directly with VCF files to view and subset by chromosome, region or samples
    * Tabix: provides indexing and fast random access to compressed VCF/BCF files (otherwise would have to decompress entire chromosomes at each step)
    * vcftools: manipulating vcfs beyond what bcftools does, to convert into plink format (.ped/.map) and from there to hapbin's required .hap/.map
* 1000 Genomes Data Preparation
    * Downloaded per-chromosome VCFs (GRCh37/hg19, Phase 3) (download_phase3.sh script, puts it in /home/vanbruggenmit/mit-ihh-pib/data/vcfs)
    * Verlified that the data is already phased (checked for 0|1 and 1|0 in chr 22)
    ```bash 
    pixi run bcftools view /home/vanbruggenmit/mit-ihh-pib/data/vcfs/ALL.chr22*.vcf.gz \
  | head -n 200 | grep -m1 '|'
    ```
    * Downloaded official 1000 Genomes sample --> population list (make_sample_lists.sh --> panels)
    * Generated custom sample lists for GRB and YRI populations (subset_vcfs.sh in subsets)
    * Subsetted VCFs by chromosome and then by population (idem)

    So: 
    
| Script               | Purpose (what it does)                                                                                          | Key inputs/params                                                                                                                                         | Outputs (what files)                                                                                                                                             | Where the files are now                                                                |
|----------------------|------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------|
| `download_phase3.sh` | Downloads the official 1000 Genomes Phase 3 **per-chromosome VCFs** (+ their tabix indexes). Creates the `vcfs/` folder if missing. | `DATA_ROOT` (default: `/home/vanbruggenmit/mit-ihh-pib/data`), `CHRS` (currently set to `22` for testing), internal `prefix`/`suffix` pointing to FTP path | `ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` and `.tbi` (for the test run on chr22)                                             | `${DATA_ROOT}/vcfs/` → `/home/vanbruggenmit/mit-ihh-pib/data/vcfs/`                    |
| `make_sample_lists.sh` | Downloads the **sample→population panel** file and builds **per-population sample lists**. Creates the `panels/` folder.           | `DATA_ROOT` (default as above). The script currently emits lists for `GBR` and `YRI` via the two `awk` lines.                                             | `integrated_call_samples_v3.20130502.ALL.panel`; `GBR.samples.txt`; `YRI.samples.txt` (one sample ID per line).                                                  | `${DATA_ROOT}/panels/` → `/home/vanbruggenmit/mit-ihh-pib/data/panels/`                |
| `subset_vcfs.sh`    | Creates **population-specific VCFs** (and tabix indexes) by subsetting the big VCF to only the samples in each population list. Creates the `subsets/` folder. | `DATA_ROOT` (default), `CHRS` (currently `22`), `POPS` (default `"GBR YRI"`). Uses `pixi run bcftools view -S … -Oz` and `pixi run tabix -p vcf …`.       | For each POP × chromosome: `<POP>.chr22.snps.vcf.gz` and `<POP>.chr22.snps.vcf.gz.tbi`.                                                                         | `${DATA_ROOT}/subsets/` → `/home/vanbruggenmit/mit-ihh-pib/data/subsets/` (by POP name) |

### Scripts 
```bash
#!/usr/bin/env bash
set -euo pipefail
DATA_ROOT=${1:-/home/vanbruggenmit/mit-ihh-pib/data}

prefix="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr"
suffix=".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

mkdir -p "${DATA_ROOT}/vcfs"

# change CHRS from "22" to "{1..22}" after testing
CHRS=${CHRS:-22}

for chr in ${CHRS}; do
  wget -c -P "${DATA_ROOT}/vcfs" "${prefix}${chr}${suffix}"
  wget -c -P "${DATA_ROOT}/vcfs" "${prefix}${chr}${suffix}.tbi"
done
```

### Pixi environment confusion 
* Encountered permission errors when trying to run/install in some cluster directories (\faststorage ...) 
* Now working with detached environments --> how to fix this? Sometimes still running into 'access denied' where can I run pixi? WHere can't i? 

### What I have done so far

I set up the project structure with separate `data/` and `scripts/` folders to keep things  organized.  
I wrote three bash scripts to automate the preparation of 1000 Genomes Phase 3 data:  

1. **Download** per-chromosome VCFs + tabix index files from the official FTP site.  
2. **Create sample lists** for specific populations (GBR, YRI) from the 1000 Genomes panel file.  
3. **Subset the VCFs** by chromosome and by population into smaller, SNP-only files.  

I started with **chromosome 22** because it’s the smallest autosome, so I could quickly test whether the scripts, tools, and indexing worked before scaling to all chromosomes.  

So far, I now have population-specific, phased, SNP-only VCFs for chr22 stored in `data/subsets/`, which are the necessary inputs for the next step: converting into hapbin’s `.hap/.legend/.sample` format for iHS analysis.

pyyaml
import yaml 
Do everything on 1 chr first (22) 

# 031025 Meeting
## What I’ve done so far

* Set up the working environment inside my repo so that pixi run uses the local pyproject.toml configuration. This way all dependencies (bcftools, samtools, pysam, etc.) are available in a reproducible way.

* Collected the required GRCh38 inputs for chromosome 22: the imputed haplotype file (.hap.gz), the corresponding legend file (.legend.gz), the human ancestor FASTA sequence for chr22, and the Decode sex-averaged recombination map (.tsv). I’m using the Decode map because it is widely applied in selection scan studies: it is high-resolution, based on a large Icelandic pedigree dataset, and provides sex-averaged genetic distances in centiMorgans, which are the format hapbin requires.

* Resolved a mismatch between the --chr argument and the FASTA header. The ancestor FASTA doesn’t use plain 22, but a long contig identifier (ANCESTOR_for_chromosome:GRCh38:22:1:50818468:1). Adjusted the script so pysam can locate the correct sequence.

* Indexed the ancestor FASTA (.fai) so pysam can quickly fetch single bases by genomic position. (Without the index, every fetch call would rescan the entire file, making the run impractically slow.)

* Successfully ran the impute2hapbin_chr22.py script (which chatgpt helped me write) with --verbose. This script reads variants from the impute2 files, looks up the ancestral allele from the FASTA, places the variant on the recombination map, and outputs hapbin-formatted .hap and .map files.

  * The run processed ~1.07M variants.

  * Of these, ~870k were kept and ~200k were skipped (due to missing ancestral alleles, non-biallelic sites, or positions not covered by the recombination map). This ratio is in line with expectations.

Generated hapbin output files (chr22.hapbin.hap and chr22.hapbin.map) that can be used as direct input to hapbin (e.g. for iHS).

## Next steps

* QC of the chr22 output

  Confirm the .hap and .map files agree in the number of variants.

  Check that physical positions are strictly increasing and that recombination map positions (cM) are non-decreasing.

* Test hapbin downstream

  Run iHS on chr22 using the generated .hap and .map to ensure hapbin parses the files correctly and produces reasonable outputs.

* Generalization of the pipeline

  Modify the conversion script so it works for all chromosomes (1–22). Ideally, the script should automatically detect the correct contig name in each ancestor FASTA file rather than requiring it to be passed manually.

* Scaling up to the whole genome

  Write a loop or job array to process all chromosomes in parallel on the cluster.

  Store hapbin outputs in a consistent directory structure for downstream analysis.

* Genome-wide selection scans

  Once hapbin inputs are available for all chromosomes, run iHS per chromosome, then standardize the results genome-wide.

  For population comparisons, generate XP-EHH scores using subsets of the samples.

  Perform QC on the distributions (e.g. proportion of |iHS| > 2, correlation with recombination rate).

# 161025: Progress report
## Summary
Converted all chromosome tsvs to .hap and .map format/ 
Completed data processing and iHS calculations for all chromosomes (1-22 + X). Discovered that chromosome X data only contains pseudoautosomal regions (PAR), which limits the project scope for studying X-linked selection.

## By now
1. Data preparation 
- VCF files: Downloaded and verified for all chromosomes from 1000 Genomes (3,202 individuals)
- Genetic maps: Created recombination maps for all chromosomes
  -  Applied 2/3 adjustment for chrX (accounts for female-only recombination)
- Phased haplotypes: Converted VCFs to .hap and .map format for hapbin

2. iHS Calculations
- Completed: iHS scan for all 23 chromosomes (chr1-22 + chrX)
- Sample size: 6,404 haplotypes (3,202 individuals × 2)
- Output: TSV files with iHS scores for ~2.86M variants per chromosome
- Status: All jobs completed successfully (except X - canceled)

3. QC and discovery of PAR issue 
Found: while monitoring the iHS hob logs, i noticed every chromosome reported 
``` Chromosomes per SNP: 6404 ```

This was **unexpected for chrX** because:
- Males have only 1 X chromosome (hemizygous)
- Expected chrX to show ~4,800 haplotypes (females: 2 × ~1,600 = 3,200; males: 1 × ~1,600 = 1,600)
- Instead, chrX showed 6,404 - same as autosomes

#### Investigation revealed:
1. **All 3,202 samples present on chrX** (should be ~3,200 if males excluded)
2. **Males show heterozygous genotypes** on chrX (impossible outside PAR regions)
3. **All 2.86M chrX variants are in PAR regions**:
   - PAR1 (positions 10,001-2,781,479): 1,368,031 variants
   - Non-PAR (positions 2,781,480-155,701,382): **0 variants**
   - PAR2 (positions 155,701,383+): 1,490,153 variants

#### Root cause:
The 1000 Genomes file specified in `config.yaml` only contains phased PAR regions:

CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz 

### SHow variant distribution across X chrom

cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib

# Count variants in each region
pixi run bcftools query -f '%POS\n' \
  /home/vanbruggenmit/mit-ihh-pib/data/grch38/raw/chrX/chrX.vcf.gz | \
  awk '{
    if ($1 <= 2781479) par1++;
    else if ($1 >= 155701383) par2++;
    else nonpar++;
  }
  END {
    print "PAR1 (0-2,781,479): " par1;
    print "Non-PAR (2,781,480-155,701,382): " nonpar;
    print "PAR2 (155,701,383+): " par2;
  }'
```

**Expected output:**
```
PAR1 (0-2,781,479): 1368031
Non-PAR (2,781,480-155,701,382): 0
PAR2 (155,701,383+): 1490153

### SHow male has heterozygous genotypes (only possible in PAR)
# HG00096 is male (confirmed from 1000G metadata)
# Check for heterozygous sites in PAR1
pixi run bcftools query -f '%CHROM\t%POS\t[%GT]\n' \
  -s HG00096 \
  -r X:10000-2781479 \
  /home/vanbruggenmit/mit-ihh-pib/data/grch38/raw/chrX/chrX.vcf.gz | \
  grep -E "0\|1|1\|0" | head -10

### Verify non-PAR region is empty 
# Try to find any variants in non-PAR region
pixi run bcftools view -H -r X:2781480-155701382 \
  /home/vanbruggenmit/mit-ihh-pib/data/grch38/raw/chrX/chrX.vcf.gz | wc -l

# Impact on project goals 
PAR regions behave like autosomes (recombine in both sexes)
Most X-linked genes are in non-PAR regions (currently missing)
Cannot study X-specific selection patterns with current data

