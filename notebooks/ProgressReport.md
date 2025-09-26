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
