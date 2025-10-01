#run with: pixi run python scripts/fetch_from_yaml.py chr22
#Minimal Python helper to read the YAML and download a test chromosome (chr22) first, validate pipeline
import sys, os, subprocess, yaml, pathlib

yml = yaml.safe_load(open("config/grch38_sources.yml"))
chrom = sys.argv[1] if len(sys.argv) > 1 else "chr22"

outdir = pathlib.Path("/home/vanbruggenmit/mit-ihh-pib/data/grch38/raw")
(outdir / chrom).mkdir(parents=True, exist_ok=True)

def grab(url, dest):
    if not dest.exists():
        subprocess.run(["wget", "-q", "-O", str(dest), url], check=True)

# sample VCF + index
vcf_url = yml["sample_vcf"][chrom]
tbi_url = yml["sample_vcf_index"][chrom]
grab(vcf_url, outdir/chrom/f"{chrom}.vcf.gz")
grab(tbi_url, outdir/chrom/f"{chrom}.vcf.gz.tbi")

# mask (per-chr FASTA.gz)
mask_url = yml["mask"][chrom]
grab(mask_url, outdir/chrom/f"{chrom}.mask.fasta.gz")

# ancestral set is a tarball; fetch once at top-level
anc_url = yml["ancestral_fa"]
grab(anc_url, outdir/"ancestral_GRCh38.tar.gz")

print("Downloaded:", vcf_url, tbi_url, mask_url, anc_url)
