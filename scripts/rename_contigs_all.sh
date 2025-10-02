#!/usr/bin/env bash
set -euo pipefail

# =========================
# Rename VCF contigs (remove 'chr' prefix) per mapping file,
# reindex (tbi or csi), and swap into place.
# Works via `pixi run bcftools ...` from repo root.
# =========================

# --- CONFIG ---
REPO_ROOT="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib"
DATA_ROOT="/home/vanbruggenmit/mit-ihh-pib/data/grch38/raw"
MAP="/home/vanbruggenmit/mit-ihh-pib/data/grch38/contig_map_chr_to_plain.tsv"
CHROMS=( {1..22} X )

# Set to "true" to force .tbi; bcftools will pick often .csi for large refs
FORCE_TBI=false

# --- UX niceties ---
log()   { printf "[%s] %s\n" "$(date +'%F %T')" "$*" >&2; }
die()   { printf "ERROR: %s\n" "$*" >&2; exit 1; }
trap 'echo "ERROR: Script failed at line $LINENO" >&2' ERR

# --- Pre-flight ---
[[ -d "$REPO_ROOT" ]] || die "Repo root not found: $REPO_ROOT"
[[ -d "$DATA_ROOT" ]] || die "Data root not found: $DATA_ROOT"
[[ -f "$MAP"       ]] || die "Mapping file not found: $MAP"

cd "$REPO_ROOT"
log "Using repo root: $REPO_ROOT"
log "Using data root: $DATA_ROOT"
log "Using mapping:   $MAP"

process_chr() {
  local C="$1"
  local in="$DATA_ROOT/chr${C}/chr${C}.vcf.gz"
  local out="$DATA_ROOT/chr${C}/chr${C}.noprefix.vcf.gz"
  local tmp="${out}.tmp.$$"

  if [[ ! -f "$in" ]]; then
    log "SKIP chr${C}: missing $in"
    return 0
  fi

  log ">>> Renaming contigs in chr${C} ..."
  # 1) annotate to a temporary file
  pixi run bcftools annotate --rename-chrs "$MAP" -Oz -o "$tmp" "$in"

  # 2) index the temp output
  if $FORCE_TBI; then
    pixi run bcftools index -f --tbi "$tmp"
  else
    pixi run bcftools index -f "$tmp"
  fi

  # 3) figure out which index was created
  local idx_ext=""
  if   [[ -f "${tmp}.tbi" ]]; then idx_ext="tbi"
  elif [[ -f "${tmp}.csi" ]]; then idx_ext="csi"
  else
    rm -f "$tmp"
    die "No index was created for $tmp"
  fi

  # 4) keep a backup of the original and its index (if present)
  local backup="${in/.vcf.gz/.withchr.vcf.gz}"
  mv "$in" "$backup"
  [[ -f "${in}.tbi" ]] && mv "${in}.tbi" "${in}.tbi.bak" || true
  [[ -f "${in}.csi" ]] && mv "${in}.csi" "${in}.csi.bak" || true

  # 5) move temp into final names atomically
  mv "$tmp" "$in"
  mv "${tmp}.${idx_ext}" "${in}.${idx_ext}"

  # 6) quick sanity: show a contig header line
  pixi run bcftools view -h "$in" | grep -m1 '^##contig=' || true

  log "Done chr${C}"
}

for C in "${CHROMS[@]}"; do
  process_chr "$C"
done

log "All done MITTTT."
