#!/usr/bin/env bash
# make_hg38_keep_lists.sh — robust repo-root detection for pixi
# Usage:
#   ./people/vanbruggenmit/mit-ihh-pib/scripts/make_hg38_keep_lists.sh \
#     --vcf data/grch38/raw/chr22/chr22.vcf.gz \
#     --panel data/grch38/panels/integrated_call_samples_v3.20130502.ALL.panel \
#     --outdir data/grch38/panels \
#     --pops "GBR YRI" \
#     [--repo-root /home/vanbruggenmit/mit-ihh-pib]

set -euo pipefail

REPO_ROOT=""
VCF=""
PANEL=""
OUTDIR=""
POPS="GBR YRI"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf)       VCF="$2"; shift 2 ;;
    --panel)     PANEL="$2"; shift 2 ;;
    --outdir)    OUTDIR="$2"; shift 2 ;;
    --pops)      POPS="$2"; shift 2 ;;
    --repo-root) REPO_ROOT="$2"; shift 2 ;;
    -h|--help)
      grep '^#' "$0" | sed 's/^# \{0,1\}//'
      exit 0
      ;;
    *) echo "Unknown arg: $1" >&2; exit 1 ;;
  endesac
done

# --- find repo root that contains pyproject.toml or pixi.toml ---
find_root() {
  local start="$1"
  while [[ "$start" != "/" && -n "$start" ]]; do
    if [[ -f "$start/pyproject.toml" || -f "$start/pixi.toml" ]]; then
      echo "$start"; return 0
    fi
    start="$(dirname "$start")"
  done
  return 1
}

if [[ -z "$REPO_ROOT" ]]; then
  SCRIPT_DIR="$(cd -- "$(dirname "$0")" && pwd -P)"
  CANDIDATES=("$PWD" "$SCRIPT_DIR" "$HOME/mit-ihh-pib")
  for c in "${CANDIDATES[@]}"; do
    if ROOT="$(find_root "$c")"; then REPO_ROOT="$ROOT"; break; fi
  done
fi

if [[ -z "$REPO_ROOT" ]]; then
  echo "Error: couldn’t locate repo root (no pyproject.toml/pixi.toml found). Use --repo-root <path>." >&2
  exit 1
fi

cd "$REPO_ROOT"
echo "[repo] Using pixi project at: $REPO_ROOT"

[[ -n "$VCF" && -n "$PANEL" && -n "$OUTDIR" ]] || {
  echo "Error: --vcf, --panel, and --outdir are required." >&2
  exit 1
}

mkdir -p "$OUTDIR"

echo "[1/3] Extracting sample IDs from VCF → vcf.hg38.samples"
pixi run bcftools query -l "$VCF" | sort > "${OUTDIR}/vcf.hg38.samples"

echo "[2/3] Building per-population sample lists from panel"
for POP in $POPS; do
  awk -v P="$POP" '$2==P{print $1}' "$PANEL" | sort > "${OUTDIR}/${POP}.panel.samples"
done

echo "[3/3] Intersecting panel lists with VCF samples"
for POP in $POPS; do
  comm -12 "${OUTDIR}/vcf.hg38.samples" "${OUTDIR}/${POP}.panel.samples" \
    > "${OUTDIR}/${POP}.hg38.keep"
  echo "  ${POP}: $(wc -l < "${OUTDIR}/${POP}.hg38.keep") samples"
done

echo "Done. Keep-lists written to: ${OUTDIR} (files: *.hg38.keep)"
