#!/usr/bin/env bash
# Run iHS analysis on a single population
#
# Usage:
#     ./run_single_population_ihs.sh <POP_CODE>
#
# Example:
#     ./run_single_population_ihs.sh YRI
#     ./run_single_population_ihs.sh CEU

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <POP_CODE>"
    echo ""
    echo "Example:"
    echo "  $0 YRI"
    echo "  $0 CEU"
    echo ""
    echo "Available populations:"
    echo "  ACB ASW BEB CDX CEU CHB CHD CHS CLM ESN"
    echo "  FIN GBR GIH GWD IBS ITU JPT KHV LWK MSL"
    echo "  MXL PEL PJL PUR STU TSI YRI"
    exit 1
fi

POP="$1"
POP="${POP^^}"  # Convert to uppercase

REPO_ROOT="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib"
cd "$REPO_ROOT"

# Setup logging
LOG_DIR="$REPO_ROOT/logs"
mkdir -p "$LOG_DIR"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
SINGLE_LOG="$LOG_DIR/population_ihs_${POP}_${TIMESTAMP}.log"
ERROR_LOG="$LOG_DIR/population_ihs_${POP}_errors_${TIMESTAMP}.log"

# Function to log messages
log() {
    echo "$1" | tee -a "$SINGLE_LOG"
}

log_error() {
    echo "ERROR: $1" | tee -a "$SINGLE_LOG" >> "$ERROR_LOG"
}

log "=========================================="
log "iHS Analysis for Population: $POP"
log "Started: $(date)"
log "Log file: $SINGLE_LOG"
log "=========================================="
echo ""

# Check if already done
RESULT_FILE="$REPO_ROOT/results/ihs_${POP}/${POP}.chrX.ihs.tsv"
if [[ -s "$RESULT_FILE" ]]; then
    log "✓ iHS results already exist for $POP"
    log "  Location: $RESULT_FILE"
    echo ""
    log "To rerun, delete the result file first:"
    log "  rm $RESULT_FILE"
    exit 0
fi

# Step 1: Check if subset exists
HAP_FILE="/home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin_${POP}/chrX_${POP}.hapbin.hap"

if [[ -s "$HAP_FILE" ]]; then
    log "✓ $POP subset already exists"
    log "  Location: $HAP_FILE"
else
    log "STEP 1: Subsetting to $POP samples..."
    log "----------------------------------------"
    if python scripts/ThisWorks/subset_population_hapbin.py "$POP" 2>&1 | tee -a "$SINGLE_LOG"; then
        log "  ✓ Subsetting successful"
    else
        log_error "Subsetting failed for $POP"
        log "Check error log: $ERROR_LOG"
        exit 1
    fi
fi

echo ""
log "STEP 2: Submitting iHS job for $POP..."
log "----------------------------------------"

if JOB_ID=$(sbatch --parsable --export=POP="$POP" scripts/ThisWorks/ihs_chrX_population.slurm 2>&1); then
    log "✓ Job submitted: $JOB_ID"
    echo ""
    log "Monitor progress:"
    log "  squeue -u $USER"
    log "  tail -f logs/ihs_chrX_${POP}_${JOB_ID}.out"
    echo ""
    log "When complete, results will be at:"
    log "  $RESULT_FILE"
    echo ""
    log "Pipeline log: $SINGLE_LOG"
else
    log_error "Failed to submit job for $POP"
    log_error "sbatch output: $JOB_ID"
    log "Check error log: $ERROR_LOG"
    exit 1
fi

echo ""
log "=========================================="
log "Job submitted successfully!"
log "=========================================="
log "Finished: $(date)"
echo ""
