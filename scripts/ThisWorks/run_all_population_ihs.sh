#!/usr/bin/env bash
# Master pipeline to run iHS analysis on ALL individual populations
#
# This pipeline will:
# 1. Subset chromosome X data to each population separately
# 2. Run iHS on each population
# 3. Enable population-specific selection signal detection for pairwise comparisons
#
# Populations analyzed:
# - African (AFR super-pop): ACB, ASW, ESN, GWD, LWK, MSL, YRI
# - European (EUR super-pop): CEU, FIN, GBR, IBS, TSI
# - East Asian (EAS super-pop): CDX, CHB, CHS, JPT, KHV
# - South Asian (SAS super-pop): BEB, GIH, ITU, PJL, STU
# - American (AMR super-pop): CLM, MXL, PEL, PUR
# - Plus: CHD (unassigned)

set -euo pipefail

REPO_ROOT="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib"
cd "$REPO_ROOT"

# Setup logging
LOG_DIR="$REPO_ROOT/logs"
mkdir -p "$LOG_DIR"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
PIPELINE_LOG="$LOG_DIR/population_ihs_pipeline_${TIMESTAMP}.log"
ERROR_LOG="$LOG_DIR/population_ihs_errors_${TIMESTAMP}.log"

# Function to log messages
log() {
    echo "$1" | tee -a "$PIPELINE_LOG"
}

log_error() {
    echo "ERROR: $1" | tee -a "$PIPELINE_LOG" >> "$ERROR_LOG"
}

log "=========================================="
log "Population-Specific iHS Pipeline"
log "Started: $(date)"
log "Log file: $PIPELINE_LOG"
log "Error log: $ERROR_LOG"
log "=========================================="
echo ""

# All populations with relate results
POPULATIONS=(
    ACB ASW BEB CDX CEU CHB CHD CHS CLM ESN
    FIN GBR GIH GWD IBS ITU JPT KHV LWK MSL
    MXL PEL PJL PUR STU TSI YRI
)

echo "=========================================="
echo "Population-Specific iHS Analysis Pipeline"
echo "=========================================="
echo ""
echo "This pipeline will run iHS on ${#POPULATIONS[@]} populations separately"
echo ""
echo "Populations to analyze:"
for pop in "${POPULATIONS[@]}"; do
    echo "  - $pop"
done
echo ""
echo "=========================================="
echo ""

# Track what's been done
declare -a ALREADY_SUBSET=()
declare -a NEED_SUBSET=()
declare -a ALREADY_RUN=()
declare -a NEED_RUN=()

# Check status for each population
echo "Checking current status..."
echo ""

for pop in "${POPULATIONS[@]}"; do
    HAP_FILE="/home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin_${pop}/chrX_${pop}.hapbin.hap"
    RESULT_FILE="$REPO_ROOT/results/ihs_${pop}/${pop}.chrX.ihs.tsv"

    if [[ -s "$HAP_FILE" ]]; then
        ALREADY_SUBSET+=("$pop")
    else
        NEED_SUBSET+=("$pop")
    fi

    if [[ -s "$RESULT_FILE" ]]; then
        ALREADY_RUN+=("$pop")
    else
        NEED_RUN+=("$pop")
    fi
done

echo "Status summary:"
echo "  Already subset: ${#ALREADY_SUBSET[@]} populations"
echo "  Need subsetting: ${#NEED_SUBSET[@]} populations"
echo "  Already have iHS results: ${#ALREADY_RUN[@]} populations"
echo "  Need iHS analysis: ${#NEED_RUN[@]} populations"
echo ""

# Show details
if [[ ${#ALREADY_RUN[@]} -gt 0 ]]; then
    echo "Populations with completed iHS:"
    for pop in "${ALREADY_RUN[@]}"; do
        echo "  ✓ $pop"
    done
    echo ""
fi

if [[ ${#NEED_RUN[@]} -eq 0 ]]; then
    echo "=========================================="
    echo "All populations already analyzed!"
    echo "=========================================="
    echo ""
    echo "Results available in: $REPO_ROOT/results/ihs_*/"
    exit 0
fi

# Ask for confirmation
echo "=========================================="
echo "Ready to proceed?"
echo "=========================================="
echo ""
echo "This will:"
echo "  1. Subset ${#NEED_SUBSET[@]} population(s) if needed"
echo "  2. Submit ${#NEED_RUN[@]} iHS job(s) to SLURM"
echo ""
read -p "Continue? (y/n) " -n 1 -r
echo ""
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Cancelled."
    exit 0
fi

echo ""
echo "=========================================="
echo "Starting pipeline..."
echo "=========================================="
echo ""

# Step 1: Subset populations that need it
declare -a SUBSET_FAILED=()
declare -a SUBSET_SUCCESS=()

if [[ ${#NEED_SUBSET[@]} -gt 0 ]]; then
    log "STEP 1: Subsetting populations..."
    log "----------------------------------------"
    echo ""

    for pop in "${NEED_SUBSET[@]}"; do
        log "Subsetting $pop..."

        # Run subsetting and capture output
        if python scripts/ThisWorks/subset_population_hapbin.py "$pop" >> "$PIPELINE_LOG" 2>&1; then
            SUBSET_SUCCESS+=("$pop")
            log "  ✓ $pop subset successful"
        else
            log_error "Failed to subset $pop"
            SUBSET_FAILED+=("$pop")
        fi
        echo ""
    done

    if [[ ${#SUBSET_FAILED[@]} -gt 0 ]]; then
        log_error "Subsetting failed for ${#SUBSET_FAILED[@]} population(s):"
        for pop in "${SUBSET_FAILED[@]}"; do
            log_error "  - $pop"
        done
        echo ""
        log "Continuing with populations that succeeded..."
        echo ""
    fi
else
    log "STEP 1: All populations already subset ✓"
    echo ""
fi

# Step 2: Submit iHS jobs
declare -a JOB_IDS=()
declare -a SUBMIT_FAILED=()

log "STEP 2: Submitting iHS jobs..."
log "----------------------------------------"
echo ""

for pop in "${NEED_RUN[@]}"; do
    # Skip populations that failed subsetting
    if [[ " ${SUBSET_FAILED[@]} " =~ " ${pop} " ]]; then
        log_error "Skipping $pop (subsetting failed)"
        continue
    fi

    # Check that subset exists
    HAP_FILE="/home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin_${pop}/chrX_${pop}.hapbin.hap"

    if [[ ! -s "$HAP_FILE" ]]; then
        log_error "HAP file missing for $pop: $HAP_FILE"
        SUBMIT_FAILED+=("$pop")
        continue
    fi

    log "Submitting iHS job for $pop..."
    if JOB_ID=$(sbatch --parsable --export=POP="$pop" scripts/ThisWorks/ihs_chrX_population.slurm 2>&1); then
        log "  ✓ Job submitted: $JOB_ID"
        JOB_IDS+=("$pop:$JOB_ID")
    else
        log_error "Failed to submit job for $pop"
        SUBMIT_FAILED+=("$pop")
    fi
    echo ""
done

echo ""
log "=========================================="
log "Pipeline Complete!"
log "=========================================="
echo ""
log "Summary:"
log "  Jobs submitted: ${#JOB_IDS[@]}"
if [[ ${#SUBSET_FAILED[@]} -gt 0 ]]; then
    log "  Subsetting failures: ${#SUBSET_FAILED[@]}"
fi
if [[ ${#SUBMIT_FAILED[@]} -gt 0 ]]; then
    log "  Job submission failures: ${#SUBMIT_FAILED[@]}"
fi
echo ""

if [[ ${#JOB_IDS[@]} -gt 0 ]]; then
    log "Jobs submitted:"
    for job_info in "${JOB_IDS[@]}"; do
        pop="${job_info%%:*}"
        job_id="${job_info##*:}"
        log "  $pop: Job $job_id"
    done
    echo ""
fi

if [[ ${#SUBSET_FAILED[@]} -gt 0 ]] || [[ ${#SUBMIT_FAILED[@]} -gt 0 ]]; then
    log "=========================================="
    log "FAILURES DETECTED"
    log "=========================================="
    echo ""

    if [[ ${#SUBSET_FAILED[@]} -gt 0 ]]; then
        log_error "Populations that failed subsetting:"
        for pop in "${SUBSET_FAILED[@]}"; do
            log_error "  - $pop"
        done
        echo ""
    fi

    if [[ ${#SUBMIT_FAILED[@]} -gt 0 ]]; then
        log_error "Populations that failed job submission:"
        for pop in "${SUBMIT_FAILED[@]}"; do
            log_error "  - $pop"
        done
        echo ""
    fi

    log "Check error log for details:"
    log "  cat $ERROR_LOG"
    echo ""
fi

log "Monitor progress:"
log "  squeue -u $USER"
echo ""

if [[ ${#JOB_IDS[@]} -gt 0 ]]; then
    log "Check individual job logs:"
    for job_info in "${JOB_IDS[@]}"; do
        pop="${job_info%%:*}"
        job_id="${job_info##*:}"
        echo "  tail -f logs/ihs_chrX_${pop}_${job_id}.out"
    done
    echo ""
fi

log "Results will be in:"
log "  $REPO_ROOT/results/ihs_<POP>/<POP>.chrX.ihs.tsv"
echo ""
log "Pipeline log: $PIPELINE_LOG"
if [[ -s "$ERROR_LOG" ]]; then
    log "Error log: $ERROR_LOG"
fi
echo ""
log "Finished: $(date)"
echo ""
