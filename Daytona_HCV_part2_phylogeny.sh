#!/bin/bash
#SBATCH --account=bphl-umbrella
#SBATCH --qos=bphl-umbrella
#SBATCH --job-name=hepatitis_phy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100gb
#SBATCH --output=hepatitis_pipeline.out
#SBATCH --error=hepatitis_pipeline.err
#SBATCH --open-mode=append
#SBATCH --time=12:00:00

# Strict mode (no -e because we handle exit codes manually for the case branches)
set -uo pipefail

module load apptainer
module load mafft

APPTAINER_CACHEDIR=./
export APPTAINER_CACHEDIR
NXF_SINGULARITY_CACHEDIR=./
export NXF_SINGULARITY_CACHEDIR

echo ""
echo "==================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] PART 2: Phylogeny analysis started"
echo "  SLURM_JOB_ID = ${SLURM_JOB_ID:-N/A}"
echo "  argument     = ${1:-genotype} (default: genotype)"
echo "==================================================================="

##### Safety check: ./output must exist (produced by Part 1)
if [ ! -d ./output ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: ./output directory not found." >&2
    echo "Please run 'sbatch Daytona_HCV_part1_nextflow.sh' first." >&2
    exit 1
fi

##### Setup logging directory for braken_phy.py
LOG_DIR=./output/logs
mkdir -p "$LOG_DIR"
TS=$(date '+%Y%m%d_%H%M%S')
export BRAKEN_LOG_TS="$TS"   # passed to braken_phy.py so its python-side log filename matches
BRAKEN_LOG="$LOG_DIR/braken_phy_${TS}.log"
BRAKEN_ERR="$LOG_DIR/braken_phy_${TS}.err"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Detailed braken_phy.py logs:"
echo "  shell-captured stdout+stderr : $BRAKEN_LOG"
echo "  shell-captured stderr only   : $BRAKEN_ERR"
echo "  python-side structured log   : $LOG_DIR/braken_phy_py_${TS}.log"

##### Phylogeny level argument parsing
PHY_ARG="${1:-genotype}"
case "$PHY_ARG" in
    genotype) PHY_LEVEL=1 ;;
    subtype)  PHY_LEVEL=2 ;;
    both)     PHY_LEVEL=3 ;;
    *)
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: invalid argument '$PHY_ARG'" >&2
        echo "The phylogeny level (genotype or subtype) should be added in commandline." >&2
        echo "  genotype : 7 HCV genotypes as references (default)" >&2
        echo "  subtype  : 61 HCV subtypes as references" >&2
        echo "  both     : both genotype and subtype phylogenies" >&2
        echo "Examples:" >&2
        echo "  sbatch Daytona_HCV_part2_phylogeny.sh" >&2
        echo "  sbatch Daytona_HCV_part2_phylogeny.sh genotype" >&2
        echo "  sbatch Daytona_HCV_part2_phylogeny.sh subtype" >&2
        echo "  sbatch Daytona_HCV_part2_phylogeny.sh both" >&2
        exit 1
        ;;
esac

##### The most dominant genotype is selected in each sample, and then a SNP-based phylogenetic tree is constructed
mkdir -p ./output/extract

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running: python -u ./braken_phy.py -p $PHY_LEVEL  (mode: $PHY_ARG)"
echo "-------------------------------------------------------------------"

# python -u  : unbuffered output so logs appear in real time
# tee        : write to dedicated braken log AND to the merged pipeline log
# 2> >(...)  : separately capture stderr to its own .err file while still propagating
# PIPESTATUS : capture python's exit code (not tee's)
python -u ./braken_phy.py -p "$PHY_LEVEL" \
    2> >(tee "$BRAKEN_ERR" >&2) \
    | tee "$BRAKEN_LOG"
BRAKEN_EXIT=${PIPESTATUS[0]}

echo "-------------------------------------------------------------------"
if [ "$BRAKEN_EXIT" -ne 0 ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: braken_phy.py exited with code $BRAKEN_EXIT" >&2
    echo "  Combined log : $BRAKEN_LOG" >&2
    echo "  Stderr only  : $BRAKEN_ERR" >&2
    echo "  Last 30 lines of stderr:" >&2
    tail -n 30 "$BRAKEN_ERR" >&2 || true
    exit "$BRAKEN_EXIT"
fi
echo "[$(date '+%Y-%m-%d %H:%M:%S')] braken_phy.py completed successfully (exit 0)"

#### move tree relevant files to the folder tree
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Moving tree files to ./output/trees/"
mkdir -p ./output/trees
mv ./output/extract/SNPs*            ./output/trees/ 2>/dev/null || echo "  WARNING: no SNPs* files to move"
mv ./output/extract/pairwise_matrix* ./output/trees/ 2>/dev/null || echo "  WARNING: no pairwise_matrix* files to move"
mv ./output/extract/mafft_*          ./output/trees/ 2>/dev/null || echo "  WARNING: no mafft_* files to move"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Archiving SLURM logs and renaming output directory"
mv ./*.out ./output/ 2>/dev/null || true
mv ./*err  ./output/ 2>/dev/null || true
dt=$(date '+%Y%m%d%H%M%S')
mv ./output ./output-"$dt"
rm -rf ./work

echo "==================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] PART 2 (Phylogeny) finished successfully"
echo "Final output directory : ./output-$dt"
echo "Merged pipeline log    : ./output-$dt/hepatitis_pipeline.out"
echo "Merged pipeline stderr : ./output-$dt/hepatitis_pipeline.err"
echo "braken_phy.py logs     : ./output-$dt/logs/braken_phy_*.{log,err}"
echo "==================================================================="
