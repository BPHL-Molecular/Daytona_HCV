#!/bin/bash
#SBATCH --account=bphl-umbrella
#SBATCH --qos=bphl-umbrella
#SBATCH --job-name=hepatitis_nf
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=300gb
#SBATCH --output=hepatitis_pipeline.out
#SBATCH --error=hepatitis_pipeline.err
#SBATCH --open-mode=truncate
#SBATCH --time=12:00:00

module load apptainer
module load nextflow

APPTAINER_CACHEDIR=./
export APPTAINER_CACHEDIR
NXF_SINGULARITY_CACHEDIR=./
export NXF_SINGULARITY_CACHEDIR

echo "==================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] PART 1: Nextflow pipeline started"
echo "  SLURM_JOB_ID = ${SLURM_JOB_ID:-N/A}"
echo "==================================================================="

##### run HCV pipeline (Nextflow only)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] run HCV pipeline"
nextflow run hcv.nf -params-file params_hcv.yaml
NF_EXIT=$?
if [ $NF_EXIT -ne 0 ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: nextflow exited with code $NF_EXIT" >&2
    exit $NF_EXIT
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Generating sum_report.txt"
sort ./output/*/report.txt | uniq > ./output/sum_report.txt
sed -i '/sampleID\tk_species|tax_ID|percent(%)|number/d' ./output/sum_report.txt
sed -i '1i sampleID\tk_species|tax_ID|percent(%)|number\treference\tstart\tend\tnum_raw_reads\tnum_clean_reads\tnum_mapped_reads\tpercent_mapped_clean_reads\tcov_bases_mapped\tpercent_genome_cov_map\tmean_depth\tmean_base_qual\tmean_map_qual' ./output/sum_report.txt

echo "==================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] PART 1 (Nextflow) finished successfully"
echo "Next step: sbatch Daytona_HCV_part2_phylogeny.sh [genotype|subtype|both]"
echo "==================================================================="
