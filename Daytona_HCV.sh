#!/bin/bash
#SBATCH --account=bphl-umbrella
#SBATCH --qos=bphl-umbrella
#SBATCH --job-name=hepatitis
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100gb
#SBATCH --output=hcv.%j.out
#SBATCH --error=hcv.%j.err
#SBATCH --time=12:00:00


module load apptainer
module load nextflow

APPTAINER_CACHEDIR=./
export APPTAINER_CACHEDIR

NXF_SINGULARITY_CACHEDIR=./
export NXF_SINGULARITY_CACHEDIR


##### run HCV pipeline
echo "run HCV pipeline"
nextflow run hcv.nf -params-file params_hcv.yaml
sort ./output/*/report.txt | uniq > ./output/sum_report.txt
sed -i '/sampleID\tk_species|percent(%)|number/d' ./output/sum_report.txt
sed -i '1i sampleID\tk_species|percent(%)|number\treference\tstart\tend\tnum_raw_reads\tnum_clean_reads\tnum_mapped_reads\tpercent_mapped_clean_reads\tcov_bases_mapped\tpercent_genome_cov_map\tmean_depth\tmean_base_qual\tmean_map_qual' ./output/sum_report.txt


#cat ./output/assemblies/*.fa > ./output/assemblies.fasta
mv ./*.out ./output
mv ./*err ./output
dt=$(date "+%Y%m%d%H%M%S")
mv ./output ./output-$dt
rm -r ./work