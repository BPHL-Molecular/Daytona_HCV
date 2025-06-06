#!/bin/bash
#SBATCH --account=bphl-umbrella
#SBATCH --qos=bphl-umbrella
#SBATCH --job-name=hepatitis
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
#SBATCH --mem=20gb
#SBATCH --output=hepatitis.%j.out
#SBATCH --error=hepatitis.%j.err
#SBATCH --time=3-00

module load apptainer
module load nextflow

APPTAINER_CACHEDIR=./
export APPTAINER_CACHEDIR

NXF_SINGULARITY_CACHEDIR=./
export NXF_SINGULARITY_CACHEDIR

if [ "$1" = "HAV" ]; then
   ##### run HAV pipeline
   echo "run HAV pipeline"
   nextflow run hav.nf -params-file params_hav.yaml

   sort ./output/*/report.txt | uniq > ./output/sum_report.txt
   sed -i '/sampleID\tk_species/d' ./output/sum_report.txt
   sed -i '1i sampleID\tk_species\tk_percent\treference\tstart\tend\tnum_raw_reads\tnum_clean_reads\tnum_mapped_reads\tpercent_mapped_clean_reads\tcov_bases_mapped\tpercent_genome_cov_map\tmean_depth\tmean_base_qual\tmean_map_qual\tassembly_length\tassembly_numN\tassembly_percent_ref_genome_cov' ./output/sum_report.txt
   
# elif [ "$1" = "HBV" ]; then
   # ##### run HBV pipeline
   # echo "run HBV pipeline"
   # nextflow run flaq_sc2_clearlabs2.nf -params-file params_clearlabs.yaml

elif [ "$1" = "HCV" ]; then
   ##### run HCV pipeline
   echo "run HCV pipeline"
   nextflow run hcv.nf -params-file params_hcv.yaml
   sort ./output/*/report.txt | uniq > ./output/sum_report.txt
   sed -i '/sampleID\tk_species/d' ./output/sum_report.txt
   sed -i '1i sampleID\tk_species\tk_percent\treference\tstart\tend\tnum_raw_reads\tnum_clean_reads\tnum_mapped_reads\tpercent_mapped_clean_reads\tcov_bases_mapped\tpercent_genome_cov_map\tmean_depth\tmean_base_qual\tmean_map_qual\tassembly_length\tassembly_numN\tassembly_percent_ref_genome_cov' ./output/sum_report.txt

else
   ###### exit with error
   echo "Virus-type parameter cannot be empty in commandline. Please select one from "HAV", "HCV"."
   echo "For example: $sbatch Daytona_hepatitis.sh HAV"
fi

cat ./output/assemblies/*.fa > ./output/assemblies.fasta
mv ./*.out ./output
mv ./*err ./output
dt=$(date "+%Y%m%d%H%M%S")
mv ./output ./output-$dt
rm -r ./work