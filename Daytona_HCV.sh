#!/bin/bash
#SBATCH --account=bphl-umbrella
#SBATCH --qos=bphl-umbrella
#SBATCH --job-name=hepatitis
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=300gb
#SBATCH --output=hepatitis.%j.out
#SBATCH --error=hepatitis.%j.err
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
sed -i '/sampleID\tk_species|tax_ID|percent(%)|number/d' ./output/sum_report.txt
sed -i '1i sampleID\tk_species|tax_ID|percent(%)|number\treference\tstart\tend\tnum_raw_reads\tnum_clean_reads\tnum_mapped_reads\tpercent_mapped_clean_reads\tcov_bases_mapped\tpercent_genome_cov_map\tmean_depth\tmean_base_qual\tmean_map_qual' ./output/sum_report.txt

##### The most dominant genotype is selected in each sample, and then a SNP-based phylogenetic tree is constructed
mkdir ./output/extract
if [ "${1:-genotype}" = "genotype" ]; then
   python ./braken_phy.py -p 1
elif [ "$1" = "subtype" ]; then
   python ./braken_phy.py -p 2
elif [ "$1" = "both" ]; then
   python ./braken_phy.py -p 3
else
   ###### exit with error
   echo "The phylogeny level (genotype or subtype) should be added in commandline. If genotype, the 7 HCV genotypes will be used as references in the phylogeny. If subtype, the 61 HCV subtypes will be used as references in the phylogeny. If both, the two phylogenies will be generated. If empty, the default is genotype."
   echo "For example:"
   echo "$sbatch Daytona_hepatitis.sh"
   echo "$sbatch Daytona_hepatitis.sh genotype"
   echo "$sbatch Daytona_hepatitis.sh subtype"
   echo "$sbatch Daytona_hepatitis.sh both"
   exit 1 # Exit with an error status
fi

#### move tree relevant files to the folder tree
mkdir ./output/trees
mv ./output/extract/SNPs* ./output/trees
mv ./output/extract/pairwise_matrix* ./output/trees
mv ./output/extract/mafft_* ./output/trees

mv ./*.out ./output
mv ./*err ./output
dt=$(date "+%Y%m%d%H%M%S")
mv ./output ./output-$dt
rm -r ./work