# Daytona_HCV        

## Introduction

This pipeline is designed for the analysis of Heptatits C virus data from Illumina paired-end sequencing. It performs HCV-species/genotype identification and abundance estimation, SNP calling, etc.

## Prerequisites
Nextflow is needed. The details of installation can be found at https://github.com/nextflow-io/nextflow. For HiPerGator users, its installation is not needed. 

Singularity/APPTAINER is needed. The details of installation can be found at https://singularity-tutorial.github.io/01-installation/. For HiPerGator users, its installation is not needed.

SLURM is needed. For HiPerGator users, its installation is not needed.

Python3 is needed. The package "pandas" should be installed by ``` pip3 install pandas ``` if not included in your python3.

The Kraken2 database PlusPF is needed. For HiPerGator users, downloading is not needed. It has been downloaded and configured in the pipeline.

## Recommended conda environment installation
   ```bash
   conda create -n HCV -c conda-forge python=3.10 pandas
   ```
   ```bash
   conda activate HCV
   ```


## How to run

1. put your data files into the directory /fastqs/. Your data file's name should look like "XZA22002292_1.fastq.gz", "XZA22002292_2.fastq.gz" 
2. open the file "params_hcv.yaml", and set the parameters absolute paths. They should be ".../.../fastqs", ".../.../output", etc. 
3. get into the top directory of the pipeline, run       
```bash
sbatch Daytona_HCV.sh
```       
## Main output
|sampleID|k_species/percent(%)/number|...|        
|:---|:---|:---|             
|xxx25002686_S1|Hepacivirus hominis/0.07/875,Hepatitis C virus genotype 4/0.07/872,Hepatitis C virus genotype 6/0.00/2|...|         


#### Note:      
If you want to get email notification when the pipeline running ends, please input your email address in the line "#SBATCH --mail-user=<EMAIL>" of Daytona_HCV.sh.  
