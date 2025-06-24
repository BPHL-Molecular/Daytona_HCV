# Daytona_hep        

## Introduction

This pipeline is designed for the analysis of Heptatits A and C data from Illumina paired-end sequencing. It performs quality control, species identification, abundance estimation, SNP calling, consensus sequence assembly, etc.

## Prerequisites
Nextflow is needed. The details of installation can be found at https://github.com/nextflow-io/nextflow. For HiPerGator users, its installation is not needed. 

Singularity/APPTAINER is needed. The details of installation can be found at https://singularity-tutorial.github.io/01-installation/. For HiPerGator users, its installation is not needed.

SLURM is needed. For HiPerGator users, its installation is not needed.

Python3 is needed. The package "pandas" should be installed by ``` pip3 install pandas ``` if not included in your python3.

The Kraken2 database minikraken2_v1_8GB is needed. For HiPerGator users, downloading is not needed. It has been downloaded and configured in the pipeline.

## Recommended conda environment installation
   ```bash
   conda create -n HEP -c conda-forge python=3.10 pandas
   ```
   ```bash
   conda activate HEP
   ```


## How to run

### For HAV data: 
1. put your data files into the directory /fastqs/hav. Your data file's name should look like  "XZA22002292_1.fastq.gz" and "XZA22002292_2.fastq.gz" 
2. open the file "params_hav.yaml", and set the parameters absolute paths. They should be ".../.../fastqs/hav", ".../.../output", etc. 
3. get to the top directory of the pipeline, run 
```bash
sbatch Daytona_hepatitis.sh HAV
```
### For HCV data: 
1. put your data files into the directory /fastqs/hcv. Your data file's name should look like "XZA22002292_1.fastq.gz", "XZA22002292_2.fastq.gz" 
2. open the file "params_hcv.yaml", and set the parameters absolute paths. They should be ".../.../fastqs/hcv", ".../.../output", etc. 
3. get into the directory of the pipeline, run 
```bash
sbatch Daytona_hepatitis.sh HCV
```
#### Note:      
If you want to get email notification when the pipeline running ends, please input your email address in the line "#SBATCH --mail-user=<EMAIL>" of Daytona_hepatitis.sh.  
