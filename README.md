# Daytona_HCV        

## Introduction

This pipeline can perform the genotype detection and phylogeny analysis of Heptatits C virus (HCV). Illumina paired-end sequencing data are required for the pipeline. 

## Prerequisites
Nextflow is needed. The details of installation can be found at https://github.com/nextflow-io/nextflow. For HiPerGator users, its installation is not needed. 

Singularity/APPTAINER is needed. The details of installation can be found at https://singularity-tutorial.github.io/01-installation/. For HiPerGator users, its installation is not needed.

SLURM is needed. For HiPerGator users, its installation is not needed.

Python3 is needed. The package "pandas" should be installed by ``` pip3 install pandas ``` if not included in your python3.

The Kraken2 database PlusPF is needed. For HiPerGator users, downloading is not needed. It has been downloaded and configured in the pipeline.

PhyTreeViz is needed. The installation can be found at https://github.com/moshi4/phyTreeViz.

## Workflow
```mermaid
gitGraph
       branch Daytona_HCV
       commit id: "Paired-end fastq data"
       
       branch QC
       checkout QC
       commit id: "Fastqc"
       commit id: "Trimmomatic"
       commit id: "bbtools"
       commit id: "multiqc"
       checkout Daytona_HCV
       merge QC
       branch SNP_calling
       checkout SNP_calling
       commit id: "Fastqc"
       commit id: "Trimmomatic"
       commit id: "bbtools"
       commit id: "multiqc"
       checkout Daytona_HCV
       merge QC 
       commit id: "Taxonomic classification" tag:"Kraken"
       branch Assembly
       checkout Assembly
       commit id: "Genome Assembly" tag:" Megahit | Skesa"
       commit id:"Read Alignment" tag:"BWA"
       commit id:"Sorting and Indexing" tag:"Samtools"
       commit id:"Consensus Assembly" tag:"Pilon"
       checkout Gentyspuds_wf
       merge Assembly
       branch phylogeny
       checkout phylogeny
       commit id: "MSA" tag:"ClustalW"
       commit id: "Phylogenetic tree" tag: "Iqtree"
       commit id: "Alignment and Tree" tag: "MegaX"
       checkout Gentyspuds_wf
       merge phylogeny
       commit id:"Local_Alignment" tag:"Blast"
       commit id: "EV_typing Classification" tag:"Kraken"
    
```



## Recommended conda environment installation
   ```bash
   conda create -n HCV -c conda-forge python=3.10
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
### 1. Summary report      
|sampleID|k_species/percent(%)/number|...|        
|:---|:---|:---|             
|xxx25002686_S1|Hepacivirus hominis/0.07/875,Hepatitis C virus genotype 4/0.07/872,Hepatitis C virus genotype 6/0.00/2|...|       
                                              
The second column of the above table indicates that 875 reads (0.07%) in the sample (xxx25002686_S1) are identified as HCV species. Among it, 872 reads (0.07%) are identified as HCV genotype 4, while 2 reads (0.00%) are identified as HCV genotype 6. Note, the reason why the two percentages are the same is due to rounding. Similarly, the same reason applies to 0 percentage.      
### 2. Variants    
|REGION|POS|REF|ALT|...|PVAL|PASS|...|        
|:---|:---|:---|:---|:---|:---|:---|:---|             
|KC248195.1|28|T|C|...|0.16|FALSE|...|
|KC248195.1|107|G|A|...|0.04|TRUE|...|                   

PASS is the result of p-value <= 0.05. If a SNP's PASS value is FALSE, it fails to pass the quality check.      
## Test data
The test data can be found in /fastqs/example_hcv/. To use them, please copy them to /fastqs first.

### Note:      
If you want to get email notification when the pipeline running ends, please input your email address in the line "#SBATCH --mail-user=<EMAIL>" of Daytona_HCV.sh.  
