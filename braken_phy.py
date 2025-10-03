import sys
import os
import argparse
current_dir = os.getcwd()
parser = argparse.ArgumentParser(description="Parameters for braken_phy.py.")
parser.add_argument("-p", "--phy", type=int, default=1, help="which phylogeny level is generated? Input 1, 2 or 3. 1 means using 7 HCV genetypes as references. 2 means using 61 HCV subtypes as references. 3 means both above.")
args = parser.parse_args()

with open("./output/sum_report.txt", 'r') as kreport:
        lines = kreport.readlines()
        for l in lines[1:]:
            l_parse = l.lstrip().rstrip().split("\t")
            sampleID = l_parse[0].lstrip()
            species_group = l_parse[1].lstrip().rstrip().split(",")
            species_items = species_group[1].lstrip().rstrip().split("|")
            tax = species_items[0].lstrip()
            taxID = species_items[1].lstrip()

            #substring1 = "Hepacivirus hominis"
            substring = "Hepatitis C"
            if(substring.lower() in tax.lower()):
                os.system(f"python ./extract_kraken_reads.py -k {current_dir}/output/{sampleID}/kraken_out/{sampleID}_kraken.out -s {current_dir}/fastqs/hcv/{sampleID}_1.fastq.gz -s2 {current_dir}/fastqs/hcv/{sampleID}_2.fastq.gz -o {current_dir}/output/extract/{sampleID}_{taxID}_1.fq -o2 {current_dir}/output/extract/{sampleID}_{taxID}_2.fq -t {taxID}")
                os.system(f"singularity exec docker://staphb/bwa:0.7.17 bwa mem {current_dir}/reference/hcv/KC248195.fasta {current_dir}/output/extract/{sampleID}_{taxID}_1.fq {current_dir}/output/extract/{sampleID}_{taxID}_2.fq > {current_dir}/output/extract/{sampleID}_{taxID}_aln.sam")
                os.system(f"singularity exec docker://staphb/samtools:1.12 samtools view -F 4 -u -h -bo {current_dir}/output/extract/{sampleID}_{taxID}_aln.bam {current_dir}/output/extract/{sampleID}_{taxID}_aln.sam")
                os.system(f"singularity exec docker://staphb/samtools:1.12 samtools sort -n -o {current_dir}/output/extract/{sampleID}_{taxID}.namesorted.bam {current_dir}/output/extract/{sampleID}_{taxID}_aln.bam")
                os.system(f"singularity exec docker://staphb/samtools:1.12 samtools fixmate -m {current_dir}/output/extract/{sampleID}_{taxID}.namesorted.bam {current_dir}/output/extract/{sampleID}_{taxID}.fixmate.bam")
                os.system(f"singularity exec docker://staphb/samtools:1.12 samtools sort -o {current_dir}/output/extract/{sampleID}_{taxID}.positionsort.bam {current_dir}/output/extract/{sampleID}_{taxID}.fixmate.bam")
                #os.system(f"singularity exec docker://staphb/samtools:1.12 samtools markdup {current_dir}/output/extract/{sampleID}_{taxID}.positionsort.bam {current_dir}/output/extract/{sampleID}_{taxID}.markdup.bam")
                os.system(f"singularity exec docker://staphb/samtools:1.12 samtools markdup -r {current_dir}/output/extract/{sampleID}_{taxID}.positionsort.bam {current_dir}/output/extract/{sampleID}_{taxID}.dedup.bam")
                os.system(f"singularity exec docker://staphb/samtools:1.12 samtools sort -o {current_dir}/output/extract/{sampleID}_{taxID}.sorted.bam {current_dir}/output/extract/{sampleID}_{taxID}.dedup.bam")
                os.system(f"singularity exec docker://staphb/samtools:1.12 samtools mpileup -A -B -d 8000 --reference {current_dir}/reference/hcv/KC248195.fasta -Q 0 {current_dir}/output/extract/{sampleID}_{taxID}.sorted.bam | ivar consensus -t 0 -m 10 -n N -p {current_dir}/output/extract/{sampleID}_{taxID}.consensus")
            
            else:
                print(f"No any Hepatitis C virus genotype is found in {sampleID}, based on the kraken output of {sampleID}.")
                print(f"Please check the issue from {sampleID}.")
                break
                
os.system(f"cat {current_dir}/output/extract/*.consensus.fa > {current_dir}/output/extract/sum_consensus.fa")
os.system(f"sed -i 's/>Consensus_/>/g; s/\.consensus_threshold_.*//g' {current_dir}/output/extract/sum_consensus.fa")
#os.system(f"mafft {current_dir}/output/extract/sum_consensus.fa > {current_dir}/output/extract/mafft_msa")
### add new sequences to the aligned reference sequences
if args.phy == 1:
    os.system(f"mafft --addfragments {current_dir}/output/extract/sum_consensus.fa {current_dir}/reference/HCV_align/hcv_7genotypes_aln > {current_dir}/output/extract/mafft_genotype_msa")
    os.system(f"singularity exec docker://staphb/snp-sites:2.3.3 snp-sites -o {current_dir}/output/extract/SNPs_genotype.fasta -c {current_dir}/output/extract/mafft_genotype_msa")
    os.system(f"singularity exec docker://staphb/snp-dists:0.6.2 snp-dists {current_dir}/output/extract/mafft_genotype_msa > {current_dir}/output/extract/pairwise_matrix_genotype.tsv")
    os.system(f"singularity exec docker://staphb/iqtree:1.6.7 iqtree -s {current_dir}/output/extract/SNPs_genotype.fasta -m MFP+ASC -nt AUTO -bb 1000 -alrt 1000 -pre {current_dir}/output/extract/SNPs_boot_genotype -ntmax 9")
    os.system(f"phytreeviz -i {current_dir}/output/extract/SNPs_boot_genotype.contree -o {current_dir}/output/extract/SNPs_boot_genotype.contree.pdf --show_confidence")
    os.system(f"phytreeviz -i {current_dir}/output/extract/SNPs_boot_genotype.contree -o {current_dir}/output/extract/SNPs_boot_genotype.contree.png --show_confidence")
elif args.phy == 2:
    os.system(f"mafft --addfragments {current_dir}/output/extract/sum_consensus.fa {current_dir}/reference/HCV_align/hcv_allsubtypes_aln > {current_dir}/output/extract/mafft_subtype_msa")
    os.system(f"singularity exec docker://staphb/snp-sites:2.3.3 snp-sites -o {current_dir}/output/extract/SNPs_subtype.fasta -c {current_dir}/output/extract/mafft_subtype_msa")
    os.system(f"singularity exec docker://staphb/snp-dists:0.6.2 snp-dists {current_dir}/output/extract/mafft_subtype_msa > {current_dir}/output/extract/pairwise_matrix_subtype.tsv")
    os.system(f"singularity exec docker://staphb/iqtree:1.6.7 iqtree -s {current_dir}/output/extract/SNPs_subtype.fasta -m MFP+ASC -nt AUTO -bb 1000 -alrt 1000 -pre {current_dir}/output/extract/SNPs_boot_subtype -ntmax 9")
    os.system(f"phytreeviz -i {current_dir}/output/extract/SNPs_boot_subtype.contree -o {current_dir}/output/extract/SNPs_boot_subtype.contree.pdf --show_confidence")
    os.system(f"phytreeviz -i {current_dir}/output/extract/SNPs_boot_subtype.contree -o {current_dir}/output/extract/SNPs_boot_subtype.contree.png --show_confidence")
elif args.phy == 3:
    os.system(f"mafft --addfragments {current_dir}/output/extract/sum_consensus.fa {current_dir}/reference/HCV_align/hcv_7genotypes_aln > {current_dir}/output/extract/mafft_genotype_msa")
    os.system(f"singularity exec docker://staphb/snp-sites:2.3.3 snp-sites -o {current_dir}/output/extract/SNPs_genotype.fasta -c {current_dir}/output/extract/mafft_genotype_msa")
    os.system(f"singularity exec docker://staphb/snp-dists:0.6.2 snp-dists {current_dir}/output/extract/mafft_genotype_msa > {current_dir}/output/extract/pairwise_matrix_genotype.tsv")
    os.system(f"singularity exec docker://staphb/iqtree:1.6.7 iqtree -s {current_dir}/output/extract/SNPs_genotype.fasta -m MFP+ASC -nt AUTO -bb 1000 -alrt 1000 -pre {current_dir}/output/extract/SNPs_boot_genotype -ntmax 9")
    os.system(f"phytreeviz -i {current_dir}/output/extract/SNPs_boot_genotype.contree -o {current_dir}/output/extract/SNPs_boot_genotype.contree.pdf --show_confidence")
    os.system(f"phytreeviz -i {current_dir}/output/extract/SNPs_boot_genotype.contree -o {current_dir}/output/extract/SNPs_boot_genotype.contree.png --show_confidence")
    os.system(f"mafft --addfragments {current_dir}/output/extract/sum_consensus.fa {current_dir}/reference/HCV_align/hcv_allsubtypes_aln > {current_dir}/output/extract/mafft_subtype_msa")
    os.system(f"singularity exec docker://staphb/snp-sites:2.3.3 snp-sites -o {current_dir}/output/extract/SNPs_subtype.fasta -c {current_dir}/output/extract/mafft_subtype_msa")
    os.system(f"singularity exec docker://staphb/snp-dists:0.6.2 snp-dists {current_dir}/output/extract/mafft_subtype_msa > {current_dir}/output/extract/pairwise_matrix_subtype.tsv")
    os.system(f"singularity exec docker://staphb/iqtree:1.6.7 iqtree -s {current_dir}/output/extract/SNPs_subtype.fasta -m MFP+ASC -nt AUTO -bb 1000 -alrt 1000 -pre {current_dir}/output/extract/SNPs_boot_subtype -ntmax 9")
    os.system(f"phytreeviz -i {current_dir}/output/extract/SNPs_boot_subtype.contree -o {current_dir}/output/extract/SNPs_boot_subtype.contree.pdf --show_confidence")
    os.system(f"phytreeviz -i {current_dir}/output/extract/SNPs_boot_subtype.contree -o {current_dir}/output/extract/SNPs_boot_subtype.contree.png --show_confidence")
else:
    print("Braken_phy.py Parameter Error!")
    print("Only int 1 or 2 or 3 is accepted for the flag --phy in the script braken_phy.py.")
    sys.exit()
