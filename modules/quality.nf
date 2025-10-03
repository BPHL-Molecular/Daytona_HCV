process quality {
   input:
      val x
   output:
      //stdout
      //path 'xfile.txt', emit: aLook
      val "${params.output}/${x}", emit: outputpath1
      //path "${params.output}/${x}_trim_2.fastq", emit: trimR2
      
   """  
   #echo ${params.input}/${x}_1.fastq.gz >> xfile.txt
   
   #mkdir -p ${params.output}/assemblies
   mkdir -p ${params.output}/variants
   #mkdir -p ${params.output}/vadr_error_reports
   mkdir -p ${params.output}/${x}
   cp ${params.input}/${x}_*.fastq.gz ${params.output}/${x}
   
   #Run fastqc on original reads
   #singularity exec --cleanenv /apps/staphb-toolkit/containers/fastqc_0.11.9.sif fastqc ${params.output}/${x}/${x}_1.fastq.gz ${params.output}/${x}/${x}_2.fastq.gz
   singularity exec --cleanenv docker://staphb/fastqc:latest fastqc ${params.output}/${x}/${x}_1.fastq.gz ${params.output}/${x}/${x}_2.fastq.gz

   mv ${params.output}/${x}/${x}_1_fastqc.html ${params.output}/${x}/${x}_1_original_fastqc.html
   mv ${params.output}/${x}/${x}_1_fastqc.zip ${params.output}/${x}/${x}_1_original_fastqc.zip
   mv ${params.output}/${x}/${x}_2_fastqc.html ${params.output}/${x}/${x}_2_original_fastqc.html
   mv ${params.output}/${x}/${x}_2_fastqc.zip ${params.output}/${x}/${x}_2_original_fastqc.zip
   
   # Run sra-human-scrubber to remove human reads
   gzip -d ${params.output}/${x}/${x}_1.fastq.gz
   gzip -d ${params.output}/${x}/${x}_2.fastq.gz
      
   singularity exec -B ${params.output}/${x}:/data docker://ncbi/sra-human-scrubber:latest /opt/scrubber/scripts/scrub.sh -r -i /data/${x}_1.fastq -o /data/${x}_1_humanclean.fastq
   singularity exec -B ${params.output}/${x}:/data docker://ncbi/sra-human-scrubber:latest /opt/scrubber/scripts/scrub.sh -r -i /data/${x}_2.fastq -o /data/${x}_2_humanclean.fastq

   gzip ${params.output}/${x}/${x}_1_humanclean.fastq
   gzip ${params.output}/${x}/${x}_2_humanclean.fastq
   
   #Run trimmomatic
   singularity exec --cleanenv docker://staphb/trimmomatic:0.39 trimmomatic PE -phred33 -trimlog ${params.output}/${x}/${x}.log ${params.output}/${x}/${x}_1_humanclean.fastq.gz ${params.output}/${x}/${x}_2_humanclean.fastq.gz ${params.output}/${x}/${x}_trim_1.fastq.gz ${params.output}/${x}/${x}_unpaired_trim_1.fastq.gz ${params.output}/${x}/${x}_trim_2.fastq.gz ${params.output}/${x}/${x}_unpaired_trim_2.fastq.gz SLIDINGWINDOW:4:30 MINLEN:75 TRAILING:20 > ${params.output}/${x}/${x}_trimstats.txt

   rm ${params.output}/${x}/${x}_unpaired_trim_*.fastq.gz
   rm ${params.output}/${x}/${x}_1.fastq ${params.output}/${x}/${x}_2.fastq

   #Run bbduk to remove Illumina adapter sequences and any PhiX contamination  
   #singularity exec --cleanenv /apps/staphb-toolkit/containers/bbtools_38.76.sif bbduk.sh in1=${params.output}/${x}/${x}_trim_1.fastq.gz in2=${params.output}/${x}/${x}_trim_2.fastq.gz out1=${params.output}/${x}/${x}_1.rmadpt.fq.gz out2=${params.output}/${x}/${x}_2.rmadpt.fq.gz ref=/bbmap/resources/adapters.fa stats=${params.output}/${x}/${x}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo
   #singularity exec --cleanenv /apps/staphb-toolkit/containers/bbtools_38.76.sif bbduk.sh in1=${params.output}/${x}/${x}_1.rmadpt.fq.gz in2=${params.output}/${x}/${x}_2.rmadpt.fq.gz out1=${params.output}/${x}/${x}_1.fq.gz out2=${params.output}/${x}/${x}_2.fq.gz outm=${params.output}/${x}/${x}_matchedphix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=${params.output}/${x}/${x}_phixstats.txt
   
   singularity exec --cleanenv docker://staphb/bbtools:38.76 bbduk.sh in1=${params.output}/${x}/${x}_trim_1.fastq.gz in2=${params.output}/${x}/${x}_trim_2.fastq.gz out1=${params.output}/${x}/${x}_1.rmadpt.fq.gz out2=${params.output}/${x}/${x}_2.rmadpt.fq.gz ref=/bbmap/resources/adapters.fa stats=${params.output}/${x}/${x}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo
   singularity exec --cleanenv docker://staphb/bbtools:38.76 bbduk.sh in1=${params.output}/${x}/${x}_1.rmadpt.fq.gz in2=${params.output}/${x}/${x}_2.rmadpt.fq.gz out1=${params.output}/${x}/${x}_1.fq.gz out2=${params.output}/${x}/${x}_2.fq.gz outm=${params.output}/${x}/${x}_matchedphix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=${params.output}/${x}/${x}_phixstats.txt

   rm ${params.output}/${x}/${x}_trim*.fastq.gz
   rm ${params.output}/${x}/${x}*rmadpt.fq.gz
   
   #Run fastqc on clean forward and reverse reads
   #singularity exec --cleanenv /apps/staphb-toolkit/containers/fastqc_0.11.9.sif fastqc ${params.output}/${x}/${x}_1.fq.gz ${params.output}/${x}/${x}_2.fq.gz
   singularity exec --cleanenv docker://staphb/fastqc:latest fastqc ${params.output}/${x}/${x}_1.fq.gz ${params.output}/${x}/${x}_2.fq.gz

   #Rename fastqc output files
   mv ${params.output}/${x}/${x}_1_fastqc.html ${params.output}/${x}/${x}_1_clean_fastqc.html
   mv ${params.output}/${x}/${x}_1_fastqc.zip ${params.output}/${x}/${x}_1_clean_fastqc.zip
   mv ${params.output}/${x}/${x}_2_fastqc.html ${params.output}/${x}/${x}_2_clean_fastqc.html
   mv ${params.output}/${x}/${x}_2_fastqc.zip ${params.output}/${x}/${x}_2_clean_fastqc.zip
   
   #Run multiqc
   #singularity exec --cleanenv /apps/staphb-toolkit/containers/multiqc_1.8.sif multiqc ${params.output}/${x}/${x}_*_fastqc.zip -o ${params.output}/${x}
   singularity exec --cleanenv docker://staphb/multiqc:latest multiqc ${params.output}/${x}/${x}_*_fastqc.zip -o ${params.output}/${x}

   #Map reads to reference
   mkdir ${params.output}/${x}/alignment
   
   #Run kraken to check original fastq data in input folder for their species
   mkdir -p ${params.output}/${x}/kraken_out/
   #singularity exec --cleanenv docker://staphb/kraken2:2.0.8-beta kraken2 --db /kraken2-db/minikraken2_v1_8GB/ --use-names --report ${params.output}/${x}/kraken_out/${x}.report --output ${params.output}/${x}/kraken_out/${x}_kraken.out --paired ${params.input}/${x}_1.fastq.gz ${params.input}/${x}_2.fastq.gz
   #singularity exec --cleanenv docker://staphb/kraken2:latest kraken2 --db /blue/bphl-florida/share/kraken_databases/viral_dba/ --use-names --report ${params.output}/${x}/kraken_out/${x}.report --output ${params.output}/${x}/kraken_out/${x}_kraken.out --paired ${params.input}/${x}_1.fastq.gz ${params.input}/${x}_2.fastq.gz
   singularity exec --cleanenv docker://staphb/kraken2:latest kraken2 --db /blue/bphl-florida/share/kraken_bracken_database/PlusPF/ --use-names --report ${params.output}/${x}/kraken_out/${x}.report --output ${params.output}/${x}/kraken_out/${x}_kraken.out --paired ${params.input}/${x}_1.fastq.gz ${params.input}/${x}_2.fastq.gz
   
   
   """
}
