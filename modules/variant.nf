process variant {
    input:
        val mypath
    output:
        //stdout
        val mypath
        //path "pyoutputs.txt", emit: pyoutputs
        
    
    """ 
    samplename=\$(echo ${mypath} | rev | cut -d "/" -f 1 | rev)
    #Call variants
    mkdir ${mypath}/variants
    #samtools mpileup -A -d 8000 --reference ${params.reference}/nCoV-2019.reference.fasta -B -Q 0 ${mypath}/alignment/\${samplename}.primertrim.sorted.bam | ivar variants -r ${params.reference}/nCoV-2019.reference.fasta -m 10 -p ${mypath}/variants/\${samplename}.variants -q 20 -t 0.25 -g ${params.reference}/GCF_009858895.2_ASM985889v3_genomic.gff
    #singularity exec /apps/staphb-toolkit/containers/samtools_1.12.sif samtools mpileup -A -d 8000 --reference ${params.reference}/nCoV-2019.reference.fasta -B -Q 0 ${mypath}/alignment/\${samplename}.primertrim.sorted.bam | ivar variants -r ${params.reference}/nCoV-2019.reference.fasta -m 10 -p ${mypath}/variants/\${samplename}.variants -q 20 -t 0.25 -g ${params.reference}/GCF_009858895.2_ASM985889v3_genomic.gff
    singularity exec docker://staphb/samtools:1.12 samtools mpileup -A -d 8000 --reference ${params.reference}/KC248195.fasta -B -Q 0 ${mypath}/alignment/\${samplename}.sorted.bam | ivar variants -r ${params.reference}/KC248195.fasta -m 10 -p ${mypath}/variants/\${samplename}.variants -q 20 -t 0.25 -g ${params.reference}/KC248195.gff3

    #Generate consensus assembly
    #mkdir ${mypath}/assembly
    #samtools mpileup -A -B -d 8000 --reference ${params.reference}/nCoV-2019.reference.fasta -Q 0 ${mypath}/alignment/\${samplename}.primertrim.sorted.bam | ivar consensus -t 0 -m 10 -n N -p ${mypath}/assembly/\${samplename}.consensus
    #singularity exec /apps/staphb-toolkit/containers/samtools_1.12.sif samtools mpileup -A -B -d 8000 --reference ${params.reference}/nCoV-2019.reference.fasta -Q 0 ${mypath}/alignment/\${samplename}.primertrim.sorted.bam | ivar consensus -t 0 -m 10 -n N -p ${mypath}/assembly/\${samplename}.consensus
    #singularity exec docker://staphb/samtools:1.12 samtools mpileup -A -B -d 8000 --reference ${params.reference}/KC248195.fasta -Q 0 ${mypath}/alignment/\${samplename}.sorted.bam | ivar consensus -t 0 -m 10 -n N -p ${mypath}/assembly/\${samplename}.consensus

    """
}
