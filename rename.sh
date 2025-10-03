for f in ./fastqs/hcv/*.fastq.gz; do
    mv "$f" "${f//_001.fastq/.fastq}"
done
for f in ./fastqs/hcv/*.fastq.gz; do
    mv "$f" "${f//_S*_L001_R/_}"
done
