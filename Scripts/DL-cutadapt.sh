#!/bin/bash

primerF="TCACCCAAAGCTGRARTTCTA"
primerR="GCGGGTTGCTGGTTTCACG"

revcomp_primerF=`echo $primerF | \
tr ACGTUWSMKRYBDHVNIacgtuwsmkrybdhvni \
TGCAAWSKMYRVHDBNNtgcaawskmyrvhdbnn | rev`
revcomp_primerR=`echo $primerR | \
tr ACGTUWSMKRYBDHVNIacgtuwsmkrybdhvni \
TGCAAWSKMYRVHDBNNtgcaawskmyrvhdbnn | rev`

for i in *_R1_001.fastq.gz; do
    SAMPLE=$(echo ${i} | sed "s/_R1_\001\.fastq\.gz//") &&
    cutadapt -g "${primerF};required...${revcomp_primerR};optional" \
    -G "${primerR};required...${revcomp_primerF};optional" --discard-untrimmed \
    -m 50 -o cutadapt/${SAMPLE}_R1_001-trimmed.fastq.gz \
    -p cutadapt/${SAMPLE}_R2_001-trimmed.fastq.gz \
    ${SAMPLE}_R1_001.fastq.gz ${SAMPLE}_R2_001.fastq.gz 1> cutadapt-report.txt
done
