#!/bin/bash
# $1 is your long read assembly as fasta, $2 is fw readfile $3 is rv readfile, $4 is longreadfile, $5 is number of rounds
# This script runs iterative pilon correction rounds on your long read assembly
# You know the drill. No spaces etc in your filename
# Needs bwa, samtools, minialign and pilon

cp "$1" round1.fasta
fw="$2"
rv="$3"
lr="$4"


for round in  $(seq 1 $5)
do
    fasta=round"$round".fasta
    echo "This is round $round"

    echo "indexing reference $fasta"
    minialign -d $fasta.mai $fasta
    bwa index $fasta

    echo "mapping short reads to $fasta"
    bwa mem -M -t 16 $fasta $fw $rv > $fasta.illumina.bam
    samtools view -b -F 4 $fasta.illumina.bam > $fasta.illumina.mapped.bam
    samtools sort $fasta.illumina.mapped.bam -o $fasta.illumina.mapped.sorted.bam
    samtools index $fasta.illumina.mapped.sorted.bam

    echo "mapping long reads to $fasta"
    minialign -t 16 $fasta $lr  > $fasta.pacbio.sam
    samtools view -b -F 4 $fasta.pacbio.sam >  $fasta.pacbio.mapped.bam
    samtools sort $fasta.pacbio.mapped.bam -o $fasta.pacbio.mapped.sorted.bam
    samtools index $fasta.pacbio.mapped.sorted.bam

    echo "Pilon polishing round $round"
    pilon --genome $fasta --frags $fasta.illumina.mapped.sorted.bam --unpaired $fasta.pacbio.mapped.sorted.bam
    mv pilon.fasta round$((round+1)).fasta
done

mv round$((round+1)).fasta pilon_corrected.fasta
infoseq pilon_corrected.fasta

echo "Pilon correction in $5 rounds finished"
