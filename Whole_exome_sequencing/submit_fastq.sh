#!/bin/bash

# copy your raw fastq files into 01.raw_reads/, then use "bash submit_fastq.sh" to submit.

raw_fastq_files=(01.raw_reads/*.f*q.gz)
i=0
    while [ $i -lt ${#raw_fastq_files[@]} ]
    do
        R1_fastq=${raw_fastq_files[$i]} 
        R2_fastq=${raw_fastq_files[$i+1]}
        outname=$(echo $R1_fastq | grep -Eo "P[0-9]{1,3}(-P?[0-9]{6,8})?-(N|T|LN)")
        echo $R1_fastq
        echo $R2_fastq
        echo $outname
        sbatch -p cn-long -A zeminz_g1 --qos=zeminzcnl -J lzd -c 20 -o ./log/$outname.preprocessing.out ./fastq_to_bam.sh $R1_fastq $R2_fastq $outname
        i=`expr $i + 2`
    done
