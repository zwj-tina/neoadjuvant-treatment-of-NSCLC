#!/bin/bash

source /lustre2/zeminz_pkuhpc/zhangwenjie/bashrc_zwj

conda activate base

R1_fastq=$1
R2_fastq=$2
outname=$3
#=================== parameters ================
number_of_threads=35


#======================= locate resources====================
reference_genome=../00.resources/001.hg38/hg38.fa
known_indels=../00.resources/003.gatk_bundle/Homo_sapiens_assembly38.known_indels.vcf.gz
other_indels=../00.resources/003.gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
agilent_v6_exome=../00.resources/004.agilent_v6_references/Exome_Agilent_V6.bed
gatk=../00.resources/002.gatk-4.2.4.0/gatk
gatk_temp_dir=./gatk_temp_dir/
mkdir $gatk_temp_dir

module load java/1.8.0_101

#========================= do pre_trim_qc=====================
pre_trim_qc=true
if $pre_trim_qc
then
    echo doing pre_trim qc
    fastqc $R1_fastq -t $number_of_threads -o 02.pre_trim_qc/
    fastqc $R2_fastq -t $number_of_threads -o 02.pre_trim_qc/
    echo finished pre_trim_qc
else
    echo skipping pre_trim_qc
fi
echo -e "\n"

#======================== do adaptor trimming ================
adaptor_trimming=true
R1_fastq_trimmed="03.post_trim_fastq/${outname}_R1_trimmed.fq.gz"
R2_fastq_trimmed="03.post_trim_fastq/${outname}_R2_trimmed.fq.gz"
if $adaptor_trimming
then
    echo doing adaptor trimming
    P7_adaptors_in_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    P5_adaptors_in_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
        echo processing:
        echo $R1_fastq
        echo $R2_fastq
        echo to be saved as:
        echo $R1_fastq_trimmed
        echo $R2_fastq_trimmed
    
        cutadapt \
        -a $P7_adaptors_in_R1 \
        -A $P5_adaptors_in_R2 \
        -O 5 -m 50 --cores=$number_of_threads \
        -o $R1_fastq_trimmed \
        -p $R2_fastq_trimmed \
        $R1_fastq \
        $R2_fastq \

    echo "finished adaptor trimming"
else
    echo skipping adaptor-trimming
fi
echo -e "\n"


#============================ do post_trim_qc ==========================


post_trim_qc=true
if $post_trim_qc
then
    echo doing post_trim qc
    echo identified the following fastq files:
    fastqc $R1_fastq_trimmed -t $number_of_threads -o 04.post_trim_qc/
    fastqc $R2_fastq_trimmed -t $number_of_threads -o 04.post_trim_qc/
    echo finished post_trim_qc
else
    echo skipping pre_trim_qc
fi
echo -e "\n"

#========================== do bwa-mem =====================
# -M is used to accomodate Picard's markDuplicate

conda activate ydc_cnv
sample=$outname
bam_file="05.bam_files/${sample}.bam"
bwa_mem=true
if $bwa_mem
then
         bwa mem -t $number_of_threads -M \
             -R "@RG\tID:${sample}\tSM:${sample}\tPL:illumina\tLB:${sample}" \
             $reference_genome \
             $R1_fastq_trimmed \
             $R2_fastq_trimmed \
        | \
        samtools view -o $bam_file
    echo "finished bwa mem"
else
    echo "skipping bwa mem"
fi
echo -e "\n"

#====================== samtools sortsam  ==================
sortsam=true
sorted_bam="06.sorted_bam/${outname}_sorted.bam"
if $sortsam
then
        samtools sort -@ $number_of_threads -m 1G -o $sorted_bam $bam_file
        samtools index -@ $number_of_threads $sorted_bam
    echo "finished sortsam and added index"
else
    echo "skipping sortsam"
fi
echo -e "\n"

#====================== gatk mark duplicates ====================
markduplicates=true
marked_bam="07.duplicates_marked_bam/${outname}_duplicates_marked.bam"
duplicates_metric="07.duplicates_marked_bam/${outname}_duplicates_metrics.txt"
if $markduplicates
then

        $gatk --java-options "-Djava.io.tmpdir=./gatk_temp_dir/" MarkDuplicates \
        -I $sorted_bam \
        -O $marked_bam \
        -M $duplicates_metric
    echo "finished markduplicates"
else
    echo "skipped markduplicates"
fi

#==================== gatk baserecalibration ==================
base_recalibrate=true
bqsr_bam="08.bqsr_bam/${outname}_bqsr.bam"
recal_table=$outname
recal_table="${recal_table}_recal.table"
if $base_recalibrate 
then
        samtools index -@ $number_of_threads $marked_bam
        $gatk --java-options "-Djava.io.tmpdir=./gatk_temp_dir/" BaseRecalibrator \
            --reference $reference_genome \
            --input $marked_bam \
            --output $recal_table \
            -L $agilent_v6_exome \
            -ip 50 \
            --use-original-qualities \
            --known-sites $known_indels \
            --known-sites $other_indels \

        echo "finished base recalibration, moving into applyBQSR"

        $gatk --java-options "-Djava.io.tmpdir=./gatk_temp_dir" ApplyBQSR \
            -R $reference_genome \
            -I $marked_bam \
            -O $bqsr_bam \
            -L $agilent_v6_exome \
            -ip 50 \
            --bqsr-recal-file $recal_table 

    echo "finished base recalibrate"
else
    echo "skipped base recalibrate"
fi

#=============  END =================
echo "whole program finished"










