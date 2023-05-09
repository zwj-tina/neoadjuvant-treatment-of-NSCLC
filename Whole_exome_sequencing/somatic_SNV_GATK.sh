#!/bin/bash

normal_bam=$1
tumor_bam=$2
sample=$3
LN=$4
#==================== parameters ================
number_of_threads=35

#================= Define References =====================

bcftools=/lustre1/zeminz_pkuhpc/01.bin/bcftools-1.9/bcftools
panel_of_normals_1000g=../00.resources/008.panel_of_normals_1000g/somatic_hg38_1000g_pon.hg38.vcf.gz
reference_genome=../00.resources/001.hg38/hg38.fa
gatk=../00.resources/002.gatk-4.2.4.0/gatk
af_only_gnomad_vcf=../00.resources/005.af-only/somatic-hg38_af-only-gnomad.hg38.vcf.gz
interval_list=../00.resources/006.interval_list/Exome_Agilent_V6.interval_list
agilent_v6_exome=../00.resources/004.agilent_v6_references/Exome_Agilent_V6.bed
gatk_temp_dir=./gatk_temp_dir/
annovar=../00.resources/007.ANNOVAR/table_annovar.pl
convert_to_avinput=../00.resources/007.ANNOVAR/convert2annovar.pl
mkdir $gatk_temp_dir

module load java/1.8.0_101
#==================== Using Mutect2 to call variants ==================
output_name="$(echo ${tumor_bam} | cut -c13- | rev | cut -c35- | rev )${LN}"
output_vcf="02.raw_variants/${output_name}_raw.vcf.gz"
output_F1R2="02.F1R2/${output_name}_F1R2.tar.gz"

call_variants=true
if $call_variants
then
        tumor_bam_RG="${sample}${LN}"
        normal_bam_RG="${sample}N"
        echo $normal_bam
        $gatk Mutect2 -I $normal_bam -normal $normal_bam_RG \
        -I $tumor_bam -tumor $tumor_bam_RG \
        -O $output_vcf \
        -A AlleleFraction \
        --reference $reference_genome \
        -L $agilent_v6_exome \
        -ip 50 \
        -pon $panel_of_normals_1000g \
        --f1r2-tar-gz $output_F1R2 \
        --germline-resource $af_only_gnomad_vcf \
        --genotype-germline-sites true \
        --genotype-pon-sites true \
        --af-of-alleles-not-in-resource  1e-6 \
        --native-pair-hmm-threads $number_of_threads
fi
#========================= doing GetpileupSummaries ======================== 

GetpileupSummaries=true
if $GetpileupSummaries
then
        normal_pileups_table="03.pileup_summary/${output_name}_normal_pileups.table"
        tumor_pileups_table="03.pileup_summary/${output_name}_tumor_pileups.table"
        $gatk GetPileupSummaries \
        -I $tumor_bam \
        -V $af_only_gnomad_vcf \
        -L $interval_list  \
        -O $tumor_pileups_table
        $gatk GetPileupSummaries \
        -I $normal_bam \
        -V $af_only_gnomad_vcf \
        -L $interval_list  \
        -O $normal_pileups_table
fi
#=============== Calculate Contamination  ================
contamination=true
if $contamination
then
    output_table="04.contamination_data/${sample}${LN}_tumor_contamination.table"
    $gatk CalculateContamination \
    -I $tumor_pileups_table \
    -matched $normal_pileups_table \
    -O $output_table
else
    echo "skip contamination"
fi
#================== Learn F1R2 model =====================

Learn_F1R2_model=true
if $Learn_F1R2_model
then
        output_F1R2_artifact=$(echo $output_F1R2 | cut -c9- | rev | cut -c13- | rev)
        output_F1R2_artifact="02.F1R2/${output_F1R2_artifact}{LN}_artifact.tar.gz"
        $gatk LearnReadOrientationModel \
        -I $output_F1R2 \
        -O $output_F1R2_artifact
fi

#======================== Using FilterMutectCalls to mark false positives =======================
FilterMutectCalls=true
if $FilterMutectCalls
then
        false_marked_vcf="05.false_marked_variants/${output_name}_false_marked.vcf.gz"
        Contamination_table=${output_table}
        F1R2_model=${output_F1R2_artifact}
        $gatk FilterMutectCalls \
        -V $output_vcf \
        -O $false_marked_vcf \
        -R $reference_genome \
        --contamination-table $Contamination_table \
        -ob-priors $F1R2_model \
        -ip 50 \
        --f-score-beta 1.0
fi
#======================== Use bcftools to filter vcf files and index them ==============
false_marked_vcf="05.false_marked_variants/${output_name}_false_marked.vcf.gz"
filtering=true
if $filtering
then
        filtered_vcf=$(echo $false_marked_vcf | cut -c26- | rev | cut -c8- | rev)
        filtered_vcf="05.filtered_variants/${filtered_vcf}_filtered.vcf.gz"
        $bcftools view -I $false_marked_vcf -o $filtered_vcf --threads $number_of_threads -f "PASS" -O z
        $bcftools index $filtered_vcf
fi


#======================== convert VCF to avinput ========================

Convertion=true
if $Convertion
then
        avinput=$(echo $filtered_vcf | cut -c22- | rev | cut -c8- | rev)
        avinput="06.avinput/${avinput}.avinput"
        echo $avinput
        echo $filtered_vcf
        perl $convert_to_avinput --withfreq --allsample --format vcf4 $filtered_vcf > $avinput
fi

#======================== Annotate the VCF files =========================

Annotation=true
if $Annotation
then
        output_csv="07.annotated_filtered_variants/${output_name}_annotated"
        perl $annovar \
         --buildver hg38 \
         --out $output_csv \
         --csvout \
         --protocol refGene,gnomad211_exome,clinvar_20210501,exac03,icgc28  \
         --operation g,f,f,f,f \
         --nastring . \
         --remove \
        $avinput ../00.resources/007.ANNOVAR/humandb
fi
#======================= Add VAF into variant files =========================

Add_VAF=true
if $Add_VAF
then
    echo $filtered_vcf
        output_added_VAF_csv="09.VAF_added.variants/${output_name}_added_VAF.csv"
        echo $sample$LN
        echo $output_added_VAF_csv
        echo "${output_csv}.hg38_multianno.csv"
        $bcftools query -f "[%AF,]\n" $filtered_vcf | sed -e 1i\N_AF,T_AF, | paste - "${output_csv}.hg38_multianno.csv" -d ""  > $output_added_VAF_csv
fi
echo "whole program done"
