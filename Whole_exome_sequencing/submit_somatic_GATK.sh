#!/bin/bash
#================= make directories ==============
mkdir 02.raw_variants
mkdir 03.pileup_summary
mkdir 04.contamination_data
mkdir 05.false_marked_variants
mkdir 06.avinput
mkdir 07.annotated_filtered_variants
mkdir 02.F1R2
mkdir 05.filtered_variants
mkdir 09.VAF_added.variants
#================== submit program =================
bqsr_bam=(08.bqsr_bam/*.bam)
i=0
for ((i = 0; i < ${#bqsr_bam[@]}; i++))
    do
        normal_bam=${bqsr_bam[$i]}
        normal_name=$(echo ${bqsr_bam[$i]} | cut -c13- | rev | cut -c35- | rev )
        normal_bam_splited=($(echo $normal_name | tr "-" "\n"))
        normal_bam_type=${normal_bam_splited[2]}
        normal_bam_sample=${normal_bam_splited[1]}
        if [[ $normal_bam_type == N ]]
        then
            echo $normal_name
        for ((j = i-10; j < ${#bqsr_bam[@]}; j++))
        do
            next_bam=${bqsr_bam[$j]}
            next_name=$(echo ${bqsr_bam[$j]} | cut -c13- | rev | cut -c35- | rev )
            next_bam_splited=($(echo $next_bam| tr "-" "\n"))
            next_bam_type=${next_bam_splited[2]}
            next_bam_sample=${next_bam_splited[1]}
            if [ $next_bam_sample = $normal_bam_sample ] 
            then
                if [ $next_bam_type = T ] && [[ ${next_bam_splited[3]} != "core" ]] && [[ ${next_bam_splited[3]} != "edge" ]]
                then 
                echo $next_name
              sbatch -p cn-long  -A zeminz_g1 --qos=zeminzcnl -J ydc -c 20 -o ./log/$next_name.Mutect2.out -e ./log/$next_name.Mutect2.err ./somatic_SNV_GATK.sh $normal_bam $next_bam $normal_bam_sample T
                elif [[ $next_bam_type != N ]]
                then
                    if [[ ${next_bam_splited[3]} != P* ]] 
                    then
                    bam_type="$next_bam_type-${next_bam_splited[3]}"
                    else
                    bam_type=$next_bam_type
                    fi
                    echo $next_name
                    echo $bam_type
               sbatch -p cn-long  -A zeminz_g1 --qos=zeminzcnl -J ydc -c 20 -o ./log/$next_name.Mutect2.out -e ./log/$next_name.Mutect2.err ./somatic_SNV_GATK.sh $normal_bam $next_bam $normal_bam_sample $bam_type
                fi
            fi  
        done
        else
            continue
        fi
    done


