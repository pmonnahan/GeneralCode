#!/bin/bash -e
#SBATCH -J HC.LUZ
#SBATCH -o /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/HC.LUZ.out
#SBATCH -e /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/HC.LUZ.err
#SBATCH -p nbi-long
#SBATCH -c 4
#SBATCH -t 14-00:00
#SBATCH --mem=100000
source jre-1.8.0_45
source GATK-3.6.0
source samtools-1.3

samtools index /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/LUZ1.sort.RG.bam
samtools idxstats /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/LUZ1.sort.RG.bam > /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/LUZ1.sort.RG.idxstats
samtools flagstat /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/LUZ1.sort.RG.bam > /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/LUZ1.sort.RG.flagstat
samtools index /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/LUZ2.sort.RG.bam
samtools idxstats /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/LUZ2.sort.RG.bam > /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/LUZ2.sort.RG.idxstats
samtools flagstat /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/LUZ2.sort.RG.bam > /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/LUZ2.sort.RG.flagstat
samtools index /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/LUZ3.sort.RG.bam
samtools idxstats /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/LUZ3.sort.RG.bam > /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/LUZ3.sort.RG.idxstats
samtools flagstat /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/LUZ3.sort.RG.bam > /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/LUZ3.sort.RG.flagstat
srun java -XX:ParallelGCThreads=2 -Xmx100g -jar /nbi/software/testing/GATK/3.6.0/src/GenomeAnalysisTK.jar -T HaplotypeCaller -I /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/LUZ1.sort.RG.bam --min_base_quality_score 10 --min_mapping_quality_score 20 -rf DuplicateRead -rf BadMate -rf BadCigar -ERC BP_RESOLUTION -R /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/ChirsutaRef/chi_v1.fa -o /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/LUZ1_HC_g.vcf.gz -ploidy 50 --pcr_indel_model NONE -nct 4
srun java -XX:ParallelGCThreads=2 -Xmx100g -jar /nbi/software/testing/GATK/3.6.0/src/GenomeAnalysisTK.jar -T HaplotypeCaller -I /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/LUZ2.sort.RG.bam --min_base_quality_score 10 --min_mapping_quality_score 20 -rf DuplicateRead -rf BadMate -rf BadCigar -ERC BP_RESOLUTION -R /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/ChirsutaRef/chi_v1.fa -o /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/LUZ2_HC_g.vcf.gz -ploidy 50 --pcr_indel_model NONE -nct 4
srun java -XX:ParallelGCThreads=2 -Xmx100g -jar /nbi/software/testing/GATK/3.6.0/src/GenomeAnalysisTK.jar -T HaplotypeCaller -I /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/LUZ3.sort.RG.bam --min_base_quality_score 10 --min_mapping_quality_score 20 -rf DuplicateRead -rf BadMate -rf BadCigar -ERC BP_RESOLUTION -R /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/ChirsutaRef/chi_v1.fa -o /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/LUZ3_HC_g.vcf.gz -ploidy 50 --pcr_indel_model NONE -nct 4