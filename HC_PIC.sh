#!/bin/bash -e
#SBATCH -J HC.PIC
#SBATCH -o /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/HC.PIC.out
#SBATCH -e /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/HC.PIC.err
#SBATCH -p nbi-long
#SBATCH -c 4
#SBATCH -t 14-00:00
#SBATCH --mem=100000
source jre-1.8.0_45
source GATK-3.6.0

srun java -XX:ParallelGCThreads=2 -Xmx100g -jar /nbi/software/testing/GATK/3.6.0/src/GenomeAnalysisTK.jar -T HaplotypeCaller -I /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/PIC1.sort.RG.bam --min_base_quality_score 10 --min_mapping_quality_score 20 -rf DuplicateRead -rf BadMate -rf BadCigar -ERC BP_RESOLUTION -R /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/ChirsutaRef/chi_v1.fa -o /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/PIC1_HC_g.vcf.gz -ploidy 104 --pcr_indel_model NONE -nct 4 --max_alternate_alleles 2
srun java -XX:ParallelGCThreads=2 -Xmx100g -jar /nbi/software/testing/GATK/3.6.0/src/GenomeAnalysisTK.jar -T HaplotypeCaller -I /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/PIC2.sort.RG.bam --min_base_quality_score 10 --min_mapping_quality_score 20 -rf DuplicateRead -rf BadMate -rf BadCigar -ERC BP_RESOLUTION -R /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/ChirsutaRef/chi_v1.fa -o /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/PIC2_HC_g.vcf.gz -ploidy 104 --pcr_indel_model NONE -nct 4 --max_alternate_alleles 2
srun java -XX:ParallelGCThreads=2 -Xmx100g -jar /nbi/software/testing/GATK/3.6.0/src/GenomeAnalysisTK.jar -T HaplotypeCaller -I /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/PIC3.sort.RG.bam --min_base_quality_score 10 --min_mapping_quality_score 20 -rf DuplicateRead -rf BadMate -rf BadCigar -ERC BP_RESOLUTION -R /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/ChirsutaRef/chi_v1.fa -o /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/PIC3_HC_g.vcf.gz -ploidy 104 --pcr_indel_model NONE -nct 4 --max_alternate_alleles 2