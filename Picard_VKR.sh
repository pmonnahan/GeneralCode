#!/bin/bash -e
#SBATCH -J HC.PIC
#SBATCH -o /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/HC.LUZ.out
#SBATCH -e /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/HC.LUZ.err
#SBATCH -p nbi-long
#SBATCH -c 4
#SBATCH -t 14-00:00
#SBATCH --mem=100000
source jre-1.8.0_45
source GATK-3.6.0
source picard-1.134
picard AddOrReplaceReadGroups I=/nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/VKR1.sort.bam O=/nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/VKR1.sort.RG.bam RGPL=illumina RGSM=VKR1 RGLB=VKR1 RGPU=VKR1
picard AddOrReplaceReadGroups I=/nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/VKR2.sort.bam O=/nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/VKR2.sort.RG.bam RGPL=illumina RGSM=VKR2 RGLB=VKR2 RGPU=VKR2
picard AddOrReplaceReadGroups I=/nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/VKR3.sort.bam O=/nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/VKR3.sort.RG.bam RGPL=illumina RGSM=VKR3 RGLB=VKR3 RGPU=VKR3