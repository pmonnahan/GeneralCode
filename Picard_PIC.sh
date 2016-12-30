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
picard AddOrReplaceReadGroups I=/nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/PIC1.sort.bam O=/nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/PIC1.sort.RG.bam RGPL=illumina RGSM=PIC1 RGLB=PIC1 RGPU=PIC1
picard AddOrReplaceReadGroups I=/nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/PIC2.sort.bam O=/nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/PIC2.sort.RG.bam RGPL=illumina RGSM=PIC2 RGLB=PIC2 RGPU=PIC2
picard AddOrReplaceReadGroups I=/nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/PIC3.sort.bam O=/nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/PIC3.sort.RG.bam RGPL=illumina RGSM=PIC3 RGLB=PIC3 RGPU=PIC3