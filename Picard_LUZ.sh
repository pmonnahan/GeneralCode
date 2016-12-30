#!/bin/bash -e
#SBATCH -J HC.LUZ
#SBATCH -o /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/HC.LUZ.out
#SBATCH -e /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/HC.LUZ.err
#SBATCH -p nbi-long
#SBATCH -c 4
#SBATCH -t 14-00:00
#SBATCH --mem=100000
source picard-1.134
picard AddOrReplaceReadGroups I=/nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/LUZ1.sort.bam O=/nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/LUZ1.sort.RG.bam RGPL=illumina RGSM=LUZ1 RGLB=LUZ1 RGPU=LUZ1
picard AddOrReplaceReadGroups I=/nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/LUZ2.sort.bam O=/nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/LUZ2.sort.RG.bam RGPL=illumina RGSM=LUZ2 RGLB=LUZ2 RGPU=LUZ2
picard AddOrReplaceReadGroups I=/nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/LUZ3.sort.bam O=/nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/LUZ3.sort.RG.bam RGPL=illumina RGSM=LUZ3 RGLB=LUZ3 RGPU=LUZ3