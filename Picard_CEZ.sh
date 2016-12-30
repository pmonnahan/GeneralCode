#SBATCH -J HC.CEZ
#SBATCH -o /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/HC.CEZ.out
#SBATCH -e /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/HC.CEZ.err
#SBATCH -p nbi-long
#SBATCH -c 4
#SBATCH -t 14-00:00
#SBATCH --mem=100000
source picard-1.134
picard AddOrReplaceReadGroups I=/nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/CEZ1.sort.bam O=/nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/CEZ1.sort.RG.bam RGPL=illumina RGSM=CEZ1 RGLB=CEZ1 RGPU=CEZ1
picard AddOrReplaceReadGroups I=/nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/CEZ2.sort.bam O=/nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/CEZ2.sort.RG.bam RGPL=illumina RGSM=CEZ2 RGLB=CEZ2 RGPU=CEZ2
picard AddOrReplaceReadGroups I=/nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/CEZ3.sort.bam O=/nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/aligned/CEZ3.sort.RG.bam RGPL=illumina RGSM=CEZ3 RGLB=CEZ3 RGPU=CEZ3