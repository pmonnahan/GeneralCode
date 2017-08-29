#!/bin/bash
#SBATCH -J DipPhase.scaf7
#SBATCH -e DipPhase.scaf7.err
#SBATCH -p nbi-medium
#SBATCH -n 10
#SBATCH -t 2-00:00
#SBATCH --mem=125000

source beagle-4.1
java -Xmx125g -jar $BEAGLE gt=/nbi/Research-Groups/JIC/Levi-Yant/300/ploidySplit/dips-withCRO/DM_HM_BP_BI_scaf_7_PASS.Dips.DP8.vcf nthreads=10 impute=false out=/nbi/Research-Groups/JIC/Levi-Yant/300/ploidySplit/dips-withCRO/Phased/DM_HM_BP_BI_scaf_7_PASS.Dips.BeagPhased
