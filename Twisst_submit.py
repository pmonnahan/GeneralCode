#!/usr/bin/env python35
import os
import subprocess

aa = ["VEL", "DRA", "VID"], ["VEL", "LAC", "VID"], ["VEL", "TZI", "VID"], ["HOC", "DRA", "VID"], ["HOC", "LAC", "VID"], ["HOC", "TZI", "VID"], ["SPI", "DRA", "VID"], ["SPI", "LAC", "VID"], ["SPI", "TZI", "VID"], ["TKO", "DRA", "VID"], ["TKO", "LAC", "VID"], ["TKO", "TZI", "VID"], ["TRE", "DRA", "VID"], ["TRE", "LAC", "VID"], ["TRE", "TZI", "VID"], [""], ["VEL", "KOW", "MIE"], ["VEL", "STE", "MIE"], ["VEL", "TBG", "MIE"], ["HOC", "KOW", "MIE"], ["HOC", "STE", "MIE"], ["HOC", "TBG", "MIE"], ["SPI", "KOW", "MIE"], ["SPI", "STE", "MIE"], ["SPI", "TBG", "MIE"], ["TKO", "KOW", "MIE"], ["TKO", "STE", "MIE"], ["TKO", "TBG", "MIE"], ["TRE", "KOW", "MIE"], ["TRE", "STE", "MIE"], ["TRE", "TBG", "MIE"]
joblist = []
for pops in aa:

    name = ".".join(pops)


    shfile1 = open('Twisst.' + name + '.sh', 'w')
    shfile1.write('#!/bin/bash\n' +
                  '#SBATCH -J Twisst.' + name + '.sh' + '\n' +
                  '#SBATCH -e /nbi/group-data/JIC/Research-Groups/Levi-Yant/Patrick/OandE/Twisst' + name + '.err' + '\n' +
                  '#SBATCH -o /nbi/group-data/JIC/Research-Groups/Levi-Yant/Patrick/OandE/Twisst' + name + '.out' + '\n' +
                  '#SBATCH -p nbi-medium\n' +
                  '#SBATCH -n 5\n' +
                  '#SBATCH -t 2-00:00\n' +
                  '#SBATCH --mem=16000\n' +
                  'source python-2.7.12\n' +
                  'cd /nbi/Research-Groups/JIC/Levi-Yant/Patrick/300/Twisst/twisst/\n' +
                  'export PATH=~/anaconda_ete/bin:$PATH\n' +
                  'for i in {1..8}\n' +
                  'do\n' +
                  'python run_twisst_parallel.py --threads 5 -t /nbi/Research-Groups/JIC/Levi-Yant/Patrick/300/Twisst/DM_HM_BP_BI.4dg.scf$i.DP8MIN230MAC2.LDphase.phyml_bionj.trees.gz -w /nbi/Research-Groups/JIC/Levi-Yant/Patrick/300/Twisst/DM_HM_BP_BI.4dg.scf$i.DP8MIN230MAC2.LDphase.phyml_bionj.' + name + '.weights.tsv --outputTopos /nbi/Research-Groups/JIC/Levi-Yant/Patrick/300/Twisst/DM_HM_BP_BI.4dg.scf$i.DP8MIN230MAC2.LDphase.phyml_bionj.' + name + '.topos.tsv -g ' + aa[0] + ' -g ' + aa[1] + ' -g ' + aa[2] + ' -g BEL --groupsFile /nbi/Research-Groups/JIC/Levi-Yant/Patrick/300/Twisst/Twisst.PopLabels.High.Cov.Indx.txt --method complete --abortCutoff 2000 --backupMethod threshold --thresholdTable /nbi/Research-Groups/JIC/Levi-Yant/Patrick/300/Twisst/twisst/threshold_tables/binom_threshold_wilson_CI95_0.05.tsv\n' +
                  'done\n')
    shfile1.close()

    cmd1 = ('sbatch Twisst.' + name + '.sh')
    p1 = subprocess.Popen(cmd1, shell=True)
    sts1 = os.waitpid(p1.pid, 0)[1]
    joblist.append(p1.pid)

    os.remove('Twisst.' + name + '.sh')
