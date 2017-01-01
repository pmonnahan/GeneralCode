#!/usr/bin/env python3

import os
import argparse
import subprocess

parser = argparse.ArgumentParser(description='This script uses the analysis ready reads for variant calling via the GATK tool ' +
                                 'HaplotypeCaller. The cohort1 and cohort2 option should be used to separate ploidies. A seperate ' +
                                 'SLURM shell script is created for each input file and sent to the NBI SLURM cluster. If the output ' +
                                 'directory does not exist, it is created automatically.')


parser.add_argument('-bamdir', type=str, metavar='', default='realigned/', help='REQUIRED: Full path to directory with input realigned bam files')
parser.add_argument('-R', type=str, metavar='reference_path', default='/nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/alygenomes.fasta', help='REQUIRED: Full path to reference genome')
parser.add_argument('-o', type=str, metavar='HC_dir', default='HC/', help='Realtive path to output directory [HC/]')
parser.add_argument('-mem', type=str, metavar='memory', default='16000', help='Total memory for each job (Mb) [16000]')
parser.add_argument('-time', type=str, metavar='time', default='2-12:00', help='Time for job [2-12:00]')
parser.add_argument('-print', type=str, metavar='print', default='false', help='If changed to true then shell files are printed to screen and not launched [false]')
args = parser.parse_args()

in_bam_list = []
for file in os.listdir(args.bamdir):
    if file.endswith('.bam'):
        in_bam_list.append(file)
in_bam_list.sort()

if args.o.endswith("/") is False:
    args.o += "/"

if os.path.exists(args.o + 'LoFreq/') is False:
    os.mkdir(args.o + 'LoFreq/')
outdir = args.o + 'LoFreq/'

for bam in in_bam_list:
    bam_basename = bam.replace('.bam', '')

    sh_file = open(outdir + bam_basename + '.sh', 'w')
    # write slurm shell file
    sh_file.write('#!/bin/bash -e\n' +
                  '#SBATCH -J LoFreq.' + bam_basename + '\n' +
                  '#SBATCH -o ' + outdir + bam_basename + '.LoFreq.out\n' +
                  '#SBATCH -e ' + outdir + bam_basename + '.LoFreq.err\n' +
                  '#SBATCH -p nbi-long\n' +
                  '#SBATCH -t ' + args.time + '\n' +
                  '#SBATCH --mem=' + args.mem + '\n' +
                  'source sametools-1.3.1\n' +
                  'source htslib-1.3.1\n' +
                  '/nbi/software/testing/lofreq/2.1.2/x86_64/lofreq call --ref /nbi/Research-Groups/JIC/Levi-Yant/Patrick/Camara/ChirsutaRef/chi_v1.fa --out ' + outdir + bam_basename[:-5] + '.Lofreq.vcf ' + outdir + bam_basename + '.bam')
    sh_file.close()

    # check if slurm shell file should be printed or sent to NBI SLURM
    if args.print == 'false':
        # send slurm job to NBI SLURM cluster
        cmd = ('sbatch ' + outdir + bam_basename + '.sh')
        p = subprocess.Popen(cmd, shell=True)
        sts = os.waitpid(p.pid, 0)[1]
    else:
        file = open(outdir + bam_basename + '.sh', 'r')
        data = file.read()
        print(data)
    os.remove(outdir + bam_basename + '.sh')
