#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess

# create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser()

parser.add_argument('-bamdir', type=str, metavar='fastq_ca', help='REQUIRED: Full path to directory with input bam files')
parser.add_argument('-o', type=str, metavar='output_directory', help='full path to output directory [aligned/]')
parser.add_argument('-R', type=str, metavar='reference_path', required=True, help='REQUIRED: Full path to reference genome')
parser.add_argument('-C', type=str, metavar='downgrade_mapping_quality', required=True, help='coefficient for downgrading mapping quality scores due to excessive mismatches; 0 disables functionality')
parser.add_argument('-c', type=int, metavar='num_cores', default='1', help='Number of requested cores for job [5]')
parser.add_argument('-mem', type=str, metavar='memory', default='16000', help='Total memory for each job (Mb) [16000]')
parser.add_argument('-time', type=str, metavar='time', default='0-02:00', help='Time for job [2-00:00]')
parser.add_argument('-print', type=str, metavar='print', default='false', help='If changed to true then shell files are printed to screen and not launched [false]')
args = parser.parse_args()


# gather list of input fastq files
bam_list = []
for file in os.listdir(args.bamdir):
    if file.endswith('.bam'):
        bam_list.append(file)
bam_list.sort()

# check if output directory exists, create it if necessary
if os.path.exists(args.o) is False:
    os.mkdir(args.o)


# loop through input fastq files
count = 0
for bam in bam_list:
    bamname = bam.split(".")[0]
    for chrom in range(1, 9):  # Loop over chromosomes
        # write slurm shell file
        sh_file = open(args.o + bamname + '.sh', 'w')
        sh_file.write('#!/bin/bash -e\n' +
                      '#SBATCH -J pileup.' + bamname + '\n' +
                      '#SBATCH -o ' + args.o + bamname + '.out\n' +
                      '#SBATCH -e ' + args.o + bamname + '.err\n' +
                      '#SBATCH -p nbi-long\n' +
                      '#SBATCH -c ' + str(args.c) + '\n' +
                      '#SBATCH -t ' + args.time + '\n' +
                      '#SBATCH --mem=' + args.mem + '\n' +
                      'source samtools-1.3\n' +
                      'samtools mpileup -C ' + args.C + ' -f ' + args.R + ' -r Chr' + str(chrom) + ' -o ' + args.o + bamname + ".C" + args.C + ".Chr" + str(chrom) + ".pileup" + args.bamdir + bam)
        sh_file.close()


        # check if slurm shell file should be printed or sent to NBI SLURM
        if args.print == 'false':
            # send slurm job to NBI SLURM cluster
            cmd = ('sbatch ' + args.o + bamname + '.sh')
            p = subprocess.Popen(cmd, shell=True)
            sts = os.waitpid(p.pid, 0)[1]
        else:
            file = open(args.o + bamname + '.sh', 'r')
            data = file.read()
            print(data)
        count += 1
        os.remove(args.o + bamname + '.sh')

# if appropriate, report how many slurm shell files were sent to NBI SLURM
if args.print == 'false':
    print('\nSent ' + str(count) + ' jobs to the NBI SLURM  cluster\n\nFinished!!\n\n')
