#!/usr/bin/env python3

import os, sys, argparse, subprocess

#create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser()

parser.add_argument('-bamdir', type=str, metavar='fastq_ca', default='fastq_ca/', help='REQUIRED: Full path to directory with input fastq files')
parser.add_argument('-o', type=str, metavar='aligned_dir', default='aligned/', help='Realtive path to output directory [aligned/]')
parser.add_argument('-R', type=str, metavar='reference_path', required=True, help='REQUIRED: Full path to reference genome')
parser.add_argument('-c', type=int, metavar='num_cores', default='1', help='Number of requested cores for job [5]')
parser.add_argument('-mem', type=str, metavar='memory', default='16000', help='Total memory for each job (Mb) [16000]')
parser.add_argument('-time', type=str, metavar='time', default='0-02:00', help='Time for job [2-00:00]')
parser.add_argument('-print', type=str, metavar='print', default='false', help='If changed to true then shell files are printed to screen and not launched [false]')
args = parser.parse_args()


# check if output directory exists, create it if necessary
if os.path.exists(args.o) == False:
    os.mkdir(args.o)

Pops = ["VKR","CEZ","LUZ","PIC"]

#loop through input fastq files
count = 0
for pop in Pops:
    
    #write slurm shell file
    sh_file = open(args.o+pop+'.sh','w')
    sh_file.write('#!/bin/bash -e\n'+
                  '#SBATCH -J ' + pop + '.pileup\n'+
                  '#SBATCH -o '+args.o+pop+'.pileup.out\n'+ 
                  '#SBATCH -e '+args.o+pop+'.pileup.err\n'+
                  '#SBATCH -p nbi-long\n'+
                  '#SBATCH -c '+str(args.c)+'\n'+
                  '#SBATCH -t '+args.time+'\n'
                  '#SBATCH --mem='+args.mem+'\n'+
                  'source samtools-1.3\n'+
                  'samtools mpileup ' + args.bamdir + pop + '1.sort.RG.bam ' + args.bamdir + pop + '2.sort.RG.bam ' + args.bamdir + pop + '3.sort.RG.bam ' + '-C ' + args.C + ' -o ' + args.o + pop + '.pileup')
    sh_file.close()

    #check if slurm shell file should be printed or sent to NBI SLURM 
    if args.print == 'false':
        #send slurm job to NBI SLURM cluster
        cmd = ('sbatch '+args.o+pop+'.sh')
        p = subprocess.Popen(cmd, shell=True)
        sts = os.waitpid(p.pid, 0)[1]
    else:
        file = open(args.o+pop+'.sh','r')
        data = file.read()
        print(data)
    count +=1
    os.remove(args.o+pop+'.sh')
    
#if appropriate, report how many slurm shell files were sent to NBI SLURM 
if args.print == 'false':
    print('\nSent '+str(count)+' jobs to the NBI SLURM  cluster\n\nFinished!!\n\n')
