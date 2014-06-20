# -*- coding: utf-8 -*-
"""
Created on Fri Nov 08 10:45:40 2013

@author: Jeff

run as python [0 or 1] [reps] file

if 0: codons are collected from all seqs in input file and subsampled if
reps > 1.  If reps = 1 codons are not subsampled.

if 1: GC by position is output individually for all seqs in input file. All
codons are used and reps must = 1.

"""
import sys
#sys.argv = ['NA', '0', '1', 'Escherichia_coli_07798_uid66001.draft.combined.ffn']

reps = sys.argv[2]
file_in = sys.argv[3]

from Bio import SeqIO
from Bio.SeqUtils import GC
import random

if sys.argv[1] == '1': ## analyze each cds individually
    with open('test_gc.txt', 'w') as output:
        for record in SeqIO.parse(open(sys.argv[3], 'r'), 'fasta'):
            seq = record.seq
            seq_list = []
            for l in range(0,len(seq),3):
                seq_list.append(seq[l:l + 3])
                boot = seq_list
                GC1 = ''
                GC2 = ''
                GC3 = ''
                for codon in boot:
                    GC1 = GC1 + codon[0]
                    GC2 = GC2 + codon[1]
                    GC3 = GC3 + codon[2]
            print str(GC(GC1) / 100) + '\t' + str(GC(GC2) / 100) + '\t'+str(GC(GC3) / 100)
            print >> output, GC(GC1) / 100, GC(GC2) / 100, GC(GC3) / 100

elif sys.argv[1] == '0': ## analyze all cds together
    GC1 = ''
    GC2 = ''
    GC3 = ''
    with open('test_gc.txt', 'w') as output:
        seq_list = []
        for record in SeqIO.parse(open(sys.argv[3], 'r'), 'fasta'):
            seq = record.seq
            if len(seq) % 3 == 0: ## some cds are not divisible by three, can't use these
                for l in range(0,len(seq),3):
                    seq_list.append(seq[l:l + 3])
        for i in range(0,int(reps)):
            if reps > 1:
                boot = random.sample(seq_list, int(round(len(seq_list) * 0.5)))
            else:
                boot = seq_list
            GC1 = ''
            GC2 = ''
            GC3 = ''
            for codon in boot:
                GC1 = GC1 + codon[0]
                GC2 = GC2 + codon[1]
                GC3 = GC3 + codon[2]

            print str(GC(GC1) / 100) + '\t' + str(GC(GC2) / 100) + '\t'+str(GC(GC3) / 100)
            print >> output, GC(GC1) / 100, GC(GC2) / 100, GC(GC3) / 100                