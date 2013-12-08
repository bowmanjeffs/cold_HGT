# -*- coding: utf-8 -*-
"""
Created on Fri Nov 08 10:45:40 2013

@author: Jeff

"""

import sys

reps = sys.argv[2]
file_in = sys.argv[3]

from Bio import SeqIO
from Bio.SeqUtils import GC
import random

if sys.argv[1] == '1':
    with open('test_gc.txt', 'w') as output:
        for record in SeqIO.parse(open(sys.argv[3], 'r'), 'fasta'):
            seq = record.seq
            seq_list = []
            for l in range(0,len(seq),3):
                seq_list.append(seq[l:l + 3])
            for i in range(0,int(reps)):
                boot = random.sample(seq_list, int(round(len(seq_list) * 0.5)))
                GC1 = ''
                GC2 = ''
                GC3 = ''
                for codon in boot:
                    GC1 = GC1 + codon[0]
                    GC2 = GC2 + codon[1]
                    GC3 = GC3 + codon[2]
                print str(GC(GC1))+'\t'+str(GC(GC2))+'\t'+str(GC(GC3))
                print >> output, GC(GC1) / 100, GC(GC2) / 100, GC(GC3) / 100

elif sys.argv[1] == '0': #analyze all cds togethre
    GC1 = ''
    GC2 = ''
    GC3 = ''
    with open('test_gc.txt', 'w') as output:
        for record in SeqIO.parse(open(sys.argv[3], 'r'), 'fasta'):
            seq = record.seq
            seq_list = []
            for l in range(0,len(seq),3):
                seq_list.append(seq[l:l + 3])
            for i in range(0,int(reps)):
                for codon in seq_list:
                    GC1 = GC1 + codon[0]
                    GC2 = GC2 + codon[1]
                    GC3 = GC3 + codon[2]

        print str(GC(GC1))+'\t'+str(GC(GC2))+'\t'+str(GC(GC3))
        print >> output, GC(GC1) / 100, GC(GC2) / 100, GC(GC3) / 100                