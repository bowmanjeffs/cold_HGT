# -*- coding: utf-8 -*-
"""
Created on Sat Nov 09 15:18:04 2013

@author: Jeff
"""

from Bio import SeqIO
from Bio.SeqUtils import GC
from math import sqrt
import gzip
import os

groups = {} ## key = fna, value = group

with open('select_genomes.final.groups', 'r') as group_file:
    for line in group_file:
        line = line.rstrip()
        line = line.split('\t')
        groups[line[0]] = line[1]
        

with gzip.open('select_genome_gc.txt.gz', 'w') as gc_out, open('select_genomes_large_gc.txt', 'w') as large_gc_out:
    for f in os.listdir('select_fna_cds'):
        if f.endswith('nuc.fasta'):
            print f
            all_seq = ''
            seq_GC = {}
            strain = f.rstrip('.combined.nuc.fasta')
            group = groups[strain+'.combined.fna']
            
            pfams = {} ## key = sequence name, value = pfam
            with open('/Volumes/deming/cold_HGT_rd3/select_fna_cds/'+strain+'.combined.pfam.txt', 'r') as pfam_file:
                for line in pfam_file:
                    if line.startswith('#') != True:
                        line = line.split()
                        pfams[line[2]] = line[0]
            
            for record in SeqIO.parse(open('select_fna_cds/'+f, 'r'), 'fasta'):
                seq = record.seq
                all_seq = all_seq + seq
                temp_GC = GC(seq)
                seq_GC[record.id] = temp_GC
            
            mean_GC = GC(all_seq)
            
            x_2 = []
            for key in seq_GC.keys():
                new_value = (seq_GC[key] - mean_GC) ** 2
                x_2.append(new_value)
            
            n = len(x_2)
            
            sd_GC = sqrt(sum(x_2)/(n - 1)) # equation 3.4 in Elementary Statistics
            
            keep = set()
            for key in seq_GC.keys():
                pfam = pfams[key]
                print >> gc_out, strain+'\t'+group+'\t'+key+'\t'+pfam+'\t'+str(seq_GC[key])+'\t'+str(mean_GC)+'\t'+str(sd_GC)
    
                if abs(seq_GC[key] - mean_GC) - sd_GC * 2 > 0:
                    print >> large_gc_out, strain+'\t'+group+'\t'+key+'\t'+pfam+'\t'+str(seq_GC[key])+'\t'+str(mean_GC)+'\t'+str(sd_GC)
                    keep.add(key)
                        
            with open('/Volumes/deming/cold_HGT_rd3/select_fna_cds/large_gc/'+strain+'.combined.nuc.largegc.fasta', 'w') as fasta_out:
                for record in SeqIO.parse(open('select_fna_cds/'+f, 'r'), 'fasta'):
                    if record.id in keep:
                        print >> fasta_out, '>'+record.id+'\n'+record.seq