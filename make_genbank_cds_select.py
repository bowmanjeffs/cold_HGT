# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 16:23:08 2013

@author: Jeff
"""
get = ['Psychroflexus_gondwanensis_ACAM_44.combined.fna']

#with open('select_genomes.final.groups', 'r') as group_file:
#    for line in group_file:
#        line = line.split('\t')
#        group = line[0]
#        get.append(group)

def find_pros_with_trans(seqname, seq, trans_table, min_protein_length):
    answer = []
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            trans = str(nuc[frame:].translate(trans_table))
            untrans = str(nuc) # not frame, you only want to count DNA from the true start
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start) # starting from aa_start find first instance of * and return index
                if aa_end == -1:
                    aa_end = trans_len
                if aa_end-aa_start >= min_protein_length:
                    start = frame+aa_start*3 # counting forward from left side of reverse complement!
                    end = frame+aa_end*3 # +3   
                    print strand, frame, aa_start, aa_end, 'dna start = ', start, 'dna end = ', end
                    pro=trans[aa_start:aa_end]
                    print 'pro=',pro[:10],'...',pro[-10:]
                    dna=untrans[start:end]
                    dna_trans=str(nuc[start:end].translate(trans_table))
                    print 'dna_trans=',dna_trans[:10],'...',dna_trans[-10:]
                    answer.append((seqname, strand, frame, start, end, pro, dna))
                aa_start = aa_end+1
    answer.sort()
    return answer

import re
from Bio import SeqIO
import subprocess
    
for f in get:
    name = re.sub('.fna', '', f)
    
    with open('select_fna_cds/'+name+'.pro.fasta', 'w') as output_pro, open('select_fna_cds/'+name+'.nuc.fasta', 'w') as output_nuc:
        i = 0            
        for record in SeqIO.parse('select_fna/'+f, 'fasta'):
            i = i + 1                
            cds_list = find_pros_with_trans(record.id, record.seq, 11, 50)
            with open('temp'+'_'+name+'_pro.fasta', 'w') as temp_pro:
                for seqname, strand, frame, start, end, pro, dna in cds_list:
                    print >> temp_pro, '>'+seqname+'_'+str(start)+':'+str(end)+'_strand='+str(strand)
                    print >> temp_pro, pro
    
            hmmer = subprocess.Popen('hmmscan -E 1e-5 --tblout select_fna_cds/'+name+'_'+str(i)+'.pfam.txt /volumes/deming/databases/Pfam-A.hmm temp'+'_'+name+'_pro.fasta', shell=True)
            hmmer.communicate()
                    
            keep = set()
            with open('select_fna_cds/'+name+'_'+str(i)+'.pfam.txt', 'r') as pfam_file:
                for line in pfam_file:
                    if line.startswith('#') == False:
                        line = line.split()
                        seq = line[2]
                        keep.add(seq)
                    
            for seqname, strand, frame, start, end, pro, dna in cds_list:
                if seqname+'_'+str(start)+':'+str(end)+'_strand='+str(strand) in keep:
                    print >> output_pro, '>'+seqname+'_'+str(start)+':'+str(end)+'_strand='+str(strand)
                    print >> output_pro, pro
                    print >> output_nuc, '>'+seqname+'_'+str(start)+':'+str(end)+'_strand='+str(strand)
                    print >> output_nuc, dna
                        
            subprocess.call('rm temp'+'_'+name+'_pro.fasta', shell = True)
        subprocess.call('cat select_fna_cds/'+name+'_*.pfam.txt > select_fna_cds/'+name+'.pfam.txt', shell = True)
        subprocess.call('rm select_fna_cds/'+name+'_*.pfam.txt', shell = True)                    
