# -*- coding: utf-8 -*-
"""
Created on Sun Jul 28 17:26:33 2013

@author: Jeff

v2 adds normalization according to Hao, 2004
"""

import os
import subprocess
import re
import gzip
import cPickle
from Bio import SeqIO
from joblib import Parallel, delayed

## desired nmer length

k = 5
k1 = k - 1
k2 = k - 2

fasta_names = set()
groups = {}
    
#### evaluate for 16S rRNA ####
    
names = set()
with open('combined_16S.fasta', 'w') as fasta_out:
    for f in os.listdir('/Volumes/deming/cold_HGT_rd4/select_fna'):
        if f.endswith('.combined.fna'):
            name = re.sub('.fna', '', f)
            if name not in names:
 
                ## find 16S rRNA gene
                print 'finding 16S genes'
                    
                print name
                blast = subprocess.Popen('blastn ' \
                '-task megablast ' \
                '-num_threads 8 ' \
                '-max_target_seqs 1 ' \
                '-evalue 1e-50 ' \
                '-db /volumes/deming/databases/16SMicrobial ' \
                '-outfmt 5 ' \
                '-out '+name+'_16S.xml ' \
                '-query /Volumes/deming/cold_HGT_rd4/select_fna/'+f+';' \
                'python parse_blast_xml_v5.py '+name+'_16S.xml', shell = True)
                blast.communicate()
                
                with open(name+'_16S.fasta', 'r') as fasta_in:
                    s = ''
                    for line in fasta_in:
                        line = line.rstrip('\n')
                        if line.startswith('>') == False:
                            s = s + line
                    print >> fasta_out, '>'+name
                    print >> fasta_out, s
                       
            names.add(name)
                                
mothur_commands = 'mothur "#align.seqs(candidate=combined_16S.fasta, flip=t, template=/volumes/deming/databases/core_set_aligned.imputed.fasta);' \
'screen.seqs(start=3000,end=4000);' \
'filter.seqs(trump=.,vertical=T);' \
'dist.seqs(output=square)"'

mothur = subprocess.Popen(mothur_commands, shell = True)
mothur.communicate()       

#### generate compositional vector ####

## if you haven't already done so generate all possible five letter words

import itertools
if str(k)+'mers_pro.set' not in os.listdir('.'):
    pro = ['A','R','N','D','C','E','Q','G','H','I',\
    'L','K','M','F','P','S','T','W','Y','V'] 
    
    ## k      
    
    pros = list(itertools.repeat(pro, k))
    bins = list(itertools.product(*pros))
    new_bins = []    
    for b in bins:
        b = ''.join(b)
        new_bins.append(b)
    bins = new_bins
    nmers = open(str(k)+'mers_pro.set', 'wb')                   
    cPickle.dump(bins, nmers)
    nmers.close()
    itertools.product()
    
    ## k1
    
    pros = list(itertools.repeat(pro, k1))
    k1_bins = list(itertools.product(*pros))
    k1_new_bins = []    
    for b in k1_bins:
        b = ''.join(b)
        k1_new_bins.append(b)
    k1_bins = k1_new_bins
    k1_nmers = open(str(k1)+'mers_pro.set', 'wb')                   
    cPickle.dump(k1_bins, k1_nmers)
    k1_nmers.close()
    itertools.product()    
    
    ## k2
    
    pros = list(itertools.repeat(pro, k2))
    k2_bins = list(itertools.product(*pros))
    k2_new_bins = []    
    for b in k2_bins:
        b = ''.join(b)
        k2_new_bins.append(b)
    k2_bins = k2_new_bins
    k2_nmers = open(str(k2)+'mers_pro.set', 'wb')                   
    cPickle.dump(k2_bins, k2_nmers)
    k2_nmers.close()
    itertools.product()     
    
else:
    nmers = open(str(k)+'mers_pro.set', 'rb')
    bins = cPickle.load(nmers)
    nmers.close()
    
    k1_nmers = open(str(k1)+'mers_pro.set', 'rb')
    k1_bins = cPickle.load(k1_nmers)
    k1_nmers.close()
    
    k2_nmers = open(str(k2)+'mers_pro.set', 'rb')
    k2_bins = cPickle.load(k2_nmers)
    k2_nmers.close()
                
#### search proteome for all n letter kmers, and tally occurrence
    
def calc_vector(name, bins, k1_bins, k2_bins):
    k1_found_bins = {}
    k1_used_bins = set()
    k2_found_bins = {}
    k2_used_bins = set()
    found_bins = {}
    used_bins = set()
    
    seqs = name+'.pro.fasta'
    for record in SeqIO.parse('/Volumes/deming/cold_HGT_rd4/select_fna_cds/'+seqs, 'fasta'):
        query = str(record.seq)
        
        ## k1 and k2
        
        for i in range(0,len(query)):
            kmer = query[i:i+k1]
            print name, 'k1', i
            if kmer not in k1_used_bins:
                k1_found_bins[kmer] = 1
                k1_used_bins.add(kmer)
            else:
                k1_found_bins[kmer] = k1_found_bins[kmer] + 1  
            
        for i in range(0,len(query)):
            kmer = query[i:i+k2]
            print name, 'k2', i
            if kmer not in k2_used_bins:
                k2_found_bins[kmer] = 1
                k2_used_bins.add(kmer)
            else:
                k2_found_bins[kmer] = k2_found_bins[kmer] + 1
                
        ## k
            
        for i in range(0,len(query)):
            kmer = query[i:i+k]
            print name, 'k', i
            if kmer not in used_bins:
                found_bins[kmer] = 1
                used_bins.add(kmer)
            else:
                found_bins[kmer] = found_bins[kmer] + 1
                
    ## k0
        
    norm_bins = {}
    for kmer in found_bins.keys():
        if len(kmer) == k:
            kmer_1 = kmer[0:-1]
            kmer_2 = kmer[1:]
            kmer_3 = kmer[1:-1]
            bigL = len(query)
            kmer_0 = ((k1_found_bins[kmer_1] * k1_found_bins[kmer_2])
            / float(k2_found_bins[kmer_3])) * (((bigL - k + 1) * (bigL - k + 3))
            / float((bigL - k + 2) ** 2))
            kmer_norm = (found_bins[kmer] - kmer_0) / kmer_0
            norm_bins[kmer] = kmer_norm
            print name, kmer, kmer_norm
    
    ## fill out dictionary with 0 values for unrepresented kmers
            
    for nmer in bins:
        if nmer not in used_bins:
            norm_bins[nmer] = 0
        
    with gzip.open(name+'_'+str(k)+'mer_bins.txt.gz', 'wb') as bins_out:
        for each in sorted(norm_bins.keys()):
            print >> bins_out, each+'\t'+str(norm_bins[each])
    
    subprocess.call('echo '+name+' >> log.txt', shell = True)
    
names = sorted(names)
subprocess.call('rm log.txt', shell = True)    
Parallel(n_jobs = -1, verbose = 5)(delayed(calc_vector)
(name, bins, k1_bins, k2_bins) for name in names)

final_bins = {}
for f in os.listdir('.'):
    if f.endswith('bins.txt.gz'):
        name = re.split('_'+str(k)+'mer_bins', f)
        name = name[0]
        print name
        with gzip.open(f, 'rb') as bin_file:
            temp = []
            for line in bin_file:
                line = line.rstrip('\n')
                line = line.split('\t')
                bin_id = line[0]
                bin_num = line[1]
                temp.append(bin_num)
            final_bins[name] = temp
            
print 'writing vector matrix'
        
with gzip.open(str(k)+'mer_normalized_phylogeny_vector_output.txt.gz', 'wb') as output:
    for name in sorted(final_bins.keys()):
        if name != sorted(final_bins.keys())[-1]:
            print >> output, name+'\t',
        else:
            print >> output, name+'\t'    
    for i in range(0,20 ** k):
        print i
        for name in sorted(final_bins.keys()):
            if name != sorted(final_bins.keys())[-1]:
                print >> output, str(final_bins[name][i])+'\t',
            else:
                print >> output, str(final_bins[name][i])+'\t'
                                    
cPickle.dump(final_bins, open(str(k)+'mer_normalized_phylogeny_vector_output.p', 'wb'))

#subprocess.call('R --no-save --'+str(k)+'mer_normalized_phylogeny_vector_output.txt.gz < compositional_vector_3.r', shell = True)