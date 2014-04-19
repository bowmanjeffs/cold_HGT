# -*- coding: utf-8 -*-
"""
Created on Mon Dec 09 07:43:14 2013

@author: Jeff
"""

import os
import re
import gzip
import subprocess

data = {}

with open('select_genomes_large_gc.txt', 'r') as file_in:
    for line in file_in:
        line = line.rstrip()
        line = line.split()
        data[line[2]] = line[0], line[1], line[3]
        

with gzip.open('master_all_out_clean_gc1.txt.gz', 'wb') as gc1_out, gzip.open('master_all_out_clean_gc2.txt.gz', 'wb') as gc2_out:
    for d in os.listdir('/Volumes/deming/cold_HGT_rd4/gc_amel'):
        if os.path.isdir('/Volumes/deming/cold_HGT_rd4/gc_amel/'+d):
            for f in os.listdir('/Volumes/deming/cold_HGT_rd4/gc_amel/'+d):
                if f == 'Rplots.pdf':
                    subprocess.call('rm /Volumes/deming/cold_HGT_rd4/gc_amel/'+d+'/'+f, shell = True)
                elif f.endswith('tminmax.csv'):
                    with open('/Volumes/deming/cold_HGT_rd4/gc_amel/'+d+'/'+f, 'r') as csv_in:
                        for line in csv_in:
                            if "V1" not in line:
                                line = line.rstrip()
                                line = re.sub('\"', '', line)
                                line = line.split(',')
                                cds = line[1]
                                cds = re.sub('.fasta', '', cds)
                                gc_1 = line[3]
                                gc_2 = line[4]
                                temp = data[cds]
                                print >> gc1_out, cds+'\t'+temp[0]+'\t'+temp[1]+'\t'+temp[2]+'\t'+str(gc_1)
                                print >> gc2_out, cds+'\t'+temp[0]+'\t'+temp[1]+'\t'+temp[2]+'\t'+str(gc_2)
                                print cds, temp[0], temp[1], temp[2], gc_1, gc_2