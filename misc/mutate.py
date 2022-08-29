#! /usr/bin/env python3

import sys
import random

fasta_file = sys.argv[1]
prob = float(sys.argv[2])

def mut_seq(seq):

    new_seq = ''
    for i in range(len(seq)):
        if seq[i] not in ['A', 'C', 'G', 'T']: 
            new_seq = new_seq + seq[i]
            continue

        if random.random() < prob:
            trand = random.random()
            tbase = ['A', 'C', 'G', 'T']
            tbase.remove(seq[i])

            if trand > (2.0 / 3):
                new_seq = new_seq + tbase[2]
            elif trand > (1.0 / 3):
                new_seq = new_seq + tbase[1]
            else:
                new_seq = new_seq + tbase[0]
        else:
            new_seq = new_seq + seq[i]

    return(new_seq)


temp_rid = None
temp_rseq = None        
with open(fasta_file, 'r') as hin:
    for line in hin:
        line = line.rstrip('\n')
        if line.startswith('>'):
            if temp_rid is not None:
                temp_rseq = mut_seq(temp_rseq)
                print(f'>{temp_rid}\n{temp_rseq}')

            temp_rid = line.lstrip('>')
            temp_rseq = ''
        else:
            temp_rseq = temp_rseq + line

    if temp_rid is not None:
        temp_rseq = mut_seq(temp_rseq)
        print(f'>{temp_rid}\n{temp_rseq}')

