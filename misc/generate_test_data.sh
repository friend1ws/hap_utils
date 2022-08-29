#! /bin/bash

reference=$1

seqtk subseq ${reference} name.lst > hap1_ref.fa

python3 mutate.py hap1_ref.fa 0.001 > hap2_ref.fa

wgsim hap1_ref.fa hap1_seq.fq /dev/null -d 50000 -N 10000 -1 15000 -2 15000 > /dev/null

wgsim hap2_ref.fa hap2_seq.fq /dev/null -d 50000 -N 10000 -1 15000 -2 15000 > /dev/null

