#!/usr/bin/python

from Bio import SeqIO
import sys

assert len(sys.argv) == 3  #fasta file, string to append to description of each record
seqs = list(SeqIO.parse(sys.argv[1], 'fasta'))
prefix = sys.argv[2] # eg. 'Contaminant_'
for seq in seqs:
    id = seq.id
    seq.id = prefix + id
SeqIO.write(seqs, sys.argv[1], 'fasta')



