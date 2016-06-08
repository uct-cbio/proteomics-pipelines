#!/home/ptgmat003/ve/bin/python

import sys
import os
import Bio; from Bio import SeqIO
#sys.path.append(os.getcwd())
import sequtils

if len(sys.argv)==3:
    infile  = sys.argv[1]
    outfile = sys.argv[2]

elif len(sys.argv) == 1:
    infile = sys.stdin
    outfile = sys.stdout

elif len(sys.argv)==2:
    infile = sys.argv[1]
    outfile = sys.stdout

seqs = list(SeqIO.parse(infile, 'fasta'))

sf = sequtils.sf_contigs(seqs)

SeqIO.write(sf, outfile, 'fasta')
