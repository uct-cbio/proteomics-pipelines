#!/usr/bin/env python

import sys
import os
import Bio; from Bio import SeqIO
import sequtils
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


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
alternative = sequtils.alt_starts_recs(seqs, starts=['ATG','GTG','TTG'])
SeqIO.write(alternative, outfile, 'fasta')
