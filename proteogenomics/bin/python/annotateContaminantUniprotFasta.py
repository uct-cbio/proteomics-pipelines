#!/usr/bin/env python

import Bio
from Bio import SeqIO
import sys

fasta = list(SeqIO.parse(sys.argv[1], 'fasta'))
if len(sys.argv) ==2:
    out = sys.stdout
else:
    out = sys.argv[2]

for seq in fasta:
   seq.description = "CONTAMINANT " + seq.description 
SeqIO.write(fasta, out, 'fasta')
