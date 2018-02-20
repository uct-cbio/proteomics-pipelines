#!/usr/bin/env python

import Bio
from Bio import SeqIO
import sys

filt = []

seqs = list(SeqIO.parse(sys.argv[1],'fasta'))

minlen = int(sys.argv[2])
maxlen = int(sys.argv[3])

output = sys.argv[4]

for rec in seqs:
    s = str(rec.seq)
    l = len(s)

    if ((l >= minlen) and (l <= maxlen)):
        filt.append(rec)

SeqIO.write(filt, output, 'fasta')

