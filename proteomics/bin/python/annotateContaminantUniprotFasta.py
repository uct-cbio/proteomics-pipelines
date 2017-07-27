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
    id_ = seq.id.split('|')
    new_ids = []
    for id in id_[1:]:
        new_ids.append(id+'_CONTAMINANT')
    seq.id = id_[0] + '|' + '|'.join(new_ids)
    seq.description = "CONTAMINANT " + seq.description 
SeqIO.write(fasta, out, 'fasta')
