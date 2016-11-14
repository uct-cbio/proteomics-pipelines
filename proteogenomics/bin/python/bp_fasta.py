#!/usr/bin/env python

import pandas as pd
import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import defaultdict
import sys

data = list(SeqIO.parse(sys.argv[1], 'fasta'))

sequences_dct=defaultdict(list)
for rec in data:
    s = str(rec.seq)
    i  = rec.id
    assert '|' not in i
    sequences_dct[s].append(i)
    
new_sequences=[]
count=0

for sequence in sequences_dct:
    id='denovo_sequence_{}'.format(count)
    descs= list(set(sequences_dct[sequence]))
    desc = 'sc={};scans='.format(len(descs)) + '|'.join(descs)
    rec=SeqRecord(id=id, seq=Seq(sequence), description=desc)
    new_sequences.append(rec)
    count +=1

SeqIO.write(new_sequences, sys.argv[2], 'fasta')
