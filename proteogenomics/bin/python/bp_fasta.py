#!/usr/bin/env python

import pandas as pd
import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import defaultdict
import sys
import pickle

data = list(SeqIO.parse(sys.argv[1], 'fasta'))

sequences_dct=defaultdict(set)
for rec in data:
    s = str(rec.seq)
    i  = rec.id
    assert '|' not in i
    sequences_dct[s].add(i)
    
new_sequences=[]
count=0

for sequence in sequences_dct:
    id='denovo_sequence_{}'.format(count)
    descs= list(sequences_dct[sequence])
    desc = 'sc={};scans='.format(len(descs)) + '|'.join(descs)
    rec=SeqRecord(id=id, seq=Seq(sequence), description=desc)
    new_sequences.append(rec)
    count +=1

pickle.dump(sequences_dct, open(sys.argv[2] +'.seq2scans.p','wb'))
SeqIO.write(new_sequences, sys.argv[2], 'fasta')
