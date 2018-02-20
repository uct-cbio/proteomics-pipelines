#!/usr/bin/env python

import Bio; from Bio import SeqIO
import sys
import collections
from collections import defaultdict
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

x = list(SeqIO.parse(sys.stdin, 'fasta'))

fasta_header_string = sys.argv[1]

f_dct = defaultdict(list)

for i in x:
    seq = str(i.seq)
    description = '|'.join(i.description.split())
    f_dct[seq].append(description)

count = 1

new_list = []
for i in f_dct:
    seq = i
    description = ' '.join(f_dct[i])
    id ='_'.join(fasta_header_string.split()) +'_'+str(count)
    rec = SeqRecord(seq = Seq(seq), description = description, id=id)
    new_list.append(rec)
    count += 1

SeqIO.write(new_list, sys.stdout, 'fasta')

