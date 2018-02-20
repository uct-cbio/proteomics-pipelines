#!/usr/bin/env python

import sys
import Bio
from Bio import SeqIO
from collections import defaultdict
import pickle

groups = defaultdict(list)

inpt = list(SeqIO.parse(sys.argv[1],'fasta'))

chunknumber = int(sys.argv[2])
assert chunknumber >= 1
for rec in inpt:
    id = ';'.join(rec.id.split(';')[:2])
    groups[id].append(rec)


count = 1
grouped =defaultdict(list)

for key in groups:
    group = groups[key]
    grouped[count].append(group)
    if count == chunknumber:
        count= 1
    else:
        count += 1

for key in grouped:
    group = grouped[key]
    new_group = []
    for g in group:
        for gg in g:
            new_group.append(gg)
    SeqIO.write(new_group, sys.argv[3] + '/tags.{}.fasta'.format(str(key)),'fasta' )
    
