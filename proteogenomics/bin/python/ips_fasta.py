#!/usr/bin/env python

import numpy as np
import pandas as pd
import sys
import importlib.machinery
import Bio; from Bio import SeqIO
import os
import shutil
from collections import defaultdict
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Seq import translate
import json

#loader = importlib.machinery.SourceFileLoader('config', sys.argv[1])
#config = loader.load_module()
input = sys.argv[1]
output = sys.argv[2]

fasta = list(SeqIO.parse(input,'fasta'))

fdict = {}

for rec in fasta:
    id = rec.id.split('|')[1]
    fdict[id] = rec


path = output
orfs = []

for i in fdict:
    try:
        d = int(i.split('.')[-1]) 
        if d == 0:
            orfs.append(fdict[i])
    except:
        d = None

orf_sequence = defaultdict(list)

for rec in orfs:
    orf_sequence[str(rec.seq)].append(rec.id)


export = []
mapping = {}

count=1
for seq in orf_sequence:
    id = 'orf_sequence_{}'.format(str(count))
    ids = orf_sequence[seq]
    
    for i in ids:
        mapping[i] = id

    description =';'.join(ids)
    seq = Seq(''.join(seq.split('*')))
    record = SeqRecord(id = id, description = description, seq = seq)
    export.append(record)
    count += 1

jstr = json.dumps(mapping)

w = open(path +'/id_mapping.json','w')
w.write(jstr)
w.close()

SeqIO.write(export, path +'/nr_translated_pg_orfs.fasta','fasta')







