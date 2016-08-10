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

loader = importlib.machinery.SourceFileLoader('config', sys.argv[1])
config = loader.load_module()
output = sys.argv[2]

fasta = list(SeqIO.parse(config.six_frame_translated,'fasta'))
fdict = {}
for rec in fasta:
    id = rec.id.split('|')[1]
    fdict[id] = rec

path = output + '/fasta'

try:
    shutil.rmtree(path)
    os.mkdir(path)
except:
    os.mkdir(path)

pg = pd.read_csv(config.mq_txt +'proteinGroups.txt',sep='\t')
pg =pg[(pg['Potential contaminant'] != '+') & (pg['Reverse'] != '+')]

ids = pg['Protein IDs'].tolist()

datum = [i.split(';') for i in ids]
new_ids = [item for sublist in datum for item in sublist]

orfs = []

for i in new_ids:
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







