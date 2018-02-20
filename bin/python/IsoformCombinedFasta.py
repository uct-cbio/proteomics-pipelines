#!/usr/bin/env python
import Bio
from Bio import SeqIO
import sys
from collections import defaultdict
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

proteins     = defaultdict(set)
descriptions = defaultdict(list)

db = SeqIO.parse(sys.argv[1],'fasta')

for rec in db:
    protein = str(rec.seq)
    id = rec.id.split('|')
    desc = ' '.join(rec.description.split()[1:])
    assert len(id) == 3
    id = id[1]
    proteins[protein].add(id)
    descriptions[protein].append(desc)

concatenated = []
for protein in proteins:
    ids = list(proteins[protein])
    new_id = '_'.join(ids)
    desc = '; '.join(descriptions[protein])
    rec = SeqRecord(seq = Seq(protein), id = 'compacted|{}|'.format(new_id), description = '')
    concatenated.append(rec)

SeqIO.write(concatenated, sys.argv[1].split('.fasta')[0] + '_compacted.fasta', 'fasta')

