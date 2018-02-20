#!/usr/bin/env python

import sys
import Bio; from Bio import SeqIO
from collections import defaultdict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

db = list(SeqIO.parse(sys.argv[1],'fasta'))
tag = sys.argv[2]

seqdct = defaultdict(set)
for rec in db:
    seqdct[str(rec.seq)].add(str(rec.description))

count = 1
newrecs = []

for seq in list(seqdct.keys()):
    id = 'generic|{}{}|'.format(tag, str(count)) 
    description = '; '.join(seqdct[seq])
    newrec = SeqRecord(seq=Seq(seq), id=id, description = description)
    newrecs.append(newrec)
    count += 1

SeqIO.write( newrecs, sys.argv[3] , 'fasta')
