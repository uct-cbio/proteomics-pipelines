#!/usr/bin/python

import Bio; from Bio import SeqIO
import json
from Bio import SeqRecord
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys
import collections
from collections import defaultdict
import json

assert len(sys.argv) > 1
db_lst = sys.argv[1:]
dbs = []
map_set = set()
map_dct = {}
mapping = defaultdict(set)
for i in db_lst:
    db = list(SeqIO.parse(i,'fasta'))
    for j in db:
        x = str(j.seq)
        map_set.add(x)
    dbs.append(db)
count = 0
for i in map_set:
    map_dct[i] = count
    count += 1
count  = 0
for db in dbs[:]:
    for rec in db[:]:
        t1 = str(rec.seq)
        db.remove(rec)
        t1_id = map_dct[t1]
        mapping[t1].add(t1_id)
        print 'Mapped:', count 
        count += 1
        for db_ in dbs:
            for rec_ in db_:
                t2 = str(rec_.seq)
                if ((t1[1:] in t2) or (t2[1:] in t1)):
                    t2_id = map_dct[t2]
                    mapping[t1].add(t2_id)
                    mapping[t2].add(t1_id)  
    dbs.remove(db)
print len(mapping)
map_ = defaultdict(list)
for i in mapping:
    map_[i] = list(mapping[i])
for i in map_:
    print map_[i]
    break
mapping = json.dumps(map_)
f = open('mapping.json', 'w')
f.write(mapping)
f.close()
