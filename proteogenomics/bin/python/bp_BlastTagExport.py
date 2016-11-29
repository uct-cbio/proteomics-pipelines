#!/usr/bin/env python

import pandas as pd
import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

data = pd.read_csv(sys.argv[1])
print(data.columns)
data = data[data['_hsp.sbjct.cleaned'] != data['_query.sequence']]
tags = data['blast_record.query'].apply(lambda x : x.split('scans=')[1].split('|')).tolist()
subjects = data['_hsp.sbjct.cleaned'].tolist()

tag_dict = defaultdict(set)
for val in range(len(subjects)):
    subjct = subjects[val]
    tag_dict[subjct].update(tags[val])
corrected = []
for subjct in tag_dict:
    scans = tag_dict[subjct]
    for scn in scans:
        rec = SeqRecord(id=scn, seq=Seq(subjct))
        corrected.append(rec)
SeqIO.write(corrected, sys.argv[2], 'fasta')


