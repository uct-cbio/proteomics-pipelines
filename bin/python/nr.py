#!/usr/bin/env python

import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys

df = pd.read_csv(sys.stdin)

nr_count = 0
records = []

for name, group in df.groupby('Translated sequence'):
    df.loc[group.index , 'Non redundant sequence id'] = nr_count
    id          = 'Non_redundant_protein_sequence_{}'.format(nr_count)
    seq         = Seq(name)
    description = 'Six_frame_sequence'
    rec = SeqRecord(id = id, seq = seq, description = description)
    records.append(rec)
    nr_count += 1


df.to_csv(sys.stdout)

