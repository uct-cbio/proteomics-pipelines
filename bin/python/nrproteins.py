#!/usr/bin/env python

import pandas as pd
import sys
import Bio
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

if len(sys.argv) == 3:
    infile  = sys.argv[1]
    outfile = sys.argv[2]
elif len(sys.argv) == 2:
    infile  = sys.stdin
    outfile = sys.argv[1]
elif len(sys.argv) == 1:
    infile  = sys.stdin
    outfile = sys.stout

df = pd.read_csv(infile)

df = df.drop_duplicates('Non redundant sequence id')
df = pd.Series(df['Translated sequence'].values, index=df['Non redundant sequence id']).to_dict()

recs = []
for key in df:
    print(key)
    seq = Seq(df[key])
    description = 'Six_frame_protein_sequence'
    id = 'Non_redundant_protein_sequence_{}'.format(key)
    rec = SeqRecord(id=id, seq = seq, description = description)
    print rec.format('fasta')
    recs.append(rec)
SeqIO.write(recs, outfile, 'fasta')
