#!/usr/bin/env python

import pandas as pd
import sys

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

seqs = df['Translated sequence'].tolist()
seqs = list(set(seqs))
seqs_dct = {}
nr_count = 0
for i in seqs:
    seqs_dct[i] = nr_count
    nr_count += 1
df["Non redundant sequence id"] = df["Translated sequence"].map(seqs_dct)
df.to_csv(outfile)
