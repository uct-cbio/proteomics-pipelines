#!/usr/bin/env python

import pandas as pd
import sys
import os

blast_file = sys.argv[1]
aln_cutoff = int(sys.argv[2])
blast_outfile =sys.argv[3]

df = pd.read_csv(blast_file, sep='\t')
df = df[df['_alignment_rank'].apply(int) <= aln_cutoff]

df.to_csv(blast_outfile, sep='\t')
