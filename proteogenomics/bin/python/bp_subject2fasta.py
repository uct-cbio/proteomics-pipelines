#!/usr/bin/env python

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import sys
import Bio; from Bio import SeqIO


data  = pd.read_csv(sys.argv[1], sep=None, engine='python')

def check_subject(value):
    test = value.split()
    if len(test) >  1:
        return False
    if '-' in value:
        return False
    else:
        return True

recs = []

peptides = list(set(data['_hsp.sbjct.cleaned'].values.tolist()))
count=1

for peptide in peptides:
    id='subject_export|peptide_{}|'.format(str(count))
    seq = Seq(peptide)
    rec = SeqRecord(id=id, seq=seq)
    recs.append(rec)
    count += 1

SeqIO.write(recs, sys.argv[2], 'fasta')
