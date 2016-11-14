#!/usr/bin/env python

import pandas as pd
import sys
import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np
import multiprocessing
from Bio import SeqIO
import gc

fname = sys.argv[1].split('/')[-1]
peptide_seqs=list(SeqIO.parse(sys.argv[1], 'fasta'))
peptide_seqs = [str(rec.seq) for rec in peptide_seqs]
nput = sys.argv[2]
tput = sys.argv[3] + '/' + fname  #give the output dir , the file name will be determined b y input
sc_cutoff=int(sys.argv[4])

df_filt = []
chunksize = 1000000
for chunk in pd.read_csv(nput, chunksize=chunksize):
    filt = chunk[chunk['_hsp.sbjct.cleaned'].isin(peptide_seqs)]
    df_filt.append(filt)
    del filt

df_filt = pd.concat(df_filt) # create a dataframe containing only the subject sequences required

gc.collect()

#nput = nput[(nput['_alignment_rank'] == 1)&(nput['_hsp_rank'] ==1)]
#nput = nput[nput['hsp.positives'] == nput['_query.sequence.length']]

seqcount=1

fastas = []

for peptide in peptide_seqs:
    group = df_filt[df_filt['_hsp.sbjct.cleaned'] == peptide]
    scanlist = group['blast_record.query'].apply(lambda x : x.split('scans=')[1].split('|'))
    comb_scans = set()
    for scans in scanlist:
        for scan in scans:
            comb_scans.add(scan)
    sc = len(comb_scans)
    if sc >= sc_cutoff:
        rec = SeqRecord(seq = Seq(peptide), id='subject_export|{}.{}|'.format(fname, str(seqcount)), description='sc={};scans={}'.format(str(sc), '|'.join(comb_scans)))
        fastas.append(rec)
        seqcount += 1

SeqIO.write(fastas, tput, 'fasta')
