#!/usr/bin/env python

import pandas as pd
import tagmatch
import Bio; from Bio import SeqIO
import sys

denovo = list(SeqIO.parse(sys.argv[1],'fasta'))
fragtol=float(sys.argv[2]) # frag tolerance in Da

def get_charge(val):
    return int(val.split(';')[-3][0])
def get_mz(val):
    return float(val.split(';')[-4])
def prec_mw(df):
    mw = tagmatch.mz2mw(df['mz'], df['charge'])
    return mw
def sequence_mw(val):
    pm = tagmatch.peptide_mass(val)
    return pm

def validate_sequence(df):
    query= df['Sequence']
    prec_mass_list = df['Precursor_MW']
    target = df['Sequence']
    tm= tagmatch.TagMatch(query, [prec_mass_list], target, tol= fragtol)
    return tm.validated

records = [i.format('fasta') for i in denovo]
headers = [i.description for i in denovo]
seqs    = [str(i.seq) for i in denovo]

table = pd.DataFrame()
table['Records'] = pd.Series(records)
table['Headers'] = pd.Series(headers)
table['Sequence'] = pd.Series(seqs)

table['mz'] = table['Headers'].apply(get_mz)
table['charge'] = table['Headers'].apply(get_charge)
table['Precursor_MW'] = table.apply(prec_mw, axis=1)
table['Sequence_MW'] = table['Sequence'].apply(sequence_mw)
table['Validated'] = table.apply(validate_sequence, axis=1)
table = table[table['Validated'] ==True]

filtered_records = table['Records'].tolist()
w =open(sys.argv[3],'w')
w.write('\n'.join(filtered_records))
w.close()




