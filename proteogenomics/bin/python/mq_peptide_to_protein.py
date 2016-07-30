#!/usr/bin/env python

import pandas as pd
import sys
import importlib.machinery
import Bio; from Bio import SeqIO
import sequtils
import shutil
import os

loader = importlib.machinery.SourceFileLoader('config', sys.argv[1])
config = loader.load_module()
output = sys.argv[2]

pg = pd.read_csv(config.mq_txt +'proteinGroups.txt',sep='\t')
pg =pg[(pg['Potential contaminant'] != '+') & (pg['Reverse'] != '+')]

peptides=pd.read_csv(config.mq_txt +'peptides.txt',sep='\t')
peptides=peptides[(peptides['Potential contaminant'] != '+') & (peptides['Reverse'] != '+')]

evidence=pd.read_csv(config.mq_txt +'evidence.txt',sep='\t')
evidence=evidence[(evidence['Potential contaminant'] != '+') & (evidence['Reverse'] != '+')]

strain_dct={}
samples = config.samples

for sample in samples:
    strain_dct[sample] = samples[sample]['STRAIN']

for row in peptides.iterrows():
    print(row[1])
    break

pep2prot=pd.DataFrame()

for row in pg.iterrows():
    peptide_ids = row[1]['Peptide IDs'].split(';')
    evidence_ids = row[1]['Evidence IDs'].split(';')

    row_peps = peptides[peptides['id'].apply(str).isin(peptide_ids)]
    row_evs = evidence[evidence['id'].apply(str).isin(evidence_ids)]

    for sample in samples:
        sample_peps = row_peps[row_peps['Experiment {}'.format(sample)] >= 1]['Sequence'].tolist()
        sample_evs = row_evs[row_evs['Experiment'] == sample]['Modified sequence'].tolist()
        pg.loc[row[0], "_Peptides {}".format(sample)] = '\n'.join(sample_evs)

    break

print(pg.head())



