#!/usr/bin/env python

import pandas as pd
import sys
from collections import defaultdict
from Bio import SeqIO

data = pd.read_csv(sys.argv[1])

peptide_list = list(SeqIO.parse(sys.argv[2], 'fasta'))
peptide_dict = defaultdict(set)
for _ in peptide_list:
    s = str(_.seq)
    ids = _.description.split('scans=')[1].split('|')
    peptide_dict[s].update(ids)
del peptide_list

data = data[[i for i in data.columns if not i.startswith('Untitled: ')]]
sample_cols = [i for i in data.columns if i.startswith('Sample ')]

nsafcols = []

for col in sample_cols:
    sample = col.split(' (msms)')[0]
    nsafcol = sample + ' (NSAF)'
    data[nsafcol] = data[col]/data['SequenceLength']
    data[nsafcol] = data[nsafcol]/data[nsafcol].sum()
    nsafcols.append(nsafcol)

data['Summed NSAF'] = data[nsafcols].sum(axis=1)

data = data.sort_values(by='Summed NSAF', ascending = False)
data = data.reset_index()
del data['index']

data['peptides_list'] = data['SubjectSequences'].apply( lambda x : x.split('\n'))
all_peptides = data['peptides_list'].tolist()
del data['peptides_list']

msms_set = set()
all_cum= []

for peptides in all_peptides:
    for peptide in peptides:
        msms_set.update(peptide_dict[peptide])
    cum_set_len = len(list(msms_set))
    all_cum.append(cum_set_len)

msms_set_len = len(list(msms_set))

data['CumulativeMSMSCount'] = pd.Series(all_cum)
data['CumulativeMSMSPercentage'] = data['CumulativeMSMSCount'] / msms_set_len * 100.0


prefix = sys.argv[3].split('.csv')[0]
data.to_csv(sys.argv[3])

data_99 = data[data['CumulativeMSMSPercentage'] <= 99]
records = data_99['Records'].tolist()
w= open(prefix +'_99perc_msmms_mapped.fasta','w')
w.write('\n'.join(records))
w.close()

data_95 = data[data['CumulativeMSMSPercentage'] <= 95]
records = data_95['Records'].tolist()
w= open(prefix +'_95perc_msmms_mapped.fasta','w')
w.write('\n'.join(records))
w.close()


data_90 = data[data['CumulativeMSMSPercentage'] <= 90]
records = data_90['Records'].tolist()
w= open(prefix +'_90perc_msmms_mapped.fasta','w')
w.write('\n'.join(records))
w.close()

data_80 = data[data['CumulativeMSMSPercentage'] <= 80]
records = data_80['Records'].tolist()
w= open(prefix +'_80perc_msmms_mapped.fasta','w')
w.write('\n'.join(records))
w.close()





