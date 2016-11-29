#!/usr/bin/env python

import pandas as pd
import sys
from collections import defaultdict
from Bio import SeqIO
from collections import Counter

data = pd.read_csv(sys.argv[1])

def get_organism(val):
    try:
        val = val.split('OS=')[1]
        if '=' in val:
            val = val.split('=')[0]
            val = ' '.join(val.split()[:-1])
        return val
    except:
        return val.split('\n')[0]

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
#data['Organism'] = data['Records'].apply(get_organism)
#org_counter = Counter(data['Organism'])

#data['OrganismCounts'] = data['Organism'].apply(lambda x : org_counter[x] )

data = data.sort_values( by='Summed NSAF', ascending = False )
data = data.reset_index()
del data['index']

data['MSMS_list'] = data['MSMS'].apply( lambda x : x.split('\n'))
all_msms = data['MSMS_list'].tolist()
del data['MSMS_list']

msms_set = set()
all_cum= []
for msms in all_msms:
    msms_set.update(msms)
    all_cum.append(len(msms_set))

msms_set_len = len(msms_set)

data['CumulativeMSMSCount'] = pd.Series(all_cum)
data['CumulativeMSMSPercentage'] = data['CumulativeMSMSCount'] / msms_set_len * 100.0

prefix = sys.argv[2].split('.csv')[0]
data.to_csv(sys.argv[2])

data_filt = data.drop_duplicates('CumulativeMSMSPercentage', keep='first')
filt_prefix=prefix + '_compacted'
data_filt.to_csv(filt_prefix +'.csv')
data_filt_recs = data_filt['Records'].tolist()
w  = open(filt_prefix + '.fasta', 'w')
w.write(''.join(data_filt_recs))
w.close()

#data_100 = data[data['CumulativeMSMSPercentage'] <= 100]
#records = data_100['Records'].tolist()
#w= open(prefix +'_100perc_msmms_mapped.fasta','w')
#w.write('\n'.join(records))
#w.close()

#data_99 = data[data['CumulativeMSMSPercentage'] <= 99]
#records = data_99['Records'].tolist()
#w= open(prefix +'_99perc_msmms_mapped.fasta','w')
#w.write('\n'.join(records))
#w.close()

#data_95 = data[data['CumulativeMSMSPercentage'] <= 95]
#records = data_95['Records'].tolist()
#w= open(prefix +'_95perc_msmms_mapped.fasta','w')
#w.write('\n'.join(records))
#w.close()

#data_90 = data[data['CumulativeMSMSPercentage'] <= 90]
#records = data_90['Records'].tolist()
#w= open(prefix +'_90perc_msmms_mapped.fasta','w')
#w.write('\n'.join(records))
#w.close()

#data_80 = data[data['CumulativeMSMSPercentage'] <= 80]
#records = data_80['Records'].tolist()
#w= open(prefix +'_80perc_msmms_mapped.fasta','w')
#w.write('\n'.join(records))
#w.close()





