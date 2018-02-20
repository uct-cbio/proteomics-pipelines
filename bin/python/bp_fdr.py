#!/usr/bin/env python

import pandas as pd
import sys

table = sys.argv[1]
score_column = 'hsp.score'
description_column = 'alignment.hit_def'
rev_tag = sys.argv[2]

data = pd.read_csv(table, sep=None, engine='python')
data = data.sort_values(by=score_column,  ascending=False)
data = data.reset_index()
del data['index']

data['Reverse'] = data[description_column].apply(lambda x : rev_tag in x)

print(data.head())

print(len(data[data['Reverse'] == False]))
print(len(data[data['Reverse'] == True]))


for n, g in data.groupby('hsp.score'):
    sample = data[data['hsp.score'] >= n]
    t =len( sample[sample['Reverse'] == False] )
    d = len(sample[sample['Reverse'] == True ] )

    perc = (float(d) / float(len(sample))) * 100

    data.loc[g.index, score_column + '.cutoff.fdr.percent'] = perc



data.to_csv('fdr_analysis.csv')
