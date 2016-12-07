#!/usr/bin/env python

import pandas as pd
import sys
import pickle
import sequtils
import collections; from collections import defaultdict

outpath= sys.argv[1]

mapping = pickle.load( open( outpath +'/fasta/id_mapping.p','rb'))

data = pd.read_csv(outpath  + '/fasta/nr_translated_pg_orfs.fasta.tsv', sep=None, engine='python', header=None)
data = data.reset_index()
data.rename(columns={'level_0':'seqid'}, inplace=True)

new_df = pd.DataFrame()
new_df['_mapped.id'] = pd.Series(list(mapping.keys()))
new_df['seqid'] = new_df['_mapped.id'].apply(lambda x : mapping[x])
merged = pd.merge(data, new_df)

#merged['seqid'] = merged['_mapped.id'].apply(lambda x : x.split('|')[1])

id2go= defaultdict(set)
gos = set()
def go(df):
    go = df[9]
    id = df['_mapped.id']
    try:
        go=go.split('|')
        id2go[id].update(go)
        for goterm in go:
            gos.add(goterm)
    except:
        pass
merged.apply(go, axis=1)
pickle.dump( id2go, open( outpath +'/fasta/id2go.p', 'wb'))

id2kegg= defaultdict(set)
keggs = set()
def kegg(df):
    id = df['_mapped.id']
    try:
        kegg = df[10].split('|')
        kegg = [i for i in kegg if i.startswith('KEGG: ')]
        kegg = [i.split('KEGG: ')[1].split('+')[0] for i in kegg]
        id2kegg[id].update(kegg)
        for keggterm in kegg:
            keggs.add(keggterm)
    except:
        pass
merged.apply(kegg, axis=1)
pickle.dump( id2kegg, open( outpath +'/fasta/id2kegg.p', 'wb'))

kegg_df = pd.DataFrame()
kegg_vals = list(keggs)
kegg_df['KEGG_ID'] = pd.Series(kegg_vals)
kegg_df.to_csv(outpath +'/fasta/kegg_terms.csv')

go_df = pd.DataFrame()
go_vals = list(gos)
go_df['GO_ID'] = pd.Series(go_vals)
go_df.to_csv(outpath +'/fasta/go_terms.csv')

