#!/usr/bin/env python3

import pandas as pd
import sys
import pickle
import sequtils
import collections; from collections import defaultdict
import os
import json
import mqparse

outpath= sys.argv[1]
keggid= sys.argv[2]
proteins = outpath +'/msnbase/normalized.csv'
proteingroups = pd.read_csv(proteins)

proteingroups = proteingroups[proteingroups['ORF.ids...all.strains'].notnull()]
orf2pg = defaultdict(list)
for row in proteingroups.iterrows():
    ids = row[1]['ORF.ids...all.strains'].split('\n')
    identifier = row[1]['Identifier']
    for i in  ids:
        orf2pg[i].append(identifier)

if not os.path.exists(outpath +'/gsea'):
    os.mkdir(outpath +'/gsea')

with open( outpath +'/fasta/id_mapping.json') as f:
    mapping = json.loads(f.read())
mapping2pg = defaultdict(set)
for key in mapping:
    if key in orf2pg:
        pg = orf2pg[key]
        m = mapping[key]
        mapping2pg[m].update(pg)

newmap = {}
for key in mapping2pg:
    newmap[key] = '|'.join(list(mapping2pg[key]))

data = pd.read_csv(outpath  + '/fasta/nr_translated_pg_orfs.fasta.tsv', sep='\t', engine='python', header=None)
data = data.reset_index()
data.rename(columns={'level_0':'seqid'}, inplace=True)

new_df = pd.DataFrame()
new_df['seqid'] = pd.Series(list(newmap.keys()))
new_df['_mapped.id'] = new_df['seqid'].apply(lambda x : newmap[x])
merged = pd.merge(data, new_df)
print(new_df.head())
print(data.head())
print(merged.head())
#del merged['index']
merged['seqid'] =merged['_mapped.id']
del merged['_mapped.id']
print(merged.head())
ipr = outpath +'/gsea/interproscan.tsv'
merged.to_csv(ipr, sep='\t')
mqparse.mq_txt.ips_genesets(None, ipr, proteingroups, outpath + '/gsea', keggid=keggid, id_col="Identifier")



def proteingroup_genesets(df, id_col, set_col, set_name):   
    mapping = defaultdict(list)
    for row in df.iterrows():
        for m in str(row[1][set_col]).split(';'):
                mapping['door2_'+ m].append(row[1][id_col])
    genes = []
    keys = []
    for key in mapping:
        keys.append(key)
        genes.append('|'.join(mapping[key]))    
    new_pg = pd.DataFrame()
    new_pg['GENES'] = pd.Series(genes)
    new_pg[set_name] = pd.Series(keys)
    return new_pg

new_pg =proteingroup_genesets(proteingroups, "Identifier", "DOOR2_Operons", "OPERON_ID")
new_pg.to_csv(outpath + '/gsea/operons2proteingroups.csv')
