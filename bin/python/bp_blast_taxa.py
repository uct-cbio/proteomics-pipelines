#!/usr/bin/env python3

import pandas as pd
import sys
from collections import Counter
import numpy as np

infile= sys.argv[1]
outfile = sys.argv[2]

hits = pd.read_csv(infile)
def clean_score(val):
    try:
        return float(val)
    except:
        return np.nan
hits['hsp.score'] = hits['hsp.score'].apply(clean_score)
hits = hits[hits['hsp.score'].notnull()]
hits['MatchProb'] = hits['hsp.score'] / hits['hsp.score'].sum()
hits['Accession'] = hits['_alignment.entry']
# Protein Prior Prob
protein_counts = Counter(hits['_alignment.entry'])
protein_probs = pd.Series(protein_counts, name='_alignment.entry')
protein_probs = protein_probs/protein_probs.sum()
protein_probs = pd.DataFrame(protein_probs)
protein_probs.rename(columns={'_alignment.entry':"AccessionProb"}, inplace=True)
protein_probs.sort_values('AccessionProb', ascending=False, inplace=True)
protein_probs = protein_probs.to_dict(orient='dict')

def get_org(val):
    try:
        val = val.split('OS=')[1]
        if '=' in val:
            val = val.split('=')[0].split()[:-1]
            val = ' '.join(val)
        return val
    except:
        return val

# Get Organism probs
hits['OS'] = hits['_alignment.description'].apply(get_org)
organism_counts = Counter(hits['OS'])
org_probs = pd.Series(organism_counts, name='OS')
org_probs = org_probs/org_probs.sum()
org_probs = pd.DataFrame(org_probs)
org_probs.rename(columns={'OS':"OSProb"}, inplace=True)
org_probs.sort_values('OSProb',ascending=False, inplace=True)
org_probs = org_probs.to_dict(orient='dict')

# Probablity Table

hits['OSProb'] = hits['OS'].apply(lambda x : org_probs['OSProb'][x])
hits['AccessionProb'] = hits['_alignment.entry'].apply(lambda x : protein_probs['AccessionProb'][x])

hits['CombinedProb'] = hits['OSProb'] * hits['AccessionProb'] * hits['MatchProb']
hits = hits.sort_values('CombinedProb', ascending=False)
print(hits.head())

# Inference table
matches = hits.drop_duplicates('_query.sequence', keep='first').reset_index()
del matches['index']
matches = matches.sort_values('CombinedProb',ascending=False)
print(matches.head())

hits.to_csv(outfile)
matches.to_csv(outfile +'.filtered.csv')

# inferred organisms
org_counts = Counter(matches['OS'].tolist())
inferred = pd.DataFrame(pd.Series(org_counts, name='OS'))
inferred = inferred.sort_values('OS',ascending=False)
inferred.to_csv(outfile +'.inferred_taxa.csv')

