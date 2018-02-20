#!/usr/bin/env python

import pandas as pd 
import sys
import re

table = pd.read_csv(sys.argv[1])

hsp_score_cutoff= int(sys.argv[2])

table = table[table['hsp.score'] >= hsp_score_cutoff ] 

table['_hsp.sbjct.cleaned'] = table['hsp.sbjct'].apply( lambda x : ''.join(x.split('-')))

#def check_match(val):
#    matches = re.findall(r"\w+",val)
#    longest = [len(i) for i in matches].max()
#    if longest >= min_match:
#        return True
#    else:
#        return False

#table = table[table['hsp.positives'] == table['_query.sequence.length']]
#table = table[(table['_alignment_rank'] == 1) & (table['_hsp_rank'] ==1)]

#table = table[table['hsp.match'].apply(check_match) ==True]

table.to_csv(sys.argv[3])
