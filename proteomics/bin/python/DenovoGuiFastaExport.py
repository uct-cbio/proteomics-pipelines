#!/usr/bin/env python

import multiprocessing
import pandas as pd
import sys
import os
from io import StringIO
import re
import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


file = sys.argv[1] 

threads=int(sys.argv[2])

def clean_pep(val):
    regex = re.compile('[^a-zA-Z]')
    rval = regex.sub('', val)
    return rval

def process_res(val):
    r = val
    r = r.split('\n')
    header = r[0]
    
    exclude =[ '# No solutions found.', 
            '# too few peaks...' , 
            '# Could not process spectrum...']
    if header != '': 
        cols = r[1]
        if cols not in exclude:

            t = StringIO('\n'.join(r[1:]))
            table = pd.read_csv(t, sep='\t')
            table['_Header'] = header
            table['_Peptide'] = table['Sequence'].apply(clean_pep)
            table['_Length'] = table['_Peptide'].apply(lambda x : len(x))
            table['_ID'] = table['_Header'].apply(lambda x: ','.join(x.split())) + ',Rank:' + table['#Index'].apply(str) 
            return table

def export_fasta(table):
    rec = SeqRecord(id = table['_ID'], seq = Seq(table['_Peptide']))
    return rec

tables = []

file =open(file).read()


results = file.split('>>')
pool = multiprocessing.Pool(threads)
res = pd.concat(pool.map(process_res, results))
tables.append(res)

combined = pd.concat(tables).reset_index()
del combined['index']

fasta = combined.apply(export_fasta, axis=1).tolist()

combined.to_csv(sys.argv[1] +'_peptide_export.csv')
SeqIO.write(fasta, sys.argv[1] + '_peptide_export.fasta', 'fasta')

