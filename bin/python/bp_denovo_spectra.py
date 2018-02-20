#!/usr/bin/env python

import sys
import Bio; from Bio import SeqIO
import pandas as pd

fasta = SeqIO.to_dict(list(SeqIO.parse(sys.argv[1], 'fasta')))

blast_table = pd.read_csv(sys.argv[2], sep='\t')

sample_df = pd.DataFrame()
sample_df['Sample'] = pd.Series(list(fasta.keys()))
sample_df['Sequence'] = pd.Series(list(fasta.values()))

sample_df.rename(columns={0 : 'Sample', 1 : 'Sequence' }, inplace = True )

def return_seq(val):
    return str(val.seq)

sample_df['Sequence'] = sample_df['Sequence'].apply(return_seq)

left_cols = sample_df.columns
right_cols = blast_table.columns

merged = pd.merge(sample_df, blast_table, how='left', left_on='Sequence', right_on="_query.sequence" )

merged.to_csv(sys.argv[3], sep='\t')


