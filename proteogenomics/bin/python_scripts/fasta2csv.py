#!/usr/bin/python

import pandas
import sys
import os
import Bio; from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import translate
import sequtils
import pandas as pd

if len(sys.argv)==3:
    infile  = sys.argv[1]
    outfile = sys.argv[2]

elif len(sys.argv) == 1:
    infile = sys.stdin
    outfile = sys.stdout

elif len(sys.argv)==2:
    infile = sys.argv[1]
    outfile = sys.stdout

def translate_start(seq, starts =['ATG','GTG','TTG']): # list of Met-initiator codons
    if seq[:3] in starts:
        translated = 'M' + str(translate(Seq(seq), table =11))[1:]
    else:
        translated = str(translate(Seq(seq), table = 11))
    if translated.endswith('*'):
        translated = translated[:-1]
    return translated


def fasta2dataframe(fasta, min_trans_len =20):
    new_df = pd.DataFrame()
    count = 0
    ids = []
    descriptions = []
    seqs         = []
    for rec in fasta:
        id = rec.id
        description = rec.description
        sequence = str(rec.seq)
        ids.append(id)
        descriptions.append(description)
        seqs.append(sequence)
    new_df['Id'] = pd.Series(ids)
    new_df['Description'] = pd.Series(descriptions)
    new_df['Sequence']    = pd.Series(seqs)
    new_df['Contig'] = new_df['Id'].apply(lambda x: '|'.join(x.split('|')[:-1]))
    new_df['ORF strand'] = new_df['Id'].apply(lambda x: x.split('|')[-1].split('(')[1].split(')')[0])
    new_df['ORF start'] = new_df['Id'].apply(lambda x: x.split('|')[-1].split('(')[1].split(')')[1].split(':')[0])
    new_df['ORF end'] = new_df['Id'].apply(lambda x: x.split('|')[-1].split('(')[1].split(')')[1].split(':')[1])
    new_df['ORF Number'] = new_df['Description'].apply(lambda x: x.split('_')[-1].split('.')[0])
    new_df['ORF Variant'] = new_df['Description'].apply(lambda x: x.split('_')[-1].split('.')[1])
    new_df['Translated sequence']  = new_df['Sequence'].apply(translate_start)
    new_df['Translated sequence length'] = new_df['Translated sequence'].apply(len)
    new_df = new_df[new_df['Translated sequence length'] >= min_trans_len]
    return new_df

seqs = list(SeqIO.parse(infile, 'fasta'))
table = fasta2dataframe(seqs)
table.to_csv(outfile)
