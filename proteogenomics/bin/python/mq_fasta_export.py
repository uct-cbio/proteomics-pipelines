#!/usr/bin/env python

import pandas as pd
import numpy as np
import pandas as pd
import sys
import importlib.machinery
import Bio; from Bio import SeqIO
import os
import shutil
from collections import defaultdict
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Seq import translate
import json

loader = importlib.machinery.SourceFileLoader('config', sys.argv[1])
config = loader.load_module()
output = sys.argv[2]

combined_fasta=[]

def export_fasta(df):
    id=df['ORF_ids']
    seq= ''.join(df['ORF_translation'].split('*'))
    seqrecord = SeqRecord(seq = Seq(seq), id = id)
    return seqrecord

for strain in config.strains:
    datum=output + '/{}/{}_mapped_peptides.csv'.format(strain,strain)
    table = pd.read_csv(datum)
    table = table[['ORF_ids', 'ORF_translation']].drop_duplicates()
    recs = table.apply(export_fasta, axis=1).tolist()
    combined_fasta += recs
try:
    os.mkdir(output + '/fasta')
except:
    pass

SeqIO.write(combined_fasta , output + '/fasta/combined_translated.fasta', 'fasta')
