#!/usr/bin/env python3

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
import pickle
import yaml

config = yaml.load(open(sys.argv[1]))


output = sys.argv[2]
os.mkdir(output + '/fasta')

combined_fasta=[]

def export_fasta(df):
    id=df['ORF_id']
    seq= ''.join(df['ORF_translation'].split('*'))
    seqrecord = SeqRecord(seq = Seq(seq), id = id)
    return seqrecord

for strain in config['strains']:
    datum=output + '/strains/{}/{}_mapped_peptides.p'.format(strain,strain)
    table = pickle.load(open(datum,'rb'))
    table = table[['ORF_id', 'ORF_translation']].drop_duplicates()
    recs = table.apply(export_fasta, axis=1).tolist()
    combined_fasta += recs

SeqIO.write(combined_fasta , output + '/fasta/combined_translated.fasta', 'fasta')
