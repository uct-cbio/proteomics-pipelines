#!/usr/bin/env python3
import json
import pandas as pd
import sys
import importlib.machinery
import Bio; from Bio import SeqIO
import sequtils
import shutil
import os
from collections import defaultdict
import pickle
import yaml
import pgfunctions
from Bio import SeqIO
from collections import Counter
import pgfunctions
import functools
import numpy as np
import mqparse

config = yaml.load(open(sys.argv[1]).read(),  Loader=yaml.Loader)

output = sys.argv[2]

annotation_folder = output +'/annotations/'

output_folder = output + '/stats/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)

for reference in config['reference']:
    fpath = annotation_folder + '/export/{}/Combined_Annotated_ORFs.csv'.format(reference)
    tpath = annotation_folder + '/export/{}/protein2kegg.csv'.format(reference)
    if os.path.exists(tpath):
        continue
    d = pd.read_csv(fpath)
    proteins = d.drop_duplicates(['BLASTP'], keep='first')
    del proteins['ORF_id']
    proteins = mqparse.mq_txt.leading_protein_ko(None, proteins, 'BLASTP')
    rename = {}
    for col in proteins.columns:
        pref = 'Leading Protein '
        if col.startswith(pref):
            target = col.split(pref)[1]
            rename[col] = target
    proteins.rename(columns=rename, inplace=True)
    proteins.to_csv(tpath)
