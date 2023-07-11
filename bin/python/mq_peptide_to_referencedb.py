#!/usr/bin/env python3

import pandas as pd
import importlib.machinery
import sys
import os
import collections
from collections import defaultdict
import json
import sequtils
import shutil
import Bio; from Bio import SeqIO
import pickle
import yaml
from yaml import Loader

config = yaml.load(open(sys.argv[1]), Loader=Loader)
output = sys.argv[2]

os.mkdir(output +'/mapping')

def refmap(proteome_id, tax_id, output):
    proteome = list(SeqIO.parse(output+ '/uniprot/{}/{}.fasta'.format(proteome_id, proteome_id),'fasta'))
    peptides = pd.read_csv(config['mq_txt'] +'/peptides.txt',sep='\t')
    peptides = peptides[(peptides['Reverse'].isnull()) & (peptides['Potential contaminant'].isnull())]

    mapped = sequtils.peptides2proteome(proteome, peptides['Sequence'].tolist(), threads=config['threads'])

    outpath = output + '/mapping/{}_peptides.json'.format(proteome_id)
    jstr = json.dumps(mapped.pepdict)

    with open(outpath, 'w') as w:
        w.write(jstr)

for strain in config['reference']:
    proteome_id = config['reference'][strain]['proteome_id']
    taxon_id = config['reference'][strain]['taxon_id']
    refmap(proteome_id, taxon_id, output)
