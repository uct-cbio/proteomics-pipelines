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

config = yaml.load(open(sys.argv[1]))
output = sys.argv[2]

os.mkdir(output +'/mapping')

proteome = list(SeqIO.parse(output+ '/uniprot/{}/{}_{}.fasta'.format(config['reference_proteome_id'], config['reference_proteome_id'], config['reference_taxid']),'fasta'))
peptides = pd.read_csv(config['mq_txt'] +'/peptides.txt',sep='\t')
peptides = peptides[(peptides['Reverse'].isnull()) & (peptides['Potential contaminant'].isnull())]

mapped = sequtils.peptides2proteome(proteome, peptides['Sequence'].tolist(), threads=config['threads'])

outpath = output + '/mapping/{}_peptides.p'.format(config['reference_proteome_id'])
pickle.dump( mapped.pepdict, open( outpath, "wb" ) )

