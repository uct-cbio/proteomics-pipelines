#!/usr/bin/env python3

import pandas as pd
import sys
import multiprocessing
import numpy as np
import pandas as pd
import sys
import importlib.machinery
import Bio; from Bio import SeqIO
import sequtils
import shutil
import algo
import os
import subprocess
from collections import defaultdict
import json
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Seq import translate
import sys
from collections import Counter
from io import StringIO
import uniprot
import pickle
from io import StringIO
import yaml
import blast
from Bio.Blast import NCBIXML
import gff3

config = yaml.load(open(sys.argv[1]).read(), Loader=yaml.Loader)
output = sys.argv[2]

query_fasta = list(SeqIO.parse(output +'/strains/all_mapped_trans_orfs.fasta','fasta'))

strain_samples = defaultdict(list)

design = pd.read_csv(config['design'])

if 'rename' in design.columns:
    scol = 'rename'
else:
    scol = 'sample'
assert scol in design.columns
design = design[design['exclude'] != '+' ] 

rename_dict = sequtils.rename_design(design)

sample2strain = design.set_index(scol)['Strain'].to_dict()
samples = list(sample2strain.keys())

print(samples)
print(rename_dict)

for sample in samples:
    st = sample2strain[sample]
    strain_samples[st].append(rename_dict[sample])

peptide_sequence  = SeqIO.to_dict(list(SeqIO.parse(output + '/blast/peptides2orfs/peptides.fasta','fasta')))
# Per strain analysis
peptides = pd.read_csv(config['mq_txt'] + '/peptides.txt', sep='\t',engine='python')
peptides = peptides[(peptides['Potential contaminant'].isnull()) & (peptides['Reverse'].isnull())]
all_peptides = peptides['Sequence'].tolist()
exp_cols = [i for i in peptides.columns if i.startswith('Experiment')]


for ref in config['reference']:
    assembly_id = config['reference'][ref]['assembly_id']
    gff3.peptide_blast_gff3(assembly_id, output, peptide_sequence, config, strain_samples, peptides)
    gff3.orf_blast_gff3(config, assembly_id, query_fasta , output)
