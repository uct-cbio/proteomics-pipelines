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

config = yaml.load(open(sys.argv[1]).read())
output = sys.argv[2]

query_fasta = list(SeqIO.parse(output +'/fasta/nr_translated_pg_orfs.fasta','fasta'))

blast_name= config['reference_genome'].split('/')[-1].split('.fasta')[0]

blast_results = output + '/blast/orfs2genome/{}.xml'.format(blast_name)

result_handle = open(blast_results)

blast_records = NCBIXML.parse(result_handle)

recdict = gff3.tblastn2gff3(blast_records=blast_records, query_fasta=query_fasta, eval_cutoff=200000, subject_cutoff=1)

strains = defaultdict(list)

for key in recdict:
    s = key.split('_')[0]
    strains[s] += recdict[key]
    

for strain in config['strains']:
    strain_rows = ['##gff-version 3\n##Index-subfeatures 1\n']
    strain_rows += strains[strain]

    gff3 = '\n'.join(strain_rows)
    
    if not os.path.exists(output + '/jbrowse/' + strain):
        os.mkdir(output + '/jbrowse/' + strain)     
    # Write the peptide features 
    
    w =open(output+ '/jbrowse/' + strain + '/{}_reference_orfs.gff3'.format(strain), 'w')
    w.write(gff3)
    w.close()
