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
    

quant = pd.read_csv(output+ '/diff/msnbase/imputed.csv')

id_dict = {}

for row in quant.iterrows():
    try: 
        orfs = row[1]['ORF.ids...all.strains'].split('\n')
    except:
        continue
    orfs = [ i.split('|')[1] for i in orfs ] 
    for orf in orfs:
        if not orf in id_dict:
            id_dict[orf] = defaultdict()
    for col in row[1].index:
        if col.startswith('iBAQ.'):
            s = col.split('.')[1]
            intensity = row[1][col]
            for orf in orfs:
                if not s in id_dict[orf]:
                    id_dict[orf][s] = intensity         
                elif not intensity > id_dict[orf][s]:
                    id_dict[orf][s] = intensity         
orf_coords = {}

for strain in config['strains']:
    
    orf_coords[strain] = defaultdict(list)

    strain_rows = ['##gff-version 3\n##Index-subfeatures 1\n']
    strain_rows += strains[strain]

    gff3 = '\n'.join(strain_rows)
    
    if not os.path.exists(output + '/jbrowse/' + strain):
        os.mkdir(output + '/jbrowse/' + strain)     
    # Write the peptide features 
    
    w =open(output+ '/jbrowse/' + strain + '/{}_reference_orfs.gff3'.format(strain), 'w')
    w.write(gff3)
    w.close()
    
    for row in strains[strain]:
        _ = row.split('\t')
        if _[8].startswith('ID='):
            start = int(_[3])
            end = int(_[4])
            coords = [i for i in range(start, end + 1)]
            Id = _[8].split(';Name=')[1]
            orf_coords[strain][Id] += coords

    samps = defaultdict(list)
    for Id in orf_coords[strain]:
        coords = orf_coords[strain][Id]
        if Id in id_dict: 
            for samp in id_dict[Id]:
                intensity = id_dict[Id][samp]
                for coord in coords:
                    tup = (coord, intensity)
                    samps[samp].append(tup)
    
    for samp in samps:
        tups = samps[samp]
        ints = {}
        for tup in tups:
            if tup[0] in ints:
                ints[tup[0]] += tup[1]
            else:
                ints[tup[0]] = tup[1]
        vals = [(k, v) for k, v in ints.items()]
        vals = sorted(vals, key=lambda x: x[0])
        vals=  [ ' '.join([str(v[0]), str(v[1])]) for v in vals ] 
        header = 'variableStep chrom=Chromosome\n' + '\n'.join(vals)
        w =open(output+ '/jbrowse/' + strain + '/{}_{}_reference_orfs.wig'.format(strain,samp), 'w')
        w.write(header)
        w.close()




