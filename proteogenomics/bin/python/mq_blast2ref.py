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
import subprocess

config = yaml.load(open(sys.argv[1]))


output = sys.argv[2]

ref = output + '/uniprot/{}/{}_{}.fasta'.format(config['reference_proteome_id'], config['reference_proteome_id'], config['reference_taxid'])
new=output+'/blast/{}/{}_{}.fasta'.format(config['reference_proteome_id'], config['reference_proteome_id'],config['reference_taxid'])

newfolder=output+'/blast/{}/'.format(config['reference_proteome_id'])

ref_id = config['reference_proteome_id']

os.mkdir(output +'/blast/{}'.format(config['reference_proteome_id']))

cmd="cp {} {} && cd {} && makeblastdb -in {} -dbtype 'prot' -out {}".format(ref, newfolder, newfolder, new, ref_id )

process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
process.wait()
assert process.returncode == 0
                 
query=output + '/fasta/nr_translated_pg_orfs.fasta'
outfmt=5
out=output + '/blast/{}.xml'.format(config['reference_proteome_id'])
db = newfolder + '/{}'.format(config['reference_proteome_id'])
num_threads=10

cmd="blastp -query {} -outfmt {} -out {} -db {} -max_target_seqs 1 -max_hsps 1 -num_threads {}".format(query, outfmt, out, db, num_threads)
process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
process.wait()
assert process.returncode == 0

cmd="blast_XML_to_csv.py {} {} {} {}".format(out, query, out +'.csv', 500)
process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
process.wait()
assert process.returncode == 0

data = pd.read_csv(out+'.csv')

data = data[(data['_alignment_rank']==1) & (data['_hsp_rank']==1)]

mp = {}

def get_mapping(df):
    global mp
    ids = df['blast_record.query'].split()[1].split(';')
    evalue= df['hsp.expect']
    for i in ids:
        i = i.split('|')[1]
        if evalue < 0.0001:
            mp[i] = df['_alignment.entry']
data.apply(get_mapping, axis=1)
mp = json.dumps(mp)

w= open(output +'/blast/{}_mapping.json'.format(config['reference_proteome_id']),'w')
w.write(mp)
w.close()






