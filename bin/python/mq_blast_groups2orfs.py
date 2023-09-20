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

config = yaml.load(open(sys.argv[1]), Loader=yaml.Loader)
output = os.path.abspath( sys.argv[2])

query_fasta  = os.path.abspath(config['search_fasta'])
target_fasta = output + '/strains/all_mapped_trans_orfs.fasta'

txt_path = os.path.abspath(config['mq_txt'])
pg = pd.read_csv(txt_path +'/proteinGroups.txt', sep='\t')
pg_ids = []
for p in pg['Protein IDs'].tolist():
    pids = p.split(';')
    for pp in pids:
        pg_ids.append(pp)
pg_ids = set(pg_ids)

queries = list(SeqIO.parse(query_fasta,'fasta'))
new_queries = []

for q in queries:
    ID = q.id
    if ID in pg_ids:
        new_queries.append(q)
    elif ID.split('|')[1] in pg_ids:
        new_queries.append(q)


newfolder = output +'/blast/groups2orfs'
new_query_path = output +'/blast/groups2orfs/filtered.fasta'
with open(new_query_path,'w') as w:
    SeqIO.write(new_queries, w, 'fasta')

cmd="cp {} {} && cd {} && makeblastdb -in {} -dbtype 'prot' -out {}".format(target_fasta, newfolder, newfolder, 'all_mapped_trans_orfs.fasta', 'mapped_orfs' )
process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
process.wait()
assert process.returncode == 0
                 
#query=output + '/fasta/nr.fasta'

outfmt=5

out=output + '/blast/groups2orfs/mapped_groups.xml'

db = newfolder + '/{}'.format('mapped_orfs')

num_threads=5

cmd="blastp -query {} -outfmt {} -out {} -db {} -max_target_seqs 500 -max_hsps 1 -num_threads {} -evalue 0.0001".format(new_query_path, outfmt, out, db, num_threads)
process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

process.wait()

assert process.returncode == 0

cmd="blast_XML_to_csv.py {} {} {} {}".format(out, query_fasta, out +'.csv', 500)

process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

process.wait()

assert process.returncode == 0

data = pd.read_csv(out+'.csv')


#data = data[(data['_alignment_rank']==1) & (data['_hsp_rank']==1)]

mp = defaultdict(list)

def get_mapping(df):
    ids = df['blast_record.query'].split()[0].split(';')
    #ids = df['blast_record.query'].split()[0].split(';')
    evalue= df['hsp.expect']
    for i in ids:
        print(i)
        i = i.split('|')[1]#.split('.')[0]
        #if evalue < 0.0001:
        if evalue == 0:
            mapped = df['_alignment.entry']
            if not mapped in mp[i]:
                mp[i].append(mapped)

data.apply(get_mapping, axis=1)

mp = json.dumps(mp)

w= open(output +'/blast/groups2orfs/mapped_orfs_mapping.json','w')
w.write(mp)
w.close()



