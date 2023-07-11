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

def _(proteome_id):
    ref = output + '/uniprot/{}/{}.fasta'.format(proteome_id, proteome_id)
    new=output+'/blast/orfs2proteins/{}/{}.fasta'.format(proteome_id, proteome_id)
    newfolder=output+'/blast/orfs2proteins/{}/'.format(proteome_id)

    ref_id = proteome_id

    os.mkdir(output +'/blast/orfs2proteins/{}'.format(proteome_id))

    cmd="cp {} {} && cd {} && makeblastdb -in {} -dbtype 'prot' -out {}".format(ref, newfolder, newfolder, new, ref_id )

    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    assert process.returncode == 0
                     
    #query=output + '/fasta/nr.fasta'
    
    query=output + '/strains/all_mapped_trans_orfs.fasta'
    outfmt=5

    out=output + '/blast/orfs2proteins/{}.xml'.format(proteome_id)

    db = newfolder + '/{}'.format(proteome_id)

    num_threads=5

    cmd="blastp -query {} -outfmt {} -out {} -db {} -max_target_seqs 500 -max_hsps 1 -num_threads {} -evalue 0.0001".format(query, outfmt, out, db, num_threads)

    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

    process.wait()

    assert process.returncode == 0

    cmd="blast_XML_to_csv.py {} {} {} {}".format(out, query, out +'.csv', 500)

    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

    process.wait()

    assert process.returncode == 0

    data = pd.read_csv(out+'.csv')

    data = data[(data['_alignment_rank']==1) & (data['_hsp_rank']==1)]

    mp = defaultdict(list)

    def get_mapping(df):
        ids = df['blast_record.query'].split()[0].split(';')
        #ids = df['blast_record.query'].split()[0].split(';')
        evalue= df['hsp.expect']
        for i in ids:
            print(i)
            i = i.split('|')[1]#.split('.')[0]
            if evalue < 0.0001:
                mapped = df['_alignment.entry']
                if not mapped in mp[i]:
                    mp[i].append(mapped)

    data.apply(get_mapping, axis=1)

    mp = json.dumps(mp)

    w= open(output +'/blast/orfs2proteins/{}_mapping.json'.format(proteome_id),'w')
    w.write(mp)
    w.close()


for strain in config['reference']:
    proteome_id = config['reference'][strain]['proteome_id']

    _(proteome_id)

