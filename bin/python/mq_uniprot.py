#!/usr/bin/env python3

import shutil
import subprocess
import importlib.machinery 
import sys
import os
import subprocess
import pickle
import yaml
from collections import defaultdict
import pandas as pd
import json

#loader = importlib.machinery.SourceFileLoader('config', sys.argv[1])
#config = loader.load_module()

config = yaml.load(open(sys.argv[1]).read(), Loader=yaml.CLoader)
#print(config)

output = sys.argv[2]

try:
    os.makedirs(output +'/uniprot')
except:
    shutil.rmtree(output +'/uniprot')
    os.mkdir(output +'/uniprot')

output = output +'/uniprot'


proteome = config['reference_proteome_id']
taxid = config['reference_taxid']
path = output +'/{}'.format(proteome)
try:
    os.makedirs(path)
except:
    pass

c1='ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README'
c2='ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/{}/{}_{}.fasta.gz'.format(proteome, proteome, taxid)
c3='ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/{}/{}_{}.gene2acc.gz'.format(proteome, proteome, taxid)
c4='ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/{}/{}_{}.idmapping.gz'.format(proteome, proteome, taxid)
c4='ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/{}/{}_{}_DNA.fasta.gz'.format(proteome, proteome, taxid)
c5='ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/{}/{}_{}_DNA.miss.gz'.format(proteome, proteome,  taxid)
c6='ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/{}/{}_{}_additional.fasta.gz'.format(proteome, proteome,  taxid)   
c7='ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/{}/{}_{}.idmapping.gz'.format(proteome, proteome,  taxid)   

# Unzip everything
c8="cd {} && gunzip * && rm -rf *.gz".format(path)

requests = [c1, c2, c3, c4, c5, c6, c7]
for request in requests:
    command = "cd {} && wget {}".format(path, request)
    print(command)
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    process.wait()

command = c8
process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
process.wait()

# Create idmapping json
idmapping = pd.read_csv(output + '/{}/{}_{}.idmapping'.format(config['reference_proteome_id'], config['reference_proteome_id'], config['reference_taxid']), sep='\t', header=None)
mapping_dct = {}
def mapping_dict(df):
    global mapping_dct
    protid =df[0]
    name=df[1]
    value=df[2]
    if not protid in mapping_dct:
        mapping_dct[protid] = defaultdict(list)
    mapping_dct[protid][name].append(value)
idmapping.apply(mapping_dict, axis=1)
idmap = json.dumps(mapping_dct)
w = open( output + '/{}/{}_{}.idmapping.json'.format(config['reference_proteome_id'], config['reference_proteome_id'], config['reference_taxid']), 'w' )
w.write(idmap)

# Get features
c1="UniProtRetrieve.sh {} {}.gff gff".format(config['reference_taxid'], config['reference_taxid'])
path = output +'/features'
try:
    os.mkdir(path)
except:
    pass

requests = [c1]
for request in requests:
    command = "cd {} && {}".format(path,  request)
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    process.wait()
