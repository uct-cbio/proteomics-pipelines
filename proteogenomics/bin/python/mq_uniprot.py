#!/usr/bin/env python

import subprocess
import importlib.machinery 
import sys
import os
import subprocess



loader = importlib.machinery.SourceFileLoader('config', sys.argv[1])
config = loader.load_module()
output = sys.argv[2]

try:
    os.mkdir(output +'/uniprot')
except:
    pass

output = output +'/uniprot'


proteome = config.reference_proteome_id
taxid = config.reference_taxid
path = output +'/{}'.format(proteome)
try:
    os.mkdir(path)
except:
    pass

c1='ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README'
c2='ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/{}_{}.fasta.gz'.format(proteome, taxid)
c3='ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/{}_{}.gene2acc.gz'.format(proteome, taxid)
c4='ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/{}_{}.idmapping.gz'.format(proteome, taxid)
c4='ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/{}_{}_DNA.fasta.gz'.format(proteome, taxid)
c5='ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/{}_{}_DNA.miss.gz'.format(proteome, taxid)
c6='ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/{}_{}_additional.fasta.gz'.format(proteome, taxid)   

# Unzip everything
c7="cd {} && gunzip * && rm -rf *.gz".format(path)

requests = [c1, c2, c3, c4, c5, c6]
for request in requests:
    command = "cd {} && wget {}".format(path, request)
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    process.wait()

command = c7
process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
process.wait()


# Get features
c1="UniProtRetrieve.sh {} {}.gff gff".format(config.reference_taxid, config.reference_taxid)
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
