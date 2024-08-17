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
import urllib.parse

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


def download(proteome, taxid, output):
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
    idmapping = pd.read_csv(output + '/{}/{}_{}.idmapping'.format(proteome, proteome, taxid), sep='\t', header=None)
    mapping_dct = {}
    def mapping_dict(df):
        #global mapping_dct
        protid =df[0]
        name=df[1]
        value=df[2]
        if not protid in mapping_dct:
            mapping_dct[protid] = defaultdict(list)
        mapping_dct[protid][name].append(value)
    idmapping.apply(mapping_dict, axis=1)
    idmap = json.dumps(mapping_dct)
    w = open( output + '/{}/{}_{}.idmapping.json'.format(proteome, proteome, taxid), 'w' )
    w.write(idmap)

    # Get features
    c1="UniProtRetrieve.sh {} {}.gff gff".format(taxid, taxid)
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

def download2(proteome_id, output):
    #url='https://www.uniprot.org/uniprot/?query=proteome:{}&format=fasta&include=yes&force=true'.format(proteome_id)
    url='https://rest.uniprot.org/uniprotkb/stream?compressed=true&download=true&format=fasta&includeIsoform=true&query=%28{}%29'.format(proteome_id)
    print(url)
    path = output +'/{}'.format(proteome_id)
    try:
        os.makedirs(path)
    except:
        pass
    #for request in requests:
    
    command = "cd {} && curl '{}' > {}.fasta.gz  && gunzip *.gz && echo '{}' > urls.txt".format(path, url, proteome_id, url)
    
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    process.wait()
    assert process.returncode == 0

    #url='https://www.uniprot.org/uniprot/?query=proteome:{}&format=gff&include=yes&force=true'.format(proteome_id)
    url='https://rest.uniprot.org/uniprotkb/stream?compressed=true&download=true&format=gff&includeIsorform=true&query=%28{}%29'.format(proteome_id)
    
    command = "cd {} && curl '{}' > {}.gff3.gz && gunzip *.gz  && echo '{}' >> urls.txt".format(path, url, proteome_id, url)
    
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    process.wait()
    assert process.returncode == 0
    
    #url='https://www.uniprot.org/uniprot/?query=proteome:UP000236349&include=yes&format=tab&force=true&columns=id,entry%20name,reviewed,protein%20names,genes,organism,length,database(KEGG),comment(PATHWAY),genes(OLN),genes(ORF),database(GeneID)'
    url='https://rest.uniprot.org/uniprotkb/stream?compressed=true&includeIsoform=true&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength%2Ccc_alternative_products%2Ccc_mass_spectrometry%2Cgene_orf%2Cgene_oln%2Cgene_primary%2Corganism_id%2Cxref_kegg%2Ccc_pathway&format=tsv&query=%28{}%29+AND+%28reviewed%3Afalse%29'.format(proteome_id)
    command = "cd {} && curl '{}' > {}.tab.gz && gunzip *.gz && echo '{}' >> urls.txt".format(path, url, proteome_id, url)
    
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    process.wait()
    assert process.returncode == 0
    return
    
output = output +'/uniprot'

for ref in config['reference']:
    proteome_id = config['reference'][ref]['proteome_id']
    download2(proteome_id, output)

    

