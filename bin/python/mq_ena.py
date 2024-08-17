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

config = yaml.load(open(sys.argv[1]).read(), Loader=yaml.CLoader)

output = sys.argv[2]

try:
    os.makedirs(output +'/ena')
except:
    shutil.rmtree(output +'/ena')
    os.mkdir(output +'/ena')


def download(genome, output):
    path = output +'/{}'.format(genome)
    try:
        os.makedirs(path)
    except:
        pass
    c1 = "wget 'https://www.ebi.ac.uk/ena/browser/api/embl/{}?download=true&gzip=true' -O {}.embl.gz && gunzip {}.embl.gz".format(genome, genome, genome)
    c2 = "wget 'https://www.ebi.ac.uk/ena/browser/api/fasta/{}?download=true&gzip=true' -O {}.fasta.gz && gunzip {}.fasta.gz".format(genome, genome, genome)
    requests = [c1]
    for request in requests:
        command = "cd {} && {}".format(path, request)
        print(command)
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        process.wait()
        assert process.returncode == 0
    c = 'cd {} && seqret -sformat embl -sequence {}.embl -feature -osformat fasta -supper -osname {} -offormat gff3 -ofname {}.gff3 -auto'.format(path,genome,genome,genome)
    process = subprocess.Popen(c, shell=True, stdout=subprocess.PIPE)
    process.wait()
    assert process.returncode == 0

output = output +'/ena'
for ref in config['reference']:
    assembly_id = config['reference'][ref]['assembly_id']
    download(assembly_id, output)

    

