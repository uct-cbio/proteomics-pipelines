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

blastfile = config['reference_genome']
blastname = blastfile.split('/')[-1].split('.')[0]
blastdir = output +'/blast/orfs2genome/{}'.format(blastname)
out = output + '/blast/orfs2genome/{}.xml'.format(blastname)

query=output + '/fasta/nr_translated_pg_orfs.fasta'

os.mkdir(blastdir)

cmd="cp {} {} && cd {} && makeblastdb -in {} -dbtype 'nucl' -out {}".format(blastfile, blastdir, blastdir, blastfile, blastname)
process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
process.wait()
assert process.returncode == 0

num_threads=1
cmd ="cd {} && tblastn -query {} -db {} -out={} -outfmt 5 -max_target_seqs 500 -max_hsps 1 -num_threads {} -evalue 0.0001".format(blastdir,query, blastname, out, num_threads) 

process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
process.wait()
assert process.returncode == 0
