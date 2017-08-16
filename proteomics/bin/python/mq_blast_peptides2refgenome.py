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

peptides = pd.read_csv(config['mq_txt'] + '/peptides.txt', sep='\t',engine='python')
peptides = peptides[(peptides['Potential contaminant'].isnull()) & (peptides['Reverse'].isnull())]
all_peptides = list(set(peptides['Sequence'].tolist()))

count = 1
recs = []

for peptide in all_peptides:
    recid = "PEPTIDE_{}".format(count)
    seq = Seq(peptide)
    rec = SeqRecord(id = recid, seq=seq)
    recs.append(rec)
    count += 1

pepfasta = output + '/blast/peptides2genome/peptides.fasta'

SeqIO.write(recs, pepfasta, 'fasta')


blastfile = config['reference_genome']

blastname = blastfile.split('/')[-1].split('.')[0]

blastdir = output +'/blast/peptides2genome/{}'.format(blastname)

out = output + '/blast/peptides2genome/{}.xml'.format(blastname)

os.mkdir(blastdir)

cmd="cp {} {} && cd {} && makeblastdb -in {} -dbtype 'nucl' -out {}".format(blastfile, blastdir, blastdir, blastfile.split('/')[-1], blastname)
process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
process.wait()

assert process.returncode == 0

num_threads=5
cmd ="cd {} && tblastn -query {} -db {} -out={} -outfmt 5 -max_target_seqs 50 -max_hsps 1 -num_threads {} -evalue 200000 -matrix PAM30 -gapopen 9 -word_size 2 -gapextend 1 -comp_based_stats 0 -window_size 15 -threshold 16 -seg 'no'".format(blastdir, pepfasta,  blastname, out, num_threads) 
#cmd ="cd {} && tblastn -query {} -db {} -out={} -outfmt 5 -max_target_seqs 50 -max_hsps 1 -num_threads {} -evalue 200000".format(blastdir, pepfasta,  blastname, out, num_threads) 

process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

process.wait()

assert process.returncode == 0
