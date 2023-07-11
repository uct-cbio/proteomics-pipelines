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
import sequtils

config = yaml.load(open(sys.argv[1]), yaml.Loader)
output =os.path.abspath( sys.argv[2])

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

pepfasta = output + '/blast/peptides2orfs/peptides.fasta'

SeqIO.write(recs, pepfasta, 'fasta')

def _(genome_path, output):
    genome = list(SeqIO.parse(genome_path, 'fasta'))

    blastname = genome_path.split('/')[-1].split('.fasta')[0]

    orfs = sequtils.sf_contigs(genome, assembly_name = blastname, table=config['translation_table'],      codons='All', peptide_length=1, translated=True )

    orffasta = output + '/blast/peptides2orfs/{}_reforfs.fasta'.format(blastname)
    SeqIO.write(orfs, orffasta, 'fasta')

    #sys.exit()

    blastdir = output +'/blast/peptides2orfs/{}'.format(blastname)

    out = output + '/blast/peptides2orfs/{}.xml'.format(blastname)

    os.mkdir(blastdir)

    cmd="cp {} {} && cd {} && makeblastdb -in {} -dbtype 'prot' -out {}".format(orffasta, blastdir, blastdir, orffasta, blastname)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    assert process.returncode == 0

    num_threads=5

    cmd ="cd {} && blastp -query {} -db {} -out={} -outfmt 5 -max_hsps 1 -num_threads {} -evalue 200000 -matrix PAM30 -gapopen 9 -word_size 2 -gapextend 1 -comp_based_stats 0 -window_size 15 -threshold 16 -seg 'no'".format(blastdir, pepfasta,  blastname, out, num_threads) 

    #cmd ="cd {} && tblastn -query {} -db {} -out={} -outfmt 5 -max_target_seqs 50 -max_hsps 1 -num_threads {} -evalue 200000".format(blastdir, pepfasta,  blastname, out, num_threads) 

    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

    process.wait()

    assert process.returncode == 0


for strain in config['reference']:
    assembly_id = config['reference'][strain]['assembly_id']
    path = output + '/ena/{}/{}.fasta'.format(assembly_id, assembly_id)
    _(path, output)
