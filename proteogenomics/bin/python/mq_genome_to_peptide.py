#!/usr/bin/env python

import pandas as pd
import sys
import importlib.machinery
import Bio; from Bio import SeqIO
import sequtils
import shutil
import os
from collections import defaultdict

loader = importlib.machinery.SourceFileLoader('config', sys.argv[1])
config = loader.load_module()

output = sys.argv[2]

peptides = pd.read_csv(config.mq_txt + 'peptides.txt', sep='\t',engine='python')

peptides = peptides[(peptides['Potential contaminant'] != '+') & (peptides['Reverse'] != '+')]

strain_samples = defaultdict(list)

samples = config.samples

for sample in samples:
    st = samples[sample]['STRAIN']
    strain_samples[st].append(sample)

exp_cols = [i for i in peptides.columns if i.startswith('Experiment')]

#print(peptides[exp_cols])
#peptides = peptides[:20] # DEV
#peptide_sequences = peptides['Sequence'].tolist()
#assert 4 == 5
# Specific peptides for each strain

for strain in config.strains:
    samples = strain_samples[strain]
    sample_columns = []

    for sample in samples:
        sample_columns.append('Experiment {}'.format(sample))

    strain_filt = peptides[peptides[sample_columns].sum(axis=1) >=1]
    strain_peptides = strain_filt['Sequence'].tolist()
    paths = config.strains[strain]
    genome = list(SeqIO.parse(paths['sf_genome'],'fasta'))
    g2p = sequtils.peptides2genome(genome, assembly_name = str(strain), translation_table=config.translation_table, peptides_list=strain_peptides, threads=config.threads)
    speps = g2p.peptides
    strainpath=output +'/' + strain
    try:
        shutil.rmtree(strainpath)
        os.mkdir(strainpath)
    except:
        os.mkdir(strainpath)    
    speps.to_csv(strainpath + '/' + '{}_mapped_peptides.csv'.format(str(strain)))
    


