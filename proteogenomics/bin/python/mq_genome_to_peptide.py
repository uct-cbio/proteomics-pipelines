#!/usr/bin/env python

import pandas as pd
import sys
import importlib.machinery
import Bio; from Bio import SeqIO
import sequtils
import shutil
import os

loader = importlib.machinery.SourceFileLoader('config', sys.argv[1])
config = loader.load_module()

output = sys.argv[2]

peptides = pd.read_csv(config.mq_txt + 'peptides.txt', sep='\t',engine='python')

peptides = peptides[(peptides['Potential contaminant'] != '+') & (peptides['Reverse'] != '+')]

#peptides = peptides[:20] # DEV

peptide_sequences = peptides['Sequence'].tolist()

# Specific peptides for each strain
for strain in config.strains:
    paths = config.strains[strain]
    genome = list(SeqIO.parse(paths['sf_genome'],'fasta'))
    g2p = sequtils.peptides2genome(genome, assembly_name = str(strain), translation_table=config.translation_table, peptides_list=peptide_sequences, threads=config.threads)
    speps = g2p.peptides
    strainpath=output +'/' + strain
    try:
        shutil.rmtree(strainpath)
        os.mkdir(strainpath)
    except:
        os.mkdir(strainpath)
    
    speps.to_csv(strainpath + '/' + '{}_mapped_peptides.csv'.format(str(strain)))
    


