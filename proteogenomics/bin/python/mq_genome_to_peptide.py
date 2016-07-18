#!/usr/bin/env python

import pandas as pd
import sys
import importlib.machinery
import Bio; from Bio import SeqIO
import sequtils

loader = importlib.machinery.SourceFileLoader('config', sys.argv[1])
config = loader.load_module()

peptides = pd.read_csv(config.mq_txt + 'peptides.txt', sep='\t',engine='python')
peptides = peptides[(peptides['Potential contaminant'] != '+') & (peptides['Reverse'] != '+')]

peptides = peptides[:100] # DEV

peptide_sequences = peptides['Sequence'].tolist()

for strain in config.strains:
    paths = config.strains[strain]
    
    genome = list(SeqIO.parse(paths['sf_genome'],'fasta'))
    orfs = sequtils.gssp(genome, assembly_name = strain, translation_table=config.translation_table, peptides_list=peptide_sequences)
    print(orfs.ORF_set_count())
    
    print()
    break

