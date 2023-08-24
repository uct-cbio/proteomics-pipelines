#!/usr/bin/env python3

import pandas as pd
import sys
import importlib.machinery
import Bio; from Bio import SeqIO
import sequtils
import shutil
import os
from collections import defaultdict
import pickle
import yaml
import mqparse


config = yaml.load(open(sys.argv[1]), Loader=yaml.Loader)

output = sys.argv[2]



output = os.path.abspath(output)
design = os.path.abspath(config['design'])

design = pd.read_csv(design)
design = design[design['exclude'] != '+']
assert "Strain" in design.columns

os.mkdir(output +'/strains')

peptides = pd.read_csv(config['mq_txt'] + '/peptides.txt', sep='\t',engine='python')
evidence = pd.read_csv(config['mq_txt'] + '/evidence.txt', sep='\t',engine='python')
Evidence = mqparse.Evidence(evidence)
print("Loaded Evidence")

peptides = peptides[(peptides['Potential contaminant'].isnull()) & (peptides['Reverse'].isnull())]

all_peptides = peptides['Sequence'].tolist()

strain_samples = defaultdict(list)

samples = design['sample'].tolist()

for row in design.iterrows():
    sample  = row[1]['sample']
    strain  = row[1]['Strain']
    #for sample in samples:
    #st = samples[sample]['STRAIN']
    strain_samples[strain].append(sample)

exp_cols = [i for i in peptides.columns if i.startswith('Experiment')]
#peptides = peptides[:20] # DEV
#peptide_sequences = peptides['Sequence'].tolist()
#assert 4 == 5
# Specific peptides for each strain

def check_identified(peptide, identified):
    if peptide in identified:
        return '+'
    else:
        return '-'

all_orfs = []
all_trans_orfs = []
for strain in config['strains']:
    samples = strain_samples[strain]
    sample_columns = []

    for sample in samples:
        sample_columns.append('Experiment {}'.format(sample))
    strain_filt = peptides[peptides[sample_columns].sum(axis=1) >=1]
    strain_peptides = strain_filt['Sequence'].tolist()
    paths = config['strains'][strain]
    strain_acetylated = Evidence.export_peptides(samples,['_(ac)'])
    strain_m_ox = Evidence.export_peptides(samples,['M(ox)'])
    print('Strain acetylated', len(strain_acetylated))
    print('Strain M(ox)', len(strain_m_ox))

    if paths['assembly'] != None:
        print(strain)
        genome = list(SeqIO.parse(paths['assembly'],'fasta'))
        
        strainpath=output +'/strains/' + strain
        try:
            shutil.rmtree(strainpath)
            os.mkdir(strainpath)
        except:
            os.mkdir(strainpath)    

        g2p = sequtils.peptides2genome(genome, assembly_name = str(strain), translation_table=config['translation_table'], peptides_list=all_peptides, outdir=strainpath, threads=config['threads'])
        
        speps = g2p.peptides
        
        speps['Strain_identified'] = speps['Peptide_sequence'].apply(lambda x : check_identified(x, strain_peptides))
        speps['Strain_Nterm_Acetylated'] = speps['Peptide_sequence'].apply(lambda x : check_identified(x, strain_acetylated))
        speps['Strain_M_Oxidation'] = speps['Peptide_sequence'].apply(lambda x : check_identified(x, strain_m_ox))
        
        #outpath=strainpath + '/' + '{}_mapped_peptides.p'.format(str(strain))   
        
        speps.to_csv(strainpath + '/' + '{}_mapped_peptides.csv'.format(str(strain)))
        SeqIO.write(g2p.mapped_orfs, strainpath + '/' + '{}_mapped_orfs.fasta'.format(str(strain)), 'fasta')
        SeqIO.write(g2p.mapped_trans_orfs, strainpath + '/' + '{}_mapped_transated_orfs.fasta'.format(str(strain)), 'fasta')
        all_orfs += g2p.mapped_orfs
        all_trans_orfs += g2p.mapped_trans_orfs
        #pickle.dump( speps, open( outpath, "wb" ) )
    

SeqIO.write(all_orfs, output + '/strains/all_mapped_orfs.fasta', 'fasta')
SeqIO.write(all_trans_orfs, output + '/strains/all_mapped_trans_orfs.fasta', 'fasta')
