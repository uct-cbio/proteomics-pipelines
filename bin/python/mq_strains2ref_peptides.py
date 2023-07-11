#!/usr/bin/env python3

import pandas as pd
import sys
import multiprocessing
import numpy as np
import pandas as pd
import sys
import importlib.machinery
import Bio; from Bio import SeqIO
import sequtils
import shutil
import algo
import os
import subprocess
from collections import defaultdict
import json
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Seq import translate
import sys
from collections import Counter
from io import StringIO
import uniprot
import pickle
from io import StringIO
import yaml
import blast
from Bio.Blast import NCBIXML
import gff3

config = yaml.load(open(sys.argv[1]).read(), Loader = yaml.Loader)
output = sys.argv[2]
strain_samples = defaultdict(list)

design = pd.read_csv(config['design'])

if 'rename' in design.columns:
    scol = 'rename'
else:
    scol = 'sample'
assert scol in design.columns
design = design[design['exclude'] != '+' ] 
sample2strain = design.set_index(scol)['Strain'].to_dict()
print(sample2strain)
samples = list(sample2strain.keys())

for sample in samples:
    st = sample2strain[sample]
    strain_samples[st].append(sample)

peptide_sequence  = SeqIO.to_dict(list(SeqIO.parse(output + '/blast/peptides2orfs/peptides.fasta','fasta')))

# Per strain analysis
peptides = pd.read_csv(config['mq_txt'] + '/peptides.txt', sep='\t',engine='python')
peptides = peptides[(peptides['Potential contaminant'].isnull()) & (peptides['Reverse'].isnull())]
all_peptides = peptides['Sequence'].tolist()
exp_cols = [i for i in peptides.columns if i.startswith('Experiment')]
#blast_name= config['reference_genome'].split('/')[-1].split('.fasta')[0]

def _(blast_name):
    orf_sequence  = SeqIO.to_dict(list(SeqIO.parse(output + '/blast/peptides2orfs/{}_reforfs.fasta'.format(blast_name),'fasta')))
    blast_results = output + '/blast/peptides2orfs/{}.xml'.format(blast_name)
    result_handle = open(blast_results)
    blast_records = NCBIXML.parse(result_handle)
    print("PARSED BLAST RESULTS FOR REFERENCE GENOMES")
    mapdict = {}
    recdict = defaultdict(list)
    recdict = gff3.orfblastp2gff3(blast_records=blast_records, query_fasta=peptide_sequence, orf_fasta=orf_sequence,  eval_cutoff=200000, subject_cutoff=1)

    ref_genome =list(SeqIO.parse(output + '/ena/{}/{}.fasta'.format(blast_name, blast_name),'fasta'))
    contigs = gff3.gff3_contig_export(ref_genome)


    for strain in config['strains']:
        samples = strain_samples[strain]
        sample_columns = []
        
        for sample in samples:
            sample_columns.append('Experiment {}'.format(sample))
        
        strain_filt = peptides[peptides[sample_columns].sum(axis=1) >=1]
        strain_peptides = list(set(strain_filt['Sequence'].tolist()))
        strain_rows = ['##gff-version 3\n##Index-subfeatures 1\n']
        #strain_rows += contigs
        for peptide in strain_peptides:
            if peptide in recdict:
                results = recdict[peptide]
                strain_rows += results
        gff_str = '\n'.join(strain_rows)

        if not os.path.exists(output + '/jbrowse/' + strain):
            os.mkdir(output + '/jbrowse/' + strain)     

        # Write the peptide features 
        w =open(output+ '/jbrowse/' + strain + '/{}_{}_peptides.gff3'.format(strain,blast_name), 'w')
        w.write(gff_str)
        w.close()

for ref in config['reference']:
    assembly_id = config['reference'][ref]['assembly_id']
    _(assembly_id)
