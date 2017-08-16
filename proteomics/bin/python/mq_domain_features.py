#!/usr/bin/env python

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

loader = importlib.machinery.SourceFileLoader('config', sys.argv[1])
config = loader.load_module()
output = sys.argv[2]

idmapping = pickle.load(open(output+'/fasta/id_mapping.p','rb'))
gff3 = sequtils.gff3(output +'/fasta/nr_translated_pg_orfs.fasta.gff3')
gff3.expand_table(idmapping)
print(gff3.table['attributes'].values)
gff3.table['strain'] = gff3.table['_mapped.id'].apply( lambda x : x.split('_')[0])

def gff3_domain_export(df):
    global gdict
    global domain_columns
    global strain

    strain = df['seqid'].split('_')[0]
    reference_sequence = df['seqid'].split('_')[1]
    contig = gdict[reference_sequence]

    source = df['source']
    type = df['type']

    feature_start = (int(df['start'])-1) * 3
    feature_end = int(df['end']) * 3
    ORF_id = df['_mapped.id'].split('|')[1]

    strand = df['_mapped.id'].split('|')[-1].split(')')[0].split('(')[1]
    start = int(df['_mapped.id'].split('|')[-1].split(')')[1].split(':')[0])
    end = int(df['_mapped.id'].split('|')[-1].split(')')[1].split(':')[1])

    score = df['score']
    phase = df['phase']
    
    try:
        Ontology_term =  df['attributes'].split('Ontology_term=')[1].split(';')[0]
    except:
        Ontology_term = 'None'
    try:
        dbxref =  df['attributes'].split('Dbxref=')[1].split(';')[0]
    except:
        dbxref = 'None'
    try:
        sigdesc =  df['attributes'].split('signature_desc=')[1].split(';')[0]
    except:
        sigdesc = 'None'
    try:
        name =  df['attributes'].split('Name=')[1].split(';')[0]
    except:
        name = 'None'
    
    if type == 'polypeptide':
        name=ORF_id

    #attributes1 = 'ID=domain-{};Name={};Note={};Ontology_term={};dbxref={}'.format(ORF_id, name, sigdesc, Ontology_term, dbxref)
    attributes2 = 'Name={};Note={};Ontology_term={};dbxref={}'.format(name, sigdesc, Ontology_term, dbxref)
    
    if strand == '+':
        gstart = start + feature_start
        gend   = start + feature_end - 1 
        nucs = str(contig.seq)[gstart-1:gend]
        assert len(str(nucs)) % 3 == 0
    
    if strand == '-':
        gend = end - feature_start
        gstart = end - feature_end + 1
        nucs = Seq(str(contig.seq)[gstart-1:gend]).reverse_complement()
        assert len(str(nucs)) % 3 == 0

    #row = reference_sequence +'\t' + 'match' +'\t' + type + '\t' + str(gstart) + '\t' + str(gend) + '\t' + score + '\t' + strand + '\t' + phase + '\t' + attributes1
    #domain_columns.append(row)
    
    row = reference_sequence +'\t' + source +'\t' + type + '\t' + str(gstart) + '\t' + str(gend) + '\t' + score + '\t' + strand + '\t' + phase + '\t' + attributes2
    domain_columns.append(row)

def gff3_orf_export(df):
    global orf_columns
    reference_sequence = df['ORF_id'].split('_')[1]
    source = 'ProteogenomicsPipeline'
    type = 'open_reading_frame'
    recno = 'recno_' + df['ORF_id'].split('|')[1].split('recno_')[1]
    ORF_id = df['ORF_id'].split('|')[1]

    ORF_sequence = df['ORF_sequence']
    ORF_translation = df['ORF_translation']
    ORF_preceding_codon = df['ORF_preceding_codons']
    ORF_Cterm_fragment = df['ORF_Cterm_fragment']
    ORF_Nterm_fragment = df['ORF_Nterm_fragment']
    Assembly=df['Assembly']
    Most_Upstream_Inferred_Start = str(int(df['Most_Upstream_Inferred_Start']))
    Most_Upstream_Inferred_Codon = df['Most_Upstream_Inferred_Codon']
    Most_Upstream_Inferred_Translation = df['Most_Upstream_Inferred_Translation']

    strand = df['ORF_id'].split('|')[-1].split(')')[0].split('(')[1]
    start = int(df['ORF_id'].split('|')[-1].split(')')[1].split(':')[0])
    end = int(df['ORF_id'].split('|')[-1].split(')')[1].split(':')[1])

    score ='.'
    phase = '.'
    #attributes = 'Name={};Note=ORF id: {}%0AAssembly: {}%0AMost upstream inferred codon: {}%0AMost upstream inferred start: {}%0APreceding codon: {}%0AN-term fragment: {}%0AC-term fragment: {}'.format(recno, ORF_id, Assembly, Most_Upstream_Inferred_Codon, Most_Upstream_Inferred_Start, ORF_preceding_codon, ORF_Nterm_fragment, ORF_Cterm_fragment)
    
    attributes = 'Name={};Note=ORF id {}'.format(recno, ORF_id)
    
    row = reference_sequence + '\t' + source + '\t' + type + '\t' + str(start) + '\t' + str(end) + '\t' + score + '\t' + strand + '\t' + phase + '\t' + attributes
    orf_columns.append(row)

def gff3_contig_export(recs):
    global contig_columns
    global strain
    for rec in recs:
        id = rec.id
        reference_sequence = id.split('|')[0]
        size = id.split('|')[1]
        start = '1'
        end = str(len(str(rec.seq)))
        strand = '.'
        source = 'Assembly'
        type='contig'
        score ='.'
        strand='.'
        phase='.'
        attributes='Name={}'.format(reference_sequence)
        row = reference_sequence + '\t' + source + '\t' + type + '\t' + start + '\t' + end + '\t' + score + '\t' + strand + '\t' + phase + '\t' + attributes
        contig_columns.append(row)


for strain in config.strains:
    contig_columns = []
    domain_columns = []
    
    table = gff3.table
    table = table[table['strain'] == strain]

    genome = list(SeqIO.parse(config.strains[strain]['sf_genome'],'fasta'))
    gff3_contig_export(genome)

    gdict = {}
    for rec in genome:
        id = rec.id.split('|')[0]
        gdict[id] = rec
    
    table.apply(gff3_domain_export,axis=1)
    
    gff3_list = ['##gff-version 3\n##Index-subfeatures 1\n']
    
    gff3_list += contig_columns
    #gff3_list += orf_columns
    gff3_list  += domain_columns
    
    gf = '\n'.join(gff3_list)
    
    w =open(output+ '/' + strain + '/{}_domains.gff3'.format(strain), 'w')
    w.write(gf)
    w.close()


    
