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
import collections
import urllib

config = yaml.load(open(sys.argv[1]).read(), Loader=yaml.Loader)
output = sys.argv[2]



    #ref_gff3 = open('/root/features.gff3').readlines()
    #ref_gff3 = [i for i in ref_gff3 if not i.startswith('#')]
    #ref_gff3 = '\n'.join(ref_gff3)
    #ref_gff3 = 'Note='.join(ref_gff3.split('description='))
    #ref_gff3 = 'Note='.join(ref_gff3.split('protein_id='))

def reference(val):
    reference_sequence = '|'.join(val.split('|')[0].split('_')[-2:])
    return reference_sequence

def gff3_peptide_export(df):
    global gdict
    global peptide_columns
    reference_sequence = reference(df['ORF_id'])
    #reference_sequence = '|'.join(df['ORF_id'].split('|')[0].split('_')[-2:])
    contig = gdict[reference_sequence]
    
    source = 'MaxQuant'
    type = 'polypeptide'

    Peptide_starts = [(int(i)-1) * 3 for i in str(df['Peptide_starts']).split(';')]

    ORF_translation = df['ORF_translation']
    Peptide_sequence = df['Peptide_sequence']
    ORF_id = df['ORF_id'].split('|')[1]
    Peptide_id = ORF_id + '_' + Peptide_sequence
    #recno = 'recno_' + df['ORF_id'].split('|')[1].split('recno_')[1]

    ORF_sequence = df['ORF_sequence']
    Peptide_first_codon = df['Peptide_first_codon']

    strand = df['ORF_id'].split('|')[-1].split(')')[0].split('(')[1]
    start = int(df['ORF_id'].split('|')[-1].split(')')[1].split(':')[0])
    end = int(df['ORF_id'].split('|')[-1].split(')')[1].split(':')[1])

    Peptide_ends = [ (i + (len(Peptide_sequence) * 3)) for i in Peptide_starts] 
    
    score ='.'
    phase = '.'
    Peptide_tryptic_nterm = df['Peptide_tryptic_nterm']
    Peptide_tryptic_cterm = df['Peptide_tryptic_cterm']
    Peptide_previous_codon = df['Peptide_previous_codon']
    Peptide_first_codon = df['Peptide_first_codon']
    Peptide_amino_acid_before = df['Peptide_amino_acid_before']
    Peptide_amino_acid_first = df['Peptide_amino_acid_first']
    Peptide_amino_acid_last = df['Peptide_amino_acid_last']
    Peptide_amino_acid_after = df['Peptide_amino_acid_after']
    Peptide_inferred_translated_sequence_count = df['Peptide_inferred_translated_sequence_count']
    Peptide_inferred_translated_sequence_specific = df['Peptide_inferred_translated_sequence_specific']
    Peptide_genome_ORF_count = df['Peptide_genome_ORF_count']
    Strain_identified = df['Strain_identified']
    #mapped_up = ' '.join(pep2entry[Peptide_sequence]) 
    
    # mapping a peptide to a reference sequence should happen here
    #try:
    #    mapped_ref = ' '.join(pep2reference[Peptide_sequence]) 
    #except:
    
    mapped_ref= 'None'

    #attributes1 = 'ID=peptide-{};Name={};Note=ORF id {},%0APeptide sequence {}'.format(ORF_id, ORF_id, ORF_id, Peptide_sequence)
    note='ORF id: {},\nPeptide sequence: {}'.format(ORF_id, Peptide_sequence)
    note =urllib.parse.quote(note)
    attributes1 = 'ID={};Name={};Note={}'.format(Peptide_id, Peptide_sequence, note)
    #attributes1 = 'ID={};Name={};Note=ORF id {},%0APeptide sequence {},%0ATryptic N-terminal {},%0ATryptic C-terminal {},%0APrevious codon {},%0AFirst codon {},%0AAmino acid before {},%0AFirst amino acid {},%0ALast amino acid {},%0AAmino acid after {},%0ASpecific {},%0APeptide ORF count {},%0AIdentified in strain {},%0AMapped reference proteins {}'.format(Peptide_id, Peptide_sequence, ORF_id, Peptide_sequence, Peptide_tryptic_nterm, Peptide_tryptic_cterm, Peptide_previous_codon, Peptide_first_codon, Peptide_amino_acid_before, Peptide_amino_acid_first, Peptide_amino_acid_last, Peptide_amino_acid_after, Peptide_inferred_translated_sequence_specific, Peptide_genome_ORF_count, Strain_identified, mapped_ref)
    #attributes1 =urllib.parse.quote(attributes1)
    attributes2 = 'Parent={};Name={}'.format(ORF_id, Peptide_id)

    #attributes2 =urllib.parse.quote(attributes2)
    #attributes2 = 'Parent={};Name={};Note=ORF id {},%0APeptide sequence {},%0ATryptic N-terminal {},%0ATryptic C-terminal {},%0APrevious codon {},%0AFirst codon {},%0AAmino acid before {},%0AFirst amino acid {},%0ALast amino acid {},%0AAmino acid after {},%0ASpecific {},%0APeptide ORF count {},%0AIdentified in strain {},%0AMapped reference proteins {}'.format(ORF_id, Peptide_sequence, ORF_id, Peptide_sequence, Peptide_tryptic_nterm, Peptide_tryptic_cterm, Peptide_previous_codon, Peptide_first_codon, Peptide_amino_acid_before, Peptide_amino_acid_first, Peptide_amino_acid_last, Peptide_amino_acid_after, Peptide_inferred_translated_sequence_specific, Peptide_genome_ORF_count, Strain_identified, mapped_ref)
    #attributes2 = 'Name={};Note=ORF id {},%0APeptide sequence {},%0ATryptic N-terminal {},%0ATryptic C-terminal {},%0APrevious codon {},%0AFirst codon {},%0AAmino acid before {},%0AFirst amino acid {},%0ALast amino acid {},%0AAmino acid after {},%0ASpecific {},%0APeptide ORF count {},%0AIdentified in strain {},%0AMapped reference proteins (Proteome ID {}) {}'.format(Peptide_sequence, ORF_id, Peptide_sequence, Peptide_tryptic_nterm, Peptide_tryptic_cterm, Peptide_previous_codon, Peptide_first_codon, Peptide_amino_acid_before, Peptide_amino_acid_first, Peptide_amino_acid_last, Peptide_amino_acid_after, Peptide_inferred_translated_sequence_specific, Peptide_genome_ORF_count, Strain_identified, config['reference_proteome_id'], mapped_ref)
    #attributes1 = 'ID=peptide-{};Name={};Note=ORF id {};Peptide sequence {};Mapped reference proteins ({}) {}, Mapped taxon proteins (TaxID {}) {}'.format(ORF_id, Peptide_sequence, ORF_id, Peptide_sequence, config.reference_proteome_id, mapped_ref, config.group_taxid, mapped_up)
    
    for i in range(len(Peptide_starts)):
        pstart = Peptide_starts[i]
        pend   = Peptide_ends[i]
        
        if strand == '+':
            gstart = start + pstart
            gend   = start + pend - 1 
            nucs = str(contig.seq)[gstart-1:gend]
            assert len(str(nucs)) % 3 == 0
            trans = str(translate(Seq(nucs), table=11, cds=False)) 
            assert trans[1:] == Peptide_sequence[1:]
        
        if strand == '-':
            gend = end - pstart
            gstart = end - pend + 1
            nucs = Seq(str(contig.seq)[gstart-1:gend]).reverse_complement()
            assert len(str(nucs)) % 3 == 0
            trans = str(translate(nucs, table=11, cds=False)) 
            assert trans[1:] == Peptide_sequence[1:]

        row = reference_sequence +'\t' + 'match' +'\t' + type  + '\t' + str(gstart) + '\t' + str(gend) + '\t' + score + '\t' + strand + '\t' + phase + '\t' + attributes1
        peptide_columns.append(row)
        
        #row = reference_sequence +'\t' + source +'\t' + 'open_reading_frame' + '\t' + str(gstart) + '\t' + str(gend) + '\t' + score + '\t' + strand + '\t' + phase + '\t' + attributes2
        #peptide_columns.append(row)


def gff3_orf_export(df):
    global orf_columns
    #reference_sequence = df['ORF_id'].split('_')[-3]
    reference_sequence = reference(df['ORF_id'])
    #print(reference_sequence)
    source = 'ProteogenomicsPipeline'
    type = 'open_reading_frame'
    #recno = 'recno_' + df['ORF_id'].split('|')[1].split('recno_')[1]
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
    
    attributes = 'ID={};Name={};Note=ORF id {}'.format(ORF_id, ORF_id, ORF_id)
    
    row = reference_sequence + '\t' + source + '\t' + type + '\t' + str(start) + '\t' + str(end) + '\t' + score + '\t' + strand + '\t' + phase + '\t' + attributes
    orf_columns.append(row)

def gff3_contig_export(recs):
    global contig_columns
    global strain
    for rec in recs:
        id = rec.id
        #print(id)
        reference_sequence = id
        #size = id.split('|')[1]
        size = len(str(rec.seq))
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
        seq_reg="##sequence-region {} {} {}".format(reference_sequence, start, end)

        contig_columns.append(seq_reg)
        #contig_columns.append(row)

for strain in config['strains']:
    contig_columns = []
    orf_columns = []
    peptide_columns = []
    
    peptides = pd.read_csv(output+ '/strains/' + strain + '/{}_mapped_peptides.csv'.format(strain))
    strain_identified = peptides[(peptides['Strain_identified']=='+') ]
    genome = list(SeqIO.parse(config['strains'][strain]['assembly'],'fasta'))
    
    gff3_contig_export(genome)

    gdict = {}
    for rec in genome:
        gdict[rec.id] = rec
    
    #strain_identified = strain_identified.head()
    orfs_identified = strain_identified.drop_duplicates('ORF_id')
    orfs_identified.apply(gff3_orf_export,axis=1)
    strain_identified.apply(gff3_peptide_export, axis=1)
    gff3_list = ['##gff-version 3']
    
    gff3_list += contig_columns
    gff3_list += orf_columns
    gff3_list  += peptide_columns
    
    gff3 = '\n'.join(gff3_list)
    ##if strain in config['reference_strains']:
    #    gff3 = gff3 +'\n' + ref_gff3

    if not os.path.exists(output + '/jbrowse/' + strain):
        os.mkdir(output + '/jbrowse/' + strain) 
    w =open(output+ '/jbrowse/' + strain + '/{}_features.gff3'.format(strain), 'w')
    w.write(gff3)
    w.close()


    
