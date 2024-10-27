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
import uuid
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

def gff3_peptide_export(df, annotations, acetylations):
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
    ORF_id_full = df['ORF_id']
    try: 
        annotation = annotations[(ORF_id_full, Peptide_sequence)]
    except:
        annotation = 'globally non-specific'
    
    annotation = annotation.lower()


    acetylation = str(acetylations[Peptide_sequence]).lower()
    
    # standardize the format
    if acetylation == '-':
        acetylation = 'false'
    elif acetylation == '+':
        acetylation = 'true'
    assert acetylation in ['true', 'false']
    
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
    note='{}'.format(Peptide_sequence)
    note =urllib.parse.quote(note)
    attributes1 = 'ID={};Name={};Note={};start_type={};nterm_acetylated={}'.format(Peptide_id, Peptide_sequence, note, annotation, acetylation)

    #attributes1 = 'ID={};Name={};Note=ORF id {},%0APeptide sequence {},%0ATryptic N-terminal {},%0ATryptic C-terminal {},%0APrevious codon {},%0AFirst codon {},%0AAmino acid before {},%0AFirst amino acid {},%0ALast amino acid {},%0AAmino acid after {},%0ASpecific {},%0APeptide ORF count {},%0AIdentified in strain {},%0AMapped reference proteins {}'.format(Peptide_id, Peptide_sequence, ORF_id, Peptide_sequence, Peptide_tryptic_nterm, Peptide_tryptic_cterm, Peptide_previous_codon, Peptide_first_codon, Peptide_amino_acid_before, Peptide_amino_acid_first, Peptide_amino_acid_last, Peptide_amino_acid_after, Peptide_inferred_translated_sequence_specific, Peptide_genome_ORF_count, Strain_identified, mapped_ref)
    #attributes1 =urllib.parse.quote(attributes1)

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

def gff3_proteoform_export(df):
    global proteoform_columns
    global proteoform_domains
    global domains
    
    reference_sequence = reference(df['VarId'])
    source = 'ProteogenomicsPipeline'
    type = 'proteoform'
    proteoform_id = df['VarId'].split('|')[1]
    
    # proteform orf information

    orf_strand = df['ORF_id'].split('|')[-1].split(')')[0].split('(')[1]
    orf_start = int(df['ORF_id'].split('|')[-1].split(')')[1].split(':')[0])
    orf_end = int(df['ORF_id'].split('|')[-1].split(')')[1].split(':')[1])
    
    # proteoform information
    proteoform_name = '_'.join(proteoform_id.split('.None.'))
    strand = df['VarId'].split('|')[-1].split(')')[0].split('(')[1]
    start = int(df['VarId'].split('|')[-1].split(')')[1].split(':')[0])
    end = int(df['VarId'].split('|')[-1].split(')')[1].split(':')[1])

    score ='.'
    phase = '.'
    
    attributes = 'ID={};Name={};Note=Proteoform id {}'.format(proteoform_id  + '_orf', proteoform_name , proteoform_id)
    
    row = reference_sequence + '\t' + source + '\t' + type + '_orf' + '\t' + str(orf_start) + '\t' + str(orf_end) + '\t' + score + '\t' + orf_strand + '\t' + phase + '\t' + attributes
    
    proteoform_columns.append(row)
    
    attributes = 'Parent={};ID={};Name={};Note=Proteoform id {}'.format(proteoform_id + '_orf', proteoform_id, proteoform_name, proteoform_id)
    
    row = reference_sequence + '\t' + source + '\t' + type + '\t' + str(start) + '\t' + str(end) + '\t' + score + '\t' + strand + '\t' + phase + '\t' + attributes
    proteoform_columns.append(row)

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

def gff3_domain_export(df, domain_columns):
    global gdict
    global strain

    reference_sequence = reference(df['ProteinAccession'])
    
    accession_parents  = []

    try:
        accession2 = df['_ProteinAccession2'].split('|')[1]
    except:
        accession2 = None
    try:
        excluded = df['Excluded'] 
    except:
        excluded = False

    contig = gdict[reference_sequence]

    source = df['Analysis']
    type = 'domain'
    feature_start = (int(df['Start'])-1) * 3
    feature_end = int(df['Stop']) * 3
    ORF_id = df['ProteinAccession'].split('|')[1]

    strand = df['ProteinAccession'].split('|')[-1].split(')')[0].split('(')[1]
    start = int(df['ProteinAccession'].split('|')[-1].split(')')[1].split(':')[0])
    end = int(df['ProteinAccession'].split('|')[-1].split(')')[1].split(':')[1])

    phase = '.'
    score='.'
    try:
        Ontology_term =  df['GoAnnotations']
    except:
        Ontology_term = 'None'
    try:
        dbxref =  df['InterProAnnotationsAccession']
    except:
        dbxref = 'None'
    try:
        ipr =  df['InterProAnnotationsDescription']
    except:
        ipr = 'None'
    try:
        pathways =  df['PathwaysAnnotations']
    except:
        pathways = 'None'
    try:
        sigdesc =  df['SignatureDescription']
    except:
        sigdesc = 'None'
    try:
        name =  df['SignatureAccession']
    except:
        name = 'None'
    
    if type == 'polypeptide':
        name=ORF_id
    ID=ORF_id + '_'  + uuid.uuid4().hex 

    attributes2 ='ID={};Name={};Note={};go={};dbxref={};pathways={}'.format(ID, urllib.parse.quote(sigdesc), urllib.parse.quote(name), urllib.parse.quote(Ontology_term), dbxref, urllib.parse.quote(pathways))
    if accession2 is not None:
        attributes2 = "Parent={};excluded={};".format(accession2 + '_orf', excluded) + attributes2
        #parent = accession2 + '_domains' 
        #if not parent in accession_parents:
        #    accession_parents.append(parent)
        #    parent_row = reference_sequence +'\t' + source +'\t' + "domains_group" + '\t' + str(gstart) + '\t' + str(gend) + '\t' + score + '\t' + strand + '\t' + phase + '\t' + attributes2
    #domain_columns.append(row)

        
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
    
    row = reference_sequence +'\t' + source +'\t' + type + '\t' + str(gstart) + '\t' + str(gend) + '\t' + score + '\t' + strand + '\t' + phase + '\t' + attributes2
    domain_columns.append(row)

def export_tracks(strain, group, strain_domains, peptide_annotations, peptides, filt_peptides, peptide_acetylations,  exclude_orfs=False, exclude_domains=False):
    '''strain strain, group group ie. experimental group, peptide_annotations ,peptides, peptide_acetlations dict, exclude_orfs True/False''' 
    proteoforms = peptide_annotations[peptide_annotations['VarId'].notnull()].drop_duplicates('VarId')
    proteoforms.apply(gff3_proteoform_export, axis=1)

    peptide_annotations_dict = peptide_annotations.set_index(['ORF_id','PeptideSequence'])['StartType'].to_dict() 
    
    gff3_contig_export(genome)
    
    proteoform_domains = get_var_domains(strain_domains, peptide_annotations)
    ################
    # All Features #
    ################

    if exclude_domains == False:
        strain_domains.apply( lambda row : gff3_domain_export(row, domain_columns) , axis=1)
    
    proteoform_domains.apply( lambda row : gff3_domain_export(row, proteoform_domain_columns) , axis=1)
    
    if exclude_orfs == False:
        orfs_identified = peptides.drop_duplicates('ORF_id') # by peptides from any strain sample
        orfs_identified.apply(gff3_orf_export,axis=1)
    
    # peptides from the group
    filt_peptides.apply( lambda row : gff3_peptide_export(row, peptide_annotations_dict, peptide_acetylations) , axis=1)
    gff3_list = ['##gff-version 3']
    
    gff3_list += contig_columns
    if exclude_orfs == False:
        gff3_list += orf_columns
    gff3_list  += proteoform_columns
    gff3_list  += peptide_columns
    gff3_list  += domain_columns

    gff3 = '\n'.join(gff3_list)
    ##if strain in config['reference_strains']:
    #    gff3 = gff3 +'\n' + ref_gff3
    


    if not os.path.exists(output + '/jbrowse/' + strain):
        os.mkdir(output + '/jbrowse/' + strain) 
    w =open(output+ '/jbrowse/' + strain + '/{}_features.gff3'.format(group), 'w')
    w.write(gff3)
    w.close()


    ###########
    # domains #
    ###########
    if exclude_domains == False:
        gff3_list = ['##gff-version 3']
        
        gff3_list += contig_columns
        #gff3_list += orf_columns
        #gff3_list  += peptide_columns
        gff3_list  += domain_columns
        
        gff3 = '\n'.join(gff3_list)
        ##if strain in config['reference_strains']:
        #    gff3 = gff3 +'\n' + ref_gff3
        


        if not os.path.exists(output + '/jbrowse/' + strain):
            os.mkdir(output + '/jbrowse/' + strain) 
        w =open(output+ '/jbrowse/' + strain + '/{}_domains.gff3'.format(group), 'w')
        w.write(gff3)
        w.close()
    
    ##############
    # ORFS       #
    ##############
    if exclude_orfs == False:
        gff3_list = ['##gff-version 3']
        
        gff3_list += contig_columns
        gff3_list += orf_columns
        #gff3_list  += peptide_columns
        #gff3_list  += domain_columns
        
        gff3 = '\n'.join(gff3_list)
        ##if strain in config['reference_strains']:
        #    gff3 = gff3 +'\n' + ref_gff3
        


        if not os.path.exists(output + '/jbrowse/' + strain):
            os.mkdir(output + '/jbrowse/' + strain) 
        w =open(output+ '/jbrowse/' + strain + '/{}_orfs.gff3'.format(group), 'w')
        w.write(gff3)
        w.close()

    ##################
    # PEPTIDES       #
    ##################
    gff3_list = ['##gff-version 3']
    
    gff3_list += contig_columns
    #gff3_list += orf_columns
    gff3_list  += peptide_columns
    #gff3_list  += domain_columns
    
    gff3 = '\n'.join(gff3_list)
    ##if strain in config['reference_strains']:
    #    gff3 = gff3 +'\n' + ref_gff3

    if not os.path.exists(output + '/jbrowse/' + strain):
        os.mkdir(output + '/jbrowse/' + strain) 
    w =open(output+ '/jbrowse/' + strain + '/{}_peptides.gff3'.format(group), 'w')
    w.write(gff3)
    w.close()
    
    #####################
    # PROTEOFORMS       #
    #####################
    gff3_list = ['##gff-version 3']
    
    gff3_list += contig_columns
    
    gff3_list  += proteoform_columns
    gff3_list  += proteoform_domain_columns
    
    gff3 = '\n'.join(gff3_list)
    ##if strain in config['reference_strains']:
    #    gff3 = gff3 +'\n' + ref_gff3

    if not os.path.exists(output + '/jbrowse/' + strain):
        os.mkdir(output + '/jbrowse/' + strain) 
    w =open(output+ '/jbrowse/' + strain + '/{}_proteoforms.gff3'.format(group), 'w')
    w.write(gff3)
    w.close()


domains = sequtils.gff3(output +'/fasta/proteins.fasta.tsv').table 
with open(output + '/annotations/global_group_acetyl.json') as f:
    global_acetyl = json.loads(f.read())

with open(output + '/annotations/global_group_identified.json') as f:
    global_identified = json.loads(f.read())

#domains = pd.read_csv(output + '/fasta/proteins.fasta.tsv',sep='\t')

#print(domains.columns)
#domains = domains[domains['InterProAnnotationsAccession'].notnull()]
##domains = domains.drop_duplicates(['ProteinAccession','Start','Stop'])

def get_var_domains(domains, annotations):
    varids = annotations['VarId'].dropna().tolist()
    domains_new = []
    parents = []
    for varid in varids:
        annot = annotations[annotations['VarId'] == varid].to_dict(orient='records')[0]
        parent_id = varid.split('.None.')[0]
        if not parent_id in parents:
            parents.append(parent_id)
            

        orf_domains = domains[domains['ProteinAccession'].apply(lambda x : x.split('|')[1] == parent_id.split('|')[1])]
        try:
            orf = orf_domains['ProteinAccession'].tolist()[0]
        except:
            continue
        aa_start = annot['StartPosition']  / 3
        orf_domains['Excluded'] = orf_domains['Start'] < aa_start
        #filt_domains = orf_domains[(orf_domains['Start'] >= aa_start)  ]
        orf_domains['_ProteinAccession2'] = varid
        domains_new.append(orf_domains)
    domains_new = pd.concat(domains_new)
    domains_new['Excluded'] = domains_new['Excluded'].apply(lambda x : str(x).lower())
    return domains_new

for strain in config['strains']:
    
    ref = config['strains'][strain]['reference'][0]
    peptide_annotations =pd.read_csv(output +'/annotations/{}_{}_peptide_annotations_groups.csv'.format(ref, strain))
    peptide_annotations = peptide_annotations[peptide_annotations['Ref_id'].isnull()]
    
    
    samples = [] 
    
    for column in peptide_annotations.columns:
        if column.startswith('Identified '):
            sample = column.split('Identified ')[1]
            samples.append(sample)
    
    assert len(samples ) == len(list(set(samples))) # make sure that sample names are unique
    contig_columns = []
    orf_columns = []
    peptide_columns = []
    domain_columns = []
    proteoform_columns = []
    proteoform_domain_columns = []

    strain_domains = domains[domains['ProteinAccession'].apply(lambda x : x.startswith(strain))]
    peptides = pd.read_csv(output+ '/strains/' + strain + '/{}_mapped_peptides.csv'.format(strain))
    
    strain_identified = peptides[(peptides['Strain_identified']=='+') ]
    genome = list(SeqIO.parse(config['strains'][strain]['assembly'],'fasta'))
    gdict = {}
    
    for rec in genome:
        gdict[rec.id] = rec
    peptide_acetylations = peptides.set_index('Peptide_sequence')['Strain_Nterm_Acetylated'].to_dict() 
    export_tracks(strain, strain, strain_domains, peptide_annotations, peptides,  strain_identified, peptide_acetylations, exclude_orfs=False)
    

    for sample in samples:
        contig_columns = []
        orf_columns = []
        peptide_columns = []
        domain_columns = []
        proteoform_columns = []
        proteoform_domain_columns = []
        
        group = sample
        

        group_identified = strain_identified[strain_identified['Peptide_sequence'].isin(set(global_identified[group]))]
        group_annotations = peptide_annotations[peptide_annotations['PeptideSequence'].isin(set(global_identified[group]))]
        
        group_identified['Acetylated'] = group_identified['Peptide_sequence'].apply(lambda x : x in set(global_acetyl[group]))
        
        group_acetylations = group_identified.set_index('Peptide_sequence')['Acetylated'].to_dict()

        

        export_tracks(strain, group, strain_domains, group_annotations, peptides,  group_identified, group_acetylations, exclude_orfs=True)
