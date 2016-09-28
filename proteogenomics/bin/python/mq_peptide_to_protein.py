#!/usr/bin/env python

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


loader = importlib.machinery.SourceFileLoader('config', sys.argv[1])
config = loader.load_module()
output = sys.argv[2]

pg = pd.read_csv(config.mq_txt +'proteinGroups.txt',sep='\t')
pg =pg[(pg['Potential contaminant'] != '+') & (pg['Reverse'] != '+')]

peptides=pd.read_csv(config.mq_txt +'peptides.txt',sep='\t')
peptides=peptides[(peptides['Potential contaminant'] != '+') & (peptides['Reverse'] != '+')]

evidence=pd.read_csv(config.mq_txt +'evidence.txt',sep='\t')
evidence=evidence[(evidence['Potential contaminant'] != '+') & (evidence['Reverse'] != '+')]

reference_peptides=json.loads(open(output +'/mapping/{}_peptides.json'.format(config.reference_proteome_id)).read())

ref_entry = reference_peptides[list(reference_peptides.keys())[0]][0].split('|')[1]


taxon_peptides=json.loads(open(output +'/mapping/pep2entry.json'.format(config.reference_proteome_id)).read())
taxon_data=json.loads(open(output +'/mapping/entrydata.json'.format(config.reference_proteome_id)).read())

reference_features = sequtils.gff3(output +'/uniprot/features/{}.gff'.format(config.reference_taxid))

samples = config.samples

try:
    shutil.rmtree(output+'/clustalw')
    os.mkdir(output+'/clustalw')
    os.mkdir(output+'/clustalw/orf')
    os.mkdir(output+'/clustalw/prot')
except: 
    os.mkdir(output+'/clustalw')
    os.mkdir(output+'/clustalw/orf')
    os.mkdir(output+'/clustalw/prot')

global_non_specific_peptides=[]
global_specific_peptides=[]

strain_map={}
strain_sets={}

for strain in config.strains:
    st =  pd.read_csv(output +'/{}/{}_mapped_peptides.csv'.format(strain, strain)) 
    peptide_map = sequtils.mapping2peptides(st, config.translation_table)
    
    global_non_specific_peptides += peptide_map.non_specific()
    global_specific_peptides += peptide_map.specific()
    
    strain_map[strain]  = peptide_map
    strain_sets[strain] = set(peptide_map.mapping['Peptide_sequence'].tolist())

global_specific_peptides = set(global_specific_peptides) - set(global_non_specific_peptides)
global_non_tryptic_nterm = global_specific_peptides
global_non_tryptic_cterm = global_specific_peptides
global_non_atg_starts    = global_specific_peptides

for strain in config.strains:
    peptide_map = strain_map[strain]
    st = peptide_map.mapping
    strain_non_tryptic_nterm = set(st[st['Peptide_tryptic_nterm']=='False']['Peptide_sequence'])
    strain_non_tryptic_cterm = set(st[st['Peptide_tryptic_cterm']=='False']['Peptide_sequence'])
    strain_non_atg_starts    = peptide_map.non_atg_m()
    global_non_tryptic_nterm = global_non_tryptic_nterm & strain_non_tryptic_nterm
    global_non_tryptic_cterm = global_non_tryptic_cterm & strain_non_tryptic_cterm
    global_non_atg_starts    = global_non_atg_starts    & strain_non_atg_starts

# Annotation status
annotated_peptides = set(reference_peptides.keys())
novel_peptides = global_specific_peptides - annotated_peptides

# Get strain_exclusive
strain_sets_exclusive={}
for strain in strain_sets:
    tmp = strain_sets[strain] 
    for comparison_strain in strain_sets:
        if comparison_strain != strain:
            comp = strain_sets[comparison_strain]
            tmp = tmp - comp  # remove all peptides in the another strain
    strain_sets_exclusive[strain] = tmp - annotated_peptides # make sure all are novel

def peptide_list_blast(peptides, targets, features):
    count=1
    blasted = []
    
    feature_dict = {}
    blast_dict = {}

    for p in peptides:
        feature_dict[p] = {}
        blast_dict[p] = {}

        for target in targets:
            pairwise_blast = []
            pairwise_features = []

            id = 'peptide_{}'.format(strain, str(count))
            seq = Seq(p)
            rec = SeqRecord(id = id, seq = seq)
            out = sequtils.pairwise_blast(rec, target, output)
            
            if len(out.differences) > 0:
                header = '***Results for peptide {} in record {}***'.format(p, target.id) 
                pairwise_blast.append(header)
                pairwise_blast.append(out.results)
                pairwise_blast.append('\nIdentified polymorphism positions: {}'.format(';'.join([str(i) for i in out.differences])))
                
                overlap = out.feature_overlap(features)
                
                if overlap != None: 
                    pairwise_features.append(header)
                    pairwise_features.append(overlap)
                count += 1

                pairwise_blast = '\n\n'.join(pairwise_blast)
                pairwise_features = '\n\n'.join(pairwise_features)
                
                feature_dict[p][target.id] = pairwise_features
                blast_dict[p][target.id] = pairwise_blast

    return feature_dict, blast_dict

# Get the mapping of samples to strain
strain_dct={}
for sample in samples:
    strain_dct[sample] = samples[sample]['STRAIN']

# 1048
pg = pg[pg['id']==934] #mutant
#pg = pg[pg['id'] ==  2356] 

for row in pg.iterrows():
    print(row[0])
    peptide_ids = row[1]['Peptide IDs'].split(';')  # Get the ids of peptides in the row
    evidence_ids = row[1]['Evidence IDs'].split(';') # Get the ids of envidence.txt entries
    row_peps = peptides[peptides['id'].apply(str).isin(peptide_ids)]
    row_evs = evidence[evidence['id'].apply(str).isin(evidence_ids)]

    all_row_peptides = row_peps['Sequence'].tolist()
    row_specific  = set(all_row_peptides) - set(global_non_specific_peptides)  
    row_novel     = row_specific & novel_peptides
    row_annotated = row_specific & annotated_peptides
    row_non_tryptic_cterm = row_specific & global_non_tryptic_cterm  
    row_non_tryptic_nterm = row_specific & global_non_tryptic_nterm  
    row_non_atg_starts = row_specific & global_non_atg_starts 
    strain_peps=defaultdict(list)
    group_evs=defaultdict(list)
    
    for sample in samples:
        strain=samples[sample]['STRAIN']
        group=samples[sample]['STRAIN_GROUP']
        sample_peps = list(set(row_peps[row_peps['Experiment {}'.format(sample)]>=1]['Sequence'].tolist()))
        sample_evs = list(set(row_evs[row_evs['Experiment'] == sample]['Modified sequence'].tolist()))
        pg.loc[row[0], "_all.peptides.sample.{}".format(sample)] = '\n'.join(list(sample_evs))
        strain_peps[strain] += sample_peps
        group_evs[group] += sample_evs
    
    for group in group_evs:
        pg.loc[row[0],"_all.peptides.group.{}".format(group)]='\n'.join(list(set(group_evs[group])))
    all_orfs_fasta=[]
    all_prot_fasta=[]
    fasta_holder = {}
    id_holder={}
    frameshifts =[]
    
    strains_exclusive=[]
         
    mapped_reference = set()
    for peptide in row_annotated:
        mapped_reference.update(reference_peptides[peptide]) 
    
    refps = uniprot.map2fasta(list(mapped_reference), taxon_data)
    
    assert len(refps) == len(mapped_reference)

    reftrie  =  sequtils.list_trie_upper(refps, all_row_peptides)
    
    variant_dict = {}
    
    min_eval_ref_blast = config.eval_cutoff
    best_eval_strain = None

    row_variant_features_ref, row_variant_blast_ref = peptide_list_blast(row_specific, refps, reference_features)
    row_variant_features_ref_list = []
    row_variant_blast_ref_list = []
    
    for pep in row_variant_blast_ref:
        row_variant_blast_ref_list += list(row_variant_blast_ref[pep].values())
        if pep in row_variant_features_ref:
            row_variant_features_ref_list += list(row_variant_features_ref[pep].values()) 
    
    mapped_orfs = []
    
    for strain in strain_peps:
        peptide_mapping = strain_map[strain]
        st_peps = set(strain_peps[strain])        
        
        specific_strain_peps = st_peps & row_specific
        genome_unmapped = specific_strain_peps - strain_sets[strain]
        strain_exclusive_novel = specific_strain_peps & strain_sets_exclusive[strain]
        x,y = peptide_mapping.return_fasta(row_specific)
        all_orfs_fasta += x
        all_prot_fasta += y
        trie = sequtils.list_trie_upper(y, st_peps)
        
        fasta_holder['{}_prots'.format(strain)] = y
        fasta_holder['{}_nucs'.format(strain)] = x
        id_holder[strain] =[rec.id for rec in y]
        
        mapped_orfs += y

        for comparison_strain in strain_peps:
            if comparison_strain != strain:
                comp_peps=set(strain_peps[comparison_strain])
                strain_exclusive_novel = strain_exclusive_novel - comp_peps
        
        if len(strain_exclusive_novel) > 0:
            strains_exclusive.append(strain)
        
        pg.loc[row[0], "_all.peptides.strain.{}".format(strain)]='\n'.join(st_peps)
        pg.loc[row[0], "_specific.peptides.strain.{}".format(strain)]='\n'.join(specific_strain_peps)
        pg.loc[row[0], "_exclusive.peptides.strain.{}".format(strain)]='\n'.join(strain_exclusive_novel)        
        variant_dict[strain] = strain_exclusive_novel
        blasted = []
        features = []
        count=1
        
        for p in strain_exclusive_novel:
            if p in row_variant_blast_ref:
                blasted += list(row_variant_blast_ref[p].values())
                if p in row_variant_features_ref:
                    features += list(row_variant_features_ref[p].values())

        pg.loc[row[0], "_exclusive.peptides.mapped.reference.blast.strain.{}".format(strain)]='\n'.join(blasted) 
        pg.loc[row[0], "_identified.polymorphism.mapped.reference.feature.overlap.strain.{}".format(strain)]='\n'.join(features) 

        pg.loc[row[0], "_unmapped.peptides.strain.{}".format(strain)]='\n'.join(genome_unmapped) 
        pg.loc[row[0], "_translated.orfs.strain.{}".format(strain)] = '\n'.join(trie)
        
        ref_blast = sequtils.reference_mapping_blast(y, refps, output)
        
        if ref_blast[0] < min_eval_ref_blast:
            min_eval_ref_blast = ref_blast[0]
            best_eval_strain = strain

        pg.loc[row[0], "_translated.orfs.mapped.reference.best.blast.evalue.strain.{}".format(strain)] = ref_blast[0]
        pg.loc[row[0], "_translated.orfs.mapped.reference.best.blast.match.{}".format(strain)] = ref_blast[1]

    for key in id_holder:
        ids = id_holder[key]
        if len(ids) > 1:
            frameshifts.append(key) 
     
    pg.loc[row[0],"_frameshift"] = '\n'.join(list(frameshifts))
    pg.loc[row[0],"_exclusive.peptide.strains"] = '\n'.join(list(strains_exclusive))
    pg.loc[row[0],"_combined.specific.peptides"] ='\n'.join(list(row_specific))
    pg.loc[row[0],"_combined.specific.annotated.peptides"] ='\n'.join(row_annotated)
    pg.loc[row[0],"_combined.specific.novel.peptides"] ='\n'.join(row_novel)
    
    pg.loc[row[0],'_reference.proteins.mapped']= '\n'.join(reftrie)
    pg.loc[row[0],'_reference.proteins.mapped.count'] = len(mapped_reference)
    pg.loc[row[0],'_reference.entries.mapped'] = ';'.join([i.split('|')[1] for i in mapped_reference])
    
    pg.loc[row[0], "_reference.variant.peptides.blast"]='\n'.join(row_variant_blast_ref_list)
    pg.loc[row[0], "_reference.variant.peptides.features"]='\n'.join(row_variant_features_ref_list)
    pg.loc[row[0],'_combined.non.tryptic.nterm.peptides'] = '\n'.join(row_non_tryptic_nterm)
    pg.loc[row[0],'_combined.non.tryptic.cterm.peptides'] = '\n'.join(row_non_tryptic_cterm)
    pg.loc[row[0],'_combined.non.atg.Met.start.peptides'] = '\n'.join(row_non_atg_starts)
    
    mapped  = uniprot.peptidesmapped(row_specific, config.reference_taxonomic_lineage, taxon_peptides, taxon_data)
    
    mapped_proteome2= mapped.pep2proteome(config.proteome2)
    mapped_taxa = mapped.best
    max_pep = mapped.best_count 
    spec_fastas = mapped.best_fasta
    
    p2 = []
    for m in  mapped_proteome2:
        d = "Entry: {}; Gene names: {}; Peptides: {}".format(m[0],m[1], m[2])
        p2.append(d)
    pg.loc[row[0],'_{}.entries.mapped'.format(config.proteome2)] = '\n'.join(p2)
    


    best = None
    rmatch = True
    if len(mapped_reference) > 0:
        max_ref = mapped.count_best_mapped([i.split('|')[1] for i in mapped_reference])
        pg.loc[row[0],'_reference.proteome.best.peptide.count'] = max_ref
        
        icds = sequtils.icds_blast(refps, output)
        if len(icds) > 0: 
            pg.loc[row[0],'_reference.proteome.mapped.non.alligned'] = '\n'.join(icds)
        
        if max_pep != None:
            if (max_ref < max_pep) and (len(icds) > 0):
               best = spec_fastas[0].id
               rmatch = False

            elif len(mapped_reference) == 1:
               best = list(mapped_reference)[0].split('|')[1]

            else:
                best=None
                rmatch = False
         
    

    pg.loc[row[0],'_taxon.best.matches'] = mapped_taxa 
    pg.loc[row[0],'_taxon.best.peptide.count'] = max_pep
    
    if spec_fastas != None:
        pg.loc[row[0],'_taxon.best.match.fasta'] = '\n'.join(sequtils.list_trie_upper([spec_fastas[0]],all_row_peptides)) 
        all_prot_fasta += spec_fastas

    if len(all_prot_fasta) > 0:
        all_prot_fasta = sequtils.remove_asterisk(all_prot_fasta)
        prot_aln,prot_dend=sequtils.clustalw(output+'/clustalw/prot/{}.fasta'.format(row[0]), all_prot_fasta) 
        
        pg.loc[row[0],'_clustalw.proteins'] = str(prot_aln)
        pg.loc[row[0],'_clustalw.proteins.newick'] = prot_dend
    
        print(prot_dend)
        prot_muscle=sequtils.muscle(all_prot_fasta)
        pg.loc[row[0],'_muscle.proteins'] = str(prot_muscle)

    if len(all_orfs_fasta) > 0:
        nuc_aln, nuc_dend = sequtils.clustalw(output+'/clustalw/orf/{}.fasta'.format(row[0]), all_orfs_fasta)
        pg.loc[row[0],'_clustalw.orfs'] = str(nuc_aln)
        pg.loc[row[0],'_clustalw.orfs.newick'] = nuc_dend 
        
        print(nuc_dend)
        nuc_muscle=sequtils.muscle(all_orfs_fasta)
        pg.loc[row[0],'_muscle.orfs'] = str(nuc_muscle)
    
    identifier = []
    

    if best_eval_strain  != None:
        pg.loc[row[0],'_strain.mapped.reference.blast.best.match.strain'] = best_eval_strain
        if rmatch == True:
            pg.loc[row[0],'_protein.group.reference.match'] = '+'

    if best != None:
        data = taxon_data[best]
        for key in data:
            pg.loc[row[0],key] = data[key]
        identifier.append(data['Gene names'])
    
    identifier.append('(protein group {})'.format(str(row[0])))

    pg.loc[row[0], 'Identifier'] = ' '.join(identifier)

pg.to_csv(output+'/combined.csv')



