#!/usr/bin/env python

import pandas as pd
import sys
import importlib.machinery
import Bio; from Bio import SeqIO
import sequtils
import shutil
import os
from collections import defaultdict
import json
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio import Phylo
import contextlib
import sys
from io import StringIO

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

strain_dct={}
samples = config.samples

try:
    shutil.rmtree(output+'/clustalw')
    os.mkdir(output+'/clustalw')
    os.mkdir(output+'/clustalw/orf')
    os.mkdir(output+'/clustalw/prot')
except:
    os.mkdir(output+'/clustalw/orf')
    os.mkdir(output+'/clustalw/prot')

def get_ns(map):
    ns =map[map['Peptide_distinct_translated_ORF_specfic'] =='-']['Peptide_sequence'].tolist()
    ns = list(set(ns))
    return ns

def get_s(map):
    s =map[map['Peptide_distinct_translated_ORF_specfic'] =='+']['Peptide_sequence'].tolist()
    s = list(set(s))
    return s

# Get the mapping of the genomic locations of identified peptides for each strain
non_specific_peptides=[]
specific_peptides=[]
strain_map={}
strain_sets={}
for strain in config.strains:
    st =  pd.read_csv(output +'/{}/{}_mapped_peptides.csv'.format(strain, strain))
    non_specific_peptides += get_ns(st)
    specific_peptides += get_s(st)
    strain_map[strain] = st
    strain_sets[strain] = set(st['Peptide_sequence'].tolist())


# Get the mapping of genomics locations of identified peptides for reference
reference_map = pd.read_csv(output +'/reference/{}_mapped_peptides.csv'.format(config.reference_taxid))
reference_nontryptic_nterm=set(reference_map[reference_map['Peptide_tryptic_nterm']=='False']['Peptide_sequence'])
reference_nontryptic_cterm=set(reference_map[reference_map['Peptide_tryptic_cterm']=='False']['Peptide_sequence'])

non_specific_peptides += get_ns(reference_map)
reference_set = set(reference_map['Peptide_sequence'].tolist())
specific_peptides += get_s(reference_map)

# Combined specific set (exclide all sequences that are non-specific in any group
global_specific_peptides = set(specific_peptides) - set(non_specific_peptides)
annotated_peptides = set(reference_peptides.keys())
exclusive_peptides = annotated_peptides - reference_set # in reference proteome but not reference genome ie. annotated frameshift proteins from swissprot
novel_peptides = global_specific_peptides - annotated_peptides

# Get strain_exclusive
strain_sets_exclusive={}
for strain in strain_sets:
    tmp = strain_sets[strain] - reference_set  # remove all peptides in reference
    for comparison_strain in strain_sets:
        if comparison_strain != strain:
            comp = strain_sets[comparison_strain]
            tmp = tmp - comp  # remove all peptides in another strain
    strain_sets_exclusive[strain] = tmp - annotated_peptides # make sure all are novel


# Get the mapping of samples to strain
for sample in samples:
    strain_dct[sample] = samples[sample]['STRAIN']

pep2prot=pd.DataFrame()

pg = pg[pg['id']==934]

def return_fasta(map, specific):
    orfs_fasta=[]
    prots_fasta=[]
    specific_map=map[map['Peptide_sequence'].isin(specific)].copy()
    
    orf_id=specific_map['ORF_ids'].tolist()[0]
    orf_nucs=specific_map['ORF_sequence'].tolist()[0]
    orf_trans=specific_map['ORF_translation'].tolist()[0]
    nucs=SeqRecord(id=orf_id,seq=Seq(orf_nucs))
    trans=SeqRecord(id=orf_id,seq=Seq(orf_trans))
    orfs_fasta.append(nucs)
    prots_fasta.append(trans)
    return orfs_fasta, prots_fasta

@contextlib.contextmanager
def stdout_redirect(where):
    sys.stdout = where
    try:
        yield where
    finally:
        sys.stdout = sys.__stdout__

def clustalw(output_file, fasta):
    SeqIO.write(fasta, output_file, 'fasta')
    cline = ClustalwCommandline("clustalw2", infile=output_file)
    stdout, stderr = cline()
    aln = output_file.split('.fasta')[0] 
    align = AlignIO.read(aln +'.aln', "clustal")
    lines=[]
    for record in align:
        ln="{} - {}".format(record.seq, record.id)
        lines.append(ln)
    lines = '\n'.join(lines)
    tree = Phylo.read(aln +".dnd", "newick") 
    with stdout_redirect(StringIO()) as new_stdout:
        Phylo.draw_ascii(tree)
    new_stdout.seek(0)
    tree = new_stdout.read()
    return lines, tree

for row in pg.iterrows():
    peptide_ids = row[1]['Peptide IDs'].split(';')  # Get the ids of peptides in the row
    evidence_ids = row[1]['Evidence IDs'].split(';') # Get the ids of envidence.txt entries
    row_peps = peptides[peptides['id'].apply(str).isin(peptide_ids)]
    row_evs = evidence[evidence['id'].apply(str).isin(evidence_ids)]
    # Divide identified peptides into categories
    row_specific  = set(row_peps['Sequence'].tolist()) - set(non_specific_peptides)  # All peptides that are not non-specific in any genome
    row_exclusive = row_specific & exclusive_peptides # All annotated peptides that are not mapped to the reference genome (excl non_specific)
    row_novel     = row_specific & novel_peptides
    row_annotated = row_specific & annotated_peptides
    row_nontryptic_cterm = row_specific & reference_nontryptic_cterm  
    row_nontryptic_nterm = row_specific & reference_nontryptic_nterm  

    strain_peps=defaultdict(list)
    group_evs=defaultdict(list)
    for sample in samples:
        strain=samples[sample]['STRAIN']
        group=samples[sample]['GROUP']
        sample_peps = list(set(row_peps[row_peps['Experiment {}'.format(sample)]>=1]['Sequence'].tolist()))
        sample_evs = list(set(row_evs[row_evs['Experiment'] == sample]['Modified sequence'].tolist()))
        pg.loc[row[0], "_all.peptides.sample.{}".format(sample)] = '\n'.join(list(sample_evs))
        strain_peps[strain] += sample_peps
        group_evs[group] += sample_evs
    for group in group_evs:
        pg.loc[row[0],"_all.peptides.group.{}".format(group)]='\n'.join(list(set(group_evs[group])))
    all_orfs_fasta=[]
    all_prot_fasta=[]
    
    # Add reference fasta sequences

    x,y=return_fasta(reference_map, row_specific)
    all_orfs_fasta += x
    all_prot_fasta += y

    strains_exclusive=[]
    for strain in strain_peps:
        peptide_mapping = strain_map[strain]
        st_peps = set(strain_peps[strain])        
        specific_strain_peps = st_peps & row_specific
        genome_unmapped = specific_strain_peps - strain_sets[strain]
        strain_exclusive_novel = specific_strain_peps & strain_sets_exclusive[strain]
        
        x,y=return_fasta(peptide_mapping, specific_strain_peps)
        all_orfs_fasta += x
        all_prot_fasta += y
        
        for comparison_strain in strain_peps:
            if comparison_strain != strain:
                comp_peps=set(strain_peps[comparison_strain])
                strain_exclusive_novel = strain_exclusive_novel - comp_peps
        if len(strain_exclusive_novel) > 0:
            strains_exclusive.append(strain)
        pg.loc[row[0], "_all.peptides.strain.{}".format(strain)]='\n'.join(st_peps)
        pg.loc[row[0], "_exclusive.gssp.peptides.strain.{}".format(strain)]='\n'.join(strain_exclusive_novel)
        pg.loc[row[0], "_genome.unmapped.peptides.strain.{}".format(strain)]='\n'.join(genome_unmapped) 
    
    pg.loc[row[0],"_exclusive.peptide.strains"] = '\n'.join(list(strains_exclusive))
    pg.loc[row[0],"_combined.specific.peptides"] ='\n'.join(list(row_specific))
    pg.loc[row[0],"_combined.specific.annotated.peptides"] ='\n'.join(row_annotated)
    pg.loc[row[0],"_combined.specific.novel.peptides"] ='\n'.join(row_novel)
    mapped_reference = set()
    for peptide in row_annotated:
        mapped_reference.update(reference_peptides[peptide])
    pg.loc[row[0],'_reference.proteins.mapped']='\n'.join(mapped_reference)
    pg.loc[row[0],'_reference.proteins.mapped.count'] = len(mapped_reference)
    pg.loc[row[0],'_combined.specific.non.tryptic.nterm'] = '\n'.join(row_nontryptic_nterm)
    pg.loc[row[0],'_combined.specific.non.tryptic.cterm'] = '\n'.join(row_nontryptic_cterm)
    prot_aln, prot_dend=clustalw(output+'/clustalw/prot/{}.fasta'.format(row[0]), all_prot_fasta)
    nuc_aln, nuc_dend = clustalw(output+'/clustalw/orf/{}.fasta'.format(row[0]), all_orfs_fasta)
    
    pg.loc[row[0],'_clustalw.orfs'] = nuc_aln
    pg.loc[row[0],'_clustalw.orfs.newick'] = nuc_dend
     
    pg.loc[row[0],'_clustalw.proteins'] = prot_aln
    pg.loc[row[0],'_clustalw.proteins.newick'] = prot_dend


newcols = [i for i in pg.columns if i.startswith('_')]
pg[newcols].to_csv(output+'/combined.csv')
