#!/usr/bin/env python

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
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from Bio import Phylo
import contextlib
import sys
from collections import Counter
from io import StringIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# Create two sequence files
from functools import lru_cache

def pairwise_blast(query, target, temp_folder):
    SeqIO.write(query, temp_folder +"/query.fasta", "fasta")
    SeqIO.write(target, temp_folder +"/target.fasta", "fasta")
    output = NcbiblastpCommandline(query=temp_folder +"/query.fasta", subject=temp_folder +"/target.fasta", outfmt=5)()[0]
    blast_result_record = NCBIXML.read(StringIO(output))    
    results = []
    # Print some information on the result
    for alignment in blast_result_record.alignments:
        for hsp in alignment.hsps:
            #print(dir(hsp))
            results.append('****Alignment****')
            results.append('sequence:         ' + str(alignment.title))
            results.append('length:           ' + str(alignment.length))
            results.append('e value:          ' + str(hsp.expect))
            results.append('hsp query:        ' + str(hsp.query))
            results.append('hsp match:        ' + str(hsp.match))
            results.append('hsp subject:      ' + str(hsp.sbjct))
            results.append('hsp align_length: ' + str(hsp.align_length))
            #results.append('hsp bits:         ' + str(hsp.bits))
            #results.append('hsp frame:        ' + str(hsp.frame))
            results.append('hsp gaps:         ' + str(hsp.gaps))
            #results.append('hsp identities:   ' + str(hsp.identities))
            results.append('hsp sbjct_start:  ' + str(hsp.sbjct_start))
            results.append('hsp sbjct_end:    ' + str(hsp.sbjct_end))
            results.append('hsp score:        ' + str(hsp.score))
            break 
    results = '\n'.join(results)
    os.remove(temp_folder +'/query.fasta')
    os.remove(temp_folder +'/target.fasta')
    return results

def return_ranked_best_taxa(peptides, reference_lineage, taxon_peptides, taxon_data, reference_list):
    mapped = []
    mapped_fastas = {}
    for peptide in peptides:
        if peptide in taxon_peptides:
            entries = taxon_peptides[peptide]
            mapped += entries
    if len(mapped) > 0:
        counted = Counter(mapped)
        biggest = max(counted.values())
        overlap_max = defaultdict(list)
        for key in counted:
            if counted[key] == biggest:
                dict = taxon_data[key]
                val = '{} - {} Taxid: {}'.format(dict['Entry name'] ,dict['Organism'],dict['Organism ID'])
                
                seq = Seq(dict['Sequence'])
                rec = SeqRecord(seq = seq, id = dict['Entry name'], description = dict['Organism'])
                mapped_fastas[dict['Entry name']] = rec 
                lineage = dict['Taxonomic lineage IDs (all)'].split(', ')
                lineage  = [int(i) for i in lineage]
                
                overlap  = set(lineage) & set(reference_lineage)        
                overlap_max[len(overlap)].append(val)
        new_biggest = max(overlap_max.keys())
        newmapped = '\n'.join(overlap_max[new_biggest])

        spec_fastas = []
        for datum in overlap_max[new_biggest]:
            d = datum.split()[0]
            fs = mapped_fastas[d]
            spec_fastas.append(fs)
        
        max_ref = []
        for ref in reference_list:
            r = ref.split('|')[1]
            refc = counted[r]
            max_ref.append(refc)
        if len(max_ref) > 0:
            max_ref = max(max_ref)
        else:
            max_ref = 0
        return newmapped, biggest, max_ref , spec_fastas




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

def map2fasta(ids, mapping):
    recs = []
    taxonomy = None
    for id in ids:
        header = id
        id = id.split('|')[1]
        data = mapping[id]
        seq = data['Sequence']
        organism = data['Organism']
        name = data['Protein names']
        description = "{} {}".format(organism, name)
        rec = SeqRecord(seq = Seq(seq), id = header, description = description)
        recs.append(rec)
    return recs

def check_for_met(codons):
    met = False
    for codon  in codons.split(';'):
        t = translate(Seq(codon), table = config.translation_table, cds=False)
        if t == 'M':
            met = True
    return met

def non_atg_m(map):
    mstarts = map[map['Peptide_amino_acid_first'] =='M']
    mstarts['MetCodon'] = mstarts['Peptide_first_codon'].apply(check_for_met)
    nonatgm = mstarts[mstarts['MetCodon'] == False]
    return nonatgm['Peptide_sequence'].tolist()

# Get the mapping of genomics locations of identified peptides for reference
reference_map = pd.read_csv(output +'/reference/{}_mapped_peptides.csv'.format(config.reference_taxid))
reference_nontryptic_nterm=set(reference_map[reference_map['Peptide_tryptic_nterm']=='False']['Peptide_sequence'])
reference_nontryptic_cterm=set(reference_map[reference_map['Peptide_tryptic_cterm']=='False']['Peptide_sequence'])
reference_non_atg_starts = set(non_atg_m(reference_map))

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

    strain_nontryptic_nterm=set(st[st['Peptide_tryptic_nterm']=='False']['Peptide_sequence'])
    strain_nontryptic_cterm=set(st[st['Peptide_tryptic_cterm']=='False']['Peptide_sequence'])
    strain_non_atg_starts = set(non_atg_m(st))

    reference_nontryptic_nterm = reference_nontryptic_nterm & strain_nontryptic_nterm
    reference_nontryptic_cterm = reference_nontryptic_cterm & strain_nontryptic_cterm
    reference_non_atg_starts = reference_non_atg_starts & strain_non_atg_starts


non_specific_peptides += get_ns(reference_map)
reference_set = set(reference_map['Peptide_sequence'].tolist())
specific_peptides += get_s(reference_map)

# Combined specific set (exclide all sequences that are non-specific in any group
global_specific_peptides = set(specific_peptides) - set(non_specific_peptides)
annotated_peptides = set(reference_peptides.keys())
exclusive_peptides = annotated_peptides - reference_set # in reference proteome but not reference genome ie. annotated frameshift proteins from swissprot
novel_peptides = global_specific_peptides - annotated_peptides
global_non_atg_starts = reference_non_atg_starts
global_nontryptic_cterm = reference_nontryptic_cterm
global_nontryptic_nterm = reference_nontryptic_nterm



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

#pg = pg[pg['id']==934]

def return_fasta(map, specific):
    orfs_fasta=[]
    prots_fasta=[]
    
    identical_translated= False

    specific_map=map[map['Peptide_sequence'].isin(specific)].copy() 
    specific_map = specific_map.drop_duplicates(['ORF_ids'])

    orf_id_list =specific_map['ORF_ids'].tolist()
    orf_nucs_list =specific_map['ORF_sequence'].tolist()
    orf_trans_list =specific_map['ORF_translation'].tolist()
    
    if len(set(orf_trans_list)) == 1:
        identical_translated=True

    for item in range(len(orf_id_list)):
        
        orf_id = orf_id_list[item]
        orf_nucs = orf_nucs_list[item]
        orf_trans = orf_trans_list[item]

        nucs=SeqRecord(id='|'.join(orf_id.split('|')[1:]),seq=Seq(orf_nucs))
        trans=SeqRecord(id='|'.join(orf_id.split('|')[1:]),seq=Seq(orf_trans))
        
        orfs_fasta.append(nucs)
        prots_fasta.append(trans)
    
    if identical_translated == True: # In case the translated ORF sequence is identical
        prot_ids = []
        for rec in prots_fasta:
            prot_ids.append(rec.id)
        prots_fasta = [SeqRecord(id = ';'.join(prot_ids), seq=rec.seq)]
    
    return orfs_fasta, prots_fasta

@contextlib.contextmanager
def stdout_redirect(where):
    sys.stdout = where
    try:
        yield where
    finally:
        sys.stdout = sys.__stdout__

def clustalw(output_file, fasta):    
    new_fasta = fasta.copy()
    for rec in new_fasta:
        new_ids = []
        ids = rec.id.split(';')
        for id in ids:
            new_ids.append(id.split('|')[0])
        new_ids = ';'.join(new_ids)
        rec.id = new_ids
    SeqIO.write(new_fasta, output_file, 'fasta')
    cline = ClustalwCommandline("clustalw2", infile=output_file)
    stdout, stderr = cline()
    aln = output_file.split('.fasta')[0] 
    align = AlignIO.read(aln +'.aln', "clustal")
    #align_str = open(aln +'.aln').read()
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
    return str(align), tree

def list_trie_upper(fasta, peptide_list):
    new_fasta = fasta.copy()
    peptide_list = list(peptide_list)
    upper = []
    pep_trie = algo.trie_graph(peptide_list)
    for rec in new_fasta:
        id = rec.id
        temp_seq = str(rec.seq)
        description= rec.description
        temp_seq = algo.trie_upper(pep_trie, temp_seq)
        newrec = SeqRecord(seq = Seq(temp_seq), id = id, description = description)
        upper.append(newrec.format('fasta'))
    return upper

def muscle(fasta):
    new_fasta = fasta.copy()
    for rec in new_fasta:
        new_ids = []
        ids = rec.id.split(';')
        for id in ids:
            new_ids.append(id.split('|')[0])
        new_ids = ';'.join(new_ids)
        rec.id = new_ids

    cline = MuscleCommandline(clwstrict=True)
    child = subprocess.Popen(str(cline),
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    universal_newlines=True,
    shell=(sys.platform!="win32"))
    SeqIO.write(new_fasta, child.stdin, "fasta")
    child.stdin.close()
    align = AlignIO.read(child.stdout, "clustal")
    lines=[]
    for record in align:
        ln="{} - {}".format(record.seq, record.id)
        lines.append(ln)
    lines = '\n'.join(lines)
    return lines

for row in pg.iterrows():
    print(row[0])
    peptide_ids = row[1]['Peptide IDs'].split(';')  # Get the ids of peptides in the row
    evidence_ids = row[1]['Evidence IDs'].split(';') # Get the ids of envidence.txt entries
    row_peps = peptides[peptides['id'].apply(str).isin(peptide_ids)]
    row_evs = evidence[evidence['id'].apply(str).isin(evidence_ids)]
    # Divide identified peptides into categories
    
    all_row_peptides = row_peps['Sequence'].tolist()
    row_specific  = set(all_row_peptides) - set(non_specific_peptides)  # All peptides that are not non-specific in any genome
    row_exclusive = row_specific & exclusive_peptides # All annotated peptides that are not mapped to the reference genome (excl non_specific)
    row_novel     = row_specific & novel_peptides
    row_annotated = row_specific & annotated_peptides
    row_nontryptic_cterm = row_specific & global_nontryptic_cterm  
    row_nontryptic_nterm = row_specific & global_nontryptic_nterm  
    row_non_atg_starts = row_specific & global_non_atg_starts 

    strain_peps=defaultdict(list)
    group_evs=defaultdict(list)
    for sample in samples:
        strain=samples[sample]['STRAIN']
        group=samples[sample]['GROUP']
        sample_peps = list(set(row_peps[row_peps['Experiment {}'.format(sample)]>=1]['Sequence'].tolist()))
        sample_evs = list(set(row_evs[row_evs['Experiment'] == sample]['Modified sequence'].tolist()))
        #pg.loc[row[0], "_all.peptides.sample.{}".format(sample)] = '\n'.join(list(sample_evs))
        strain_peps[strain] += sample_peps
        group_evs[group] += sample_evs
    #for group in group_evs:
        #pg.loc[row[0],"_all.peptides.group.{}".format(group)]='\n'.join(list(set(group_evs[group])))
    all_orfs_fasta=[]
    all_prot_fasta=[]
    fasta_holder = {}
    id_holder={}
    frameshifts =[]
    # Add reference fasta sequences
    x,y=return_fasta(reference_map, row_specific)
    all_orfs_fasta += x
    all_prot_fasta += y
    fasta_holder['Reference_prots'] = y
    fasta_holder['Reference_nucs'] = x
    id_holder['Reference'] = [rec.id for rec in y]
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
        trie = list_trie_upper(y, st_peps)
        fasta_holder['{}_prots'.format(strain)] = y
        fasta_holder['{}_nucs'.format(strain)] = x
        id_holder[strain] =[rec.id for rec in y]
        for comparison_strain in strain_peps:
            if comparison_strain != strain:
                comp_peps=set(strain_peps[comparison_strain])
                strain_exclusive_novel = strain_exclusive_novel - comp_peps
        if len(strain_exclusive_novel) > 0:
            strains_exclusive.append(strain)
        pg.loc[row[0], "_all.peptides.strain.{}".format(strain)]='\n'.join(st_peps)
        pg.loc[row[0], "_exclusive.peptides.strain.{}".format(strain)]='\n'.join(strain_exclusive_novel)
        blasted = []
        count=1
        for p in strain_exclusive_novel:
            ref_recs = fasta_holder['Reference_prots']
            id = '{}_exclusive_peptide_{}'.format(strain, str(count))
            seq = Seq(p)
            rec = SeqRecord(id = id, seq = seq)
            out = pairwise_blast(rec, ref_recs, output)
            blasted.append(out)
            count += 1
        
        pg.loc[row[0], "_exclusive.peptides.blast.strain.{}".format(strain)]='\n'.join(blasted) 
        pg.loc[row[0], "_unmapped.peptides.strain.{}".format(strain)]='\n'.join(genome_unmapped) 
        pg.loc[row[0], "_translated.orfs.strain.{}".format(strain)] = '\n'.join(trie)
    for key in id_holder:
        ids = id_holder[key]
        if len(ids) > 1:
            frameshifts.append(key) 
    ref_recs = fasta_holder['Reference_prots']
    trie = list_trie_upper(ref_recs, all_row_peptides)
    pg.loc[row[0],"_translated.orfs.reference"] ='\n'.join(trie)
    pg.loc[row[0],"_frameshift"] = '\n'.join(list(frameshifts))
    pg.loc[row[0],"_annotated.frameshift.peptides"] = '\n'.join(list(row_exclusive))
    pg.loc[row[0],"_exclusive.peptide.strains"] = '\n'.join(list(strains_exclusive))
    pg.loc[row[0],"_combined.specific.peptides"] ='\n'.join(list(row_specific))
    pg.loc[row[0],"_combined.specific.annotated.peptides"] ='\n'.join(row_annotated)
    pg.loc[row[0],"_combined.specific.novel.peptides"] ='\n'.join(row_novel)
    mapped_reference = set()
    for peptide in row_annotated:
        mapped_reference.update(reference_peptides[peptide])
    
    refps =  map2fasta(list(mapped_reference), taxon_data)
    refps =  list_trie_upper(refps, all_row_peptides)
    
    pg.loc[row[0],'_reference.proteins.mapped']= '\n'.join(refps)
    pg.loc[row[0],'_reference.proteins.mapped.count'] = len(mapped_reference)
    pg.loc[row[0],'_reference.entries.mapped'] = ';'.join([i.split('|')[1] for i in mapped_reference])
    
    pg.loc[row[0],'_combined.non.tryptic.nterm.peptides'] = '\n'.join(row_nontryptic_nterm)
    pg.loc[row[0],'_combined.non.tryptic.cterm.peptides'] = '\n'.join(row_nontryptic_cterm)
    pg.loc[row[0],'_combined.non.atg.Met.start.peptides'] = '\n'.join(row_non_atg_starts)
    
    try:
        mapped_taxa, max_pep, max_ref, spec_fastas= return_ranked_best_taxa(row_specific, config.reference_taxonomic_lineage, taxon_peptides, taxon_data, mapped_reference)
        pg.loc[row[0],'_reference.proteome.best.peptide.count'] = max_ref
        pg.loc[row[0],'_taxon.best.matches'] = mapped_taxa 
        pg.loc[row[0],'_taxon.best.peptide.count'] = max_pep
        pg.loc[row[0],'_taxon.best.match.fasta'] = '\n'.join(list_trie_upper([spec_fastas[0]],all_row_peptides)) 
        all_prot_fasta += spec_fastas
    except:
        pass

    if len(all_prot_fasta) > 0:
        prot_aln,prot_dend=clustalw(output+'/clustalw/prot/{}.fasta'.format(row[0]), all_prot_fasta) 
        prot_muscle=muscle(all_prot_fasta)
        pg.loc[row[0],'_clustalw.proteins'] = prot_aln
        pg.loc[row[0],'_clustalw.proteins.newick'] = prot_dend
        pg.loc[row[0],'_muscle.proteins'] = prot_muscle
        print(prot_dend)
    if len(all_orfs_fasta) > 0:
        nuc_aln, nuc_dend = clustalw(output+'/clustalw/orf/{}.fasta'.format(row[0]), all_orfs_fasta)
        nuc_muscle=muscle(all_orfs_fasta)
        pg.loc[row[0],'_clustalw.orfs'] = nuc_aln
        pg.loc[row[0],'_clustalw.orfs.newick'] = nuc_dend
        pg.loc[row[0],'_muscle.orfs'] = nuc_muscle
        print(nuc_dend)

newcols = [i for i in pg.columns if i.startswith('_')]
pg.to_csv(output+'/combined.csv')

master = pg
#print(master.stack())

# Export Frameshifts
fs_anno = master[((master['_frameshift'].notnull()) & (master["_combined.specific.novel.peptides"]=='') & (master['_reference.proteins.mapped.count']==1))]

fs_novel = master[~((master['_frameshift'].notnull()) & (master["_combined.specific.novel.peptides"].apply(str) =='') & (master['_reference.proteins.mapped.count']==1))]

fs_anno.to_csv(output+'/validated_ICDSs.csv')
fs_novel.to_csv(output+'/novel_ICDSs.csv')



