#!/usr/bin/env python4

import pprint
from functools import reduce
import pandas as pd
import scipy
import pandas as pd
import numpy as np
import os 
import shutil
import scipy.stats
import scikit_posthocs as ph
import numpy as np
import Bio
import scipy.stats
import scikit_posthocs as ph
import numpy as np
import matplotlib_venn
from matplotlib_venn import venn3, venn2
import matplotlib.pyplot as plt
from pylab import *
import matplotlib.pyplot as plt
from pylab import *
#import skbio
import pandas as pd
import Bio
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from collections import defaultdict
from collections import Counter
import os
import matplotlib.pyplot as plt
from pylab import *
import pandas as pd
#from IPython.display import Image
#from IPython.display import HTML
from matplotlib_venn import venn3, venn2
from natsort import realsorted
import operator

# Get a representative sequence for each ORF cluster, using shortest length as a heuristic
def get_orf_cluster(table):
    table['ORF_length'] = table['ORF_sequence'].str.len()
    table = table.sort_values(['ORF_length'])
    table = table.drop_duplicates(['Orthologous_ORFs','Peptide_sequence'], keep='first')
    return table

def tss(table):
    potential_starts = ['ATG','GTG','TTG']
    if (table['Peptide_first_codon'] in ['GTG' or 'TTG']) \
    and (table['Peptide_amino_acid_first'] == 'M') and (table['Strain_identified'] =='+' ):
        return 'Methionine translation of non-ATG'
    elif (table['Peptide_previous_codon'] in potential_starts) \
    and not (table['Peptide_first_codon'] in potential_starts) and (table['Strain_identified'] =='+' ):
        return 'MAP cleaved start codon'
    elif (table['Peptide_first_codon'] in potential_starts) and (table['Peptide_tryptic_nterm']== 'False') and (table['Strain_identified'] =='+' ):
        return 'Non-tryptic start codon'
    elif (table['Peptide_first_codon'] in potential_starts) and (table['Peptide_previous_codon'] in potential_starts) and (table['Strain_identified'] =='+' ):
        return 'Ambiguous: peptide previous codon and peptide first codon are start codons'
    elif (table['Peptide_first_codon'] in potential_starts) and (table['Peptide_tryptic_nterm']== 'True') and (table['Strain_Nterm_Acetylated'] == '+') and (table['Strain_identified'] =='+' ):
        return 'Tryptic start codon {} with n-term acetylation'.format(table['Peptide_first_codon'])
    elif (table['Strain_Nterm_Acetylated'] == '+') and (table['Strain_identified'] =='+' ):
        return "Other n-term acetylation"

def start_codons(table):
    if table['StartType'] == 'MAP cleaved start codon':
        return table['Peptide_previous_codon']
    elif table['StartType'] == 'Non-tryptic start codon':
        return table['Peptide_first_codon']
    elif table['StartType'] == 'Methionine translation of non-ATG':
        return table['Peptide_first_codon']
    elif table['StartType'] is None:
        return None
    else:
        print(table['StartType'])
        assert 4==5
        
def evaluate_starts(start, offset=0):
    
    try:
        start = int(start)
        start -= offset
        return int(start)
    except:
        return "ambiguous"
        
def start_position(table):
    if table['StartType'] == 'MAP cleaved start codon':
        return evaluate_starts(table['Peptide_starts'], offset=-1)
    elif table['StartType'] == 'Non-tryptic start codon':
        return evaluate_starts(table['Peptide_starts'])
    elif table['StartType'] == 'Methionine translation of non-ATG':
        return evaluate_starts(table['Peptide_starts'])
    elif table['StartType'] is None:
        return None
    else:
        print(table['StartType'])
        assert 4==5
        
def nterm_acetylation_MAP(table):
    if table['Strain_Nterm_Acetylated'] != '+' and table['StartType'] == 'MAP cleaved start codon':
        return 'Nterm not acetylated MAP TSS peptides'
    elif table['Strain_Nterm_Acetylated'] == '+' and table['StartType'] == 'MAP cleaved start codon':
        return 'Nterm acetylated MAP TSS peptides'
    
def nterm_acetylation_iMet(table):
    if table['Strain_Nterm_Acetylated'] != '+' and table['StartType'] != 'MAP cleaved start codon':
        return 'Nterm not acetylated iMet TSS peptides'
    elif table['Strain_Nterm_Acetylated'] == '+' and table['StartType'] != 'MAP cleaved start codon':
        return 'Nterm acetylated iMet TSS peptides'

def nterm_acetylation_all(table):
    if table['Strain_Nterm_Acetylated'] != '+':
        return 'Nterm not acetylated TSS peptides'
    elif table['Strain_Nterm_Acetylated'] == '+':
        return 'Nterm acetylated TSS peptides'


def start_site_count_map(table):
    start_map = {}
    for name, group in table.groupby('ORF_id'):
        group = group[group['StartType'].apply(str) != 'None' ]
        group = group[group['StartPosition'].apply(str) != 'None' ]
        unique_starts = group.drop_duplicates('StartPosition')
        start_map[name]=len(unique_starts)
    return start_map

def orf_protein_map(table):
    protein_map = {}
    for name, group in table.groupby('ORF_id'):
        sequences = group['Peptide_sequence']
        reference_ids = []
        for sequence in sequences:
            if sequence in peptide2protein:
                id_list = peptide2protein[sequence]
                reference_ids += id_list
        reference_ids = list(set(reference_ids))
        #group = group[group['StartType'].apply(str) != 'None' ]
        #unique_starts = group.drop_duplicates('Peptide_starts')
        protein_map[name]=reference_ids
    return protein_map

def start_codon_analysis(df, strain='Strain'):
    table = df.copy()
    print(table.head()) 
    analysis = pd.DataFrame()
    table['StartType'] = table.apply(tss, axis=1)
    #temp_table = table[table['Peptide_starts'].str.contains(';') == True ]
    
    #print(temp_table[:1].stack())
    #return
    table = table[table['Strain_identified'] == '+']
    table = table[table['Peptide_translated_ORF_cluster_specific'] == '+']
    #table = table[table['Peptide_translated_ORF_cluster_count'] == 1]
    
    #all_orf_table = table.copy() # keep the paralogs for the start position count analysis
    table = get_orf_cluster(table)

    
    #table['Peptide_mapped_reference_proteins'] = table['Peptide_sequence'].map(peptide2protein)
    
    table = table[table['StartType'].notnull()]
    
   
    # characterize TSS peptides counting each peptide only once
    nr_table = table.drop_duplicates('Peptide_sequence', keep='first')
    #assert len(table.drop_duplicates('Peptide_sequence', keep='first')) == len(table.drop_duplicates(subset=['Peptide_sequence','StartType'], keep='first'))
    identified = Counter(nr_table['StartType'])

    
    total = 0
    for key in list(identified.keys()):
        if not key is None:
            analysis.loc[strain, key] = int(identified[key])
            analysis[key] = analysis[key].astype('int')
            total += identified[key]
    analysis.loc[strain,'Total TSS peptides'] = int(total)
    analysis['Total TSS peptides'] = analysis['Total TSS peptides'].astype(int)
    analysis.loc[strain,'MAP percentage'] = np.round(identified['MAP cleaved start codon']/total * 100,2)
    table['StartCodon'] = table.apply(start_codons, axis=1)
    table['StartPosition'] = table.apply(start_position, axis=1)
    codon_table = table.drop_duplicates(['StartPosition', 'ORF_id'])
    codons = Counter(codon_table['StartCodon'])
    del codons[None] 
    for key in list(codons.keys()):
        analysis.loc[strain, 'Start codon count ' + key] = int(codons[key])
        analysis.loc[strain, 'Start codon count ' + key + ' (%)'] = np.round((codons[key] / total ) * 100 ,2)
    
    # counting acetylated start peptides
    table['AcetylatedStartAll'] =  table.apply(nterm_acetylation_all, axis=1)
    acetyl_all = Counter(table['AcetylatedStartAll'])
    del acetyl_all[None]
    table['AcetylatedStartMAP'] =  table.apply(nterm_acetylation_MAP, axis=1)
    acetyl_map = Counter(table['AcetylatedStartMAP'])
    del acetyl_map[None]
    table['AcetylatedStartMet'] =  table.apply(nterm_acetylation_iMet, axis=1)
    acetyl_imet = Counter(table['AcetylatedStartMet'])
    del acetyl_imet[None]
    print(acetyl_imet)
    
    int_keys = []
    for key in list(acetyl_all.keys()):
        analysis.loc[strain, key] = int(acetyl_all[key])
        int_keys.append(key)
    for key in list(acetyl_map.keys()):
        analysis.loc[strain, key] = int(acetyl_map[key])
        int_keys.append(key)
    for key in list(acetyl_imet.keys()):
        analysis.loc[strain, key] = int(acetyl_imet[key])
        int_keys.append(key)
    for key in int_keys:
        analysis[key] = analysis[key].astype(int)

    start_site_map = start_site_count_map(codon_table)
    codon_table['ORF_start_count'] = codon_table['ORF_id'].map(start_site_map)
    
    #protein_map = orf_protein_map(table)
    #print(protein_map)
    #table['ORF_mapped_reference_proteins'] = table['ORF_id'].map(protein_map)
    #table['ORF_mapped_reference_proteins'] = table['ORF_mapped_reference_proteins'].apply(lambda x : ';'.join(x))
    orf_table = codon_table.drop_duplicates('ORF_id') # got to here, are we using TSS counts for orfs even belonging to ORF clusters? or maybe take the shorters ORF of a given cluster? Or ORFs that only occur once?
    start_counts = Counter(orf_table['ORF_start_count'])
    for key in list(start_counts.keys()):
        analysis.loc[strain, 'ORFs with ' + str(key) + ' TSSs'] = int(start_counts[key])
    return analysis

