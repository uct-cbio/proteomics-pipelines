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
import pgfunctions

config = yaml.load(open(sys.argv[1]).read(),  Loader=yaml.Loader)

output = sys.argv[2]

peptides = pd.read_csv(config['mq_txt'] + '/peptides.txt', sep='\t',engine='python')

#evidence = pd.read_csv(config['mq_txt'] + '/evidence.txt', sep='\t',engine='python')


peptides = peptides[(peptides['Potential contaminant'].isnull()) & (peptides['Reverse'].isnull())]

peptides_set = set(peptides['Sequence'].tolist())

#print(config)
#for strain in config['strains']:
#    proteogenomics_out = config['outdir'] + '/strains/{}/{}_mapped_peptides.csv'.format(strain, strain)
#    pg = pd.read_csv(proteogenomics_out)
#    print(pg.head())
#    print(pg.columns)
#    print(strain)
#
#quit()


summary_list = []
reference_names = []
column_order  = []

strains = []


for reference in config['reference']:
    
    column_order  = []

    reference_names.append(reference)
    data = pd.read_csv(output +'/{}_combined.csv'.format(reference), sep='\t')
    #print(data.columns.tolist())
    stats = output +'/stats/{}'.format(reference)
    counts = {}
    try:
        os.makedirs(stats)
    except:
        shutil.rmtree(stats)
        os.makedirs(stats)


    def collist(df, col):
        peps = set()
        for pepset in df[col].dropna().tolist():
            for pep in pepset.split('\n'):
                peps.add(pep)
        return peps
        
    novel_peptides = collist(data,"Specific novel peptides - all strains")
    specific_peptides = collist(data,"Specific peptides - all strains")
    
    strain_all_peptides = defaultdict(set)
    combined_identified_peptides = set()
    exclusive_peptide_set = set()
    non_genomic_peptide_set = set() 
    frameshift_orf_pg_set = set()
    specific_protein_groups = set()
    paralogous_specific_protein_groups = set()
    all_protein_groups = set()

    for col in data.columns:
        rex = "Reference BLAST strain "
        if col.startswith(rex):
            strain = col.split(rex)[1]

        rex = 'All peptides strain '
        if col.startswith(rex):
            strain = col.split(rex)[1]
            if not strain in strains:
                strains.append(strain)
            
            


            strain_all_pg_list = [str(i) for i in list(data[col].dropna().index)] 
            coln = ('identified protein groups', strain)
            column_order.append(coln)
            counts[coln] = len(strain_all_pg_list)
            all_protein_groups.update(set(strain_all_pg_list))
            

            peptides = data[col].dropna().tolist()
            for pepset in peptides:
                for pep in pepset.split('\n'):
                    strain_all_peptides[strain].add(pep)
            
            st_peptides = strain_all_peptides[strain]
            
            combined_identified_peptides.update(st_peptides)

            coln = ('identified peptides', strain)
            column_order.append(coln)
            counts[coln] = len(st_peptides)
            
            w = open( stats +'/all.peptides.strain.{}.txt'.format(strain),'w')
            w.write('\n'.join(st_peptides))
            w.close()
            
            st_novel =  st_peptides & novel_peptides
            coln = ('novel peptides', strain)
            column_order.append(coln)
            counts[coln] = len(st_novel)

            w = open( stats +'/specific.novel.peptides.strain.{}.txt'.format(strain),'w')
            w.write('\n'.join(st_novel))
            w.close()
             
        
    


        rex = 'Specific peptides strain '
        if col.startswith(rex):
            strain = col.split(rex)[1]
            strain_pg = data[data[col].notnull()] 
            
            strain_pg_list = [str(i) for i in list(strain_pg.index)] 
            specific_protein_groups.update(set(strain_pg_list))

            w = open( stats +'/protein.groups.strain.{}.txt'.format(strain),'w')
            w.write('\n'.join(strain_pg_list))
            w.close()
            counts[('specific protein groups',strain)] = len(strain_pg)
            
            strain_pg_paralogous=strain_pg[strain_pg["Row paralogous specific peptides"].notnull()]
            strain_pg_paralogous_list  = [str(i) for i in list(strain_pg_paralogous.index)] 
            paralogous_specific_protein_groups.update(set(strain_pg_paralogous_list))

            counts[('paralogous specific protein groups',strain)] = len(strain_pg_paralogous)
            
            strain_pg_novel  = strain_pg[strain_pg["Reference BLAST strain "+strain].isnull()]
            counts[('novel protein groups',strain)] = len(strain_pg_novel)
            
            strain_pg_annotated  = strain_pg[strain_pg["Reference BLAST strain "+strain].notnull()]
            counts[('annotated protein groups',strain)] = len(strain_pg_annotated)

            peptides = collist(data, col)
            w = open( stats +'/specific.peptides.strain.{}.txt'.format(strain),'w')
            w.write('\n'.join(peptides))
            w.close()
            
            counts[('specific peptides', strain)] = len(peptides)
            
        rex = "Exclusive peptides strain "
        if col.startswith(rex):
            strain = col.split(rex)[1]
            peptides = collist(data, col)
    
            w = open( stats +'/exclusive.peptides.strain.{}.txt'.format(strain),'w')
            w.write('\n'.join(peptides))
            w.close()
            counts[('exclusive novel peptides',strain)] = len(peptides)
            
            exclusive_peptide_set.update(peptides)

        rex = "Non-genomic peptides strain " 
        if col.startswith(rex):
            strain = col.split(rex)[1]
            peptides = collist(data, col)

            w = open( stats +'/unmapped.peptides.strain.{}.txt'.format(strain),'w')
            w.write('\n'.join(peptides))
            w.close()
            counts[('non-genomic peptides',strain)] = len(peptides)
            
            non_genomic_peptide_set.update(peptides)
        
        rex='Frameshift validated strain '
        if col.startswith(rex):
            strain_pg = data[data[col].notnull()] 
            
            strain_pg_list = [str(i) for i in list(data[col].dropna().index)] 
            
            frameshift_orf_pg_set.update(set(strain_pg_list))

            strain = col.split(rex)[1]

            w = open( stats +'/frameshift.validated.strain.{}.txt'.format(strain),'w')
            w.write('\n'.join(strain_pg_list))
            w.close()
            counts[('frameshift validated ORFs',strain)] = len(strain_pg_list)
            

    coln = ('novel peptides', 'combined')
    column_order.append(coln)
    counts[coln] = len(novel_peptides)
    
    coln = ('identified peptides', 'combined')
    column_order.append(coln)
    counts[coln] = len(combined_identified_peptides)
    
    coln = ('specific peptides', 'combined')
    column_order.append(coln)
    counts[coln] = len(specific_peptides)
    
    coln = ('exclusive novel peptides', 'combined')
    column_order.append(coln)
    counts[coln] = len(exclusive_peptide_set)
    
    coln = ('non-genomic peptides', 'combined')
    column_order.append(coln)
    counts[coln] = len(non_genomic_peptide_set)
            
    counts[('paralogous specific protein groups','combined')]=len(list(paralogous_specific_protein_groups))
    counts[('specific protein groups','combined')] = len(list(specific_protein_groups))
    counts[('identified protein groups','combined')] = len(list(all_protein_groups))

    ref = data[data["Reference BLAST - all strains"].notnull()]
    #ref = data[data["Reference proteins mapped - all strains"].notnull()]
    

    counts[('annotated protein groups','combined')] = len(ref)
    ref_unmapped = data[~data["Reference BLAST - all strains"].notnull()]
    #ref_unmapped = data[~data["Reference proteins mapped - all strains"].notnull()]
    counts[('novel protein groups','combined')] = len(ref_unmapped)
    
    counts[('frameshift validated ORFs','combined')] = len(frameshift_orf_pg_set)
    summary_list.append(counts)
    #print(counts)

    strain_df = pd.DataFrame([counts], index=[reference]).transpose()

    unmapped_peptides = peptides_set - combined_identified_peptides
    
    w = open( stats +'/unmapped.peptides.{}.txt'.format(reference),'w')
    w.write('\n'.join(unmapped_peptides))
    w.close()
    


summary_df = pd.DataFrame(summary_list, index= reference_names).transpose()

###################
# identifications #
###################

summary_columns = ['identified peptides',
                   'specific peptides',
                   'non-genomic peptides',
                   'identified protein groups',
                   'specific protein groups',
                   'paralogous specific protein groups']
df_indx = []

for column in summary_columns:
    for strain in strains:
        indx = (column, strain)
        df_indx.append(indx)
    indx = (column, 'combined')
    df_indx.append(indx)
#print(config)
identification_df = summary_df.loc[df_indx,:]

tuples = list( identification_df.index)
index = pd.MultiIndex.from_tuples(tuples)
identification_df.index = index
identifications = identification_df[list(config['reference'].keys())[0]]
identifications = identifications.unstack(level=-1)
identifications = identifications.loc[summary_columns]
print("\n1) Identifications")
print(identifications)
identifications.to_csv( output + '/stats/identifications.csv')

###############
# annotations #
###############
print("\n2) Annotations")
df_indx=[]
annotation_columns = ['novel peptides',
                      'exclusive novel peptides',
                      'annotated protein groups',
                      'novel protein groups']

#                      'Annotated proteome ICDS putative']
#        'Annotated proteome ICDS validated',
#        'Strain ICDS putative',
#        'Strain ICDS validated',
#        'Strain ICDS sequencing error / programmed frameshift']


for column in annotation_columns:
    for strain in strains:
        indx = (column, strain)
        df_indx.append(indx)
    indx = (column, 'combined')
    df_indx.append(indx)
annotation_df = summary_df.loc[df_indx,:]

tuples = list( annotation_df.index)

index = pd.MultiIndex.from_tuples(tuples)
annotation_df.index = index
print(annotation_df)
annotation_df.to_csv( output + '/stats/annotations.csv')

################
# tss peptides #
################
print("\n3) TSS peptides")
strain_analyses = []
for strain in config['strains']:
    st_peptides = pd.read_csv(output + '/strains/{}/{}_mapped_peptides.csv'.format(strain, strain))
    
    analysis = pgfunctions.start_codon_analysis(st_peptides , strain)
    strain_analyses.append(analysis)

combined = pd.concat(strain_analyses)
tss_count_cols = []
tss_other_cols = []
for col in combined:
    if col.startswith("ORFs with "):
        combined[col]= combined[col].fillna(0)
        combined[col]= combined[col].astype(int)
        tss_count_cols.append(col)
    else:
        tss_other_cols.append(col)
tss_count_cols.sort(reverse=False, key=lambda x : int(x.split("ORFs with ")[1].split()[0]))
combined_cols = tss_other_cols + tss_count_cols
combined = combined[combined_cols]

print(combined)
combined.to_csv(output+ '/stats/tss.csv')
